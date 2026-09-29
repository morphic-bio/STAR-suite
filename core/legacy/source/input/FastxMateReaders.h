#ifndef CODE_input_FastxMateReaders
#define CODE_input_FastxMateReaders

// Per-mate FASTX reader threads for STAR's standard chunk reader.
//
// Each mate's stream (the FIFO or file STAR already opens for it) is parsed
// on its own thread, in order from the start. The parsed fields go into
// record batches of kFastxMateBatchRecords records, counted from the start of
// each lane file, so record i of lane L pairs across mates by count alone: the
// readers exchange no offsets or positions. A mapping thread holding STAR's
// input lock pairs the records and writes the same chunk text as the
// single-threaded loop in ReadAlignChunk_processChunks.cpp, so everything
// downstream sees identical input.
//
// The parse repeats the single-threaded loop's istream calls, in the same
// order, on each mate's stream, through a stream buffer that gives back the
// reader's compute permit while it waits for input bytes. A reader never
// holds a permit while it waits for input or for queue space.

#include "input/BgzfRangeReader.h"  // BgzfWorkPermitHooks

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <istream>
#include <memory>
#include <mutex>
#include <streambuf>
#include <string>
#include <thread>
#include <vector>

namespace star {
namespace input {

static constexpr uint32_t kFastxMateBatchRecords = 2048;
static constexpr size_t kFastxMateBatchesInFlight = 32;
static constexpr size_t kFastxMateReadBufferBytes = size_t(1) << 20;
static constexpr uint32_t kFastxMateMaxMates = 3;

// STAR's line limits (IncludeDefine.h), passed in so this module does not
// depend on STAR's global headers.
struct FastxMateLimits {
    long long nameSeqLineMax = 0;  // DEF_readNameSeqLengthMax
    long long seqLineMax = 0;      // DEF_readSeqLengthMax
};

// One record of one mate. Offsets index the batch arena. seq/qual hold the
// exact bytes the single-threaded loop copies for that line, including the
// trailing newline (FASTA: the joined sequence lines plus the newline).
struct FastxMateRecord {
    uint32_t idOff = 0, idLen = 0;        // FASTQ mate 0: ID token; FASTA: own token
    uint32_t extraOff = 0, extraLen = 0;  // FASTQ header text after the ID token
    uint32_t seqOff = 0, seqLen = 0;
    uint32_t qualOff = 0, qualLen = 0;    // FASTQ only
    char filter = 'N';                    // FASTQ mate 0: Illumina filter flag
    char format = '@';                    // '@' FASTQ, '>' FASTA
};

enum class FastxMateEnd : uint8_t {
    None,   // batch full; more records of the same lane follow
    Lane,   // a FILE marker follows: nextLane
    Input,  // end of this mate's input
    Error   // unreadable record start (mate 0) or reader failure
};

struct FastxMateBatch {
    int lane = 0;
    std::vector<FastxMateRecord> records;
    std::vector<char> arena;
    FastxMateEnd end = FastxMateEnd::None;
    int nextLane = 0;           // End::Lane
    int endChar = -1;           // End::Input: character peeked at the record start
    bool endSilent = false;     // End::Input: the stream had already failed before the record start
    uint64_t laneRecords = 0;   // records of this lane read so far, this batch included
    std::string errorWord;      // End::Error, mate 0: first word of the bad line
    std::string errorRest;      // End::Error, mate 0: rest of the bad line
    std::string errorText;      // End::Error: reader failure message

    void reset(int laneIn) {
        lane = laneIn;
        records.clear();
        arena.clear();
        end = FastxMateEnd::None;
        nextLane = 0;
        endChar = -1;
        endSilent = false;
        laneRecords = 0;
        errorWord.clear();
        errorRest.clear();
        errorText.clear();
    }
    uint32_t append(const char* data, size_t size) {
        const uint32_t at = static_cast<uint32_t>(arena.size());
        arena.insert(arena.end(), data, data + size);
        return at;
    }
    const char* at(uint32_t offset) const { return arena.data() + offset; }
};

// What the chunk filler needs from STAR, so the harness can drive the same
// code with its own counters and callbacks.
struct FastxChunkFillContext {
    unsigned long long chunkInSizeBytes = 0;  // P.chunkInSizeBytes: stop once mate 0 or 1 reaches it
    unsigned long long chunkArrayBytes = 0;   // P.chunkInSizeBytesArray: buffer size per mate
    unsigned long long readMapNumber = static_cast<unsigned long long>(-1);
    bool fastqReadIdNumber = false;  // P.outSAMreadIDnumber
    bool fastaReadIdNumber = false;  // P.outSAMreadID == "Number"
    int thread = 0;                  // for the end-of-stream log line
    unsigned long long* iReadAll = nullptr;
    int* readFilesIndex = nullptr;
    // The single-threaded loop's Log.out lines.
    std::function<void(int lane)> onLaneStart;
    std::function<void(const std::string& line)> onLog;
    std::function<void(const std::string& text)> onWarning;
    // Mate 0 has a record start that is neither a read nor a FILE marker.
    // STAR reports it with the single-threaded loop's message; must not return.
    std::function<void(const std::string& word, const std::string& rest)> onBadRecordStart;
    // Any other fatal input error; must not return in STAR.
    std::function<void(const std::string& text)> onFatal;
};

class FastxMateReader {
public:
    struct Stats {
        uint64_t records = 0;
        uint64_t batches = 0;
        uint64_t bytes = 0;              // input bytes read from the stream
        uint64_t parseNs = 0;            // parsing, input waits excluded
        uint64_t inputWaitNs = 0;        // blocked reading the stream
        uint64_t freeWaitNs = 0;         // blocked on a full queue (back-pressure)
        uint64_t consumerWaitNs = 0;     // the chunk filler waited for this mate
        uint64_t permitAcquires = 0;
        uint64_t permitWaitNs = 0;
    };

    FastxMateReader(uint32_t mate, std::istream* source, int initialLane,
                    const BgzfWorkPermitHooks& hooks, const FastxMateLimits& limits);
    ~FastxMateReader();
    FastxMateReader(const FastxMateReader&) = delete;
    FastxMateReader& operator=(const FastxMateReader&) = delete;

    void start();
    void requestStop();
    void join();

    // Consumer side (called under STAR's input lock). popReady returns
    // nullptr only when the reader stopped without a final batch.
    FastxMateBatch* popReady();
    void recycle(FastxMateBatch* batch);

    // Valid after join().
    Stats stats() const;

private:
    class InputBuf : public std::streambuf {
    public:
        InputBuf(FastxMateReader* owner, std::streambuf* source);
    protected:
        int_type underflow() override;
    private:
        FastxMateReader* owner_;
        std::streambuf* source_;
        std::vector<char> buffer_;
    };

    void run();
    bool fillBatch(FastxMateBatch* batch, std::istream& in);
    bool parseFastq(FastxMateBatch* batch, std::istream& in, FastxMateRecord* record,
                    const std::string& token);
    void parseFasta(FastxMateBatch* batch, std::istream& in, FastxMateRecord* record,
                    const std::string& token);
    bool appendLine(FastxMateBatch* batch, std::istream& in, uint32_t* offset, uint32_t* length);
    void markLane(FastxMateBatch* batch, int nextLane);
    FastxMateBatch* takeFree();
    void pushReady(FastxMateBatch* batch);
    void acquirePermit();
    void releasePermit();
    // Reports the queue to the permit allocator (hooks.observe); mutex_ held.
    void observeLocked();

    uint32_t mate_;
    std::istream* source_;
    int lane_;
    uint64_t laneRecords_ = 0;
    BgzfWorkPermitHooks hooks_;
    FastxMateLimits limits_;
    std::vector<char> line_;

    std::thread thread_;
    std::atomic<bool> stop_{false};
    bool started_ = false;

    mutable std::mutex mutex_;
    std::condition_variable readyCv_, freeCv_;
    std::deque<FastxMateBatch*> ready_;
    std::vector<FastxMateBatch*> free_;
    std::vector<std::unique_ptr<FastxMateBatch>> owned_;
    bool finished_ = false;         // the reader thread has exited
    bool filling_ = false;          // the reader is filling a batch
    bool consumerWaiting_ = false;  // the chunk filler is waiting for a batch

    // Reader-thread state.
    FastxMateBatch* current_ = nullptr;
    bool holding_ = false;
    uint64_t holdWaitNs_ = 0;
    uint64_t holdStartNs_ = 0;
    size_t holdArenaStart_ = 0;

    Stats stats_;
};

class FastxMateReaderGroup {
public:
    // streams[m] is mate m's input stream; it must stay open until
    // stopAndJoin() returns. initialLane is P.readFilesIndex at open.
    FastxMateReaderGroup(const std::vector<std::istream*>& streams, int initialLane,
                         const BgzfWorkPermitHooks& hooks, const FastxMateLimits& limits);
    ~FastxMateReaderGroup();
    FastxMateReaderGroup(const FastxMateReaderGroup&) = delete;
    FastxMateReaderGroup& operator=(const FastxMateReaderGroup&) = delete;

    uint32_t mates() const { return static_cast<uint32_t>(readers_.size()); }

    // Starts the reader threads on first use. A stopped group stays stopped
    // and reports the end of input.
    void ensureStarted();
    // Stops, wakes and joins the reader threads. Safe to call repeatedly.
    void stopAndJoin();
    bool started() const { return started_; }

    // Fills one chunk (chunkIn[m], totals[m]; totals has at least two
    // elements and starts at zero) exactly as the single-threaded loop does.
    // Call only under STAR's input lock.
    void fillChunk(char* const* chunkIn, unsigned long long* totals,
                   FastxChunkFillContext& context);

    std::string summary() const;

private:
    struct Cursor {
        FastxMateBatch* batch = nullptr;
        size_t index = 0;
        bool dead = false;  // the reader stopped without a final batch
    };
    enum class PairStatus { Pair, Lane, End, Error };

    PairStatus nextPair(FastxChunkFillContext& context, int* laneOut, int* endCharOut);
    void normalize(uint32_t mate);
    void drainLane(uint32_t mate);
    bool hasRecord(uint32_t mate) const;
    FastxMateEnd endOf(uint32_t mate) const;
    std::string laneCount(uint32_t mate) const;
    bool appendPair(char* const* chunkIn, unsigned long long* totals,
                    FastxChunkFillContext& context);

    std::vector<std::unique_ptr<FastxMateReader>> readers_;
    std::vector<Cursor> cursors_;
    bool started_ = false;
    bool stopped_ = false;
    bool endLogged_ = false;
    bool extraWarned_ = false;
    uint64_t pairs_ = 0;
    uint64_t truncatedLanes_ = 0;
};

}  // namespace input
}  // namespace star

#endif
