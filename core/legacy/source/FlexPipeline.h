#ifndef H_FlexPipeline
#define H_FlexPipeline

#include "IncludeDefine.h"
#include <atomic>
#include <condition_variable>
#include <cstring>
#include <memory>
#include <mutex>
#include <thread>
#include <string>
#include <vector>
#include <zlib.h>
#include "FlexGdna.h"
#include "input/BgzfStarAdapter.h"

static constexpr uint32_t kFlexPipeNameMax = 512;
static constexpr uint32_t kFlexPipeSeqMax = DEF_readSeqLengthMax + 1;
static constexpr uint32_t kFlexPipeCbMatchIndMax = 4;

struct ReadPacket {
    char name[kFlexPipeNameMax];
    char seq[2][kFlexPipeSeqMax];
    char qual[2][kFlexPipeSeqMax];
    uint32_t readLen[2];
    uint64_t iReadAll;
    uint64_t laneOrdinal;
    uint8_t  laneId;
    uint32_t readFilesIndex;
    char     readFilter;
    bool     eof;

    ReadPacket() : readLen{0,0}, iReadAll(0), laneOrdinal(UINT64_MAX), laneId(0), readFilesIndex(0), readFilter('Y'), eof(false) {
        name[0] = '\0';
        seq[0][0] = '\0'; seq[1][0] = '\0';
        qual[0][0] = '\0'; qual[1][0] = '\0';
    }
};

static constexpr uint32_t kFlexPipeBarcodeSeqMax = 64;

struct DecisionPacket {
    enum Verdict : uint8_t { KEEP = 0, DENY = 1 };

    uint64_t iReadAll;
    uint32_t readFilesIndex;
    Verdict  verdict;
    uint16_t geneIdx15;
    uint8_t  cacheClass;
    FlexGdnaRegion probeRegion;
    const char *denyReason;
    bool     eof;

    // Raw R1 (barcode read) for deferred CB/UMI extraction by Solo consumer
    char     barcodeSeq[kFlexPipeBarcodeSeqMax];
    char     barcodeQual[kFlexPipeBarcodeSeqMax];
    uint32_t barcodeLen;
    char     readName[kFlexPipeNameMax];

    // Raw 8-byte sample tag region from R2 for deferred sample detection
    char     sampleTagSeq[8];
    uint32_t sampleTagLen;

    DecisionPacket() : iReadAll(0), readFilesIndex(0),
                       verdict(DENY), geneIdx15(0), cacheClass(0), probeRegion(FlexGdnaUnknown),
                       denyReason(nullptr), eof(false),
                       barcodeLen(0), sampleTagLen(0) {
        barcodeSeq[0] = '\0'; barcodeQual[0] = '\0';
        readName[0] = '\0';
        std::memset(sampleTagSeq, 0, sizeof(sampleTagSeq));
    }
};

struct EnrichedPacket : ReadPacket {
    int      cbMatch;
    uint32_t cbMatchInd[kFlexPipeCbMatchIndMax];
    uint32_t cbMatchIndN;
    uint32_t umiB;
    uint8_t  detectedSampleToken;
    uint16_t hashScreenSampleIdx;
    uint16_t residualAnchorGeneIdx15;

    EnrichedPacket() : ReadPacket(), cbMatch(-1), cbMatchIndN(0),
                       umiB(0), detectedSampleToken(0), hashScreenSampleIdx(0),
                       residualAnchorGeneIdx15(0) {
        std::memset(cbMatchInd, 0, sizeof(cbMatchInd));
    }
};

template <typename T>
class BoundedQueue {
public:
    explicit BoundedQueue(size_t capacity = 256)
        : capacity_(capacity), buf_(capacity), head_(0), tail_(0), count_(0) {}

    void push(T&& item) {
        std::unique_lock<std::mutex> lock(mu_);
        cvFull_.wait(lock, [this]{ return count_ < capacity_ || closed_; });
        if (closed_) return;
        buf_[tail_] = std::move(item);
        tail_ = (tail_ + 1) % capacity_;
        ++count_;
        cvEmpty_.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mu_);
        cvEmpty_.wait(lock, [this]{ return count_ > 0 || closed_; });
        if (count_ == 0) return false;
        item = std::move(buf_[head_]);
        head_ = (head_ + 1) % capacity_;
        --count_;
        cvFull_.notify_one();
        return true;
    }

    // Non-blocking variants for producers that must not wait on a consumer
    // that may not exist yet (fully-fused Flex; see enqueueForAlign). On a
    // full queue try_push leaves `item` untouched and returns false; on a
    // closed queue it drops the item exactly as push() does and returns true.
    bool try_push(T& item) {
        std::lock_guard<std::mutex> lock(mu_);
        if (closed_) return true;
        if (count_ >= capacity_) return false;
        buf_[tail_] = std::move(item);
        tail_ = (tail_ + 1) % capacity_;
        ++count_;
        cvEmpty_.notify_one();
        return true;
    }

    bool try_pop(T& item) {
        std::lock_guard<std::mutex> lock(mu_);
        if (count_ == 0) return false;
        item = std::move(buf_[head_]);
        head_ = (head_ + 1) % capacity_;
        --count_;
        cvFull_.notify_one();
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> lock(mu_);
        closed_ = true;
        cvEmpty_.notify_all();
        cvFull_.notify_all();
    }

    size_t size() {
        std::lock_guard<std::mutex> lock(mu_);
        return count_;
    }

private:
    size_t capacity_;
    std::vector<T> buf_;
    size_t head_;
    size_t tail_;
    size_t count_;
    bool closed_ = false;
    std::mutex mu_;
    std::condition_variable cvFull_;
    std::condition_variable cvEmpty_;
};

struct FlexPipelineCounters {
    std::atomic<uint64_t> readsTotal{0};
    std::atomic<uint64_t> triageKeep{0};
    std::atomic<uint64_t> triageDeny{0};
    std::atomic<uint64_t> triageSampleReject{0};
    std::atomic<uint64_t> triageMiss{0};
    // Reads a fused producer aligned itself because alignQ was full.
    std::atomic<uint64_t> alignHelped{0};
    std::atomic<uint64_t> perLaneReads[64];

    FlexPipelineCounters() {
        for (size_t lane = 0; lane < 64; ++lane) {
            perLaneReads[lane].store(0, std::memory_order_relaxed);
        }
    }
};

struct LaneFiles {
    std::string r2path;
    std::string r1path;
};

struct FlexCbqRangeTask {
    int laneId = 0;
    uint64_t firstRecord = 0;
    uint64_t recordCount = 0;
    uint64_t globalFirst = 0;
};

struct FlexBgzfLane {
    bool range = false;
    std::shared_ptr<star::input::BgzfStarAdapter> adapter;
};

struct FlexBgzfRangeTask {
    int laneId = 0;
};

// One record of a FASTQ batch: byte offsets into the batch arena. The five
// strings are stored NUL-terminated so the consumer can hand plain char* to the
// existing barcode and screen code without copying again.
struct FlexFastqRecordRef {
    // offName/offSeq0/offQual0 index the mate-0 arena, offSeq1/offQual1 the
    // mate-1 arena: the two mates are read by different threads.
    uint32_t offName = 0, offSeq0 = 0, offQual0 = 0, offSeq1 = 0, offQual1 = 0;
    uint32_t nameLen = 0, len0 = 0, len1 = 0;
};

// One mate-1 record as read by the barcode-mate reader.
struct FlexFastqMateRecordRef {
    uint32_t offSeq = 0, offQual = 0, len = 0;
};

// A run of mate-1 records. The barcode mate is read on its own thread and
// handed over in chunks of the same size as the mate-0 batches, so chunk k
// covers the same record range as batch k: counting records is the whole
// contract between the two readers, no offsets are exchanged.
struct FlexFastqMateChunk {
    std::vector<FlexFastqMateRecordRef> recs;
    std::vector<char> data;
    bool eof = false;

    void reset() { recs.clear(); data.clear(); eof = false; }
    uint32_t append(const char *text, uint32_t length) {
        const uint32_t at = static_cast<uint32_t>(data.size());
        data.insert(data.end(), text, text + length);
        data.push_back('\0');
        return at;
    }
};

// A batch of FASTQ records read from one lane. The lane reader fills these and
// hands them to the fused workers, so decompression and the per-read work
// (sample tag, hash screen, CB/UMI, recording) no longer share a thread. With a
// single delivered file pair there is one lane, and without this the reader was
// the only busy thread while the rest of the pool idled.
struct FlexFastqBatch {
    int laneId = -1;
    uint64_t globalFirst = 0;
    uint64_t laneFirst = 0;
    std::vector<FlexFastqRecordRef> recs;
    std::vector<char> data;    // mate 0: names, sequences, qualities
    std::vector<char> data1;   // mate 1: sequences and qualities

    void reset(int lane) {
        laneId = lane;
        globalFirst = 0;
        laneFirst = 0;
        recs.clear();
        data.clear();
        data1.clear();
    }
    // Appends a NUL-terminated copy and returns its offset in the arena.
    uint32_t append(const char *text, uint32_t length) {
        const uint32_t at = static_cast<uint32_t>(data.size());
        data.insert(data.end(), text, text + length);
        data.push_back('\0');
        return at;
    }
    const char *at(uint32_t offset) const { return data.data() + offset; }
    char *at(uint32_t offset) { return data.data() + offset; }
    char *at1(uint32_t offset) { return data1.data() + offset; }
};

static constexpr uint32_t kFlexFastqBatchRecords = 2048;
static constexpr size_t kFlexFastqBatchPool = 64;
// Barcode-mate chunks in flight per lane reader.
static constexpr size_t kFlexFastqMateChunks = 8;

struct FlexPipelineState {
    BoundedQueue<ReadPacket> readerQ;
    std::vector<BoundedQueue<DecisionPacket>*> soloQ;
    BoundedQueue<EnrichedPacket> alignQ;
    FlexPipelineCounters counters;

    std::atomic<uint64_t> iReadAllGlobal{0};
    std::atomic<int> readersFinished{0};
    std::atomic<int> triageFinished{0};
    std::atomic<bool> pipelineDone{false};
    int nLanes = 0;
    int nSolo = 0;
    int nTriage = 1;

    // Lane work-stealing: atomic counter for dynamic lane claim
    std::atomic<int> nextLaneIdx{0};
    std::vector<LaneFiles> laneFiles;
    // Track how many fused threads are still active (reading or aligning)
    int nFusedThreads = 0;
    std::vector<FlexCbqRangeTask> cbqRangeTasks;
    std::atomic<int> nextCbqRangeIdx{0};
    std::vector<FlexBgzfLane> bgzfLanes;
    std::vector<FlexBgzfRangeTask> bgzfRangeTasks;
    std::atomic<int> nextBgzfRangeIdx{0};
    int bgzfReaderWorkers = 0;
    bool bgzfRangeActive = false;
    // Set once before worker launch. MAP covers bounded residual-alignment
    // batches; FEATURE covers bounded BGZF inflate work items.
    bool dynamicPermitsEnabled = false;

    // FASTQ lane readers hand batches here; every fused thread consumes them.
    BoundedQueue<FlexFastqBatch*> fastqReadyQ{kFlexFastqBatchPool};
    BoundedQueue<FlexFastqBatch*> fastqFreeQ{kFlexFastqBatchPool};
    std::atomic<int> fastqReadersDone{0};

    std::atomic<bool> inputFailed{false};
    std::mutex inputErrorMutex;
    std::string inputError;

    ~FlexPipelineState() {
        for (auto* q : soloQ) delete q;
        FlexFastqBatch *batch = nullptr;
        while (fastqReadyQ.try_pop(batch)) delete batch;
        while (fastqFreeQ.try_pop(batch)) delete batch;
    }

    void init(int lanes, int soloConsumers, int triageThreads = 1, size_t queueCapacity = 256) {
        nLanes = lanes;
        nSolo = soloConsumers;
        nTriage = triageThreads;
        for (int i = 0; i < nSolo; ++i)
            soloQ.push_back(new BoundedQueue<DecisionPacket>(queueCapacity));
        for (size_t i = 0; i < kFlexFastqBatchPool; ++i) {
            FlexFastqBatch *batch = new FlexFastqBatch();
            batch->recs.reserve(kFlexFastqBatchRecords);
            if (!fastqFreeQ.try_push(batch)) delete batch;
        }
    }

    // Atomically claim the next unprocessed lane. Returns -1 if all lanes claimed.
    int claimNextLane() {
        int lane = nextLaneIdx.fetch_add(1, std::memory_order_relaxed);
        return (lane < nLanes) ? lane : -1;
    }

    bool claimNextCbqRange(FlexCbqRangeTask *task) {
        int index = nextCbqRangeIdx.fetch_add(1, std::memory_order_relaxed);
        if (index < 0 || index >= static_cast<int>(cbqRangeTasks.size())) {
            return false;
        }
        if (task != nullptr) {
            *task = cbqRangeTasks[static_cast<size_t>(index)];
        }
        return true;
    }

    bool claimNextBgzfRange(FlexBgzfRangeTask *task) {
        int index = nextBgzfRangeIdx.fetch_add(1, std::memory_order_relaxed);
        if (index < 0 || index >= static_cast<int>(bgzfRangeTasks.size())) {
            return false;
        }
        if (task != nullptr) {
            *task = bgzfRangeTasks[static_cast<size_t>(index)];
        }
        return true;
    }

    bool laneUsesBgzfRange(int lane) const {
        return lane >= 0 && lane < static_cast<int>(bgzfLanes.size()) &&
               bgzfLanes[static_cast<size_t>(lane)].range;
    }

    void failInput(const std::string& message) {
        inputFailed.store(true, std::memory_order_relaxed);
        std::lock_guard<std::mutex> lock(inputErrorMutex);
        if (inputError.empty()) {
            inputError = message;
        }
    }
};

class ReadAlign;
class SoloReadFeature;
class SoloReadBarcode;
class Parameters;
class Stats;

struct FlexLaneReaderArgs {
    FlexPipelineState* state;
    Parameters *P;
    gzFile gzR2;
    gzFile gzR1;
    int laneId;
    SoloReadFeature *readFeat;  // non-null in fully-fused mode
    Stats *stats;               // non-null in fully-fused mode
    ReadAlign *RA;              // non-null for role-switch to alignment worker
    int threadId;               // logical thread index (0..nFusedThreads-1)
    SoloReadBarcode *readBar;   // thread-owned CB evidence retained after join
};

struct FlexSoloConsumerArgs {
    FlexPipelineState *state;
    Parameters *P;
    int consumerId;
    SoloReadFeature *readFeat;
    Stats *stats;
    SoloReadBarcode *readBar;   // consumer-owned CB evidence retained after join
};

struct FlexAlignWorkerArgs {
    FlexPipelineState *state;
    ReadAlign *RA;
};

struct FlexTriageArgs {
    FlexPipelineState *state;
    Parameters *P;
};

struct FlexStatsReporterArgs {
    FlexPipelineState *state;
    Parameters *P;
};

void *flexLaneReaderThread(void *arg);
void *flexLaneReaderRouterThread(void *arg);
void *flexLaneReaderFullThread(void *arg);
void *flexTriageThread(void *arg);
void *flexSoloConsumerThread(void *arg);
void *flexAlignWorkerThread(void *arg);
void *flexStatsReporterThread(void *arg);
bool flexPrepareCbqRangeTasks(FlexPipelineState *state, Parameters &P,
                              int nWorkers, std::string *reason);
bool flexPrepareBgzfRangeTasks(FlexPipelineState *state, Parameters &P,
                               int nWorkers, std::string *reason,
                               bool *fatalError);

#endif
