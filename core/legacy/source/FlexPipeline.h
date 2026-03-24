#ifndef H_FlexPipeline
#define H_FlexPipeline

#include "IncludeDefine.h"
#include <atomic>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <vector>
#include <zlib.h>

static constexpr uint32_t kFlexPipeNameMax = 512;
static constexpr uint32_t kFlexPipeSeqMax = DEF_readSeqLengthMax + 1;
static constexpr uint32_t kFlexPipeCbMatchIndMax = 4;

struct ReadPacket {
    char name[kFlexPipeNameMax];
    char seq[2][kFlexPipeSeqMax];
    char qual[2][kFlexPipeSeqMax];
    uint32_t readLen[2];
    uint64_t iReadAll;
    uint8_t  laneId;
    uint32_t readFilesIndex;
    char     readFilter;
    bool     eof;

    ReadPacket() : readLen{0,0}, iReadAll(0), laneId(0), readFilesIndex(0), readFilter('Y'), eof(false) {
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
                       verdict(DENY), geneIdx15(0), cacheClass(0), denyReason(nullptr), eof(false),
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

    EnrichedPacket() : ReadPacket(), cbMatch(-1), cbMatchIndN(0),
                       umiB(0), detectedSampleToken(0), hashScreenSampleIdx(0) {
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
    std::atomic<uint64_t> triageMiss{0};
    uint64_t perLaneReads[64];

    FlexPipelineCounters() { std::memset(perLaneReads, 0, sizeof(perLaneReads)); }
};

struct LaneFiles {
    std::string r2path;
    std::string r1path;
};

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

    ~FlexPipelineState() {
        for (auto* q : soloQ) delete q;
    }

    void init(int lanes, int soloConsumers, int triageThreads = 1, size_t queueCapacity = 256) {
        nLanes = lanes;
        nSolo = soloConsumers;
        nTriage = triageThreads;
        for (int i = 0; i < nSolo; ++i)
            soloQ.push_back(new BoundedQueue<DecisionPacket>(queueCapacity));
    }

    // Atomically claim the next unprocessed lane. Returns -1 if all lanes claimed.
    int claimNextLane() {
        int lane = nextLaneIdx.fetch_add(1, std::memory_order_relaxed);
        return (lane < nLanes) ? lane : -1;
    }
};

class ReadAlign;
class SoloReadFeature;
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
};

struct FlexSoloConsumerArgs {
    FlexPipelineState *state;
    Parameters *P;
    int consumerId;
    SoloReadFeature *readFeat;
    Stats *stats;
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

#endif
