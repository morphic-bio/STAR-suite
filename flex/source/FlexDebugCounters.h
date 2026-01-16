#ifndef CODE_FlexDebugCounters
#define CODE_FlexDebugCounters

// Thread-local debug counters for Flex pipeline
// These replace global atomics to avoid cache line bouncing in the hot mapping path
// Enable with -DDEBUG_FLEX_COUNTERS at compile time

#include <cstdint>
#include <vector>
#include <mutex>

#ifdef DEBUG_FLEX_COUNTERS

struct FlexDebugCounters {
    // SoloReadFeature_record.cpp counters (probe/genomic resolution)
    uint64_t probeAlignCount = 0;
    uint64_t genomicAlignCount = 0;
    uint64_t probeResolvedCount = 0;
    uint64_t genomicResolvedCount = 0;
    uint64_t resolverDropped = 0;
    uint64_t probeMissingIdx = 0;
    uint64_t resolverDropProbeDisagree = 0;
    uint64_t resolverDropGenomicDisagree = 0;
    uint64_t resolverDropMixed = 0;
    uint64_t resolverDropNoCandidates = 0;
    uint64_t resolverKeepProbe = 0;
    uint64_t resolverKeepGenomic = 0;
    uint64_t genomicAlignWithProbeGenes = 0;
    uint64_t genomicOnlyReadsWithProbeGenes = 0;
    uint64_t genomicOnlyProbeGeneCount = 0;
    uint64_t genomicDroppedMapq = 0;
    uint64_t genomicDroppedNM = 0;
    uint64_t genomicDroppedMMrate = 0;
    uint64_t probeDroppedNM = 0;
    uint64_t probeDroppedMMrate = 0;
    
    // SoloReadBarcode_getCBandUMI.cpp counters (CB correction)
    uint64_t inlineCbExact = 0;
    uint64_t inlineCbCorrected = 0;
    uint64_t inlineCbNRescued = 0;
    uint64_t inlineCbRejected = 0;
    uint64_t inlineCbWithN = 0;
    
    // BAMoutput.cpp counters (sample detection)
    uint64_t debugSampleDetectionCount = 0;
    uint64_t recordsWithoutSample = 0;
    uint64_t ambiguousCbCount = 0;
    
    // Add this counter's values to another
    void addTo(FlexDebugCounters &other) const {
        other.probeAlignCount += probeAlignCount;
        other.genomicAlignCount += genomicAlignCount;
        other.probeResolvedCount += probeResolvedCount;
        other.genomicResolvedCount += genomicResolvedCount;
        other.resolverDropped += resolverDropped;
        other.probeMissingIdx += probeMissingIdx;
        other.resolverDropProbeDisagree += resolverDropProbeDisagree;
        other.resolverDropGenomicDisagree += resolverDropGenomicDisagree;
        other.resolverDropMixed += resolverDropMixed;
        other.resolverDropNoCandidates += resolverDropNoCandidates;
        other.resolverKeepProbe += resolverKeepProbe;
        other.resolverKeepGenomic += resolverKeepGenomic;
        other.genomicAlignWithProbeGenes += genomicAlignWithProbeGenes;
        other.genomicOnlyReadsWithProbeGenes += genomicOnlyReadsWithProbeGenes;
        other.genomicOnlyProbeGeneCount += genomicOnlyProbeGeneCount;
        other.genomicDroppedMapq += genomicDroppedMapq;
        other.genomicDroppedNM += genomicDroppedNM;
        other.genomicDroppedMMrate += genomicDroppedMMrate;
        other.probeDroppedNM += probeDroppedNM;
        other.probeDroppedMMrate += probeDroppedMMrate;
        other.inlineCbExact += inlineCbExact;
        other.inlineCbCorrected += inlineCbCorrected;
        other.inlineCbNRescued += inlineCbNRescued;
        other.inlineCbRejected += inlineCbRejected;
        other.inlineCbWithN += inlineCbWithN;
        other.debugSampleDetectionCount += debugSampleDetectionCount;
        other.recordsWithoutSample += recordsWithoutSample;
        other.ambiguousCbCount += ambiguousCbCount;
    }
};

// Thread-local counter instance - each thread gets its own copy
extern thread_local FlexDebugCounters g_flexCounters;

// Global registry of all thread counter instances for summation
extern std::vector<FlexDebugCounters*> g_allFlexCounters;
extern std::mutex g_flexCountersMutex;

// Register current thread's counters (call once per thread at start)
void flexCountersRegisterThread();

// Sum all thread counters and return total
FlexDebugCounters flexCountersSumAll();

// Convenience macros for incrementing counters
#define FLEX_COUNT_INC(field) (g_flexCounters.field++)
#define FLEX_COUNT_ADD(field, n) (g_flexCounters.field += (n))

#else // !DEBUG_FLEX_COUNTERS

// No-op macros when counters are disabled
#define FLEX_COUNT_INC(field) ((void)0)
#define FLEX_COUNT_ADD(field, n) ((void)0)

// Stub functions that do nothing
inline void flexCountersRegisterThread() {}

#endif // DEBUG_FLEX_COUNTERS

#endif // CODE_FlexDebugCounters
