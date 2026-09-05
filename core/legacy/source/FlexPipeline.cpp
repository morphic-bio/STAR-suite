#include "FlexPipeline.h"
#include "FlexHashScreen.h"
#include "SoloReadBarcode.h"
#include "SoloReadFeature.h"
#include "SoloRead.h"
#include "SoloFeatureTypes.h"
#include "SampleDetector.h"
#include "SequenceFuns.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include "input/CbqInputModule.h"
#include "input/BgzfStarAdapter.h"

#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <chrono>

namespace {

static constexpr int kGzBufSize = 1 << 20;
static constexpr size_t kFlexAlignPermitBatch = 64;

ThreadControl::PermitHookContext kFlexBgzfPermitContext{
    ThreadControl::PermitDomain::FEATURE};

uint64_t flexBgzfPermitAcquire(void *context) {
    const ThreadControl::PermitHookContext *permit =
        static_cast<const ThreadControl::PermitHookContext *>(context);
    const ThreadControl::PermitDomain domain = permit == nullptr
        ? ThreadControl::PermitDomain::FEATURE : permit->domain;
    return g_threadChunks.mapPermitAcquireForDomain(domain);
}

void flexBgzfPermitRelease(void *context,
                           uint64_t waitNs,
                           uint64_t workUnits,
                           uint64_t workBytes,
                           uint64_t workNs) {
    const ThreadControl::PermitHookContext *permit =
        static_cast<const ThreadControl::PermitHookContext *>(context);
    const ThreadControl::PermitDomain domain = permit == nullptr
        ? ThreadControl::PermitDomain::FEATURE : permit->domain;
    g_threadChunks.mapPermitReleaseForDomain(
        domain, waitNs, workUnits, workBytes, workNs);
}

bool gzReadLine(gzFile gz, char *buf, int maxLen, uint32_t *lengthOut = nullptr) {
    char *ret = gzgets(gz, buf, maxLen);
    if (ret == Z_NULL) return false;
    int len = static_cast<int>(strlen(buf));
    if (len > 0 && buf[len - 1] == '\n') {
        buf[--len] = '\0';
    }
    if (lengthOut != nullptr) {
        *lengthOut = static_cast<uint32_t>(len);
    }
    return true;
}

bool gzConsumeLine(gzFile gz, char *buf, int maxLen) {
    return gzgets(gz, buf, maxLen) != Z_NULL;
}

// FASTQ qualities normally have exactly the sequence length.  Check the known
// newline position first and retain the general scan as a compatibility
// fallback for malformed or non-canonical records.
bool gzReadQualityLine(gzFile gz, char *buf, int maxLen, uint32_t expectedLength) {
    if (expectedLength < static_cast<uint32_t>(maxLen)) {
        buf[expectedLength] = '\0';
    }
    if (gzgets(gz, buf, maxLen) == Z_NULL) {
        return false;
    }
    if (expectedLength < static_cast<uint32_t>(maxLen) &&
        buf[expectedLength] == '\n') {
        buf[expectedLength] = '\0';
        return true;
    }

    const size_t length = std::strlen(buf);
    if (length > 0 && buf[length - 1] == '\n') {
        buf[length - 1] = '\0';
    }
    return true;
}

size_t fastqReadNameLength(const char *name, size_t length) {
    size_t end = length;
    const void *space = std::memchr(name, ' ', end);
    if (space != nullptr) {
        end = static_cast<const char *>(space) - name;
    }
    const void *tab = std::memchr(name, '\t', end);
    if (tab != nullptr) {
        end = static_cast<const char *>(tab) - name;
    }
    return end;
}

} // namespace

void *flexLaneReaderThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    const int laneId = ctx->laneId;
    gzFile gzR2 = ctx->gzR2;
    gzFile gzR1 = ctx->gzR1;

    char lineBuf[kFlexPipeSeqMax + 256];

    while (true) {
        ReadPacket pkt;
        pkt.laneId = static_cast<uint8_t>(laneId);
        pkt.readFilesIndex = static_cast<uint32_t>(laneId);
        pkt.eof = false;

        // R2: @name, seq, +, qual
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        {
            // Parse name line: @name [iReadAll] [readFilter] [readFilesIndex]
            const char *src = lineBuf;
            if (*src == '@') ++src;
            size_t nameLen = 0;
            while (src[nameLen] && src[nameLen] != ' ' && src[nameLen] != '\t')
                ++nameLen;
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(pkt.name, src, nameLen);
            pkt.name[nameLen] = '\0';
            // Parse readFilter from Illumina header if present
            const char *rest = src + nameLen;
            while (*rest == ' ' || *rest == '\t') ++rest;
            // iReadAll field (skip)
            while (*rest && *rest != ' ' && *rest != '\t') ++rest;
            while (*rest == ' ' || *rest == '\t') ++rest;
            // readFilter field
            if (*rest) pkt.readFilter = *rest;
        }

        if (!gzReadLine(gzR2, pkt.seq[0], kFlexPipeSeqMax, &pkt.readLen[0])) break;

        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break; // + line
        if (!gzReadLine(gzR2, pkt.qual[0], kFlexPipeSeqMax)) break;

        // R1: @name, seq, +, qual
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, pkt.seq[1], kFlexPipeSeqMax, &pkt.readLen[1])) break;
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, pkt.qual[1], kFlexPipeSeqMax)) break;

        pkt.iReadAll = st->iReadAllGlobal.fetch_add(1);

        st->counters.perLaneReads[laneId]++;
        st->readerQ.push(std::move(pkt));
    }

    gzclose(gzR2);
    gzclose(gzR1);

    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nLanes) {
        st->readerQ.close();
    }

    return nullptr;
}

void *flexLaneReaderRouterThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    const int laneId = ctx->laneId;
    gzFile gzR2 = ctx->gzR2;
    gzFile gzR1 = ctx->gzR1;

    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;

    char lineBuf[kFlexPipeSeqMax + 256];
    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];

    while (true) {
        // R2: @name, seq, +, qual
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        {
            const char *src = lineBuf;
            if (*src == '@') ++src;
            size_t nameLen = 0;
            while (src[nameLen] && src[nameLen] != ' ' && src[nameLen] != '\t')
                ++nameLen;
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(name, src, nameLen);
            name[nameLen] = '\0';
        }
        uint32_t readLen0 = 0;
        if (!gzReadLine(gzR2, seq0, kFlexPipeSeqMax, &readLen0)) break;
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR2, qual0, kFlexPipeSeqMax)) break;

        // R1: @name, seq, +, qual
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        uint32_t readLen1 = 0;
        if (!gzReadLine(gzR1, seq1, kFlexPipeSeqMax, &readLen1)) break;
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, qual1, kFlexPipeSeqMax)) break;

        uint64_t iReadAll = st->iReadAllGlobal.fetch_add(1);
        st->counters.perLaneReads[laneId]++;

        FlexHashScreenDecision decision = cache.classifyReadH0Offset0(seq0, readLen0);

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {
            DecisionPacket dp;
            dp.iReadAll = iReadAll;
            dp.readFilesIndex = static_cast<uint32_t>(laneId);
            dp.eof = false;

            const uint32_t bcLen = std::min(readLen1, kFlexPipeBarcodeSeqMax - 1);
            std::memcpy(dp.barcodeSeq, seq1, bcLen);
            dp.barcodeSeq[bcLen] = '\0';
            std::memcpy(dp.barcodeQual, qual1, bcLen);
            dp.barcodeQual[bcLen] = '\0';
            dp.barcodeLen = bcLen;
            std::memcpy(dp.readName, name, kFlexPipeNameMax);

            if (sampleTagOffset + 8 <= readLen0) {
                std::memcpy(dp.sampleTagSeq, seq0 + sampleTagOffset, 8);
                dp.sampleTagLen = 8;
            } else {
                dp.sampleTagLen = 0;
            }

            if (decision.action == FlexHashScreenDecision::Keep) {
                dp.verdict = DecisionPacket::KEEP;
                dp.geneIdx15 = decision.geneIdx15;
                dp.cacheClass = decision.cacheClass;
                dp.probeRegion = decision.probeRegion;
                dp.denyReason = nullptr;
                st->counters.triageKeep.fetch_add(1);
            } else {
                dp.verdict = DecisionPacket::DENY;
                dp.geneIdx15 = 0;
                dp.cacheClass = 0;
                dp.denyReason = "NEG_PROBE_AMBIG";
                st->counters.triageDeny.fetch_add(1);
            }

            int shard = static_cast<int>(iReadAll % static_cast<uint64_t>(st->nSolo));
            st->soloQ[shard]->push(std::move(dp));
        } else {
            EnrichedPacket ep;
            std::memcpy(ep.name, name, kFlexPipeNameMax);
            std::memcpy(ep.seq[0], seq0, readLen0 + 1);
            std::memcpy(ep.seq[1], seq1, readLen1 + 1);
            std::memcpy(ep.qual[0], qual0, readLen0 + 1);
            std::memcpy(ep.qual[1], qual1, readLen1 + 1);
            ep.readLen[0] = readLen0;
            ep.readLen[1] = readLen1;
            ep.iReadAll = iReadAll;
            ep.laneId = static_cast<uint8_t>(laneId);
            ep.readFilesIndex = static_cast<uint32_t>(laneId);
            ep.readFilter = 'Y';
            ep.eof = false;
            ep.cbMatch = -1;
            ep.cbMatchIndN = 0;
            ep.umiB = 0;
            ep.detectedSampleToken = 0xFF;
            ep.hashScreenSampleIdx = 0;

            st->counters.triageMiss.fetch_add(1);
            st->alignQ.push(std::move(ep));
        }

        st->counters.readsTotal.fetch_add(1);
    }

    gzclose(gzR2);
    gzclose(gzR1);

    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nLanes) {
        for (int i = 0; i < st->nSolo; ++i)
            st->soloQ[i]->close();
        st->alignQ.close();
    }

    return nullptr;
}

// Use the same detector for the early rejection and downstream sample assignment.
// When sample tagging is not configured, every read remains eligible. When it is
// configured, an unmatched/short tag cannot contribute to a per-sample matrix and
// therefore never needs residual alignment.
static bool detectConfiguredSampleTag(const char *readSeq, uint32_t readLen,
                                      const ParametersSolo &pSolo,
                                      SampleDetector *sampleDet,
                                      bool sampleDetReady,
                                      uint8_t &detectedSampleToken) {
    detectedSampleToken = 0xFF;
    if (!sampleDetReady || sampleDet == nullptr) {
        return true;
    }
    auto detectAt = [&](int64_t offset) -> bool {
        if (offset < 0 || static_cast<uint64_t>(offset) + 8 > readLen) {
            return false;
        }
        uint8_t packedTag[4];
        nuclPackBAM(const_cast<char *>(readSeq + offset),
                    reinterpret_cast<char *>(packedTag), 8);
        const uint32_t sampleIdx = sampleDet->detectSampleFromPackedTag(packedTag);
        if (sampleIdx == 0) {
            return false;
        }
        detectedSampleToken = static_cast<uint8_t>(sampleIdx & 0x1Fu);
        return true;
    };

    const int64_t primary = static_cast<int64_t>(pSolo.sampleProbeOffset);
    if (detectAt(primary)) {
        return true;
    }
    if (!pSolo.sampleStrictMatch && pSolo.sampleSearchNearby) {
        static const int deltas[] = {-1, 1, -2, 2};
        for (int delta : deltas) {
            if (detectAt(primary + delta)) {
                return true;
            }
        }
    }
    return false;
}

// Process one lane: read all FASTQ records, hash screen, inline Solo for hits, push misses to alignQ.
// All per-thread state (localBar, sampleDet, scratch buffers) is owned by the caller
// and reused across lane claims — no per-lane allocation.
// Diagnostic only: STAR_FLEX_HASH_H0_ONLY=1 makes the fused path consult the
// exact-match (H0) tier alone, skipping the H1/deny lookup that normally runs
// on an H0 miss. Output changes (H1 denies become misses), so this is for
// timing where the hash cost lies, never for a production run.
static bool flexHashH0OnlyDiagnostic()
{
    static const bool enabled = [] {
        const char *value = std::getenv("STAR_FLEX_HASH_H0_ONLY");
        return value != nullptr && value[0] == '1';
    }();
    return enabled;
}

// Per-thread tallies for the shared pipeline counters. Every fused thread used
// to bump three or four of these atomics per read, and they all sit on one
// cache line, so 32 threads spent their time bouncing that line between cores.
// The only consumers are a ten-second progress line and the final totals, so
// each thread counts locally and publishes every kFlushEvery reads and on exit
// (the destructor flushes, so early returns still publish exact totals).
struct FlexLocalCounters {
    static const uint64_t kFlushEvery = 65536;
    FlexPipelineCounters &shared;
    int laneId;
    uint64_t total = 0, keep = 0, deny = 0, sampleReject = 0, miss = 0, lane = 0, sinceFlush = 0;

    FlexLocalCounters(FlexPipelineCounters &counters, int lane_)
        : shared(counters), laneId(lane_) {}
    ~FlexLocalCounters() { flush(); }

    void countRead() {
        ++total;
        if (++sinceFlush >= kFlushEvery) flush();
    }
    void flush() {
        if (total) shared.readsTotal.fetch_add(total, std::memory_order_relaxed);
        if (keep)  shared.triageKeep.fetch_add(keep, std::memory_order_relaxed);
        if (deny)  shared.triageDeny.fetch_add(deny, std::memory_order_relaxed);
        if (sampleReject)
            shared.triageSampleReject.fetch_add(sampleReject, std::memory_order_relaxed);
        if (miss)  shared.triageMiss.fetch_add(miss, std::memory_order_relaxed);
        if (lane && laneId >= 0 && laneId < 64)
            shared.perLaneReads[laneId].fetch_add(lane, std::memory_order_relaxed);
        total = keep = deny = sampleReject = miss = lane = sinceFlush = 0;
    }
};

// Hand a hash miss to the aligner. In fully-fused mode every thread is a
// producer until it runs out of input, so nobody drains alignQ while reading
// is in progress: a blocking push on a full queue would wait forever for a
// consumer that does not exist yet (fused mode with --flexNoAlign 0 hung as
// soon as misses exceeded the queue capacity). A producer that finds the
// queue full therefore aligns queued packets itself until there is room.
// Consumer capacity grows exactly with queue pressure, no thread is reserved,
// and the path works at --runThreadN 1. Callers without an aligner (the
// staged modes, which spawn dedicated alignment workers) keep the blocking
// push.
static void enqueueForAlign(FlexPipelineState *st, EnrichedPacket &&ep, ReadAlign *RA)
{
    if (RA == nullptr) {
        st->alignQ.push(std::move(ep));
        return;
    }
    while (!st->alignQ.try_push(ep)) {
        EnrichedPacket queued;
        if (st->alignQ.try_pop(queued)) {
            uint64_t aligned = 0;
            uint64_t workBytes = 0;
            uint64_t waitNs = 0;
            if (st->dynamicPermitsEnabled) {
                waitNs = g_threadChunks.mapPermitAcquireForDomain(
                    ThreadControl::PermitDomain::MAP);
            }
            const std::chrono::steady_clock::time_point workStart =
                std::chrono::steady_clock::now();
            do {
                workBytes += queued.readLen[0] + queued.readLen[1];
                RA->oneReadFromPacket(queued);
                ++aligned;
            } while (st->dynamicPermitsEnabled &&
                     aligned < kFlexAlignPermitBatch &&
                     st->alignQ.try_pop(queued));
            if (st->dynamicPermitsEnabled) {
                const uint64_t workNs = static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - workStart).count());
                g_threadChunks.mapPermitReleaseForDomain(
                    ThreadControl::PermitDomain::MAP,
                    waitNs, aligned, workBytes, workNs);
            }
            st->counters.alignHelped.fetch_add(aligned, std::memory_order_relaxed);
        } else {
            std::this_thread::yield();
        }
    }
}

// Consume one already-claimed packet plus a short immediately available tail
// under one MAP lease. Batching keeps scheduler overhead small while retaining
// millisecond-scale handoff between BGZF inflation and residual alignment.
static void alignAvailableBatch(FlexPipelineState *st,
                                EnrichedPacket &packet,
                                ReadAlign *RA) {
    if (!st->dynamicPermitsEnabled) {
        RA->oneReadFromPacket(packet);
        return;
    }

    const uint64_t waitNs = g_threadChunks.mapPermitAcquireForDomain(
        ThreadControl::PermitDomain::MAP);
    const std::chrono::steady_clock::time_point workStart =
        std::chrono::steady_clock::now();
    uint64_t aligned = 0;
    uint64_t workBytes = 0;
    do {
        workBytes += packet.readLen[0] + packet.readLen[1];
        RA->oneReadFromPacket(packet);
        ++aligned;
    } while (aligned < kFlexAlignPermitBatch && st->alignQ.try_pop(packet));
    const uint64_t workNs = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - workStart).count());
    g_threadChunks.mapPermitReleaseForDomain(
        ThreadControl::PermitDomain::MAP,
        waitNs, aligned, workBytes, workNs);
}

static uint64_t processOneLane(
    FlexPipelineState *st, Parameters &P, int laneId,
    gzFile gzR2, gzFile gzR1,
    SoloReadFeature *readFeat, Stats *stats,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign = false,
    ReadAlign *RA = nullptr)
{
    auto &cache = FlexHashScreenCache::instance();
    FlexLocalCounters tally(st->counters, laneId);
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    char lineBuf[kFlexPipeSeqMax + 256];
    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];
    const std::string readNameExtra;
    uint64_t nReads = 0;

    while (true) {
        uint32_t headerLength = 0;
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf), &headerLength)) break;
        {
            const char *src = lineBuf;
            size_t available = headerLength;
            if (available > 0 && *src == '@') {
                ++src;
                --available;
            }
            size_t nameLen = fastqReadNameLength(src, available);
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(name, src, nameLen);
            name[nameLen] = '\0';
        }
        uint32_t readLen0 = 0;
        if (!gzReadLine(gzR2, seq0, kFlexPipeSeqMax, &readLen0)) break;
        if (!gzConsumeLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        if (noAlign) {
            if (!gzConsumeLine(gzR2, lineBuf, kFlexPipeSeqMax)) break;
        } else if (!gzReadQualityLine(gzR2, qual0, kFlexPipeSeqMax, readLen0)) {
            break;
        }

        if (!gzConsumeLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        uint32_t readLen1 = 0;
        if (!gzReadLine(gzR1, seq1, kFlexPipeSeqMax, &readLen1)) break;
        if (!gzConsumeLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadQualityLine(gzR1, qual1, kFlexPipeSeqMax, readLen1)) break;

        uint64_t iReadAll = st->iReadAllGlobal.fetch_add(1);
        ++tally.lane;
        nReads++;

        uint8_t detectedSampleToken = 0xFF;
        const bool sampleOK = detectConfiguredSampleTag(
            seq0, readLen0, P.pSolo, sampleDet, sampleDetReady,
            detectedSampleToken);

        FlexHashScreenDecision decision;
        if (sampleOK) {
            decision = flexHashH0OnlyDiagnostic()
                    ? cache.classifyReadH0Offset0(seq0, readLen0)
                    : cache.classifyReadH0H1Offset0(seq0, readLen0);
        } else {
            // The tag is outside the configured sample universe, including
            // all accepted variants. Alignment cannot make this read eligible
            // for a per-sample output, so reject it before residual alignment.
            decision.action = FlexHashScreenDecision::Deny;
        }

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {

            char *readSeqPtrs[2]  = { dummySeq, seq1 };
            char *readQualPtrs[2] = { dummyQual, qual1 };
            uint64 readLens[2]    = { 0, readLen1 };
            localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                                  static_cast<uint32_t>(laneId), name);

            localBar.detectedSampleToken = detectedSampleToken;

            if (decision.action == FlexHashScreenDecision::Keep) {
                record_flex_hash_screen_keep(readFeat, localBar, iReadAll,
                                             decision.geneIdx15, decision.cacheClass,
                                             decision.probeRegion);
                stats->hashScreenKeep++;
                if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
                ++tally.keep;
            } else {
                record_flex_hash_screen_deny(readFeat, localBar, iReadAll,
                                             sampleOK ? "NEG_PROBE_AMBIG" : "UNMATCHED_TAG");
                if (sampleOK) {
                    stats->hashScreenDeny++;
                    ++tally.deny;
                } else {
                    stats->hashScreenSampleReject++;
                    ++tally.sampleReject;
                }
            }
        } else {
            ++tally.miss;
            if (!noAlign) {
                EnrichedPacket ep;
                std::memcpy(ep.name, name, kFlexPipeNameMax);
                std::memcpy(ep.seq[0], seq0, readLen0 + 1);
                std::memcpy(ep.seq[1], seq1, readLen1 + 1);
                std::memcpy(ep.qual[0], qual0, readLen0 + 1);
                std::memcpy(ep.qual[1], qual1, readLen1 + 1);
                ep.readLen[0] = readLen0;
                ep.readLen[1] = readLen1;
                ep.iReadAll = iReadAll;
                ep.laneId = static_cast<uint8_t>(laneId);
                ep.readFilesIndex = static_cast<uint32_t>(laneId);
                ep.readFilter = 'Y';
                ep.eof = false;
                ep.cbMatch = -1;
                ep.cbMatchIndN = 0;
                ep.umiB = 0;
                ep.detectedSampleToken = detectedSampleToken;
                ep.hashScreenSampleIdx = 0;
                enqueueForAlign(st, std::move(ep), RA);
            }
        }

        tally.countRead();
    }

    gzclose(gzR2);
    gzclose(gzR1);

    return nReads;
}

static uint64_t processOneBgzfRange(
    FlexPipelineState *st, Parameters &P,
    const FlexBgzfRangeTask &task,
    SoloReadFeature *readFeat, Stats *stats,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign, ReadAlign *RA)
{
    const FlexBgzfLane &lanePlan = st->bgzfLanes[static_cast<size_t>(task.laneId)];
    std::string inputError;
    FlexLocalCounters tally(st->counters, task.laneId);
    if (lanePlan.adapter == nullptr) {
        st->failInput("Flex BGZF lane adapter is missing for lane " +
                      std::to_string(task.laneId));
        return 0;
    }

    auto &cache = FlexHashScreenCache::instance();
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};
    const std::string readNameExtra;

    uint64_t nReads = 0;
    const size_t bgzfRecordBatchSize = 256;
    std::vector<star::input::BgzfStarRecord> recordBatch(bgzfRecordBatchSize);
    star::input::BgzfStarBatchLease recordLease;
    while (!st->inputFailed.load(std::memory_order_relaxed)) {
        size_t recordsReturned = 0;
        const star::input::InputStatus status =
            lanePlan.adapter->next_records(recordBatch.data(), recordBatch.size(),
                                           &recordsReturned, &inputError,
                                           &recordLease);
        if (status == star::input::InputStatus::End) {
            break;
        }
        if (status == star::input::InputStatus::Error) {
            std::ostringstream message;
            message << "Flex BGZF range failed in lane " << task.laneId
                    << ": " << inputError;
            st->failInput(message.str());
            break;
        }
        const uint64_t batchGlobalFirst = st->iReadAllGlobal.fetch_add(
            recordsReturned, std::memory_order_relaxed);
        for (size_t batchIndex = 0; batchIndex < recordsReturned; ++batchIndex) {
        star::input::BgzfStarRecord &record = recordBatch[batchIndex];
        const size_t readLen0Size = record.mates[0].sequenceLength;
        const size_t readLen1Size = record.mates[1].sequenceLength;
        if (readLen0Size >= kFlexPipeSeqMax || readLen1Size >= kFlexPipeSeqMax ||
            record.mates[0].qualityLength >= kFlexPipeSeqMax ||
            record.mates[1].qualityLength >= kFlexPipeSeqMax) {
            std::ostringstream message;
            message << "Flex BGZF record " << record.read_ordinal
                    << " in lane " << task.laneId
                    << " exceeds STAR's maximum read length";
            st->failInput(message.str());
            break;
        }
        const uint32_t readLen0 = static_cast<uint32_t>(readLen0Size);
        const uint32_t readLen1 = static_cast<uint32_t>(readLen1Size);
        const char *seq0 = record.mates[0].sequence;
        const char *seq1 = record.mates[1].sequence;
        const char *qual0 = record.mates[0].quality_data();
        const char *qual1 = record.mates[1].quality_data();
        const size_t nameLength = std::min(static_cast<size_t>(record.read_name_length),
                                           static_cast<size_t>(kFlexPipeNameMax - 1));
        const char *name = record.mates[0].name_data();

        const uint64_t iReadAll = batchGlobalFirst + batchIndex;
        ++tally.lane;
        ++nReads;

        uint8_t detectedSampleToken = 0xFF;
        const bool sampleOK = detectConfiguredSampleTag(
            seq0, readLen0, P.pSolo, sampleDet, sampleDetReady,
            detectedSampleToken);

        FlexHashScreenDecision decision;
        if (sampleOK) {
            decision = flexHashH0OnlyDiagnostic()
                    ? cache.classifyReadH0Offset0(seq0, readLen0)
                    : cache.classifyReadH0H1Offset0(seq0, readLen0);
        } else {
            decision.action = FlexHashScreenDecision::Deny;
        }

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {
            char *readSeqPtrs[2] = {
                dummySeq, const_cast<char *>(seq1)
            };
            char *readQualPtrs[2] = {
                dummyQual, const_cast<char *>(qual1)
            };
            uint64 readLens[2] = {0, readLen1};
            localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                                 static_cast<uint32_t>(task.laneId), name, nameLength);

            localBar.detectedSampleToken = detectedSampleToken;

            if (decision.action == FlexHashScreenDecision::Keep) {
                record_flex_hash_screen_keep(readFeat, localBar, iReadAll,
                                             decision.geneIdx15, decision.cacheClass,
                                             decision.probeRegion);
                stats->hashScreenKeep++;
                if (localBar.cbMatch < 0) {
                    stats->hashScreenKeepNoBarcode++;
                }
                ++tally.keep;
            } else {
                record_flex_hash_screen_deny(readFeat, localBar, iReadAll,
                                             sampleOK ? "NEG_PROBE_AMBIG" : "UNMATCHED_TAG");
                if (sampleOK) {
                    stats->hashScreenDeny++;
                    ++tally.deny;
                } else {
                    stats->hashScreenSampleReject++;
                    ++tally.sampleReject;
                }
            }
        } else {
            ++tally.miss;
            if (!noAlign) {
                EnrichedPacket packet;
                std::memcpy(packet.name, name, nameLength);
                packet.name[nameLength] = '\0';
                std::memcpy(packet.seq[0], seq0, readLen0);
                std::memcpy(packet.seq[1], seq1, readLen1);
                std::memcpy(packet.qual[0], qual0, readLen0);
                std::memcpy(packet.qual[1], qual1, readLen1);
                packet.seq[0][readLen0] = '\0';
                packet.seq[1][readLen1] = '\0';
                packet.qual[0][readLen0] = '\0';
                packet.qual[1][readLen1] = '\0';
                packet.readLen[0] = readLen0;
                packet.readLen[1] = readLen1;
                packet.iReadAll = iReadAll;
                packet.laneId = static_cast<uint8_t>(task.laneId);
                packet.readFilesIndex = static_cast<uint32_t>(task.laneId);
                packet.readFilter = record.read_filter;
                packet.eof = false;
                packet.cbMatch = -1;
                packet.cbMatchIndN = 0;
                packet.umiB = 0;
                packet.detectedSampleToken = detectedSampleToken;
                packet.hashScreenSampleIdx = 0;
                enqueueForAlign(st, std::move(packet), RA);
            }
        }
        tally.countRead();
        }
    }
    return nReads;
}

static uint64_t processCbqModuleRecords(
    FlexPipelineState *st, Parameters &P, int laneId,
    const std::string &cbqPath,
    star::input::CbqInputModule &module,
    SoloReadFeature *readFeat, Stats *stats,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign,
    bool deterministicReadIds,
    uint64_t globalFirst,
    ReadAlign *RA)
{
    auto &cache = FlexHashScreenCache::instance();
    FlexLocalCounters tally(st->counters, laneId);

    std::string inputError;
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};
    const std::string readNameExtra;

    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];

    uint64_t nReads = 0;
    star::input::CbqReadBatchView batch;
    while (true) {
        inputError.clear();
        const star::input::InputStatus status = module.next_batch(&batch, &inputError);
        if (status == star::input::InputStatus::End) break;
        if (status == star::input::InputStatus::Error) {
            P.inOut->logMain << "ERROR: Flex CBQ pipeline failed while reading " << cbqPath
                             << ": " << inputError << "\n" << std::flush;
            break;
        }

        for (uint32_t i = 0; i < batch.record_count; ++i) {
            const star::input::CbqReadView &record = batch.records[i];
            if (record.segment_count < 2) {
                P.inOut->logMain << "ERROR: Flex CBQ pipeline expected paired reads in "
                                 << cbqPath << "\n" << std::flush;
                module.close();
                return nReads;
            }

            size_t readLen0Size = 0;
            size_t readLen1Size = 0;
            if (!star::input::materialize_cbq_segment_sequence_to_buffer(
                    record.segments[0], seq0, sizeof(seq0), &readLen0Size, &inputError) ||
                !star::input::materialize_cbq_segment_sequence_to_buffer(
                    record.segments[1], seq1, sizeof(seq1), &readLen1Size, &inputError)) {
                P.inOut->logMain << "ERROR: Flex CBQ sequence materialization failed in "
                                 << cbqPath << ": " << inputError << "\n" << std::flush;
                module.close();
                return nReads;
            }

            if (record.segments[0].quality.size != readLen0Size ||
                record.segments[1].quality.size != readLen1Size ||
                (readLen0Size != 0 && record.segments[0].quality.data == nullptr) ||
                (readLen1Size != 0 && record.segments[1].quality.data == nullptr)) {
                P.inOut->logMain << "ERROR: Flex CBQ quality length does not match sequence in "
                                 << cbqPath << "\n" << std::flush;
                module.close();
                return nReads;
            }
            const uint32_t readLen0 = static_cast<uint32_t>(readLen0Size);
            const uint32_t readLen1 = static_cast<uint32_t>(readLen1Size);
            const char *qual0 = record.segments[0].quality.data != nullptr
                ? record.segments[0].quality.data : dummyQual;
            const char *qual1 = record.segments[1].quality.data != nullptr
                ? record.segments[1].quality.data : dummyQual;
            const char *name = record.read_name.data != nullptr ? record.read_name.data : "";
            const size_t nameLength = std::min(record.read_name.size,
                                               static_cast<size_t>(kFlexPipeNameMax - 1));

            const uint64_t localOrdinal = nReads;
            uint64_t iReadAll = deterministicReadIds
                ? globalFirst + localOrdinal
                : st->iReadAllGlobal.fetch_add(1);
            ++tally.lane;
            nReads++;

            uint8_t detectedSampleToken = 0xFF;
            const bool sampleOK = detectConfiguredSampleTag(
                seq0, readLen0, P.pSolo, sampleDet, sampleDetReady,
                detectedSampleToken);

            FlexHashScreenDecision decision;
            if (sampleOK) {
                decision = flexHashH0OnlyDiagnostic()
                    ? cache.classifyReadH0Offset0(seq0, readLen0)
                    : cache.classifyReadH0H1Offset0(seq0, readLen0);
            } else {
                decision.action = FlexHashScreenDecision::Deny;
            }

            if (decision.action == FlexHashScreenDecision::Keep ||
                decision.action == FlexHashScreenDecision::Deny) {

                char *readSeqPtrs[2]  = { dummySeq, seq1 };
                char *readQualPtrs[2] = { dummyQual, const_cast<char *>(qual1) };
                uint64 readLens[2]    = { 0, readLen1 };

                localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                                      static_cast<uint32_t>(laneId), name, nameLength);

                localBar.detectedSampleToken = detectedSampleToken;

                if (decision.action == FlexHashScreenDecision::Keep) {
                    record_flex_hash_screen_keep(readFeat, localBar, iReadAll,
                                                 decision.geneIdx15, decision.cacheClass,
                                                 decision.probeRegion);
                    stats->hashScreenKeep++;
                    if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
                    ++tally.keep;
                } else {
                    record_flex_hash_screen_deny(readFeat, localBar, iReadAll,
                                                 sampleOK ? "NEG_PROBE_AMBIG" : "UNMATCHED_TAG");
                    if (sampleOK) {
                        stats->hashScreenDeny++;
                        ++tally.deny;
                    } else {
                        stats->hashScreenSampleReject++;
                        ++tally.sampleReject;
                    }
                }
            } else {
                ++tally.miss;
                if (!noAlign) {
                    EnrichedPacket ep;
                    std::memcpy(ep.name, name, nameLength);
                    ep.name[nameLength] = '\0';
                    std::memcpy(ep.seq[0], seq0, readLen0);
                    std::memcpy(ep.seq[1], seq1, readLen1);
                    std::memcpy(ep.qual[0], qual0, readLen0);
                    std::memcpy(ep.qual[1], qual1, readLen1);
                    ep.seq[0][readLen0] = '\0';
                    ep.seq[1][readLen1] = '\0';
                    ep.qual[0][readLen0] = '\0';
                    ep.qual[1][readLen1] = '\0';
                    ep.readLen[0] = readLen0;
                    ep.readLen[1] = readLen1;
                    ep.iReadAll = iReadAll;
                    ep.laneId = static_cast<uint8_t>(laneId);
                    ep.readFilesIndex = static_cast<uint32_t>(laneId);
                    ep.readFilter = record.read_filter;
                    ep.eof = false;
                    ep.cbMatch = -1;
                    ep.cbMatchIndN = 0;
                    ep.umiB = 0;
                    ep.detectedSampleToken = detectedSampleToken;
                    ep.hashScreenSampleIdx = 0;
                    enqueueForAlign(st, std::move(ep), RA);
                }
            }

            tally.countRead();
        }
    }

    return nReads;
}

static star::input::InputSourcePlan makeSingleCbqLanePlan(const std::string &cbqPath) {
    std::vector<std::vector<std::string>> readFiles(1);
    readFiles[0].push_back(cbqPath);
    return star::input::make_cbq_input_source_plan(readFiles, std::vector<std::string>(), 2);
}

static uint64_t processOneCbqLane(
    FlexPipelineState *st, Parameters &P, int laneId,
    const std::string &cbqPath,
    SoloReadFeature *readFeat, Stats *stats,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign = false,
    ReadAlign *RA = nullptr)
{
    const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
    star::input::CbqInputModule module;
    std::string inputError;
    if (!module.configure(plan, &inputError) || !module.open(&inputError)) {
        P.inOut->logMain << "ERROR: Flex CBQ pipeline could not open " << cbqPath
                         << ": " << inputError << "\n" << std::flush;
        return 0;
    }
    const uint64_t nReads = processCbqModuleRecords(
        st, P, laneId, cbqPath, module, readFeat, stats,
        sampleDet, sampleDetReady, localBar, noAlign, false, 0, RA);
    module.close();
    return nReads;
}

static uint64_t processOneCbqRange(
    FlexPipelineState *st, Parameters &P,
    const FlexCbqRangeTask &task,
    const std::string &cbqPath,
    SoloReadFeature *readFeat, Stats *stats,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign, ReadAlign *RA)
{
    const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
    star::input::CbqInputModule module;
    std::string inputError;
    if (!module.configure(plan, &inputError) ||
        !module.open_range(0, task.firstRecord, task.recordCount, &inputError)) {
        P.inOut->logMain << "ERROR: Flex CBQ range pipeline could not open "
                         << cbqPath << " lane=" << task.laneId
                         << " first=" << task.firstRecord
                         << " count=" << task.recordCount
                         << ": " << inputError << "\n" << std::flush;
        return 0;
    }
    const uint64_t nReads = processCbqModuleRecords(
        st, P, task.laneId, cbqPath, module, readFeat, stats,
        sampleDet, sampleDetReady, localBar, noAlign, true, task.globalFirst, RA);
    module.close();
    return nReads;
}

bool flexPrepareCbqRangeTasks(FlexPipelineState *state, Parameters &P,
                              int nWorkers, std::string *reason) {
    (void)P;
    if (state == nullptr || nWorkers <= 0) {
        if (reason != nullptr) {
            *reason = "invalid range task planner inputs";
        }
        return false;
    }
    state->cbqRangeTasks.clear();
    state->nextCbqRangeIdx.store(0, std::memory_order_relaxed);

    std::vector<uint64_t> laneCounts(static_cast<size_t>(state->nLanes), 0);
    uint64_t totalRecords = 0;
    for (int lane = 0; lane < state->nLanes; ++lane) {
        const std::string &cbqPath = state->laneFiles[static_cast<size_t>(lane)].r2path;
        const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
        star::input::CbqInputModule module;
        std::string inputError;
        if (!module.configure(plan, &inputError) ||
            !module.open_range(0, 0, std::numeric_limits<uint64_t>::max(), &inputError)) {
            if (reason != nullptr) {
                *reason = "lane " + std::to_string(lane) + ": " + inputError;
            }
            return false;
        }
        laneCounts[static_cast<size_t>(lane)] = module.current_lane_record_count();
        totalRecords += laneCounts[static_cast<size_t>(lane)];
        module.close();
    }

    if (totalRecords == 0) {
        if (reason != nullptr) {
            *reason = "no CBQ records";
        }
        return false;
    }

    const uint64_t chunkSize =
        (totalRecords + static_cast<uint64_t>(nWorkers) - 1U) /
        static_cast<uint64_t>(nWorkers);
    for (int worker = 0; worker < nWorkers; ++worker) {
        const uint64_t globalFirst = static_cast<uint64_t>(worker) * chunkSize;
        if (globalFirst >= totalRecords) {
            break;
        }
        const uint64_t globalEnd = std::min(totalRecords, globalFirst + chunkSize);
        uint64_t laneGlobalFirst = 0;
        for (int lane = 0; lane < state->nLanes; ++lane) {
            const uint64_t laneCount = laneCounts[static_cast<size_t>(lane)];
            const uint64_t laneGlobalEnd = laneGlobalFirst + laneCount;
            if (laneCount > 0 && globalFirst < laneGlobalEnd && globalEnd > laneGlobalFirst) {
                const uint64_t overlapFirst = std::max(globalFirst, laneGlobalFirst);
                const uint64_t overlapEnd = std::min(globalEnd, laneGlobalEnd);
                FlexCbqRangeTask task;
                task.laneId = lane;
                task.firstRecord = overlapFirst - laneGlobalFirst;
                task.recordCount = overlapEnd - overlapFirst;
                task.globalFirst = overlapFirst;
                state->cbqRangeTasks.push_back(task);
            }
            laneGlobalFirst = laneGlobalEnd;
            if (laneGlobalFirst >= globalEnd) {
                break;
            }
        }
    }

    if (reason != nullptr) {
        std::ostringstream out;
        out << state->cbqRangeTasks.size() << " ranges across "
            << state->nLanes << " lanes and " << totalRecords << " records";
        *reason = out.str();
    }
    return !state->cbqRangeTasks.empty();
}

bool flexPrepareBgzfRangeTasks(FlexPipelineState *state, Parameters &P,
                               int nWorkers, std::string *reason,
                               bool *fatalError) {
    if (fatalError != nullptr) {
        *fatalError = false;
    }
    if (state == nullptr || nWorkers <= 0) {
        if (reason != nullptr) {
            *reason = "invalid BGZF range task planner inputs";
        }
        if (fatalError != nullptr) {
            *fatalError = true;
        }
        return false;
    }
    state->bgzfLanes.assign(static_cast<size_t>(state->nLanes), FlexBgzfLane());
    state->bgzfRangeTasks.clear();
    state->nextBgzfRangeIdx.store(0, std::memory_order_relaxed);
    state->bgzfReaderWorkers = 0;
    state->bgzfRangeActive = false;

    const std::string mode = P.readFilesBgzfMode;
    if (mode == "off") {
        if (reason != nullptr) {
            *reason = "disabled by --readFilesBgzfMode off";
        }
        return false;
    }
    if (mode != "auto" && mode != "range") {
        if (reason != nullptr) {
            *reason = "unrecognized --readFilesBgzfMode=" + mode;
        }
        if (fatalError != nullptr) {
            *fatalError = true;
        }
        return false;
    }
    if (P.readFilesTypeN != 1 || P.readFilesNames.size() < 2) {
        if (reason != nullptr) {
            *reason = "BGZF range mode requires paired FASTQ input";
        }
        if (fatalError != nullptr) {
            *fatalError = mode == "range";
        }
        return false;
    }

    const int configuredWorkers = P.bgzfReaderThreads == 0
        ? nWorkers : P.bgzfReaderThreads;
    std::vector<int> rangeLaneIds;
    for (int lane = 0; lane < state->nLanes; ++lane) {
        if (P.readFilesNames[0].size() <= static_cast<size_t>(lane) ||
            P.readFilesNames[1].size() <= static_cast<size_t>(lane)) {
            if (reason != nullptr) {
                *reason = "FASTQ mate list is incomplete at lane " + std::to_string(lane);
            }
            if (fatalError != nullptr) {
                *fatalError = true;
            }
            return false;
        }
        const std::string &mate0Path = P.readFilesNames[0][static_cast<size_t>(lane)];
        const std::string &mate1Path = P.readFilesNames[1][static_cast<size_t>(lane)];
        star::input::BgzfDetection detections[2];
        std::string inputError;
        if (!star::input::detect_bgzf(mate0Path, &detections[0], &inputError) ||
            !star::input::detect_bgzf(mate1Path, &detections[1], &inputError)) {
            if (reason != nullptr) {
                *reason = "lane " + std::to_string(lane) + ": " + inputError;
            }
            if (fatalError != nullptr) {
                *fatalError = true;
            }
            return false;
        }
        const bool pairedBgzf = detections[0].isBgzf && detections[1].isBgzf;
        if (!pairedBgzf) {
            if (mode == "range") {
                if (reason != nullptr) {
                    *reason = "lane " + std::to_string(lane) +
                              " has a non-BGZF mate";
                }
                if (fatalError != nullptr) {
                    *fatalError = true;
                }
                return false;
            }
            continue;
        }

        state->bgzfLanes[static_cast<size_t>(lane)].range = true;
        rangeLaneIds.push_back(lane);
    }

    if (rangeLaneIds.empty()) {
        if (reason != nullptr) {
            *reason = "no lane has two BGZF mates";
        }
        return false;
    }

    const int streamCount = static_cast<int>(rangeLaneIds.size()) * 2;
    std::vector<uint64_t> streamBytes(static_cast<size_t>(streamCount), 0);
    for (size_t laneIndex = 0; laneIndex < rangeLaneIds.size(); ++laneIndex) {
        const int lane = rangeLaneIds[laneIndex];
        for (int mate = 0; mate < 2; ++mate) {
            const std::string &path =
                P.readFilesNames[static_cast<size_t>(mate)][static_cast<size_t>(lane)];
            struct stat info;
            if (::stat(path.c_str(), &info) != 0 || info.st_size < 0) {
                if (reason != nullptr) {
                    *reason = "could not determine BGZF stream size for " + path;
                }
                if (fatalError != nullptr) {
                    *fatalError = true;
                }
                return false;
            }
            streamBytes[laneIndex * 2 + static_cast<size_t>(mate)] =
                static_cast<uint64_t>(info.st_size);
        }
    }

    // Give every stream an asynchronous inflater when the budget permits,
    // then place remaining workers where they reduce compressed bytes per
    // worker the most. A zero-worker stream remains valid and inflates
    // synchronously in its fused consumer.
    std::vector<uint32_t> streamThreads(static_cast<size_t>(streamCount), 0);
    int remainingWorkers = configuredWorkers;
    if (remainingWorkers >= streamCount) {
        std::fill(streamThreads.begin(), streamThreads.end(), 1U);
        remainingWorkers -= streamCount;
    }
    while (remainingWorkers-- > 0) {
        size_t best = 0;
        long double bestBytesPerWorker = -1.0L;
        for (size_t stream = 0; stream < streamBytes.size(); ++stream) {
            const long double bytesPerWorker =
                static_cast<long double>(streamBytes[stream]) /
                static_cast<long double>(streamThreads[stream] + 1U);
            if (bytesPerWorker > bestBytesPerWorker) {
                best = stream;
                bestBytesPerWorker = bytesPerWorker;
            }
        }
        ++streamThreads[best];
    }

    int streamIndex = 0;
    uint32_t mateWorkers[2] = {0, 0};
    for (size_t laneIndex = 0; laneIndex < rangeLaneIds.size(); ++laneIndex) {
        const int lane = rangeLaneIds[laneIndex];
        star::input::BgzfStarAdapterOptions options;
        options.lane_index = static_cast<uint32_t>(lane);
        options.mate0_reader_threads = streamThreads[static_cast<size_t>(streamIndex++)];
        options.mate1_reader_threads = streamThreads[static_cast<size_t>(streamIndex++)];
        options.store_mate0_quality = P.pSolo.flexNoAlign == 0;
        // Fused Flex pairs mates by ordered record ordinal and equal counts.
        // Avoid copying and comparing the barcode-read name in production.
        options.validate_read_names = false;
        mateWorkers[0] += options.mate0_reader_threads;
        mateWorkers[1] += options.mate1_reader_threads;
        options.crc_check = P.bgzfCrcCheck == 1;
        if (state->dynamicPermitsEnabled) {
            options.inflate_permit_hooks.context = &kFlexBgzfPermitContext;
            options.inflate_permit_hooks.acquire = &flexBgzfPermitAcquire;
            options.inflate_permit_hooks.release = &flexBgzfPermitRelease;
        }

        FlexBgzfLane &lanePlan = state->bgzfLanes[static_cast<size_t>(lane)];
        lanePlan.adapter.reset(new star::input::BgzfStarAdapter());
        std::string inputError;
        if (!lanePlan.adapter->open(
                P.readFilesNames[0][static_cast<size_t>(lane)],
                P.readFilesNames[1][static_cast<size_t>(lane)],
                options, &inputError)) {
            if (reason != nullptr) {
                *reason = "lane " + std::to_string(lane) + ": " + inputError;
            }
            if (fatalError != nullptr) {
                *fatalError = true;
            }
            return false;
        }
    }

    // Give every fused consumer one long-lived lane claim. Multiple claims can
    // share an adapter: only ordered parse/pair is serialized, while the Flex
    // work after each returned record remains parallel.
    state->bgzfReaderWorkers = nWorkers;
    for (size_t laneIndex = 0; laneIndex < rangeLaneIds.size(); ++laneIndex) {
        FlexBgzfRangeTask task;
        task.laneId = rangeLaneIds[laneIndex];
        state->bgzfRangeTasks.push_back(task);
    }
    for (int worker = static_cast<int>(rangeLaneIds.size()); worker < nWorkers; ++worker) {
        FlexBgzfRangeTask task;
        task.laneId = rangeLaneIds[static_cast<size_t>(worker) % rangeLaneIds.size()];
        state->bgzfRangeTasks.push_back(task);
    }

    state->bgzfRangeActive = true;
    if (reason != nullptr) {
        std::ostringstream out;
        out << rangeLaneIds.size() << " BGZF lanes using "
            << configuredWorkers << " on-demand BC/BSIZE inflate readers and "
            << state->bgzfReaderWorkers << " fused consumers; R2/R1 inflater allocation "
            << mateWorkers[0] << '/' << mateWorkers[1] << "; no pre-index";
        *reason = out.str();
    }
    return true;
}

void *flexLaneReaderFullThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    SoloReadFeature *readFeat = ctx->readFeat;
    Stats *stats = ctx->stats;
    ReadAlign *RA = ctx->RA;

    const bool noAlign = (P.pSolo.flexNoAlign != 0);

    // Preallocate once per thread — reused across all lane claims. Fully-fused
    // callers retain this object after join so exact-CB priors are not lost.
    std::unique_ptr<SoloReadBarcode> fallbackBar;
    if (ctx->readBar == nullptr) {
        fallbackBar.reset(new SoloReadBarcode(P));
    }
    SoloReadBarcode &localBar = ctx->readBar != nullptr ? *ctx->readBar : *fallbackBar;
    SampleDetector *sampleDet = nullptr;
    bool sampleDetReady = false;

    auto ensureSampleDetector = [&]() {
        if (sampleDet == nullptr &&
            !P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
            !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
            sampleDet = new SampleDetector(P.pSolo);
            if (sampleDet->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
                sampleDet->loadProbes(P.pSolo.sampleProbesPath)) {
                sampleDetReady = true;
            } else {
                delete sampleDet;
                sampleDet = nullptr;
            }
        }
    };

    const bool useCbqRangeTasks = (P.readFilesTypeN == 20 && P.cbqInputActive &&
                                   !st->cbqRangeTasks.empty());
    if (useCbqRangeTasks) {
        FlexCbqRangeTask task;
        while (!st->inputFailed.load(std::memory_order_relaxed) &&
               st->claimNextCbqRange(&task)) {
            ensureSampleDetector();
            const std::string &cbqPath =
                st->laneFiles[static_cast<size_t>(task.laneId)].r2path;
            processOneCbqRange(st, P, task, cbqPath, readFeat, stats,
                               sampleDet, sampleDetReady,
                               localBar, noAlign, RA);
        }
    } else {
        if (st->bgzfRangeActive && ctx->threadId < st->bgzfReaderWorkers) {
            FlexBgzfRangeTask task;
            while (!st->inputFailed.load(std::memory_order_relaxed) &&
                   st->claimNextBgzfRange(&task)) {
                ensureSampleDetector();
                processOneBgzfRange(st, P, task, readFeat, stats,
                                    sampleDet, sampleDetReady,
                                    localBar, noAlign, RA);
            }
        }

        // Phase 1: Lane work-stealing loop — claim and process lanes until none remain
        int lane;
        while (!st->inputFailed.load(std::memory_order_relaxed) &&
               (lane = st->claimNextLane()) >= 0) {
            if (st->laneUsesBgzfRange(lane)) {
                continue;
            }
            // Lazy-init SampleDetector on first lane claim (threads that never
            // claim a lane skip the file I/O entirely)
            ensureSampleDetector();

            const std::string &r2path = st->laneFiles[lane].r2path;
            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                processOneCbqLane(st, P, lane, r2path, readFeat, stats,
                                  sampleDet, sampleDetReady, localBar, noAlign, RA);
            } else {
                const std::string &r1path = st->laneFiles[lane].r1path;
                gzFile gzR2 = gzopen(r2path.c_str(), "rb");
                gzFile gzR1 = gzopen(r1path.c_str(), "rb");
                if (!gzR2 || !gzR1) {
                    if (gzR2) gzclose(gzR2);
                    if (gzR1) gzclose(gzR1);
                    continue;
                }
                gzbuffer(gzR2, kGzBufSize);
                gzbuffer(gzR1, kGzBufSize);

                processOneLane(st, P, lane, gzR2, gzR1, readFeat, stats,
                               sampleDet, sampleDetReady, localBar, noAlign, RA);
            }
        }
    }

    // Signal that this reader is done; last one closes alignQ
    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nFusedThreads) {
        st->alignQ.close();
    }

    // Phase 2: Role switch — drain alignQ as an alignment worker (skipped in noAlign mode)
    if (!noAlign && RA != nullptr) {
        EnrichedPacket ep;
        while (st->alignQ.pop(ep)) {
            alignAvailableBatch(st, ep, RA);
        }
    }

    delete sampleDet;
    return nullptr;
}

void *flexTriageThread(void *arg) {
    FlexTriageArgs *ctx = static_cast<FlexTriageArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;

    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;

    ReadPacket rpkt;
    while (st->readerQ.pop(rpkt)) {

        FlexHashScreenDecision decision = cache.classifyReadH0Offset0(
            rpkt.seq[0], rpkt.readLen[0]);

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {
            DecisionPacket dp;
            dp.iReadAll = rpkt.iReadAll;
            dp.readFilesIndex = rpkt.readFilesIndex;
            dp.eof = false;

            // Carry raw R1 (barcode read) for deferred CB/UMI extraction
            const uint32_t bcLen = std::min(rpkt.readLen[1], kFlexPipeBarcodeSeqMax - 1);
            std::memcpy(dp.barcodeSeq, rpkt.seq[1], bcLen);
            dp.barcodeSeq[bcLen] = '\0';
            std::memcpy(dp.barcodeQual, rpkt.qual[1], bcLen);
            dp.barcodeQual[bcLen] = '\0';
            dp.barcodeLen = bcLen;
            std::memcpy(dp.readName, rpkt.name, kFlexPipeNameMax);

            // Carry 8-byte sample tag region from R2
            if (sampleTagOffset + 8 <= rpkt.readLen[0]) {
                std::memcpy(dp.sampleTagSeq, rpkt.seq[0] + sampleTagOffset, 8);
                dp.sampleTagLen = 8;
            } else {
                dp.sampleTagLen = 0;
            }

            if (decision.action == FlexHashScreenDecision::Keep) {
                dp.verdict = DecisionPacket::KEEP;
                dp.geneIdx15 = decision.geneIdx15;
                dp.cacheClass = decision.cacheClass;
                dp.probeRegion = decision.probeRegion;
                dp.denyReason = nullptr;
                st->counters.triageKeep.fetch_add(1);
            } else {
                dp.verdict = DecisionPacket::DENY;
                dp.geneIdx15 = 0;
                dp.cacheClass = 0;
                dp.denyReason = "NEG_PROBE_AMBIG";
                st->counters.triageDeny.fetch_add(1);
            }

            int shard = static_cast<int>(rpkt.iReadAll % static_cast<uint64_t>(st->nSolo));
            st->soloQ[shard]->push(std::move(dp));
        } else {
            // MISS — route to alignment queue (alignment workers do their own CB/UMI)
            EnrichedPacket ep;
            std::memcpy(ep.name, rpkt.name, kFlexPipeNameMax);
            std::memcpy(ep.seq[0], rpkt.seq[0], rpkt.readLen[0] + 1);
            std::memcpy(ep.seq[1], rpkt.seq[1], rpkt.readLen[1] + 1);
            std::memcpy(ep.qual[0], rpkt.qual[0], rpkt.readLen[0] + 1);
            std::memcpy(ep.qual[1], rpkt.qual[1], rpkt.readLen[1] + 1);
            ep.readLen[0] = rpkt.readLen[0];
            ep.readLen[1] = rpkt.readLen[1];
            ep.iReadAll = rpkt.iReadAll;
            ep.laneId = rpkt.laneId;
            ep.readFilesIndex = rpkt.readFilesIndex;
            ep.readFilter = rpkt.readFilter;
            ep.eof = false;

            ep.cbMatch = -1;
            ep.cbMatchIndN = 0;
            ep.umiB = 0;
            ep.detectedSampleToken = 0xFF;
            ep.hashScreenSampleIdx = 0;

            st->counters.triageMiss.fetch_add(1);
            st->alignQ.push(std::move(ep));
        }

        st->counters.readsTotal.fetch_add(1);
    }

    // Last triage thread to finish closes downstream queues
    int triageDone = st->triageFinished.fetch_add(1) + 1;
    if (triageDone == st->nTriage) {
        for (int i = 0; i < st->nSolo; ++i) {
            st->soloQ[i]->close();
        }
        st->alignQ.close();
    }

    return nullptr;
}

void *flexSoloConsumerThread(void *arg) {
    FlexSoloConsumerArgs *ctx = static_cast<FlexSoloConsumerArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    SoloReadFeature *readFeat = ctx->readFeat;
    Stats *stats = ctx->stats;
    const int consumerId = ctx->consumerId;

    // The staged pipeline must retain exact-CB counts just like the fully
    // fused path; ambiguous-CB posterior weights are formed after join.
    std::unique_ptr<SoloReadBarcode> fallbackBar;
    if (ctx->readBar == nullptr) {
        fallbackBar.reset(new SoloReadBarcode(P));
    }
    SoloReadBarcode &localBar = ctx->readBar != nullptr ? *ctx->readBar : *fallbackBar;

    // Solo consumer owns sample detection (deferred from triage)
    SampleDetector *sampleDet = nullptr;
    bool sampleDetReady = false;
    if (!P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
        !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
        sampleDet = new SampleDetector(P.pSolo);
        if (sampleDet->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
            sampleDet->loadProbes(P.pSolo.sampleProbesPath)) {
            sampleDetReady = true;
        } else {
            delete sampleDet;
            sampleDet = nullptr;
        }
    }

    uint8_t packedScratch[4];  // 8 bases in 4-bit BAM = 4 bytes

    // Scratch arrays for getCBandUMI — R2 not needed, provide dummy
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    DecisionPacket dp;
    while (st->soloQ[consumerId]->pop(dp)) {

        // Deferred CB/UMI extraction from raw R1 carried in packet
        char *readSeqPtrs[2]  = { dummySeq, dp.barcodeSeq };
        char *readQualPtrs[2] = { dummyQual, dp.barcodeQual };
        uint64 readLens[2]    = { 0, dp.barcodeLen };
        std::string readNameExtra;

        localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                              dp.readFilesIndex, dp.readName);

        // Deferred sample detection from the 8-byte raw tag carried in the packet
        uint8_t detectedSampleToken = 0xFF;
        if (sampleDetReady && dp.sampleTagLen == 8) {
            nuclPackBAM(dp.sampleTagSeq, reinterpret_cast<char *>(packedScratch), 8);
            uint32_t detIdx = sampleDet->detectSampleFromPackedTag(packedScratch);
            if (detIdx > 0) {
                detectedSampleToken = static_cast<uint8_t>(detIdx & 0x1Fu);
            }
        }
        localBar.detectedSampleToken = detectedSampleToken;

        if (dp.verdict == DecisionPacket::KEEP) {
            record_flex_hash_screen_keep(readFeat, localBar, dp.iReadAll,
                                         dp.geneIdx15, dp.cacheClass, dp.probeRegion);
            stats->hashScreenKeep++;
            if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
        } else {
            record_flex_hash_screen_deny(readFeat, localBar, dp.iReadAll, dp.denyReason);
            stats->hashScreenDeny++;
        }
    }

    delete sampleDet;
    return nullptr;
}

void *flexAlignWorkerThread(void *arg) {
    FlexAlignWorkerArgs *ctx = static_cast<FlexAlignWorkerArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    ReadAlign *RA = ctx->RA;

    EnrichedPacket ep;
    while (st->alignQ.pop(ep)) {
        alignAvailableBatch(st, ep, RA);
    }

    return nullptr;
}

void *flexStatsReporterThread(void *arg) {
    FlexStatsReporterArgs *ctx = static_cast<FlexStatsReporterArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;

    auto tStart = std::chrono::steady_clock::now();
    uint64_t prevTotal = 0;

    while (!st->pipelineDone.load(std::memory_order_relaxed)) {
        sleep(10);
        if (st->pipelineDone.load(std::memory_order_relaxed)) break;

        auto tNow = std::chrono::steady_clock::now();
        int elapsedSec = static_cast<int>(
            std::chrono::duration_cast<std::chrono::seconds>(tNow - tStart).count());

        uint64_t total = st->counters.readsTotal.load(std::memory_order_relaxed);
        uint64_t keep  = st->counters.triageKeep.load(std::memory_order_relaxed);
        uint64_t deny  = st->counters.triageDeny.load(std::memory_order_relaxed);
        uint64_t sampleReject =
            st->counters.triageSampleReject.load(std::memory_order_relaxed);
        uint64_t miss  = st->counters.triageMiss.load(std::memory_order_relaxed);

        size_t rqSize = st->readerQ.size();
        size_t aqSize = st->alignQ.size();

        std::ostringstream soloQStr;
        if (st->nSolo > 0) {
            for (int i = 0; i < st->nSolo; ++i) {
                if (i > 0) soloQStr << ",";
                soloQStr << st->soloQ[i]->size();
            }
        } else {
            soloQStr << "n/a";
        }

        uint64_t delta = total - prevTotal;
        prevTotal = total;

        P.inOut->logMain << "PIPELINE_STATS t=" << elapsedSec << "s"
                         << " total=" << total
                         << " keep=" << keep
                         << " deny=" << deny
                         << " sampleReject=" << sampleReject
                         << " miss=" << miss
                         << " readerQ=" << rqSize
                         << " soloQ=" << soloQStr.str()
                         << " alignQ=" << aqSize
                         << " delta10s=" << delta
                         << "\n" << std::flush;
    }

    return nullptr;
}
