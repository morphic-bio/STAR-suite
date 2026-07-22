#include "SoloReadInfoSink.h"
#include "SoloFeature.h"
#include "SoloReadInfoLoader.h"
#include "SoloMemoryProfile.h"
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <unordered_set>
#include <vector>

#ifdef DEBUG_CB_UB_PARITY
static std::unordered_set<uint32_t> buildTraceReadSetSink() {
    std::unordered_set<uint32_t> out;
    const char* env = std::getenv("STAR_DEBUG_TRACE_READS");
    if (!env || env[0]=='\0') return out;
    std::string s(env);
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try {
            out.insert(static_cast<uint32_t>(std::stoul(tok)));
        } catch (...) {
        }
    }
    return out;
}
static const std::unordered_set<uint32_t> g_traceReads = buildTraceReadSetSink();
#endif

namespace {

bool countingSinkReadToCbEnabled()
{
    static int cached = -1;
    if (cached < 0) {
        const char *env = std::getenv("STAR_DEBUG_COUNTING_SINK_READ_TO_CB");
        cached = (env != nullptr && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return cached > 0;
}

uint64_t bucketCapacityBytes(const std::unordered_map<uint32_t, std::vector<CountingSinkRow>> &buckets)
{
    uint64_t cap = 0;
    for (const auto &kv : buckets) {
        cap += kv.second.capacity() * sizeof(CountingSinkRow);
    }
    return cap;
}

} // namespace

void MinimalSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    if (rec.status == 0 || rec.featureId == (uint32_t)-1) {
        feature.recordReadInfo((uint32_t)rec.readId, 0, 0, 0);
        return;
    }
    feature.recordReadInfo((uint32_t)rec.readId, rec.cbIdx, rec.umi, rec.status);
}

void MinimalSink::finalize(SoloFeature &feature) {
    (void)feature;
}

void CountingSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    if (rec.status != 1) return;
    if (rec.featureId == (uint32_t)-1) return;
    if (feature.pSolo.cbWLsize == 0) return;
    if (rec.cbIdx >= feature.pSolo.cbWLsize) return;

#ifdef DEBUG_CB_UB_PARITY
    if (!g_traceReads.empty()) {
        uint32_t rid = (rec.readIndex != (uint32_t)-1) ? rec.readIndex : (uint32_t)rec.readId;
        if (g_traceReads.count(rid)) {
            fprintf(stderr, "[TRACE buffer] read=%u wlIdx=%u feature=%u umi=%u status=%u readIndex=%u\n",
                    rid, rec.cbIdx, rec.featureId, rec.umi, rec.status, rec.readIndex);
        }
    }
#endif

    if (countingSinkReadToCbEnabled() && rec.readIndex != (uint32_t)-1) {
        auto it = readToCb.find(rec.readIndex);
        if (it == readToCb.end()) {
            readToCb.emplace(rec.readIndex, rec.cbIdx);
        } else if (it->second != rec.cbIdx) {
            fprintf(stderr, "[ERROR] Conflicting CB for readId=%u existing=%u new=%u\n",
                    rec.readIndex, it->second, rec.cbIdx);
            std::exit(1);
        }
    }

    const uint32_t readRef = (rec.readIndex != (uint32_t)-1)
                                 ? rec.readIndex
                                 : static_cast<uint32_t>(rec.readId);
    CountingSinkRow row{rec.featureId, rec.umi, readRef};
    auto &bucket = buckets[rec.cbIdx];
    bucket.push_back(row);
    totalRecords++;
}

void CountingSink::finalize(SoloFeature &feature) {
    if (buckets.empty() || totalRecords == 0) {
        soloMemoryProfileCheckpoint(feature.P.inOut->logMain,
                                  "CountingSink::finalize_skipped_empty_buffer");
        buckets.clear();
        readToCb.clear();
        totalRecords = 0;
        return;
    }

    feature.indCB.clear();
    feature.indCBwl.assign(feature.pSolo.cbWLsize, (uint32)-1);
    feature.nCB = 0;
    feature.setRGUStride(feature.pSolo.readIndexYes[feature.featureType] ? 3u : 2u);

    std::vector<uint32_t> activeWL;
    activeWL.reserve(buckets.size());
    for (const auto &kv : buckets) {
        if (!kv.second.empty()) {
            activeWL.push_back(kv.first);
        }
    }
    std::sort(activeWL.begin(), activeWL.end());

    std::vector<uint32_t> cbWLtoICB(feature.pSolo.cbWLsize, (uint32_t)-1);
    for (uint32_t wl : activeWL) {
        cbWLtoICB[wl] = feature.nCB;
        feature.indCB.push_back(wl);
        if (wl < feature.indCBwl.size()) {
            feature.indCBwl[wl] = feature.nCB;
        }
        feature.nCB++;
    }
    if (feature.nCB == 0) {
        buckets.clear();
        readToCb.clear();
        totalRecords = 0;
        return;
    }

    std::vector<uint32_t> counts(feature.nCB, 0);
    const uint64_t totalRecs = totalRecords;
    for (uint32_t wl : activeWL) {
        const uint32_t icb = cbWLtoICB[wl];
        counts[icb] = static_cast<uint32_t>(buckets[wl].size());
    }

#ifdef DEBUG_CB_UB_PARITY
    if (feature.parityEnabled) {
        feature.dbgBufferedRecords = totalRecs;
        feature.dbgBufferedCBs = feature.nCB;
    }
#endif

    {
        std::ostringstream extra;
        extra << "counting_sink_buffered_records=" << totalRecs
              << " counting_sink_bucket_cbs=" << buckets.size()
              << " counting_sink_bucket_capacity_bytes=" << bucketCapacityBytes(buckets)
              << " counting_sink_readToCb=" << readToCb.size()
              << " counting_sink_detected_cbs=" << feature.nCB
              << " counting_sink_row_bytes=" << sizeof(CountingSinkRow)
              << " rGeneUMI_triplets_est=" << totalRecs
              << " rGeneUMI_bytes_est=" << (totalRecs * feature.getRGUStride() * sizeof(uint32_t));
        soloMemoryProfileCheckpoint(feature.P.inOut->logMain,
                                  "CountingSink::finalize_before_rGeneUMI",
                                  extra.str());
    }

    feature.rGeneUMI = new uint32[feature.getRGUStride() * totalRecs]();
    feature.rCBp = new uint32*[feature.nCB + 1];
    feature.rCBp[0] = feature.rGeneUMI;
    for (uint32_t i = 0; i < feature.nCB; ++i) {
        feature.rCBp[i + 1] = feature.rCBp[i] + feature.getRGUStride() * counts[i];
    }

    const uint32_t stride = feature.getRGUStride();
    for (uint32_t wl : activeWL) {
        const uint32_t icb = cbWLtoICB[wl];
        uint32_t *blockStart = feature.rCBp[icb];
        uint32_t *dst = blockStart;
        auto &bucket = buckets[wl];
        for (const auto &r : bucket) {
            dst[0] = r.featureId;
            dst[1] = r.umi;
            if (stride == 3) {
                dst[2] = r.readRef;
            }
            dst += stride;
        }
        feature.rCBp[icb] = blockStart;
        std::vector<CountingSinkRow>().swap(bucket);
    }
    buckets.clear();
    readToCb.clear();
    totalRecords = 0;

    if (feature.rCBn) {
        delete[] feature.rCBn;
        feature.rCBn = nullptr;
    }
    feature.rCBn = new uint32[feature.nCB];
    for (uint32_t iCB = 0; iCB < feature.nCB; ++iCB) {
        feature.rCBn[iCB] = counts[iCB];
    }

    std::ostringstream extraDone;
    extraDone << "rGeneUMI_triplets=" << totalRecs
              << " rGeneUMI_bytes=" << (totalRecs * stride * sizeof(uint32_t))
              << " rCBp_slots=" << (feature.nCB + 1);
    soloMemoryProfileCheckpoint(feature.P.inOut->logMain,
                              "CountingSink::finalize_after_fill",
                              extraDone.str());
}
