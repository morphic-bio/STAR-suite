#include "SoloReadInfoSink.h"
#include "SoloFeature.h"
#include "SoloReadInfoLoader.h"
#include <unordered_set>
#include <sstream>

#ifdef DEBUG_CB_UB_PARITY
// Optional tracing of specific readIds via STAR_DEBUG_TRACE_READS=1,2,3
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
            // ignore malformed entries
        }
    }
    return out;
}
static const std::unordered_set<uint32_t> g_traceReads = buildTraceReadSetSink();
#endif

void MinimalSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    // Match CountingSink filter: only persist records with valid status and valid feature
    // This keeps both modes in sync and preserves legacy behavior for no-feature reads
    if (rec.status != 1 || rec.featureId == (uint32_t)-1) {
        // Explicitly write sentinel for rejected/no-feature reads to ensure binary writer sees status=0
        feature.recordReadInfo((uint32_t)rec.readId, 0, 0, 0);
        return;
    }
    
    // Pass original values; recordReadInfo will translate per-backend
    feature.recordReadInfo((uint32_t)rec.readId, rec.cbIdx, rec.umi, rec.status);
}

void MinimalSink::finalize(SoloFeature &feature) {
    (void)feature;
}

void CountingSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    // Buffer counted records per WL CB for later materialization into rGeneUMI.
    // Only accept good CB/UMI and valid featureId.
    if (rec.status != 1) return;
    if (rec.featureId == (uint32_t)-1) return;
    // Guard: if no whitelist, skip safely
    if (feature.pSolo.cbWLsize==0) return;
    if (perWL.empty()) perWL.assign(feature.pSolo.cbWLsize, {});
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
    // Hard guard: a readId must not appear under different CBs in the buffered stream
    {
        static std::unordered_map<uint32_t,uint32_t> readToCb;
        uint32_t ridKey = (rec.readIndex != (uint32_t)-1) ? rec.readIndex : (uint32_t)rec.readId;
        auto it = readToCb.find(ridKey);
        if (it == readToCb.end()) {
            readToCb.emplace(ridKey, rec.cbIdx);
        } else if (it->second != rec.cbIdx) {
            fprintf(stderr, "[ERROR] Conflicting CB for readId=%u existing=%u new=%u\n",
                    ridKey, it->second, rec.cbIdx);
            std::exit(1);
        }
    }
    perWL[rec.cbIdx].push_back(rec);
    
    // Collect UR histogram for UMI correction
    if (feature.pSolo.umiCorrectionMode > 0) {
        feature.collectURHistogram(rec.readId, rec.cbIdx, rec.featureId);
    }
}

void CountingSink::finalize(SoloFeature &feature) {
    // Materialize buffered triplets into rGeneUMI/rCBp and run collapseUMIall().
    if (perWL.empty()) return;
    // Reset state for this finalize run
    feature.indCB.clear();
    feature.indCBwl.assign(feature.pSolo.cbWLsize, (uint32) -1);
    feature.nCB = 0;
    feature.setRGUStride(feature.pSolo.readIndexYes[feature.featureType] ? 3u : 2u);

    // Build detected CB list in WL order for determinism
    std::vector<uint32_t> cbWLtoICB(perWL.size(), (uint32_t)-1);
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        if (!perWL[wl].empty()) {
            cbWLtoICB[wl] = feature.nCB;
            feature.indCB.push_back(wl);
            if (wl < feature.indCBwl.size()) feature.indCBwl[wl] = feature.nCB;
            feature.nCB++;
        }
    }
    if (feature.nCB==0) return;

    // Compute per-detected-CB counts and total
    std::vector<uint32_t> counts(feature.nCB, 0);
    uint64_t totalRecs = 0;
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        uint32_t icb = cbWLtoICB[wl];
        if (icb==(uint32_t)-1) continue;
        counts[icb] = (uint32_t)perWL[wl].size();
        totalRecs += perWL[wl].size();
    }
    #ifdef DEBUG_CB_UB_PARITY
    if (feature.parityEnabled) {
        feature.dbgBufferedRecords = totalRecs;
        feature.dbgBufferedCBs = feature.nCB;
    }
    #endif

    // Allocate rGeneUMI and rCBp
    feature.rGeneUMI = new uint32[feature.getRGUStride() * totalRecs]();
    feature.rCBp = new uint32*[feature.nCB + 1];
    feature.rCBp[0] = feature.rGeneUMI;
    for (uint32_t i=0; i<feature.nCB; ++i) {
        feature.rCBp[i+1] = feature.rCBp[i] + feature.getRGUStride() * counts[i];
    }

    // Fill triplets per CB in WL order
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        uint32_t icb = cbWLtoICB[wl];
        if (icb==(uint32_t)-1) continue;
        uint32_t *blockStart = feature.rCBp[icb];
        uint32_t *dst = blockStart;
        for (const auto &r : perWL[wl]) {
#ifdef DEBUG_CB_UB_PARITY
            if (!g_traceReads.empty()) {
                uint32_t rid = (r.readIndex != (uint32_t)-1) ? r.readIndex : (uint32_t)r.readId;
                if (g_traceReads.count(rid)) {
                    fprintf(stderr, "[TRACE bucket] read=%u wlIdx=%u icb=%u feature=%u status=%u dst_offset=%td\n",
                            rid, wl, icb, r.featureId, r.status, dst - blockStart);
                }
            }
#endif
            dst[0] = r.featureId;
            dst[1] = r.umi;
            if (feature.getRGUStride()==3) dst[2] = (r.readIndex != (uint32_t)-1) ? r.readIndex : (uint32_t)r.readId;
            dst += feature.getRGUStride();
        }
        // Ensure rCBp[icb] points to the start of the block (already set)
        feature.rCBp[icb] = blockStart;
    }

    // Prepare per-CB read counts for caller
    if (feature.rCBn) { delete[] feature.rCBn; feature.rCBn = nullptr; }
    feature.rCBn = new uint32[feature.nCB];
    for (uint32_t iCB=0; iCB<feature.nCB; ++iCB) feature.rCBn[iCB] = counts[iCB];

    // Clear buffers for reuse and free memory
    for (auto &v : perWL) { std::vector<ReadInfoRecord>().swap(v); }
    perWL.clear();
    perWL.shrink_to_fit();
}
