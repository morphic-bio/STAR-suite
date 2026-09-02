#include "FlexHashCacheGenerate.h"
#include "FlexHashScreen.h"
#include "FlexGdna.h"
#include "Genome.h"
#include "Parameters.h"
#include "ProbeListIndex.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include "alignment_model.h"
#include "SampleDetector.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "klib/khash.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct GenomeProbeRow {
    char seq[51];
    uint16_t geneIdx15 = 0;
    uint32_t chrIdx = 0;
    FlexGdnaRegion probeRegion = FlexGdnaUnknown;
};

struct SeqPairHash {
    size_t operator()(const std::pair<uint64_t, uint64_t>& p) const noexcept {
        return p.first ^ (p.second + 0x9e3779b97f4a7c15ULL + (p.first << 6) + (p.first >> 2));
    }
};

/** Dedup key: H1/H2 use sampleKey==0 (global); H0 uses runtime hash-screen sample index (matches FH01SEQ1 record). */
struct DedupKey {
    uint64_t seqLo = 0;
    uint64_t seqHi = 0;
    uint16_t sampleKey = 0;
    bool operator==(const DedupKey& o) const noexcept {
        return seqLo == o.seqLo && seqHi == o.seqHi && sampleKey == o.sampleKey;
    }
};

struct DedupBucket {
    bool denyAny = false;
    bool hasKeep = false;
    bool conflict = false;
    uint16_t gene = 0;
    uint8_t cacheClass = 0;
    FlexGdnaRegion probeRegion = FlexGdnaUnknown;
};

static inline khint_t dedupKeyHashFn(DedupKey key) {
    const std::pair<uint64_t, uint64_t> p{key.seqLo, key.seqHi};
    return static_cast<khint_t>(SeqPairHash{}(p) ^
                                (static_cast<size_t>(key.sampleKey) + 0x9e3779b97f4a7c15ULL));
}

static inline int dedupKeyEq(DedupKey a, DedupKey b) {
    return a.seqLo == b.seqLo && a.seqHi == b.seqHi && a.sampleKey == b.sampleKey;
}

KHASH_INIT(flexdedup, DedupKey, DedupBucket, 1, dedupKeyHashFn, dedupKeyEq)

static char numToAcgt(char g) {
    switch (static_cast<unsigned char>(g)) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return 'N';
    }
}

static void trimToken(std::string& s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
}

static void parseTiers(const std::string& s, bool& h0, bool& h1, bool& h2) {
    h0 = h1 = h2 = false;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        trimToken(tok);
        for (auto& c : tok) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        if (tok == "H0") {
            h0 = true;
        } else if (tok == "H1") {
            h1 = true;
        } else if (tok == "H2") {
            h2 = true;
        }
    }
    if (!h0 && !h1 && !h2) {
        h0 = h1 = h2 = true;
    }
}

static std::vector<GenomeProbeRow> extractGenomeProbes(const Genome& g, const ProbeListIndex& pli, std::ostream& log) {
    std::vector<GenomeProbeRow> out;
    for (uint i = 0; i < g.nChrReal; ++i) {
        const std::string& nm = g.chrName[i];
        if (nm.size() < 4 || nm.rfind("ENSG", 0) != 0) {
            continue;
        }
        if (g.chrLength[i] != 50) {
            continue;
        }
        const uint64_t pos = g.chrStart[i];
        GenomeProbeRow row;
        for (int j = 0; j < 50; ++j) {
            row.seq[j] = numToAcgt(g.G[pos + j]);
        }
        row.seq[50] = '\0';
        const size_t pipe = nm.find('|');
        const std::string geneId = (pipe == std::string::npos) ? nm : nm.substr(0, pipe);
        row.geneIdx15 = pli.geneIndex15(geneId);
        if (row.geneIdx15 == 0) {
            continue;
        }
        row.chrIdx = i;
        row.probeRegion = FlexGdnaProbeMetadata::instance().regionForProbeId(nm);
        out.push_back(row);
    }
    log << "[HASH-CACHE-GEN] extracted " << out.size() << " probe pseudo-chromosomes (ENSG, 50bp, probe-list hit)\n";
    return out;
}

static void fillR2Layout(char* buf90, const char* var50, const char* tag8, uint32_t tagOffset) {
    for (int i = 0; i < 90; ++i) {
        buf90[i] = 'A';
    }
    std::memcpy(buf90, var50, 50);
    if (tag8 != nullptr && tagOffset + 8 <= 90) {
        std::memcpy(buf90 + tagOffset, tag8, 8);
    }
}

static void buildR1FromParams(char* buf, uint32_t len, const ParametersSolo& ps, const std::string& cbSeq) {
    for (uint32_t i = 0; i < len; ++i) {
        buf[i] = 'A';
    }
    // Solo CB/UMI starts are 1-based; synthetic R1 must match extraction exactly.
    const uint32_t cbStart0 = (ps.cbS > 0) ? (ps.cbS - 1) : 0;
    const uint32_t umiStart0 = (ps.umiS > 0) ? (ps.umiS - 1) : 0;
    const uint32_t ncopy = std::min<uint32_t>(static_cast<uint32_t>(cbSeq.size()), ps.cbL);
    for (uint32_t i = 0; i < ncopy && cbStart0 + i < len; ++i) {
        buf[cbStart0 + i] = cbSeq[i];
    }
    for (uint32_t i = 0; i < ps.umiL && umiStart0 + i < len; ++i) {
        buf[umiStart0 + i] = 'A';
    }
    if (len < DEF_readSeqLengthMax) {
        buf[len] = '\0';
    }
}

// verdict: 1=KEEP, 0=DENY, -1=DEAD (skip, don't store)
static void appendVariantRecord(std::vector<FlexHashScreenCache::Record>& out, const char* var50, uint16_t gene15,
                                uint8_t cacheClass, FlexGdnaRegion probeRegion, int verdict) {
    if (verdict < 0)
        return; // DEAD: unmapped variant, no value in caching
    FlexHashScreenCache::Record r;
    if (!FlexHashScreenCache::encodeProbeWindow(var50, 0, r.seqLo, r.seqHi)) {
        return;
    }
    r.sampleIdx = 0;
    if (verdict > 0) {
        r.resolvedGeneIdx15 = gene15;
        r.cacheClass = cacheClass;
        r.probeRegion = probeRegion;
        r.negativeCode = 0;
    } else {
        r.resolvedGeneIdx15 = 0;
        r.cacheClass = 2;
        r.negativeCode = FlexHashNegProbeAmbig;
    }
    out.push_back(r);
}

static void mergeDedupRecord(khash_t(flexdedup)* buckets, const FlexHashScreenCache::Record& rec) {
    DedupKey key;
    key.seqLo = rec.seqLo;
    key.seqHi = rec.seqHi;
    key.sampleKey = rec.sampleIdx;
    int absent = 0;
    khiter_t it = kh_put(flexdedup, buckets, key, &absent);
    if (absent) {
        kh_val(buckets, it) = DedupBucket();
    }
    DedupBucket& b = kh_val(buckets, it);
    if (rec.cacheClass == 2 || rec.resolvedGeneIdx15 == 0) {
        b.denyAny = true;
        return;
    }
    if (!b.hasKeep) {
        b.hasKeep = true;
        b.gene = static_cast<uint16_t>(rec.resolvedGeneIdx15);
        b.cacheClass = rec.cacheClass;
        b.probeRegion = rec.probeRegion;
    } else {
        if (b.gene != static_cast<uint16_t>(rec.resolvedGeneIdx15)) {
            b.conflict = true;
        }
        b.cacheClass = std::min<uint8_t>(b.cacheClass, rec.cacheClass);
        b.probeRegion = flexGdnaMergeRegion(b.probeRegion, rec.probeRegion);
    }
}

static std::vector<FlexHashScreenCache::Record> finalizeFromBuckets(khash_t(flexdedup)* buckets) {
    std::vector<FlexHashScreenCache::Record> out;
    out.reserve(kh_size(buckets));
    for (khiter_t it = kh_begin(buckets); it != kh_end(buckets); ++it) {
        if (!kh_exist(buckets, it)) {
            continue;
        }
        FlexHashScreenCache::Record r;
        const DedupKey& key = kh_key(buckets, it);
        const DedupBucket& b = kh_val(buckets, it);
        r.seqLo = key.seqLo;
        r.seqHi = key.seqHi;
        if (b.denyAny || b.conflict || !b.hasKeep) {
            r.resolvedGeneIdx15 = 0;
            r.cacheClass = 2;
            r.negativeCode = FlexHashNegProbeAmbig;
            r.sampleIdx = 0;
        } else {
            r.resolvedGeneIdx15 = b.gene;
            r.cacheClass = b.cacheClass;
            r.probeRegion = b.probeRegion;
            r.negativeCode = 0;
            // H0: sampleKey is hash-screen index per sample; H1/H2 use sampleKey==0 (global fallback in findRecord).
            r.sampleIdx = key.sampleKey;
        }
        out.push_back(r);
    }
    return out;
}

} // namespace

void runFlexHashCacheGenerate(Parameters& P, Genome& genome, Transcriptome* transcriptomeMain,
                              const libem::Transcriptome* libemTr) {
    if (transcriptomeMain == nullptr) {
        exitWithError("EXITING: --runMode hashCacheGenerate requires transcriptome (e.g. --soloFeatures Gene).\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    ProbeListIndex probeIdx;
    uint32_t deprecatedCount = 0;
    if (!probeIdx.load(P.pSolo.probeListPath, P.pSolo.removeDeprecated, &deprecatedCount)) {
        exitWithError("EXITING: hashCacheGenerate could not load --soloProbeList\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_PARAMETER, P);
    }

    FlexGdnaProbeMetadata& gdnaMetadata = FlexGdnaProbeMetadata::instance();
    if (!gdnaMetadata.ready()) {
        const std::string probeCsv = FlexGdnaProbeMetadata::discoverProbeCsv(
            P.pSolo.flexGdnaProbeSetPath, P.pSolo.probeListPath, P.pGe.gDir);
        std::string gdnaError;
        if (!probeCsv.empty()) {
            P.pSolo.flexGdnaReady =
                gdnaMetadata.load(probeCsv, P.pSolo.probeListPath, &gdnaError);
        }
        if (!P.pSolo.flexGdnaReady) {
            P.inOut->logMain
                << "[HASH-CACHE-GEN] probe-region metadata unavailable; writing backward-compatible v2 cache"
                << (gdnaError.empty() ? "" : ": " + gdnaError) << "\n";
        }
    } else {
        P.pSolo.flexGdnaReady = true;
    }

    std::vector<GenomeProbeRow> probes = extractGenomeProbes(genome, probeIdx, P.inOut->logMain);
    if (probes.empty()) {
        exitWithError("EXITING: hashCacheGenerate found no ENSG 50bp probe chromosomes with probe-list hits\n", std::cerr,
                      P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    if (P.pSolo.hashCacheParentLimit > 0 && probes.size() > P.pSolo.hashCacheParentLimit) {
        probes.resize(P.pSolo.hashCacheParentLimit);
        P.inOut->logMain << "[HASH-CACHE-GEN] capped to --hashCacheParentLimit=" << P.pSolo.hashCacheParentLimit
                         << " probes\n";
    }

    bool wantH0 = false, wantH1 = false, wantH2 = false;
    parseTiers(P.pSolo.hashCacheTiers, wantH0, wantH1, wantH2);

    std::vector<std::unique_ptr<ReadAlignChunk>> chunks;
    chunks.reserve(P.runThreadN);
    for (int t = 0; t < P.runThreadN; ++t) {
        chunks.emplace_back(new ReadAlignChunk(P, genome, transcriptomeMain, t, libemTr));
    }
    ReadAlign* ra0 = chunks[0]->RA;
    if (ra0 == nullptr || ra0->sampleDet_ == nullptr || !ra0->sampleDet_->ready()) {
        exitWithError(
            "EXITING: hashCacheGenerate requires SampleDetector (--soloSampleWhitelist and --soloSampleProbes)\n",
            std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    const std::string tag8 = ra0->sampleDet_->canonicalFor(1);
    if (tag8.size() != 8) {
        exitWithError("EXITING: hashCacheGenerate could not read 8bp canonical tag for sample index 1 from whitelist\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    if (!P.pSolo.cbWLyes || P.pSolo.cbWLstr.empty()) {
        exitWithError("EXITING: hashCacheGenerate requires a loaded cell barcode whitelist (--soloCBwhitelist)\n", std::cerr,
                      P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    if (P.pSolo.cbumiL == 0 || P.pSolo.cbumiL > DEF_readSeqLengthMax) {
        exitWithError("EXITING: hashCacheGenerate invalid solo CB+UMI length (cbumiL)\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_PARAMETER, P);
    }
    const std::string cb0 = P.pSolo.cbWLstr[0];
    std::vector<char> r1buf(P.pSolo.cbumiL + 1);
    buildR1FromParams(r1buf.data(), P.pSolo.cbumiL, P.pSolo, cb0);

    const uint32_t tagOff = P.pSolo.sampleProbeOffset;
    if (tagOff + 8 > 90) {
        exitWithError("EXITING: hashCacheGenerate requires soloSampleProbeOffset + 8 <= 90\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_PARAMETER, P);
    }

    // Per whitelist sequential row: hash-screen sample index (same chain as ReadAlign_oneRead: detectSampleIndex → token → sampleIndexForToken).
    const uint32_t nSamplesSeq = ra0->sampleDet_->sequentialSampleCount();
    if (nSamplesSeq == 0u) {
        exitWithError("EXITING: hashCacheGenerate found no samples in --soloSampleWhitelist\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_PARAMETER, P);
    }
    std::vector<uint16_t> h0RowSampleIdx(static_cast<size_t>(nSamplesSeq) + 1u, 0);
    char r2detect[90];
    std::vector<uint8_t> packedDetect(static_cast<size_t>((90u + 1u) / 2u));
    uint32_t h0DerivedCount = 0;
    for (uint32_t s = 1u; s <= nSamplesSeq; ++s) {
        const std::string tagS = ra0->sampleDet_->canonicalFor(s);
        if (tagS.size() != 8u) {
            continue;
        }
        fillR2Layout(r2detect, probes[0].seq, tagS.c_str(), tagOff);
        nuclPackBAM(r2detect, reinterpret_cast<char*>(packedDetect.data()), 90u);
        const uint32_t detIdx = ra0->sampleDet_->detectSampleIndex(packedDetect.data(), 90, false);
        if (detIdx > 0u) {
            h0RowSampleIdx[s] = SampleDetector::sampleIndexForToken(static_cast<uint8_t>(detIdx & 0x1Fu));
            if (h0RowSampleIdx[s] != 0u) {
                ++h0DerivedCount;
            }
        }
    }
    if (wantH0 && h0DerivedCount == 0u) {
        exitWithError(
            "EXITING: hashCacheGenerate could not derive any H0 hash-screen sampleIdx from whitelist tags (check "
            "--soloSampleWhitelist / --soloSampleProbes / --soloSampleProbeOffset)\n",
            std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    if (wantH0) {
        P.inOut->logMain << "[HASH-CACHE-GEN] H0 tier: " << h0DerivedCount << " whitelist row(s) with non-zero sampleIdx "
                         << "out of " << nSamplesSeq << " sequential samples (one H0 cache row per derived index)\n";
    }

    std::vector<std::vector<FlexHashScreenCache::Record>> threadRecords(static_cast<size_t>(P.runThreadN));

    P.inOut->logMain << "[HASH-CACHE-GEN] tiers H0=" << wantH0 << " H1=" << wantH1 << " H2=" << wantH2
                     << " threads=" << P.runThreadN << "\n";

    // ---- Pass 1: H0 + H1 ----
#pragma omp parallel num_threads(P.runThreadN)
    {
        const int tid = omp_get_thread_num();
        ReadAlign* RA = chunks[tid]->RA;
        std::vector<FlexHashScreenCache::Record>& local = threadRecords[static_cast<size_t>(tid)];
        char var[51];
        char r2[90];

#pragma omp for schedule(dynamic, 4)
        for (size_t pi = 0; pi < probes.size(); ++pi) {
            const GenomeProbeRow& pr = probes[pi];
            std::memcpy(var, pr.seq, 51);

            if (wantH0) {
                for (uint32_t s = 1u; s <= nSamplesSeq; ++s) {
                    const uint16_t rowS = h0RowSampleIdx[s];
                    if (rowS == 0u) {
                        continue;
                    }
                    FlexHashScreenCache::Record h0;
                    if (FlexHashScreenCache::encodeProbeWindow(var, 0, h0.seqLo, h0.seqHi)) {
                        h0.resolvedGeneIdx15 = pr.geneIdx15;
                        h0.cacheClass = 0;
                        h0.probeRegion = pr.probeRegion;
                        h0.negativeCode = 0;
                        h0.sampleIdx = rowS;
                        local.push_back(h0);
                    }
                }
            }

            if (wantH1) {
                for (int pos = 0; pos < 50; ++pos) {
                    const char refb = var[pos];
                    for (const char* p = "ACGT"; *p; ++p) {
                        const char alt = *p;
                        if (alt == refb) {
                            continue;
                        }
                        var[pos] = alt;
                        fillR2Layout(r2, var, tag8.c_str(), tagOff);
                        const int verdict =
                            RA->flexHashCacheValidateSyntheticPair(r2, 90, r1buf.data(), P.pSolo.cbumiL, pr.geneIdx15);
                        appendVariantRecord(local, var, pr.geneIdx15, 1, pr.probeRegion, verdict);
                        var[pos] = refb;
                    }
                }
            }
        }
    }

    // ---- Build H0/H1 sequence lookup set for H2 pre-check ----
    khash_t(flexdedup)* buckets = kh_init(flexdedup);
    size_t totalH01Records = 0;
    for (const auto& tr : threadRecords) {
        totalH01Records += tr.size();
    }
    if (totalH01Records > 0) {
        kh_resize(flexdedup, buckets, totalH01Records);
    }
    for (const auto& tr : threadRecords) {
        for (const auto& rec : tr) {
            mergeDedupRecord(buckets, rec);
        }
    }

    typedef std::pair<uint64_t, uint64_t> SeqPair;
    std::unordered_set<SeqPair, SeqPairHash> h01SeqSet;
    if (wantH2) {
        h01SeqSet.reserve(kh_size(buckets));
        for (khiter_t it = kh_begin(buckets); it != kh_end(buckets); ++it) {
            if (!kh_exist(buckets, it)) {
                continue;
            }
            const DedupKey& key = kh_key(buckets, it);
            h01SeqSet.emplace(key.seqLo, key.seqHi);
        }
        P.inOut->logMain << "[HASH-CACHE-GEN] H2 pre-check set: " << h01SeqSet.size()
                         << " unique H0/H1 sequences\n";
    }
    threadRecords.clear();
    threadRecords.shrink_to_fit();

    // ---- Pass 2: H2 (KEEP-only, streamed to per-thread temp files) ----
    std::vector<std::string> h2TmpPaths;
    if (wantH2) {
        std::atomic<uint64_t> h2Skipped{0};
        std::atomic<uint64_t> h2Aligned{0};
        std::atomic<uint64_t> h2Kept{0};
        h2TmpPaths.resize(static_cast<size_t>(P.runThreadN));

#pragma omp parallel num_threads(P.runThreadN)
        {
            const int tid = omp_get_thread_num();
            ReadAlign* RA = chunks[tid]->RA;
            char var[51];
            char r2[90];
            uint64_t localSkipped = 0;
            uint64_t localAligned = 0;
            uint64_t localKept = 0;

            std::string tmpPath = P.outFileNamePrefix + "h2_tmp_t" + std::to_string(tid) + ".bin";
            h2TmpPaths[static_cast<size_t>(tid)] = tmpPath;
            std::ofstream tmpOut(tmpPath, std::ios::binary);

#pragma omp for schedule(dynamic, 4)
            for (size_t pi = 0; pi < probes.size(); ++pi) {
                const GenomeProbeRow& pr = probes[pi];
                std::memcpy(var, pr.seq, 51);

                for (int p0 = 0; p0 < 50; ++p0) {
                    const char r0 = var[p0];
                    for (const char* a0 = "ACGT"; *a0; ++a0) {
                        if (*a0 == r0)
                            continue;
                        var[p0] = *a0;
                        for (int p1 = p0 + 1; p1 < 50; ++p1) {
                            const char r1b = var[p1];
                            for (const char* a1 = "ACGT"; *a1; ++a1) {
                                if (*a1 == r1b)
                                    continue;
                                var[p1] = *a1;

                                uint64_t sLo, sHi;
                                if (FlexHashScreenCache::encodeProbeWindow(var, 0, sLo, sHi) &&
                                    h01SeqSet.count(SeqPair(sLo, sHi))) {
                                    ++localSkipped;
                                    var[p1] = r1b;
                                    continue;
                                }

                                fillR2Layout(r2, var, tag8.c_str(), tagOff);
                                const int verdict2 = RA->flexHashCacheValidateSyntheticPair(
                                    r2, 90, r1buf.data(), P.pSolo.cbumiL, pr.geneIdx15);
                                ++localAligned;
                                if (verdict2 > 0) {
                                    FlexHashScreenCache::Record rec;
                                    rec.seqLo = sLo;
                                    rec.seqHi = sHi;
                                    rec.resolvedGeneIdx15 = pr.geneIdx15;
                                    rec.cacheClass = 3;
                                    rec.probeRegion = pr.probeRegion;
                                    rec.negativeCode = 0;
                                    rec.sampleIdx = 0;
                                    tmpOut.write(reinterpret_cast<const char*>(&rec), sizeof(rec));
                                    ++localKept;
                                }
                                var[p1] = r1b;
                            }
                        }
                        var[p0] = r0;
                    }
                }
            }
            tmpOut.close();
            h2Skipped += localSkipped;
            h2Aligned += localAligned;
            h2Kept += localKept;
        }
        P.inOut->logMain << "[HASH-CACHE-GEN] H2: " << h2Aligned.load() << " aligned, "
                         << h2Skipped.load() << " skipped (already in H0/H1), "
                         << h2Kept.load() << " KEEP\n";
    }

    // ---- Stream H2 temp files into dedup buckets ----
    uint64_t h2TmpRecords = 0;
    for (const auto& tmpPath : h2TmpPaths) {
        std::ifstream tmpIn(tmpPath, std::ios::binary | std::ios::ate);
        if (!tmpIn.good()) {
            continue;
        }
        std::streamoff sz = tmpIn.tellg();
        if (sz > 0) {
            h2TmpRecords += static_cast<uint64_t>(sz / static_cast<std::streamoff>(sizeof(FlexHashScreenCache::Record)));
        }
    }
    if (h2TmpRecords > 0) {
        kh_resize(flexdedup, buckets, kh_size(buckets) + static_cast<khint_t>(h2TmpRecords));
    }
    for (const auto& tmpPath : h2TmpPaths) {
        std::ifstream tmpIn(tmpPath, std::ios::binary);
        if (!tmpIn.good())
            continue;
        FlexHashScreenCache::Record rec;
        while (tmpIn.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
            mergeDedupRecord(buckets, rec);
        }
        tmpIn.close();
        std::remove(tmpPath.c_str());
    }

    std::vector<FlexHashScreenCache::Record> finalRecs = finalizeFromBuckets(buckets);
    kh_destroy(flexdedup, buckets);

    std::string err;
    if (!FlexHashScreenCache::writeHashCacheFile(
            P.pSolo.hashCacheOutput, finalRecs, &err, P.pSolo.flexGdnaReady)) {
        ostringstream e;
        e << "EXITING: hash cache write failed: " << err << "\n";
        exitWithError(e.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    P.inOut->logMain << "[HASH-CACHE-GEN] wrote " << finalRecs.size() << " records to " << P.pSolo.hashCacheOutput << "\n";
}
