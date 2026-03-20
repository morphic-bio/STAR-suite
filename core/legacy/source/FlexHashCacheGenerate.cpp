#include "FlexHashCacheGenerate.h"
#include "FlexHashScreen.h"
#include "Genome.h"
#include "Parameters.h"
#include "ProbeListIndex.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"
#include "alignment_model.h"
#include "SampleDetector.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct GenomeProbeRow {
    char seq[51];
    uint16_t geneIdx15 = 0;
    uint32_t chrIdx = 0;
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

struct DedupKeyHash {
    size_t operator()(const DedupKey& k) const noexcept {
        const std::pair<uint64_t, uint64_t> p{k.seqLo, k.seqHi};
        return SeqPairHash{}(p) ^ (static_cast<size_t>(k.sampleKey) + 0x9e3779b97f4a7c15ULL);
    }
};

struct DedupBucket {
    bool denyAny = false;
    bool hasKeep = false;
    bool conflict = false;
    uint16_t gene = 0;
    uint8_t cacheClass = 0;
};

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
    const uint32_t ncopy = std::min<uint32_t>(static_cast<uint32_t>(cbSeq.size()), ps.cbL);
    for (uint32_t i = 0; i < ncopy && ps.cbS + i < len; ++i) {
        buf[ps.cbS + i] = cbSeq[i];
    }
    for (uint32_t i = 0; i < ps.umiL && ps.umiS + i < len; ++i) {
        buf[ps.umiS + i] = 'A';
    }
    if (len < DEF_readSeqLengthMax) {
        buf[len] = '\0';
    }
}

// verdict: 1=KEEP, 0=DENY, -1=DEAD (skip, don't store)
static void appendVariantRecord(std::vector<FlexHashScreenCache::Record>& out, const char* var50, uint16_t gene15,
                                uint8_t cacheClass, int verdict) {
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
        r.negativeCode = 0;
    } else {
        r.resolvedGeneIdx15 = 0;
        r.cacheClass = 2;
        r.negativeCode = FlexHashNegProbeAmbig;
    }
    out.push_back(r);
}

static void mergeDedupBuckets(const std::vector<FlexHashScreenCache::Record>& in,
                              std::unordered_map<DedupKey, DedupBucket, DedupKeyHash>& buckets) {
    buckets.reserve(buckets.size() + in.size());
    for (const auto& rec : in) {
        DedupKey key;
        key.seqLo = rec.seqLo;
        key.seqHi = rec.seqHi;
        key.sampleKey = rec.sampleIdx;
        DedupBucket& b = buckets[key];
        if (rec.cacheClass == 2 || rec.resolvedGeneIdx15 == 0) {
            b.denyAny = true;
            continue;
        }
        if (!b.hasKeep) {
            b.hasKeep = true;
            b.gene = static_cast<uint16_t>(rec.resolvedGeneIdx15);
            b.cacheClass = rec.cacheClass;
        } else {
            if (b.gene != static_cast<uint16_t>(rec.resolvedGeneIdx15)) {
                b.conflict = true;
            }
            b.cacheClass = std::min<uint8_t>(b.cacheClass, rec.cacheClass);
        }
    }
}

static std::vector<FlexHashScreenCache::Record> finalizeFromBuckets(
    const std::unordered_map<DedupKey, DedupBucket, DedupKeyHash>& buckets) {
    std::vector<FlexHashScreenCache::Record> out;
    out.reserve(buckets.size());
    for (const auto& kv : buckets) {
        FlexHashScreenCache::Record r;
        r.seqLo = kv.first.seqLo;
        r.seqHi = kv.first.seqHi;
        const DedupBucket& b = kv.second;
        if (b.denyAny || b.conflict || !b.hasKeep) {
            r.resolvedGeneIdx15 = 0;
            r.cacheClass = 2;
            r.negativeCode = FlexHashNegProbeAmbig;
            r.sampleIdx = 0;
        } else {
            r.resolvedGeneIdx15 = b.gene;
            r.cacheClass = b.cacheClass;
            r.negativeCode = 0;
            // H0: sampleKey is hash-screen index per sample; H1/H2 use sampleKey==0 (global fallback in findRecord).
            r.sampleIdx = kv.first.sampleKey;
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
                        appendVariantRecord(local, var, pr.geneIdx15, 1, verdict);
                        var[pos] = refb;
                    }
                }
            }

            if (wantH2) {
                for (int p0 = 0; p0 < 50; ++p0) {
                    const char r0 = var[p0];
                    for (const char* a0 = "ACGT"; *a0; ++a0) {
                        if (*a0 == r0) {
                            continue;
                        }
                        var[p0] = *a0;
                        for (int p1 = p0 + 1; p1 < 50; ++p1) {
                            const char r1b = var[p1];
                            for (const char* a1 = "ACGT"; *a1; ++a1) {
                                if (*a1 == r1b) {
                                    continue;
                                }
                                var[p1] = *a1;
                                fillR2Layout(r2, var, tag8.c_str(), tagOff);
                                const int verdict2 = RA->flexHashCacheValidateSyntheticPair(r2, 90, r1buf.data(),
                                                                                           P.pSolo.cbumiL, pr.geneIdx15);
                                appendVariantRecord(local, var, pr.geneIdx15, 3, verdict2);
                                var[p1] = r1b;
                            }
                        }
                        var[p0] = r0;
                    }
                }
            }
        }
    }

    std::unordered_map<DedupKey, DedupBucket, DedupKeyHash> buckets;
    for (const auto& tr : threadRecords) {
        mergeDedupBuckets(tr, buckets);
    }
    std::vector<FlexHashScreenCache::Record> finalRecs = finalizeFromBuckets(buckets);

    std::string err;
    if (!FlexHashScreenCache::writeHashCacheFile(P.pSolo.hashCacheOutput, finalRecs, &err)) {
        ostringstream e;
        e << "EXITING: hash cache write failed: " << err << "\n";
        exitWithError(e.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    P.inOut->logMain << "[HASH-CACHE-GEN] wrote " << finalRecs.size() << " records to " << P.pSolo.hashCacheOutput << "\n";
}
