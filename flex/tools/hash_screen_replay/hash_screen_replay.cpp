// Standalone replay tool for hash screen classification.
// Loads a cache file and a dump produced by STAR_DUMP_HASH_SCREEN, runs flat
// and/or tiered lookups, and compares against the ground truth in the dump.
//
// Usage:
//   hash_screen_replay <cache.bin> <dump.bin> [--mode flat|tiered|both]
//                      [--summary] [--stats] [--diff <out.tsv>]

#include "FlexHashScreenFlat.h"
#include "FlexHashScreenTiered.h"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <set>

static const char* actionName(FlexHashScreenDecision::Action a) {
    switch (a) {
        case FlexHashScreenDecision::Disabled: return "Disabled";
        case FlexHashScreenDecision::Pass:     return "Pass";
        case FlexHashScreenDecision::Keep:     return "Keep";
        case FlexHashScreenDecision::Deny:     return "Deny";
    }
    return "?";
}

struct ActionCounts {
    uint64_t disabled = 0, pass = 0, keep = 0, deny = 0;

    void add(FlexHashScreenDecision::Action a) {
        switch (a) {
            case FlexHashScreenDecision::Disabled: ++disabled; break;
            case FlexHashScreenDecision::Pass:     ++pass; break;
            case FlexHashScreenDecision::Keep:     ++keep; break;
            case FlexHashScreenDecision::Deny:     ++deny; break;
        }
    }

    void print(const char* label, double ms = -1, size_t nReads = 0) const {
        if (ms >= 0 && nReads > 0) {
            printf("%-22s KEEP=%llu  DENY=%llu  PASS=%llu  Disabled=%llu  (%.1f ms, %.0f ns/read)\n",
                   label, (unsigned long long)keep, (unsigned long long)deny,
                   (unsigned long long)pass, (unsigned long long)disabled,
                   ms, ms * 1e6 / nReads);
        } else {
            printf("%-22s KEEP=%llu  DENY=%llu  PASS=%llu  Disabled=%llu\n",
                   label, (unsigned long long)keep, (unsigned long long)deny,
                   (unsigned long long)pass, (unsigned long long)disabled);
        }
    }
};

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s <cache.bin> <dump.bin> [options]\n"
            "  --mode flat|tiered|both   (default: both)\n"
            "  --summary                 print summary stats\n"
            "  --stats                   print tier composition from cache\n"
            "  --h2-scan                 estimate H2 coverage for PASS reads\n"
            "  --diff <out.tsv>          write per-read mismatch details\n",
            argv[0]);
        return 1;
    }

    const char* cachePath = argv[1];
    const char* dumpPath  = argv[2];
    enum { ModeFlat, ModeTiered, ModeBoth } mode = ModeBoth;
    bool showStats = false;
    bool h2Scan = false;
    const char* diffPath = nullptr;

    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            ++i;
            if (strcmp(argv[i], "flat") == 0)        mode = ModeFlat;
            else if (strcmp(argv[i], "tiered") == 0)  mode = ModeTiered;
            else if (strcmp(argv[i], "both") == 0)    mode = ModeBoth;
            else { fprintf(stderr, "Unknown mode: %s\n", argv[i]); return 1; }
        } else if (strcmp(argv[i], "--summary") == 0) {
            // Summary is always printed; flag accepted for compatibility.
        } else if (strcmp(argv[i], "--stats") == 0) {
            showStats = true;
        } else if (strcmp(argv[i], "--h2-scan") == 0) {
            h2Scan = true;
        } else if (strcmp(argv[i], "--diff") == 0 && i + 1 < argc) {
            diffPath = argv[++i];
        }
    }

    // Load cache
    std::vector<Record> allRecords;
    std::string err;
    if (!loadCacheRecords(cachePath, allRecords, &err)) {
        fprintf(stderr, "Failed to load cache: %s\n", err.c_str());
        return 1;
    }
    fprintf(stderr, "Loaded cache: %zu records\n", allRecords.size());

    FlatCache flat;
    TieredCache tiered;

    bool runFlat   = (mode == ModeFlat || mode == ModeBoth);
    bool runTiered = (mode == ModeTiered || mode == ModeBoth);

    if (runFlat)   flat.init(allRecords);
    if (runTiered) tiered.init(allRecords);

    if (showStats) {
        printf("=== Cache Tier Composition ===\n");
        printf("Total records:  %zu\n", allRecords.size());
        if (runTiered) {
            printf("H0 (class 0):   %zu  (%.1f MB, depth %u)\n",
                   tiered.h0Count(), tiered.h0Count() * 24.0 / 1e6,
                   tiered.h0Count() > 0 ? (unsigned)__builtin_clzll(1) - (unsigned)__builtin_clzll(tiered.h0Count()) + 1 : 0);
            printf("H1 (class 1):   %zu  (%.1f MB, depth %u)\n",
                   tiered.h1Count(), tiered.h1Count() * 24.0 / 1e6,
                   tiered.h1Count() > 0 ? (unsigned)__builtin_clzll(1) - (unsigned)__builtin_clzll(tiered.h1Count()) + 1 : 0);
            printf("H2 (class 3):   %zu  (%.1f MB, depth %u)\n",
                   tiered.h2Count(), tiered.h2Count() * 24.0 / 1e6,
                   tiered.h2Count() > 0 ? (unsigned)__builtin_clzll(1) - (unsigned)__builtin_clzll(tiered.h2Count()) + 1 : 0);
            printf("Deny (ambig):   %zu  (%.1f MB, depth %u)\n",
                   tiered.denyCount(), tiered.denyCount() * 24.0 / 1e6,
                   tiered.denyCount() > 0 ? (unsigned)__builtin_clzll(1) - (unsigned)__builtin_clzll(tiered.denyCount()) + 1 : 0);
            printf("Dropped:        %zu  (%.1f MB)\n",
                   tiered.droppedCount(), tiered.droppedCount() * 24.0 / 1e6);
            if (tiered.crossTierDuplicates() > 0)
                printf("WARNING: %zu cross-tier duplicates (same seq+sample across tiers)\n",
                       tiered.crossTierDuplicates());
        }
        printf("\n");
    }

    // Load dump
    std::vector<DumpRecord> dumpRecords;
    if (!loadDumpRecords(dumpPath, dumpRecords, &err)) {
        fprintf(stderr, "Failed to load dump: %s\n", err.c_str());
        return 1;
    }
    fprintf(stderr, "Loaded dump: %zu reads\n", dumpRecords.size());

    FILE* diffFp = nullptr;
    if (diffPath) {
        diffFp = fopen(diffPath, "w");
        if (!diffFp) { perror("diff fopen"); return 1; }
        fprintf(diffFp, "readIdx\treadLen\tsampleIdx\tmismatch_type\t"
                        "truth_action\tflat_action\ttiered_action\t"
                        "truth_gene\tflat_gene\ttiered_gene\t"
                        "truth_negCode\tflat_negCode\ttiered_negCode\t"
                        "truth_class\tflat_class\ttiered_class\n");
    }

    // Warm up (touch cache lines, exclude from timing)
    if (!dumpRecords.empty()) {
        if (runFlat)
            flat.classifyRead(dumpRecords[0].readSeq.c_str(),
                              dumpRecords[0].readLen, dumpRecords[0].sampleIdx);
        if (runTiered)
            tiered.classifyRead(dumpRecords[0].readSeq.c_str(),
                                dumpRecords[0].readLen, dumpRecords[0].sampleIdx);
    }

    ActionCounts truthCounts, flatCounts, tieredCounts;
    uint64_t flatVsTruthMismatch = 0;
    uint64_t tieredVsTruthMismatch = 0;
    uint64_t flatVsTieredMismatch = 0;
    uint64_t disabledSkipped = 0;

    // Disabled reads are counted but not compared against cache lookups
    // (the cache would return Pass/Keep/Deny, not Disabled — Disabled is
    // a STAR routing decision, not a cache decision).

    // Timed flat pass
    auto tFlatStart = std::chrono::high_resolution_clock::now();
    auto tFlatEnd = tFlatStart;
    if (runFlat) {
        tFlatStart = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < dumpRecords.size(); ++i) {
            const DumpRecord& dr = dumpRecords[i];
            if (dr.truth.action == FlexHashScreenDecision::Disabled) continue;
            FlexHashScreenDecision d = flat.classifyRead(
                dr.readSeq.c_str(), dr.readLen, dr.sampleIdx);
            flatCounts.add(d.action);
            if (d != dr.truth) ++flatVsTruthMismatch;
        }
        tFlatEnd = std::chrono::high_resolution_clock::now();
    }

    // Timed tiered pass
    auto tTieredStart = std::chrono::high_resolution_clock::now();
    auto tTieredEnd = tTieredStart;
    if (runTiered) {
        tTieredStart = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < dumpRecords.size(); ++i) {
            const DumpRecord& dr = dumpRecords[i];
            if (dr.truth.action == FlexHashScreenDecision::Disabled) continue;
            FlexHashScreenDecision d = tiered.classifyRead(
                dr.readSeq.c_str(), dr.readLen, dr.sampleIdx);
            tieredCounts.add(d.action);
            if (d != dr.truth) ++tieredVsTruthMismatch;
        }
        tTieredEnd = std::chrono::high_resolution_clock::now();
    }

    // Diff pass: emit per-read details for ALL mismatch types
    for (size_t i = 0; i < dumpRecords.size(); ++i) {
        const DumpRecord& dr = dumpRecords[i];
        truthCounts.add(dr.truth.action);

        if (dr.truth.action == FlexHashScreenDecision::Disabled) {
            ++disabledSkipped;
            continue;
        }

        FlexHashScreenDecision dF, dT;
        bool haveF = false, haveT = false;

        if (runFlat) {
            dF = flat.classifyRead(dr.readSeq.c_str(), dr.readLen, dr.sampleIdx);
            haveF = true;
        }
        if (runTiered) {
            dT = tiered.classifyRead(dr.readSeq.c_str(), dr.readLen, dr.sampleIdx);
            haveT = true;
        }

        if (mode == ModeBoth && haveF && haveT && dF != dT)
            ++flatVsTieredMismatch;

        if (diffFp) {
            bool anyMismatch = false;
            const char* mismatchType = "";

            if (haveF && dF != dr.truth) {
                anyMismatch = true;
                mismatchType = "flat_vs_truth";
            }
            if (haveT && dT != dr.truth) {
                anyMismatch = true;
                mismatchType = "tiered_vs_truth";
            }
            if (haveF && haveT && dF != dT) {
                anyMismatch = true;
                mismatchType = "flat_vs_tiered";
            }
            // Most specific type wins (flat_vs_tiered is the most useful)
            if (haveF && dF != dr.truth && haveT && dT != dr.truth)
                mismatchType = "both_vs_truth";
            if (haveF && haveT && dF != dT && (dF != dr.truth || dT != dr.truth))
                mismatchType = "diverged";

            if (anyMismatch) {
                fprintf(diffFp, "%zu\t%u\t%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                        i, dr.readLen, dr.sampleIdx, mismatchType,
                        actionName(dr.truth.action),
                        haveF ? actionName(dF.action) : "-",
                        haveT ? actionName(dT.action) : "-",
                        dr.truth.geneIdx15,
                        haveF ? dF.geneIdx15 : 0,
                        haveT ? dT.geneIdx15 : 0,
                        dr.truth.negativeCode,
                        haveF ? dF.negativeCode : 0u,
                        haveT ? dT.negativeCode : 0u,
                        dr.truth.cacheClass,
                        haveF ? dF.cacheClass : 0u,
                        haveT ? dT.cacheClass : 0u);
            }
        }
    }

    if (diffFp) fclose(diffFp);

    double flatMs = std::chrono::duration<double, std::milli>(tFlatEnd - tFlatStart).count();
    double tieredMs = std::chrono::duration<double, std::milli>(tTieredEnd - tTieredStart).count();
    size_t activeReads = dumpRecords.size() - disabledSkipped;

    // Print results
    printf("=== Hash Screen Replay Results ===\n");
    printf("Total reads:          %zu\n", dumpRecords.size());
    if (disabledSkipped > 0)
        printf("Disabled (skipped):   %llu\n", (unsigned long long)disabledSkipped);
    printf("Active reads:         %zu\n", activeReads);
    printf("\n");
    truthCounts.print("Ground truth:");

    if (runFlat) {
        flatCounts.print("Flat:", flatMs, activeReads);
        printf("Flat vs truth:        %llu mismatches\n",
               (unsigned long long)flatVsTruthMismatch);
    }
    if (runTiered) {
        tieredCounts.print("Tiered:", tieredMs, activeReads);
        printf("Tiered vs truth:      %llu mismatches\n",
               (unsigned long long)tieredVsTruthMismatch);
    }
    if (mode == ModeBoth) {
        printf("Flat vs tiered:       %llu mismatches\n",
               (unsigned long long)flatVsTieredMismatch);
    }

    bool ok = true;
    if (runFlat && flatVsTruthMismatch > 0) ok = false;
    if (runTiered && tieredVsTruthMismatch > 0) ok = false;
    if (mode == ModeBoth && flatVsTieredMismatch > 0) ok = false;

    printf("\nResult: %s\n", ok ? "PASS" : "FAIL");

    // ── H2 scan: estimate how many PASS reads are within Hamming distance 2
    //    of a known probe (H0 record). This measures the potential uplift from
    //    adding H2 records to the cache.
    if (h2Scan) {
        printf("\n=== H2 Coverage Scan ===\n");
        auto tH2Start = std::chrono::high_resolution_clock::now();

        // Build set of unique probe sequences from H0 records (cacheClass==0).
        struct SeqPair {
            uint64_t lo, hi;
            bool operator<(const SeqPair& o) const {
                return lo < o.lo || (lo == o.lo && hi < o.hi);
            }
        };
        std::set<SeqPair> uniqueSet;
        for (const auto& r : allRecords) {
            if (r.cacheClass == 0)
                uniqueSet.insert({r.seqLo, r.seqHi});
        }
        std::vector<SeqPair> probes(uniqueSet.begin(), uniqueSet.end());
        printf("Unique H0 probe sequences: %zu\n", probes.size());

        // Base-level Hamming distance between two 50bp packed sequences.
        auto baseHamming = [](uint64_t loA, uint64_t hiA,
                              uint64_t loB, uint64_t hiB) -> int {
            uint64_t xLo = loA ^ loB;
            uint64_t xHi = hiA ^ hiB;
            // Collapse each 2-bit pair: any nonzero bit means the base differs.
            uint64_t dLo = (xLo | (xLo >> 1)) & 0x5555555555555555ULL;
            uint64_t dHi = (xHi | (xHi >> 1)) & 0x5555555555555555ULL;
            return __builtin_popcountll(dLo) + __builtin_popcountll(dHi);
        };

        uint64_t passReads = 0;
        uint64_t passEncodeFail = 0;
        uint64_t h2Hits = 0;     // within Hamming distance 2 of a probe
        uint64_t h2Dist[4] = {}; // breakdown: dist 0, 1, 2, 3+

        for (size_t i = 0; i < dumpRecords.size(); ++i) {
            const DumpRecord& dr = dumpRecords[i];
            if (dr.truth.action != FlexHashScreenDecision::Pass) continue;
            ++passReads;

            int bestDist = 999;
            bool anyEncoded = false;

            for (size_t oi = 0; oi < kNumOffsets; ++oi) {
                const int32_t start = kProbeStartOffset + kRelativeProbeOffsets[oi];
                if (start < 0) continue;
                const uint32_t off = static_cast<uint32_t>(start);
                if (off + kCacheKmerLength > dr.readLen) continue;

                uint64_t seqLo = 0, seqHi = 0;
                if (!encodeWindow(dr.readSeq.c_str(), off, seqLo, seqHi)) continue;
                anyEncoded = true;

                for (const auto& p : probes) {
                    int d = baseHamming(seqLo, seqHi, p.lo, p.hi);
                    if (d < bestDist) {
                        bestDist = d;
                        if (d == 0) break;
                    }
                }
                if (bestDist == 0) break;
            }

            if (!anyEncoded) {
                ++passEncodeFail;
                continue;
            }

            if (bestDist <= 2) ++h2Hits;
            if (bestDist <= 3)
                ++h2Dist[bestDist < 3 ? bestDist : 3];
            else
                ++h2Dist[3];

            if ((passReads & 0xFFFFF) == 0)
                fprintf(stderr, "  H2 scan: %llu / %llu PASS reads...\r",
                        (unsigned long long)passReads,
                        (unsigned long long)(truthCounts.pass));
        }

        auto tH2End = std::chrono::high_resolution_clock::now();
        double h2Ms = std::chrono::duration<double, std::milli>(tH2End - tH2Start).count();

        printf("PASS reads scanned:     %llu\n", (unsigned long long)passReads);
        printf("  encode failures (N):  %llu\n", (unsigned long long)passEncodeFail);
        printf("  H2 candidates (d≤2):  %llu  (%.1f%% of PASS)\n",
               (unsigned long long)h2Hits,
               passReads > 0 ? 100.0 * h2Hits / passReads : 0.0);
        printf("  Distance breakdown:\n");
        printf("    d=0 (H0 miss?):     %llu\n", (unsigned long long)h2Dist[0]);
        printf("    d=1 (H1 miss?):     %llu\n", (unsigned long long)h2Dist[1]);
        printf("    d=2 (H2 candidate): %llu\n", (unsigned long long)h2Dist[2]);
        printf("    d≥3 (beyond H2):    %llu\n", (unsigned long long)h2Dist[3]);
        printf("  Scan time: %.1f s (%.0f µs/read)\n",
               h2Ms / 1000.0,
               passReads > 0 ? h2Ms * 1000.0 / passReads : 0.0);

        uint64_t potentialNewKeep = h2Dist[2];
        printf("\n  Potential H2 uplift: +%llu KEEP  (%.1f%% of current PASS → aligner)\n",
               (unsigned long long)potentialNewKeep,
               passReads > 0 ? 100.0 * potentialNewKeep / passReads : 0.0);
        printf("  With H2, aligner-bound reads: %llu → %llu  (%.1f%% reduction)\n",
               (unsigned long long)passReads,
               (unsigned long long)(passReads - potentialNewKeep),
               passReads > 0 ? 100.0 * potentialNewKeep / passReads : 0.0);
    }

    return ok ? 0 : 1;
}
