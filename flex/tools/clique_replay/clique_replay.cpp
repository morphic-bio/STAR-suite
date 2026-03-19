// Standalone replay tool for clique UMI correction.
// Reads a binary dump produced by STAR_DUMP_CLIQUE_GROUPS and runs the same
// correctClique logic, emitting summary metrics and optional per-group output.
//
// Usage:
//   clique_replay <dump.bin> [--corrections <out.tsv>] [--override-params minCount ratioThresh maxComponentSize]

#include "UMICorrector.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <string>
#include <vector>

static const uint32_t CLIQUE_DUMP_MAGIC   = 0x434C5155;
static const uint32_t CLIQUE_DUMP_VERSION = 1;

static void decodeUMI12(uint32_t packed, char *out) {
    static const char bases[] = "ACGT";
    for (int i = 11; i >= 0; --i) {
        out[i] = bases[packed & 3];
        packed >>= 2;
    }
    out[12] = '\0';
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <dump.bin> [--corrections <out.tsv>] "
                "[--override-params minCount ratioThresh maxComponentSize]\n", argv[0]);
        return 1;
    }

    const char *dumpFile = argv[1];
    const char *corrFile = nullptr;
    bool overrideParams = false;
    int overrideMinCount = 0;
    double overrideRatioThresh = 0;
    int overrideMaxComponent = 0;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--corrections") == 0 && i + 1 < argc) {
            corrFile = argv[++i];
        } else if (strcmp(argv[i], "--override-params") == 0 && i + 3 < argc) {
            overrideParams = true;
            overrideMinCount = atoi(argv[++i]);
            overrideRatioThresh = atof(argv[++i]);
            overrideMaxComponent = atoi(argv[++i]);
        }
    }

    FILE *f = fopen(dumpFile, "rb");
    if (!f) { perror("fopen"); return 1; }

    uint32_t magic, version;
    fread(&magic, 4, 1, f);
    fread(&version, 4, 1, f);
    if (magic != CLIQUE_DUMP_MAGIC || version != CLIQUE_DUMP_VERSION) {
        fprintf(stderr, "Bad magic/version: got 0x%08X v%u, expected 0x%08X v%u\n",
                magic, version, CLIQUE_DUMP_MAGIC, CLIQUE_DUMP_VERSION);
        fclose(f);
        return 1;
    }

    int32_t minCount;
    double ratioThresh;
    int32_t maxComponentSize;
    fread(&minCount, sizeof(minCount), 1, f);
    fread(&ratioThresh, sizeof(ratioThresh), 1, f);
    fread(&maxComponentSize, sizeof(maxComponentSize), 1, f);

    uint64_t nGroups;
    fread(&nGroups, sizeof(nGroups), 1, f);

    if (overrideParams) {
        fprintf(stderr, "Overriding params: minCount=%d->%d ratioThresh=%.4f->%.4f maxComponentSize=%d->%d\n",
                minCount, overrideMinCount, ratioThresh, overrideRatioThresh,
                maxComponentSize, overrideMaxComponent);
        minCount = overrideMinCount;
        ratioThresh = overrideRatioThresh;
        maxComponentSize = overrideMaxComponent;
    }

    fprintf(stderr, "Dump: %llu groups, params: minCount=%d ratioThresh=%.4f maxComponentSize=%d\n",
            (unsigned long long)nGroups, minCount, ratioThresh, maxComponentSize);

    UMIParams params(minCount, ratioThresh, maxComponentSize);

    FILE *cf = nullptr;
    if (corrFile) {
        cf = fopen(corrFile, "w");
        if (!cf) { perror("corrections fopen"); fclose(f); return 1; }
        fprintf(cf, "groupKey\tur_packed\tub_packed\tur_seq\tub_seq\n");
    }

    uint64_t totalMerges = 0, totalComponents = 0;
    uint64_t totalCapped = 0, totalBelowThreshold = 0;
    uint64_t umisBeforeTotal = 0, umisAfterTotal = 0;
    uint64_t readsBeforeTotal = 0, readsAfterTotal = 0;
    uint64_t totalEntries = 0;
    uint32_t compSizeHist[5] = {0};
    uint32_t maxCompSeen = 0;

    std::vector<UMICount> counts;

    auto t0 = std::chrono::high_resolution_clock::now();

    for (uint64_t g = 0; g < nGroups; ++g) {
        uint64_t groupKey;
        uint32_t nEntries;
        fread(&groupKey, sizeof(groupKey), 1, f);
        fread(&nEntries, sizeof(nEntries), 1, f);

        counts.clear();
        counts.reserve(nEntries);
        for (uint32_t i = 0; i < nEntries; ++i) {
            uint32_t umi24, cnt;
            fread(&umi24, 4, 1, f);
            fread(&cnt, 4, 1, f);
            counts.emplace_back(umi24, cnt);
        }
        totalEntries += nEntries;

        uint64_t groupReads = 0;
        for (const auto &c : counts) groupReads += c.readCount;
        readsBeforeTotal += groupReads;

        UMICorrectionResult result = UMICorrector::correctClique(counts, params);

        umisBeforeTotal += result.uniqueUmisInput;
        umisAfterTotal += result.uniqueUmisPostFilter - result.merges;
        readsAfterTotal += groupReads;
        totalMerges += result.merges;
        totalComponents += result.components;
        totalCapped += result.componentsCapped;
        totalBelowThreshold += result.componentsBelowThreshold;

        for (uint32_t sz : result.componentSizes) {
            if (sz == 1) compSizeHist[0]++;
            else if (sz == 2) compSizeHist[1]++;
            else if (sz == 3) compSizeHist[2]++;
            else if (sz == 4) compSizeHist[3]++;
            else if (sz > 4) compSizeHist[4]++;
            if (sz > maxCompSeen) maxCompSeen = sz;
        }

        if (cf) {
            char urSeq[13], ubSeq[13];
            for (const auto &corr : result.urToUb) {
                if (corr.first == corr.second) continue;
                decodeUMI12(corr.first, urSeq);
                decodeUMI12(corr.second, ubSeq);
                fprintf(cf, "%llu\t%u\t%u\t%s\t%s\n",
                        (unsigned long long)groupKey, corr.first, corr.second, urSeq, ubSeq);
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    fclose(f);
    if (cf) fclose(cf);

    printf("=== Clique Replay Results ===\n");
    printf("Groups:               %llu\n", (unsigned long long)nGroups);
    printf("Total entries:        %llu\n", (unsigned long long)totalEntries);
    printf("UMIs before:          %llu\n", (unsigned long long)umisBeforeTotal);
    printf("UMIs after:           %llu\n", (unsigned long long)umisAfterTotal);
    printf("Merges:               %llu\n", (unsigned long long)totalMerges);
    printf("Components:           %llu\n", (unsigned long long)totalComponents);
    printf("Components capped:    %llu\n", (unsigned long long)totalCapped);
    printf("Below threshold:      %llu\n", (unsigned long long)totalBelowThreshold);
    printf("Reads before:         %llu\n", (unsigned long long)readsBeforeTotal);
    printf("Reads after:          %llu\n", (unsigned long long)readsAfterTotal);
    printf("Max component seen:   %u\n", maxCompSeen);
    printf("Component size hist:  [1]=%u [2]=%u [3]=%u [4]=%u [>4]=%u\n",
           compSizeHist[0], compSizeHist[1], compSizeHist[2], compSizeHist[3], compSizeHist[4]);
    printf("Elapsed:              %.1f ms\n", elapsed_ms);

    return 0;
}
