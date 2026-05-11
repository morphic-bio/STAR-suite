#include "SlamDump.h"
#include "SlamQuant.h"
#include "SlamCompat.h"
#include "SlamSolver.h"
#include "SlamQcOutput.h"
#include "SlamVarianceAnalysis.h"

#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <cctype>

struct Args {
    std::string dumpPath;
    std::string outPrefix;
    std::string maskPath;
    std::string autoTrimMode;
    std::string trimScope = "first";
    std::string strandness = "none";
    std::string weightMode = "dump";
    std::string weightFile;
    std::string weightMatch = "auto";
    std::string qcReportPrefix;
    std::string dumpOut;
    std::string dumpWeightsOut;
    std::string cbOutFile;
    std::string cbFormat = "star";
    std::array<int, 2> trim5pMate{{0, 0}};
    std::array<int, 2> trim3pMate{{0, 0}};
    double errorRate = -1.0;
    double convRate = -1.0;
    int cbOut = 0;
    uint32_t autoTrimMaxReads = 100000;
    uint32_t autoTrimMinReads = 1000;
    uint32_t autoTrimSmoothWindow = 5;
    uint32_t autoTrimSegMinLen = 3;
    uint32_t autoTrimMaxTrim = 15;
    uint32_t minCallableLength = 30;
};

static bool parseArgs(int argc, char** argv, Args* args) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                exit(2);
            }
            return argv[++i];
        };
        if (a == "--dump") args->dumpPath = next("--dump");
        else if (a == "--out") args->outPrefix = next("--out");
        else if (a == "--slamSnpMaskIn") args->maskPath = next("--slamSnpMaskIn");
        else if (a == "--trim5p") {
            int v = std::stoi(next("--trim5p"));
            args->trim5pMate[0] = v;
            args->trim5pMate[1] = v;
        }
        else if (a == "--trim3p") {
            int v = std::stoi(next("--trim3p"));
            args->trim3pMate[0] = v;
            args->trim3pMate[1] = v;
        }
        else if (a == "--autoTrim") args->autoTrimMode = next("--autoTrim");
        else if (a == "--trimScope") args->trimScope = next("--trimScope");
        else if (a == "--strandness") args->strandness = next("--strandness");
        else if (a == "--slamWeightMode") args->weightMode = next("--slamWeightMode");
        else if (a == "--slamWeightFile") args->weightFile = next("--slamWeightFile");
        else if (a == "--slamWeightMatch") args->weightMatch = next("--slamWeightMatch");
        else if (a == "--slamQcReport") args->qcReportPrefix = next("--slamQcReport");
        else if (a == "--dumpOut") args->dumpOut = next("--dumpOut");
        else if (a == "--dumpWeightsOut") args->dumpWeightsOut = next("--dumpWeightsOut");
        else if (a == "--slamCbOut") args->cbOut = std::stoi(next("--slamCbOut"));
        else if (a == "--slamCbOutFile") args->cbOutFile = next("--slamCbOutFile");
        else if (a == "--slamCbFormat") args->cbFormat = next("--slamCbFormat");
        else if (a == "--errorRate") args->errorRate = std::stod(next("--errorRate"));
        else if (a == "--convRate") args->convRate = std::stod(next("--convRate"));
        else if (a == "--autoTrimDetectionReads") args->autoTrimMaxReads = static_cast<uint32_t>(std::stoul(next("--autoTrimDetectionReads")));
        else if (a == "--autoTrimMinReads") args->autoTrimMinReads = static_cast<uint32_t>(std::stoul(next("--autoTrimMinReads")));
        else if (a == "--autoTrimSmoothWindow") args->autoTrimSmoothWindow = static_cast<uint32_t>(std::stoul(next("--autoTrimSmoothWindow")));
        else if (a == "--autoTrimSegMinLen") args->autoTrimSegMinLen = static_cast<uint32_t>(std::stoul(next("--autoTrimSegMinLen")));
        else if (a == "--autoTrimMaxTrim") args->autoTrimMaxTrim = static_cast<uint32_t>(std::stoul(next("--autoTrimMaxTrim")));
        else if (a == "--slamMinCallableLength" || a == "--minCallableLength") {
            args->minCallableLength = static_cast<uint32_t>(std::stoul(next(a.c_str())));
        }
        else {
            std::cerr << "Unknown arg: " << a << "\n";
            return false;
        }
    }
    if (args->dumpPath.empty() || args->outPrefix.empty()) {
        std::cerr << "Usage: --dump <path> --out <prefix> [--slamSnpMaskIn <bed.gz>] [--trim5p N --trim3p N] "
                     "[--slamWeightMode dump|alignments|uniform] [--slamWeightFile <path>] [--slamWeightMatch auto|order|key] "
                     "[--dumpOut <path>] [--dumpWeightsOut <path>] "
                     "[--slamCbOut 0|1] [--slamCbOutFile <path>] [--slamCbFormat star|ezbakr] "
                     "[--slamMinCallableLength N]\n";
        return false;
    }
    if (args->cbOut != 0 && args->cbOut != 1) {
        std::cerr << "--slamCbOut must be 0 or 1\n";
        return false;
    }
    if (args->cbFormat.empty() || args->cbFormat == "-") {
        args->cbFormat = "star";
    }
    std::transform(args->cbFormat.begin(), args->cbFormat.end(), args->cbFormat.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (args->cbFormat != "star" && args->cbFormat != "ezbakr") {
        std::cerr << "--slamCbFormat must be star or ezbakr\n";
        return false;
    }
    return true;
}

static int strandnessToInt(const std::string& s) {
    std::string x = s;
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x == "none" || x == "0") return 0;
    if (x == "sense" || x == "1") return 1;
    if (x == "antisense" || x == "2") return 2;
    return 0;
}

static int parseWeightMode(const std::string& s) {
    std::string x = s;
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x.empty() || x == "dump" || x == "stored") return 0;
    if (x == "alignments" || x == "nh" || x == "ntr" || x == "weight") return 1;
    if (x == "uniform" || x == "none" || x == "all") return 2;
    return -1;
}

static int parseWeightMatch(const std::string& s) {
    std::string x = s;
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
    if (x.empty() || x == "auto") return 0;
    if (x == "order" || x == "ordered") return 1;
    if (x == "key" || x == "keyed" || x == "map") return 2;
    return -1;
}

struct WeightKeyHash {
    size_t operator()(const SlamWeightKey& k) const {
        return static_cast<size_t>(k.h1 ^ (k.h2 + 0x9e3779b97f4a7c15ull + (k.h1 << 6) + (k.h1 >> 2)));
    }
};
struct WeightKeyEq {
    bool operator()(const SlamWeightKey& a, const SlamWeightKey& b) const {
        return a.h1 == b.h1 && a.h2 == b.h2;
    }
};

static void writeSlamOut(const std::string& outFile,
                         const std::vector<std::string>& geneIds,
                         const std::vector<std::string>& geneNames,
                         const SlamQuant& quant,
                         double errorRate,
                         double convRate) {
    std::ofstream out(outFile.c_str());
    if (!out.good()) return;
    out << "Gene\tSymbol\tReadCount\tConversions\tCoverage\tNTR\tMAP\tSigma\tLogLikelihood\n";
    SlamSolver solver(errorRate, convRate);
    const auto& genes = quant.genes();
    for (size_t i = 0; i < genes.size(); ++i) {
        const SlamGeneStats& stats = genes[i];
        if (stats.readCount <= 0.0) continue;
        SlamResult res = solver.solve(stats.histogram);
        const std::string& gid = (i < geneIds.size()) ? geneIds[i] : std::string("GENE_") + std::to_string(i);
        const std::string& gname = (i < geneNames.size() && !geneNames[i].empty()) ? geneNames[i] : gid;
        out << gid << "\t" << gname << "\t"
            << stats.readCount << "\t"
            << stats.conversions << "\t"
            << stats.coverage << "\t"
            << res.ntr << "\t"
            << res.ntr << "\t"
            << res.sigma << "\t"
            << res.log_likelihood << "\n";
    }
}

static SlamVarianceTrimResult computeTrimForReads(const std::vector<SlamBufferedRead>& reads,
                                                  uint32_t maxReads,
                                                  uint32_t minReads,
                                                  uint32_t smoothWindow,
                                                  uint32_t minSegLen,
                                                  uint32_t maxTrim) {
    uint32_t mateLen0_max = 0;
    uint32_t mateLen1_max = 0;
    for (const auto& r : reads) {
        mateLen0_max = std::max(mateLen0_max, r.readLength0);
        mateLen1_max = std::max(mateLen1_max, r.readLength1);
    }
    if (mateLen0_max == 0 && mateLen1_max == 0) {
        mateLen0_max = 100;
    }
    const uint32_t concatLen = mateLen0_max + mateLen1_max;

    SlamVarianceAnalyzer analyzer(maxReads, minReads, smoothWindow, minSegLen, maxTrim);
    const bool sep = mateLen1_max > 0;
    analyzer.setSeparateMateHistograms(sep);

    for (const auto& r : reads) {
        if (!analyzer.recordRead()) break;
        const uint32_t split = r.readLength0;
        for (const auto& p : r.positions) {
            bool isT = false;
            bool isTc = false;
            if (!r.isMinus) {
                isT = (p.refBase == 3);
                isTc = (p.refBase == 3 && p.readBase == 1);
            } else {
                isT = (p.refBase == 0);
                isTc = (p.refBase == 0 && p.readBase == 2);
            }
            if (sep) {
                const uint32_t mateIndex = p.secondMate ? 1u : 0u;
                uint32_t mateLocalPos = p.readPos;
                if (mateIndex != 0u) {
                    mateLocalPos = (p.readPos >= split) ? (p.readPos - split) : 0u;
                }
                analyzer.recordPositionMate(mateLocalPos, static_cast<uint8_t>(mateIndex),
                                            p.qual, isT, isTc);
            } else {
                analyzer.recordPosition(p.readPos, p.qual, isT, isTc);
            }
        }
    }
    return analyzer.computeTrimUnified(concatLen, mateLen0_max, mateLen1_max, sep);
}

int main(int argc, char** argv) {
    Args args;
    if (!parseArgs(argc, argv, &args)) {
        return 2;
    }

    SlamDumpMetadata meta;
    std::vector<SlamBufferedRead> reads;
    std::string err;
    if (!readSlamDump(args.dumpPath, &meta, &reads, &err)) {
        std::cerr << "Failed to read dump: " << err << "\n";
        return 1;
    }
    if (meta.geneIds.empty()) {
        std::cerr << "Dump contains no genes\n";
        return 1;
    }
    int weightMode = parseWeightMode(args.weightMode);
    if (weightMode < 0) {
        std::cerr << "Unknown --slamWeightMode: " << args.weightMode << "\n";
        return 2;
    }
    int weightMatch = parseWeightMatch(args.weightMatch);
    if (weightMatch < 0) {
        std::cerr << "Unknown --slamWeightMatch: " << args.weightMatch << "\n";
        return 2;
    }

    double errorRate = (args.errorRate >= 0.0) ? args.errorRate : meta.errorRate;
    double convRate = (args.convRate >= 0.0) ? args.convRate : meta.convRate;
    if (errorRate <= 0.0) errorRate = 0.001;

    SlamSnpMask mask;
    SlamSnpMask* maskPtr = nullptr;
    if (!args.maskPath.empty()) {
        if (!mask.loadBedWithChrMap(args.maskPath, meta.chrNames, meta.chrStart, &err)) {
            std::cerr << "Failed to load mask: " << err << "\n";
            return 1;
        }
        maskPtr = &mask;
    }

    if (!args.weightFile.empty()) {
        SlamWeightMetadata wmeta;
        std::vector<SlamWeightRecord> records;
        if (!readSlamWeights(args.weightFile, &wmeta, &records, &err)) {
            std::cerr << "Failed to read weights: " << err << "\n";
            return 1;
        }
        if (wmeta.nReads != reads.size()) {
            std::cerr << "WARNING: weight file read count (" << wmeta.nReads
                      << ") differs from dump (" << reads.size() << ")\n";
        }
        bool tryOrder = (weightMatch == 0 || weightMatch == 1);
        bool matchedOrder = false;
        if (tryOrder && records.size() >= reads.size()) {
            matchedOrder = true;
            for (size_t i = 0; i < reads.size(); ++i) {
                SlamWeightKey key = computeSlamWeightKey(reads[i]);
                if (records[i].key.h1 != key.h1 || records[i].key.h2 != key.h2) {
                    matchedOrder = false;
                    break;
                }
            }
        }
        if (matchedOrder) {
            for (size_t i = 0; i < reads.size(); ++i) {
                reads[i].weight = records[i].weight;
            }
            std::cerr << "Using weights (order match)\n";
        } else {
            if (weightMatch == 1) {
                std::cerr << "Weight order mismatch and --slamWeightMatch=order set\n";
                return 1;
            }
            std::unordered_map<SlamWeightKey, double, WeightKeyHash, WeightKeyEq> weightMap;
            weightMap.reserve(records.size());
            for (const auto& rec : records) {
                weightMap[rec.key] = rec.weight;
            }
            size_t matched = 0;
            for (auto& r : reads) {
                SlamWeightKey key = computeSlamWeightKey(r);
                auto it = weightMap.find(key);
                if (it != weightMap.end()) {
                    r.weight = it->second;
                    ++matched;
                }
            }
            std::cerr << "Using weights (key match): matched " << matched << "/" << reads.size() << "\n";
        }
    } else {
        if (weightMode == 0) {
            std::string mode = (meta.weightMode == 1) ? "Uniform" : "Alignments";
            std::cerr << "Using dump weights (" << mode << ")\n";
        } else if (weightMode == 1) {
            std::cerr << "Using recomputed weights: Alignments (1/nTr)\n";
        } else {
            std::cerr << "Using recomputed weights: Uniform (1.0)\n";
        }
    }

    // Compute trims if requested
    std::unordered_map<uint32_t, SlamVarianceTrimResult> perFileTrim;
    SlamVarianceTrimResult globalTrim;
    if (args.autoTrimMode == "variance") {
        if (args.trimScope == "per-file") {
            std::unordered_map<uint32_t, std::vector<SlamBufferedRead>> byFile;
            for (const auto& r : reads) {
                byFile[r.fileIndex].push_back(r);
            }
            for (const auto& kv : byFile) {
                perFileTrim[kv.first] = computeTrimForReads(kv.second,
                                                           args.autoTrimMaxReads,
                                                           args.autoTrimMinReads,
                                                           args.autoTrimSmoothWindow,
                                                           args.autoTrimSegMinLen,
                                                           args.autoTrimMaxTrim);
                if (!perFileTrim[kv.first].success) {
                    std::cerr << "ERROR: --autoTrim variance requested but insufficient stdev data to compute trims "
                              << "for file_index=" << kv.first
                              << " reads_analyzed=" << perFileTrim[kv.first].readsAnalyzed
                              << " min_reads=" << args.autoTrimMinReads << "\n";
                    return 1;
                }
            }
        } else {
            globalTrim = computeTrimForReads(reads,
                                             args.autoTrimMaxReads,
                                             args.autoTrimMinReads,
                                             args.autoTrimSmoothWindow,
                                             args.autoTrimSegMinLen,
                                             args.autoTrimMaxTrim);
            if (!globalTrim.success) {
                std::cerr << "ERROR: --autoTrim variance requested but insufficient stdev data to compute trims "
                          << "(reads_analyzed=" << globalTrim.readsAnalyzed
                          << " min_reads=" << args.autoTrimMinReads << ")\n";
                return 1;
            }
            args.trim5pMate[0] = globalTrim.mates[0].trim5p;
            args.trim5pMate[1] = globalTrim.mates[1].trim5p;
            args.trim3pMate[0] = globalTrim.mates[0].trim3p;
            args.trim3pMate[1] = globalTrim.mates[1].trim3p;
        }
    }

    auto runPartition = [&](const std::vector<SlamBufferedRead>& partReads,
                            const std::array<int, 2>& trim5pVals,
                            const std::array<int, 2>& trim3pVals) {
        SlamQuant q(meta.geneIds.size(), false);
        q.enableReadBuffer(static_cast<uint64_t>(partReads.size()));
        for (const auto& r : partReads) {
            SlamBufferedRead adjusted(r);
            if (args.weightFile.empty()) {
                if (weightMode == 1) {
                    size_t n = adjusted.geneIds.size();
                    adjusted.weight = (n > 0) ? (1.0 / static_cast<double>(n)) : 0.0;
                } else if (weightMode == 2) {
                    adjusted.weight = adjusted.geneIds.empty() ? 0.0 : 1.0;
                }
            }
            q.bufferRead(std::move(adjusted));
        }
        SlamCompatConfig cfg;
        cfg.trim5p = trim5pVals;
        cfg.trim3p = trim3pVals;
        SlamCompat compat(cfg, {}, {});
        q.replayBufferedReads(&compat, maskPtr, strandnessToInt(args.strandness),
                              args.minCallableLength);
        return q;
    };

    SlamQuant merged(meta.geneIds.size(), false);
    if (args.trimScope == "per-file") {
        std::unordered_map<uint32_t, std::vector<SlamBufferedRead>> byFile;
        for (const auto& r : reads) {
            byFile[r.fileIndex].push_back(r);
        }
        for (const auto& kv : byFile) {
            std::array<int, 2> t5 = args.trim5pMate;
            std::array<int, 2> t3 = args.trim3pMate;
            auto it = perFileTrim.find(kv.first);
            if (it != perFileTrim.end() && it->second.success) {
                t5[0] = it->second.mates[0].trim5p;
                t5[1] = it->second.mates[1].trim5p;
                t3[0] = it->second.mates[0].trim3p;
                t3[1] = it->second.mates[1].trim3p;
            }
            SlamQuant qPart = runPartition(kv.second, t5, t3);
            merged.merge(qPart);
        }
    } else {
        SlamQuant q = runPartition(reads, args.trim5pMate, args.trim3pMate);
        merged.merge(q);
    }

    std::string outBase = args.outPrefix;
    writeSlamOut(outBase + "SlamQuant.out", meta.geneIds, meta.geneNames, merged, errorRate, convRate);
    if (args.cbOut != 0) {
        std::string cbPath = (args.cbOutFile.empty() || args.cbOutFile == "-")
                             ? (outBase + "SlamQuant.cB.tsv")
                             : args.cbOutFile;
        if (!merged.writeCountBinomial(meta.geneIds, meta.geneNames, cbPath, outBase, args.cbFormat)) {
            std::cerr << "Failed to write count-binomial output: " << cbPath << "\n";
            return 1;
        }
    }
    merged.writeDiagnostics(outBase + "SlamQuant.out.diagnostics");
    merged.writeTransitions(outBase + "SlamQuant.out.transitions.tsv");
    merged.writeMismatches(outBase + "SlamQuant.out.mismatches.tsv", outBase);
    merged.writeMismatchDetails(outBase + "SlamQuant.out.mismatchdetails.tsv");

    if (!args.dumpOut.empty()) {
        SlamDumpMetadata outMeta = meta;
        outMeta.errorRate = errorRate;
        outMeta.convRate = convRate;
        if (!args.weightFile.empty()) {
            outMeta.weightMode = 0;  // alignments (weights provided separately if requested)
        } else if (weightMode == 2) {
            outMeta.weightMode = 1;  // uniform
        } else if (weightMode == 1) {
            outMeta.weightMode = 0;  // alignments
        }
        SlamReadBuffer outBuf(static_cast<uint64_t>(reads.size()));
        for (const auto& r : reads) {
            SlamBufferedRead tmp(r);
            outBuf.addRead(std::move(tmp));
        }
        std::vector<const SlamReadBuffer*> bufs{&outBuf};
        std::string dumpErr;
        if (!writeSlamDump(args.dumpOut, outMeta, bufs, 0, &dumpErr)) {
            std::cerr << "Failed to write dumpOut: " << dumpErr << "\n";
            return 1;
        }
        if (!args.dumpWeightsOut.empty()) {
            std::string wErr;
            if (!writeSlamWeights(args.dumpWeightsOut, outMeta, bufs, 0, nullptr, &wErr)) {
                std::cerr << "Failed to write dumpWeightsOut: " << wErr << "\n";
                return 1;
            }
        }
    }

    if (!args.qcReportPrefix.empty()) {
        std::string jsonPath = args.qcReportPrefix + ".slam_qc.json";
        std::string htmlPath = args.qcReportPrefix + ".slam_qc.html";
        SlamVarianceTrimResult* trimPtr = nullptr;
        if (args.autoTrimMode == "variance" && globalTrim.success) {
            trimPtr = &globalTrim;
        }
        uint32_t qcML0 = 0;
        uint32_t qcML1 = 0;
        for (const auto& r : reads) {
            qcML0 = std::max(qcML0, r.readLength0);
            qcML1 = std::max(qcML1, r.readLength1);
        }
        if (qcML0 == 0 && qcML1 == 0) {
            qcML0 = 100;
        }
        writeSlamQcComprehensiveJson(merged, jsonPath, args.trim5pMate[0], args.trim3pMate[0], args.trim5pMate[1],
                                    args.trim3pMate[1], trimPtr, nullptr, nullptr, qcML0, qcML1);
        writeSlamQcComprehensiveHtml(jsonPath, htmlPath);
    }

    std::cerr << "Requant complete: " << outBase << "SlamQuant.out\n";
    return 0;
}
