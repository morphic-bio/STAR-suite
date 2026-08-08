// Merge per-shard trim-QC dumps into one FastQC-like report.
//
// The sharded bulk pipeline maps on N nodes and never assembles the reads
// anywhere, so there is no single place to run QC over them. But the collector
// is a purely additive set of histograms, so QC does not need the reads to be
// co-located -- only the counters. Each shard writes its accumulator with
// --trimQcShardOut and this merges them into the same trim_qc.json/html the
// in-process path emits.
//
// This is the same argument as the VB/EM gather: move the summary, not the data.

#include "TrimQc.h"
#include "TrimQcShard.h"
#include "TrimQcOutput.h"
#include "Stats.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

void usage() {
    std::cerr <<
      "usage: trim_qc_merge --out-prefix PREFIX SHARD.trimqc...\n"
      "       trim_qc_merge --json PATH [--html PATH] SHARD.trimqc...\n"
      "\n"
      "  --out-prefix   writes PREFIX.trim_qc.json and PREFIX.trim_qc.html\n"
      "  --json/--html  explicit output paths\n"
      "  --stage LABEL  stage label recorded in the JSON (default: trimmed)\n"
      "\n"
      "Inputs are the dumps produced by STAR --trimQcShardOut.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string out_prefix, json_path, html_path, stage = "trimmed";
    std::vector<std::string> inputs;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << "EXITING: " << a << " needs a value\n"; std::exit(2); }
            return argv[++i];
        };
        if      (a == "--out-prefix") out_prefix = next();
        else if (a == "--json")       json_path  = next();
        else if (a == "--html")       html_path  = next();
        else if (a == "--stage")      stage      = next();
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (!a.empty() && a[0] == '-') { std::cerr << "EXITING: unknown flag " << a << "\n"; return 2; }
        else inputs.push_back(a);
    }

    if (inputs.empty()) { usage(); return 2; }
    if (!out_prefix.empty()) {
        json_path = out_prefix + ".trim_qc.json";
        html_path = out_prefix + ".trim_qc.html";
    }
    if (json_path.empty()) {
        std::cerr << "EXITING: need --out-prefix or --json\n";
        return 2;
    }

    TrimQcCollector merged;
    TrimQcShardCounters counters;
    std::size_t nmerged = 0;
    for (const std::string& path : inputs) {
        TrimQcCollector qc;
        TrimQcShardCounters c;
        std::string err;
        if (!readTrimQcShard(path, qc, c, &err)) {
            std::cerr << "EXITING: " << err << "\n";
            return 1;
        }
        // merge() adopts the first enabled collector wholesale, so an
        // empty accumulator is a valid starting point.
        merged.merge(qc);
        counters.merge(c);
        ++nmerged;
    }
    if (!merged.enabled()) {
        std::cerr << "EXITING: no shard carried trim-QC data\n";
        return 1;
    }

    // writeTrimQcJson reads only the eight trim counters off Stats, and those
    // are plain sums across shards.
    Stats stats;
    stats.trimReadsProcessed      = counters.reads_processed;
    stats.trimReadsTrimmed        = counters.reads_trimmed;
    stats.trimReadsTooShort       = counters.reads_too_short;
    stats.trimBasesQualityTrimmed = counters.bases_quality_trimmed;
    stats.trimBasesAdapterTrimmed = counters.bases_adapter_trimmed;
    stats.trimPairsProcessed      = counters.pairs_processed;
    stats.trimPairsDropped        = counters.pairs_dropped;
    stats.trimPairsKept           = counters.pairs_kept;

    if (!writeTrimQcJson(merged, stats, json_path, stage)) {
        std::cerr << "EXITING: cannot write " << json_path << "\n";
        return 1;
    }
    if (!html_path.empty() && html_path != "-") {
        if (!writeTrimQcHtml(json_path, html_path)) {
            std::cerr << "EXITING: cannot write " << html_path << "\n";
            return 1;
        }
    }
    std::printf("merged %zu shard dumps: %llu reads, %zu mates -> %s\n",
                nmerged, (unsigned long long) merged.totalReads(),
                merged.mateCount(), json_path.c_str());
    return 0;
}
