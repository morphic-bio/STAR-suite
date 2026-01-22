/**
 * @file star_feature_call.cpp
 * @brief STAR-suite pipeline entry point for feature barcode extraction and calling
 * 
 * This tool provides end-to-end perturb-seq processing:
 *   1. Feature extraction: FASTQ → MEX (via pf_api)
 *   2. Feature calling: MEX → calls (via call_features)
 * 
 * Supports CR9-compatible output layout with --compat-perturb mode.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>

// C APIs from process_features
extern "C" {
#include "../../features/process_features/include/pf_api.h"
#include "../../features/process_features/include/call_features.h"
#include "../../features/process_features/include/mex_10x.h"
}

// ============================================================================
// Configuration
// ============================================================================

struct Config {
    // Mode flags
    bool extract_mode = true;        // Run FASTQ → MEX extraction
    bool call_mode = true;           // Run MEX → calls
    bool compat_perturb = false;     // CR9-compatible output layout
    bool gmm_mode = true;            // Use GMM calling (vs ratio test)
    
    // Input paths
    std::string feature_ref;         // Feature reference CSV
    std::string whitelist;           // Barcode whitelist
    std::string filtered_barcodes;   // Filtered barcodes (cells)
    std::string fastq_dir;           // FASTQ input directory
    std::string mex_dir;             // MEX input (for call-only mode)
    
    // Output
    std::string output_dir;          // Output directory
    
    // Feature extraction params
    int barcode_length = 16;
    int umi_length = 12;
    int max_hamming = 2;
    int threads = 8;
    
    // GMM calling params
    int min_umi = 3;
    
    // Ratio test params (for non-GMM mode)
    int min_counts = 2;
    double dominance_fraction = 0.8;
    int dominance_margin = 1;
};

// ============================================================================
// Utility Functions
// ============================================================================

static int mkdir_p(const char *path) {
    char tmp[4096];
    char *p = nullptr;
    size_t len;
    
    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }
    
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            mkdir(tmp, 0755);
            *p = '/';
        }
    }
    return mkdir(tmp, 0755);
}

static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s [options]\n\n", prog);
    fprintf(stderr, "STAR-suite Feature Barcode Processing Tool\n\n");
    fprintf(stderr, "Modes:\n");
    fprintf(stderr, "  Full pipeline (FASTQ → MEX → Calls):\n");
    fprintf(stderr, "    %s --feature-ref REF --whitelist WL --fastq-dir DIR --output-dir OUT\n\n", prog);
    fprintf(stderr, "  Call-only (MEX → Calls):\n");
    fprintf(stderr, "    %s --call-only --mex-dir DIR --output-dir OUT\n\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Input:\n");
    fprintf(stderr, "    --feature-ref FILE      Feature reference CSV (required for extraction)\n");
    fprintf(stderr, "    --whitelist FILE        Barcode whitelist (required for extraction)\n");
    fprintf(stderr, "    --filtered-barcodes FILE  Filtered barcodes/cells (recommended)\n");
    fprintf(stderr, "    --fastq-dir DIR         FASTQ input directory\n");
    fprintf(stderr, "    --mex-dir DIR           MEX input directory (for --call-only)\n");
    fprintf(stderr, "    --output-dir DIR        Output directory (required)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Mode:\n");
    fprintf(stderr, "    --call-only             Skip extraction, only run calling on MEX\n");
    fprintf(stderr, "    --extract-only          Skip calling, only run FASTQ → MEX extraction\n");
    fprintf(stderr, "    --compat-perturb        CR9-compatible output layout (crispr_analysis/)\n");
    fprintf(stderr, "    --ratio-test            Use ratio test instead of GMM for calling\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Extraction parameters:\n");
    fprintf(stderr, "    --barcode-length N      Barcode length (default: 16)\n");
    fprintf(stderr, "    --umi-length N          UMI length (default: 12)\n");
    fprintf(stderr, "    --max-hamming N         Max Hamming distance for matching (default: 2)\n");
    fprintf(stderr, "    --threads N             Number of threads (default: 8)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  GMM calling parameters:\n");
    fprintf(stderr, "    --min-umi N             Minimum UMI threshold (default: 3)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Ratio test parameters:\n");
    fprintf(stderr, "    --min-counts N          Minimum counts threshold (default: 2)\n");
    fprintf(stderr, "    --dominance-fraction F  Required fraction for dominant (default: 0.8)\n");
    fprintf(stderr, "    --dominance-margin N    Required margin over second-best (default: 1)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Other:\n");
    fprintf(stderr, "    --help                  Show this help message\n");
    fprintf(stderr, "    --version               Show version\n");
}

// ============================================================================
// Feature Extraction (FASTQ → MEX)
// ============================================================================

static int run_extraction(const Config &cfg) {
    fprintf(stderr, "\n=== Feature Extraction (FASTQ → MEX) ===\n");
    
    // Create config
    pf_config *pf_cfg = pf_config_create();
    if (!pf_cfg) {
        fprintf(stderr, "Error: Failed to create pf_config\n");
        return -1;
    }
    
    pf_config_set_barcode_length(pf_cfg, cfg.barcode_length);
    pf_config_set_umi_length(pf_cfg, cfg.umi_length);
    pf_config_set_max_hamming_distance(pf_cfg, cfg.max_hamming);
    pf_config_set_threads(pf_cfg, cfg.threads);
    
    // Initialize context
    pf_context *ctx = pf_init(pf_cfg);
    pf_config_destroy(pf_cfg);
    
    if (!ctx) {
        fprintf(stderr, "Error: Failed to initialize pf_context\n");
        return -1;
    }
    
    // Load references
    fprintf(stderr, "Loading feature reference: %s\n", cfg.feature_ref.c_str());
    if (pf_load_feature_ref(ctx, cfg.feature_ref.c_str()) != PF_OK) {
        fprintf(stderr, "Error: Failed to load feature reference: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return -1;
    }
    fprintf(stderr, "  Loaded %d features\n", pf_get_num_features(ctx));
    
    fprintf(stderr, "Loading whitelist: %s\n", cfg.whitelist.c_str());
    if (pf_load_whitelist(ctx, cfg.whitelist.c_str()) != PF_OK) {
        fprintf(stderr, "Error: Failed to load whitelist: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return -1;
    }
    
    if (!cfg.filtered_barcodes.empty()) {
        fprintf(stderr, "Loading filtered barcodes: %s\n", cfg.filtered_barcodes.c_str());
        if (pf_load_filtered_barcodes(ctx, cfg.filtered_barcodes.c_str()) != PF_OK) {
            fprintf(stderr, "Warning: Failed to load filtered barcodes: %s\n", pf_get_error(ctx));
            // Continue without filtered barcodes
        }
    }
    
    // Process FASTQs
    fprintf(stderr, "Processing FASTQs from: %s\n", cfg.fastq_dir.c_str());
    fprintf(stderr, "Output directory: %s\n", cfg.output_dir.c_str());
    
    pf_stats stats = {};
    pf_error err = pf_process_fastq_dir(ctx, cfg.fastq_dir.c_str(), cfg.output_dir.c_str(), &stats);
    
    if (err != PF_OK) {
        fprintf(stderr, "Error: Failed to process FASTQs: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return -1;
    }
    
    // Print stats
    fprintf(stderr, "\nExtraction Statistics:\n");
    fprintf(stderr, "  Total reads:        %zu\n", stats.total_reads);
    fprintf(stderr, "  Matched reads:      %zu (%.1f%%)\n", 
            stats.matched_reads, 100.0 * stats.matched_reads / (stats.total_reads > 0 ? stats.total_reads : 1));
    fprintf(stderr, "  Total barcodes:     %zu\n", stats.total_barcodes);
    fprintf(stderr, "  Whitelisted:        %zu\n", stats.whitelisted_barcodes);
    fprintf(stderr, "  Deduped counts:     %zu\n", stats.total_deduped_counts);
    fprintf(stderr, "  Processing time:    %.2f sec\n", stats.processing_time_sec);
    
    pf_destroy(ctx);
    return 0;
}

// ============================================================================
// Feature Calling (MEX → Calls)
// ============================================================================

#include <dirent.h>

// Check if a directory contains MEX files (matrix.mtx or matrix.mtx.gz)
static bool has_mex_files(const std::string &dir) {
    std::string mtx_path = dir + "/matrix.mtx";
    std::string mtx_gz_path = dir + "/matrix.mtx.gz";
    struct stat st;
    return (stat(mtx_path.c_str(), &st) == 0) || (stat(mtx_gz_path.c_str(), &st) == 0);
}

// Find sample subdirectories that contain MEX files
// When prefer_filtered=true (for CR9 compat), prefer <sample>/filtered/ over <sample>/
static std::vector<std::string> find_sample_dirs(const std::string &base_dir, bool prefer_filtered) {
    std::vector<std::string> sample_dirs;
    
    // First check if base_dir itself has MEX files (direct MEX input)
    if (has_mex_files(base_dir)) {
        sample_dirs.push_back(base_dir);
        return sample_dirs;
    }
    
    // Otherwise, scan for subdirectories with MEX files
    DIR *dir = opendir(base_dir.c_str());
    if (!dir) {
        return sample_dirs;
    }
    
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;  // Skip hidden and . / ..
        if (strcmp(entry->d_name, "filtered") == 0) continue;  // Skip 'filtered' at top level
        
        std::string subdir = base_dir + "/" + entry->d_name;
        struct stat st;
        if (stat(subdir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            // For CR9 compat mode, prefer filtered/ subdirectory if it exists
            std::string filtered_subdir = subdir + "/filtered";
            if (prefer_filtered && has_mex_files(filtered_subdir)) {
                sample_dirs.push_back(filtered_subdir);
            } else if (has_mex_files(subdir)) {
                sample_dirs.push_back(subdir);
            }
        }
    }
    closedir(dir);
    
    return sample_dirs;
}

// Run calling on a single MEX directory
static int run_calling_on_mex(const std::string &mex_dir, const std::string &call_output_dir,
                               const Config &cfg) {
    fprintf(stderr, "  MEX directory: %s\n", mex_dir.c_str());
    fprintf(stderr, "  Output directory: %s\n", call_output_dir.c_str());
    
    mkdir_p(call_output_dir.c_str());
    
    int ret;
    if (cfg.gmm_mode) {
        cf_gmm_config *gmm_cfg = cf_gmm_config_create();
        if (!gmm_cfg) {
            fprintf(stderr, "Error: Failed to create GMM config\n");
            return -1;
        }
        gmm_cfg->min_umi_threshold = cfg.min_umi;
        
        ret = cf_process_mex_dir_gmm(mex_dir.c_str(), call_output_dir.c_str(), gmm_cfg);
        cf_gmm_config_destroy(gmm_cfg);
    } else {
        cf_config *ratio_cfg = cf_config_create();
        if (!ratio_cfg) {
            fprintf(stderr, "Error: Failed to create ratio config\n");
            return -1;
        }
        ratio_cfg->min_deduped_counts = cfg.min_counts;
        ratio_cfg->dominance_fraction = cfg.dominance_fraction;
        ratio_cfg->dominance_margin = cfg.dominance_margin;
        
        ret = cf_process_mex_dir(mex_dir.c_str(), call_output_dir.c_str(), ratio_cfg);
        cf_config_destroy(ratio_cfg);
    }
    
    return ret;
}

static int run_calling(const Config &cfg) {
    fprintf(stderr, "\n=== Feature Calling (MEX → Calls) ===\n");
    
    // Determine base MEX directory
    std::string base_mex_dir = cfg.mex_dir;
    if (base_mex_dir.empty()) {
        // Use output from extraction - MEX files are in output_dir/<sample>/
        base_mex_dir = cfg.output_dir;
    }
    
    fprintf(stderr, "Calling mode: %s\n", cfg.gmm_mode ? "GMM (CR9-compatible)" : "Ratio test");
    if (cfg.gmm_mode) {
        fprintf(stderr, "  min_umi: %d\n", cfg.min_umi);
    } else {
        fprintf(stderr, "  min_counts: %d\n", cfg.min_counts);
        fprintf(stderr, "  dominance_fraction: %.2f\n", cfg.dominance_fraction);
        fprintf(stderr, "  dominance_margin: %d\n", cfg.dominance_margin);
    }
    
    // Find sample directories containing MEX files
    // In compat-perturb mode, prefer filtered/ subdirectories (extraction writes filtered MEX there)
    std::vector<std::string> sample_dirs = find_sample_dirs(base_mex_dir, cfg.compat_perturb);
    
    if (sample_dirs.empty()) {
        fprintf(stderr, "Error: No MEX files found in %s or its subdirectories\n", base_mex_dir.c_str());
        fprintf(stderr, "  Expected matrix.mtx (or matrix.mtx.gz) in the directory or sample subdirectories.\n");
        return -1;
    }
    
    fprintf(stderr, "Found %zu sample(s) with MEX files\n", sample_dirs.size());
    
    int total_failures = 0;
    for (const auto &mex_dir : sample_dirs) {
        // Derive sample directory (strip trailing /filtered for compat runs)
        std::string sample_dir = mex_dir;
        if (cfg.compat_perturb) {
            const std::string filtered_suffix = "/filtered";
            if (sample_dir.size() > filtered_suffix.size() &&
                sample_dir.compare(sample_dir.size() - filtered_suffix.size(),
                                   filtered_suffix.size(), filtered_suffix) == 0) {
                sample_dir.erase(sample_dir.size() - filtered_suffix.size());
            }
        }

        // Extract sample name from sample directory path
        std::string sample_name;
        size_t last_slash = sample_dir.rfind('/');
        if (last_slash != std::string::npos && last_slash < sample_dir.length() - 1) {
            sample_name = sample_dir.substr(last_slash + 1);
        } else {
            sample_name = "sample";
        }
        
        fprintf(stderr, "\nProcessing sample: %s\n", sample_name.c_str());
        
        // Determine output directory for this sample's calls
        std::string call_output_dir;
        if (cfg.compat_perturb) {
            // CR9 layout: output_dir/<sample>/crispr_analysis/
            // If mex_dir == base_mex_dir (direct MEX input), use output_dir/crispr_analysis/
            if (sample_dirs.size() == 1 && mex_dir == base_mex_dir) {
                call_output_dir = cfg.output_dir + "/crispr_analysis";
            } else {
                // For multi-sample or extracted output, put crispr_analysis under the sample dir
                call_output_dir = sample_dir + "/crispr_analysis";
            }
        } else {
            // Non-compat mode: output_dir/<sample>/
            if (sample_dirs.size() == 1 && mex_dir == base_mex_dir) {
                call_output_dir = cfg.output_dir;
            } else {
                call_output_dir = mex_dir;
            }
        }

        std::string mex_dir_for_call = mex_dir;
        if (cfg.compat_perturb) {
            std::string tenx_dir = sample_dir + "/sample_filtered_feature_bc_matrix";
            if (pf_write_mex_10x(mex_dir.c_str(),
                                 tenx_dir.c_str(),
                                 "CRISPR Guide Capture",
                                 1) != 0) {
                fprintf(stderr, "Warning: Failed to write 10x MEX in %s; using raw MEX\n", tenx_dir.c_str());
            } else {
                mex_dir_for_call = tenx_dir;
            }
        }

        int ret = run_calling_on_mex(mex_dir_for_call, call_output_dir, cfg);
        if (ret != 0) {
            fprintf(stderr, "Warning: Feature calling failed for sample %s\n", sample_name.c_str());
            total_failures++;
        }
    }
    
    if (total_failures > 0) {
        fprintf(stderr, "\nFeature calling completed with %d failure(s)\n", total_failures);
        return -1;
    }
    
    fprintf(stderr, "\nFeature calling complete for all %zu sample(s).\n", sample_dirs.size());
    return 0;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char *argv[]) {
    Config cfg;
    
    static struct option long_options[] = {
        {"feature-ref",         required_argument, 0, 'f'},
        {"whitelist",           required_argument, 0, 'w'},
        {"filtered-barcodes",   required_argument, 0, 'b'},
        {"fastq-dir",           required_argument, 0, 'q'},
        {"mex-dir",             required_argument, 0, 'm'},
        {"output-dir",          required_argument, 0, 'o'},
        {"call-only",           no_argument,       0, 'C'},
        {"extract-only",        no_argument,       0, 'E'},
        {"compat-perturb",      no_argument,       0, 'P'},
        {"ratio-test",          no_argument,       0, 'R'},
        {"barcode-length",      required_argument, 0, 1001},
        {"umi-length",          required_argument, 0, 1002},
        {"max-hamming",         required_argument, 0, 1003},
        {"threads",             required_argument, 0, 't'},
        {"min-umi",             required_argument, 0, 1004},
        {"min-counts",          required_argument, 0, 1005},
        {"dominance-fraction",  required_argument, 0, 1006},
        {"dominance-margin",    required_argument, 0, 1007},
        {"help",                no_argument,       0, 'h'},
        {"version",             no_argument,       0, 'v'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "f:w:b:q:m:o:t:CEPRhv", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'f': cfg.feature_ref = optarg; break;
            case 'w': cfg.whitelist = optarg; break;
            case 'b': cfg.filtered_barcodes = optarg; break;
            case 'q': cfg.fastq_dir = optarg; break;
            case 'm': cfg.mex_dir = optarg; break;
            case 'o': cfg.output_dir = optarg; break;
            case 'C': cfg.extract_mode = false; cfg.call_mode = true; break;
            case 'E': cfg.extract_mode = true; cfg.call_mode = false; break;
            case 'P': cfg.compat_perturb = true; break;
            case 'R': cfg.gmm_mode = false; break;
            case 't': cfg.threads = atoi(optarg); break;
            case 1001: cfg.barcode_length = atoi(optarg); break;
            case 1002: cfg.umi_length = atoi(optarg); break;
            case 1003: cfg.max_hamming = atoi(optarg); break;
            case 1004: cfg.min_umi = atoi(optarg); break;
            case 1005: cfg.min_counts = atoi(optarg); break;
            case 1006: cfg.dominance_fraction = atof(optarg); break;
            case 1007: cfg.dominance_margin = atoi(optarg); break;
            case 'v':
                printf("star_feature_call version 1.0.0\n");
                printf("Part of STAR-suite\n");
                return 0;
            case 'h':
            default:
                print_usage(argv[0]);
                return (opt == 'h') ? 0 : 1;
        }
    }
    
    // Validate required arguments
    if (cfg.output_dir.empty()) {
        fprintf(stderr, "Error: --output-dir is required\n");
        print_usage(argv[0]);
        return 1;
    }
    
    if (cfg.extract_mode) {
        if (cfg.feature_ref.empty()) {
            fprintf(stderr, "Error: --feature-ref is required for extraction\n");
            print_usage(argv[0]);
            return 1;
        }
        if (cfg.whitelist.empty()) {
            fprintf(stderr, "Error: --whitelist is required for extraction\n");
            print_usage(argv[0]);
            return 1;
        }
        if (cfg.fastq_dir.empty()) {
            fprintf(stderr, "Error: --fastq-dir is required for extraction\n");
            print_usage(argv[0]);
            return 1;
        }
    }
    
    if (cfg.call_mode && !cfg.extract_mode && cfg.mex_dir.empty()) {
        fprintf(stderr, "Error: --mex-dir is required for --call-only mode\n");
        print_usage(argv[0]);
        return 1;
    }
    
    // Validate compat-perturb constraints
    if (cfg.compat_perturb) {
        // CR9 compat prefers GMM mode - ratio test outputs non-CR9 files
        if (!cfg.gmm_mode) {
            fprintf(stderr, "Warning: --compat-perturb with --ratio-test\n");
            fprintf(stderr, "  CR9-compatible output normally uses GMM calling; results may diverge.\n");
        }

        // Warn if filtered barcodes not provided (parity not guaranteed)
        if (cfg.extract_mode && cfg.filtered_barcodes.empty()) {
            fprintf(stderr, "Warning: --compat-perturb without --filtered-barcodes\n");
            fprintf(stderr, "  CR9 parity requires filtered barcodes to match CR9's cell selection.\n");
            fprintf(stderr, "  Output will use all whitelisted barcodes, not just filtered cells.\n");
        }
    }
    
    // Create output directory
    mkdir_p(cfg.output_dir.c_str());
    
    // Print configuration
    fprintf(stderr, "=== star_feature_call ===\n");
    fprintf(stderr, "Mode: %s%s%s\n", 
            cfg.extract_mode ? "extract" : "",
            (cfg.extract_mode && cfg.call_mode) ? " + " : "",
            cfg.call_mode ? "call" : "");
    if (cfg.compat_perturb) {
        fprintf(stderr, "Output layout: CR9-compatible (crispr_analysis/)\n");
    }
    if (cfg.call_mode) {
        fprintf(stderr, "Calling mode: %s\n", cfg.gmm_mode ? "GMM (CR9-compatible)" : "Ratio test");
    }
    
    // Run extraction if requested
    if (cfg.extract_mode) {
        int ret = run_extraction(cfg);
        if (ret != 0) {
            return ret;
        }
    }
    
    // Run calling if requested
    if (cfg.call_mode) {
        int ret = run_calling(cfg);
        if (ret != 0) {
            return ret;
        }
    }
    
    fprintf(stderr, "\n=== Done ===\n");
    return 0;
}
