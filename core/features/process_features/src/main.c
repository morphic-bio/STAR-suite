#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/utils.h"
#include "../include/io.h"
#include "../include/pf_api.h"
#include <stdio.h>

static void print_usage(const char *prog){
    fprintf(stderr, "\nUsage: %s [options] <FASTQ directories or files>\n\n", prog);
    fprintf(stderr, "Required flags:\n");
    fprintf(stderr, "  -w, --whitelist <file>            10x-style barcode whitelist (one per line)\n");
    fprintf(stderr, "  -f, --featurelist <file>          CSV with 'name' and 'sequence' columns\n");
    fprintf(stderr, "  -d, --directory  <path>           Output directory; one subdir per sample\n\n");

    fprintf(stderr, "Input & Output Files:\n");
    fprintf(stderr, "      --barcode_fastqs    <list>    Comma-separated R1 FASTQ files\n");
    fprintf(stderr, "      --forward_fastqs    <list>    Comma-separated R2 FASTQ files\n");
    fprintf(stderr, "      --reverse_fastqs    <list>    Comma-separated R3 FASTQ files\n");
    fprintf(stderr, "      --barcode_fastq_pattern <str> Pattern to identify barcode FASTQs (default _R1_)\n");
    fprintf(stderr, "      --forward_fastq_pattern <str> Pattern to identify forward FASTQs (default _R2_)\n");
    fprintf(stderr, "      --reverse_fastq_pattern <str> Pattern to identify reverse FASTQs (default _R3_)\n");
    fprintf(stderr, "      --filtered_barcodes <file>    File with barcodes to process, one per line\n");
    fprintf(stderr, "  -k, --keep_existing               Skip processing if output files exist\n");
    fprintf(stderr, "  -a, --as_named                    Treat all input files as single sample\n\n");

    fprintf(stderr, "Barcode & Feature Processing:\n");
    fprintf(stderr, "  -b, --barcode_length    <int>     Length of cell barcode (default 16)\n");
    fprintf(stderr, "  -u, --umi_length        <int>     Length of UMI (default 12)\n");
    fprintf(stderr, "  -o, --feature_constant_offset <int> Global feature offset (default: auto-detect from pattern column)\n");
    fprintf(stderr, "  -B, --barcode_constant_offset <int> Start position of barcode and UMI (default 0)\n");
    fprintf(stderr, "      --limit_search      <int>     Hard bound: search only within N bases of offset (-1 = entire read).\n");
    fprintf(stderr, "                                    No global fallback outside the window is performed.\n");
    fprintf(stderr, "      --feature_limited_mode <in_window_full|in_window_simple>\n");
    fprintf(stderr, "                                    In-window matcher when --limit_search >= 0 (default: in_window_full).\n");
    fprintf(stderr, "                                    Matching is strictly confined to the search window.\n");
    fprintf(stderr, "      --feature_limited_fallback    (deprecated alias for --feature_limited_mode; accepts full/simple)\n");
    fprintf(stderr, "      --feature_prehash_max_hamming <int> Build feature prehash up to this distance (0-2, default 2)\n");
    fprintf(stderr, "      --feature_prehash_max_entries <int> Entry budget for feature prehash build (default 50000000)\n");
    fprintf(stderr, "      --feature_prehash_memory_budget <bytes> Memory budget for prehash (0 = auto-detect, default 0)\n");
    fprintf(stderr, "      --force_individual_offsets    Use per-feature offsets from pattern column (slower for large sets)\n");
    fprintf(stderr, "      --use_feature_offset_array    (alias for --force_individual_offsets)\n");
    fprintf(stderr, "      --strict-offset-check         Error (instead of warn) on heterogeneous feature offsets\n");
    fprintf(stderr, "      --use_feature_anchor_search   Find pattern anchor via strstr before matching features\n");
    fprintf(stderr, "      --require_feature_anchor_match Require anchor match (no fallback search)\n");
    fprintf(stderr, "      --feature_mode_bootstrap_reads <int> Bootstrap N reads to learn per-feature offsets\n");
    fprintf(stderr, "  -r, --reverse_complement_whitelist Reverse complement whitelist barcodes\n\n");

    fprintf(stderr, "Error Correction & Thresholds:\n");
    fprintf(stderr, "  -m, --maxHammingDistance <int>    Max Hamming distance for feature match (default 1)\n");
    fprintf(stderr, "  -s, --stringency        <int>     UMI dedup stringency (see README, default 1)\n");
    fprintf(stderr, "  -i, --min_counts        <int>     Min reads in UMI clique for counting (default 0)\n");
    fprintf(stderr, "  -M, --min_posterior     <float>   Min posterior probability for barcode rescue (default 0.975)\n");
    fprintf(stderr, "      --max_barcode_mismatches <int> Max mismatches to rescue sequence barcode (default 3)\n");
    fprintf(stderr, "      --legacy-cb-rescue       Use legacy order-dependent pending barcode rescue\n");
    fprintf(stderr, "      --feature_n         <int>     Max 'N' bases allowed in feature sequence (default 1)\n");
    fprintf(stderr, "      --barcode_n         <int>     Max 'N' bases allowed in sequence barcode (default 1)\n");
    fprintf(stderr, "      --max_reads         <long>    Max reads to process per FASTQ (0 = all)\n");
    fprintf(stderr, "      --min_prediction    <int>     Min prediction threshold for feature assignment (default 1)\n");
    fprintf(stderr, "      --min_heatmap       <int>     Min deduped count for feature in heatmap (default 0)\n\n");
    fprintf(stderr, "      --skip_qc_outputs            Skip feature histograms and heatmaps\n\n");

    fprintf(stderr, "Performance & Parallelism:\n");
    fprintf(stderr, "  -t, --threads           <int>     Max concurrent sample processes (default 8)\n");
    fprintf(stderr, "  -S, --search_threads    <int>     Threads per consumer for feature search (default 4)\n");
    fprintf(stderr, "  -c, --consumer_threads_per_set <int> Consumer threads per sample (default 1)\n");
    fprintf(stderr, "  -R, --read_buffer_lines <int>     Lines for read buffer (default 1024)\n");
    fprintf(stderr, "  -L, --average_read_length <int>   Estimated avg read length for buffer allocation (default 300)\n\n");

    fprintf(stderr, "EmptyDrops Options:\n");
    fprintf(stderr, "      --skip_empty_drops            Skip EmptyDrops cell calling\n");
    fprintf(stderr, "      --emptydrops_expected_cells <int> Expected number of cells (default: auto)\n");
    fprintf(stderr, "      --emptydrops_failure_fatal    Treat EmptyDrops failure as error\n\n");
    fprintf(stderr, "      --emptydrops_use_fdr          Use FDR gate for tail rescue (default: raw p-value)\n\n");

    fprintf(stderr, "Namespace & Compatibility:\n");
    fprintf(stderr, "      --translate_NXT               Complement positions 8 and 9 of cell barcodes at output/filter stages\n");
    fprintf(stderr, "      --source_namespace <NXT|TRU>  Namespace of filtered barcode file (required with --filtered_barcodes)\n");
    fprintf(stderr, "      --target_namespace <NXT|TRU>  Namespace of assignment output (required with --filtered_barcodes)\n");
    fprintf(stderr, "      --allow_union_whitelist       Accept mixed NXT+TRU whitelist/filtered barcode files.\n");
    fprintf(stderr, "                                    Expands filtered barcodes at ingress so both namespace forms\n");
    fprintf(stderr, "                                    are present. Legacy compat for raw 3M-february-2018.txt.\n\n");

    fprintf(stderr, "Miscellaneous:\n");
    fprintf(stderr, "  -v, --debug                       Enable verbose debug output\n");
    fprintf(stderr, "  -h, --help                        Show this help and exit\n\n");
}

static int parse_feature_limited_mode(const char *value, int *mode_out) {
    if (!value || !mode_out) return 0;
    if (strcmp(value, "in_window_full") == 0 ||
        strcmp(value, "full") == 0 || strcmp(value, "0") == 0) {
        *mode_out = 0;
        return 1;
    }
    if (strcmp(value, "in_window_simple") == 0 ||
        strcmp(value, "simple") == 0 || strcmp(value, "1") == 0) {
        *mode_out = 1;
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[])
{   
    omp_set_nested(1);
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();
    feature_arrays *features=0;
    int reverse_complement_whitelist=0;
    char *barcodeFastqFilesString=0;
    char *forwardFastqFilesString=0;
    char *reverseFastqFilesString=0;
    char *filtered_barcodes_filename=0;
    char barcode_pattern[LINE_LENGTH], forward_pattern[LINE_LENGTH], reverse_pattern[LINE_LENGTH];
    strcpy(barcode_pattern, "_R1_");
    strcpy(forward_pattern, "_R2_");
    strcpy(reverse_pattern, "_R3_");
    
    int maxHammingDistance=1;
    char directory[LINE_LENGTH];
    char whitelist_filename[4096];
    char sample_flag=1;
    char keep_existing=0;
    uint16_t stringency=1;
    uint16_t min_counts=0;
    int read_buffer_lines=READ_BUFFER_LINES;
    int average_read_length=AVERAGE_READ_LENGTH;
    int feature_constant_offset=-1;  /* sentinel: -1 means auto-detect */
    int feature_constant_offset_explicit=0;  /* track if user provided --feature_constant_offset */
    int barcode_constant_offset=0;
    double min_posterior=MIN_POSTERIOR;
    int use_feature_offset_array_cli=0;
    int use_feature_anchor_search_cli=0;
    int require_feature_anchor_match_cli=0;
    int feature_mode_bootstrap_reads_cli=0;
    int feature_limited_fallback_mode_cli=0;
    int strict_offset_check_cli=0;
    int autodetect_chemistry_cli=0;
    int autodetect_chemistry_reads_cli=10000;
    int autodetect_chemistry_min_hits_cli=50;
    int allow_union_whitelist_cli=0;
    pf_namespace_t source_namespace_cli = PF_NS_UNKNOWN;
    pf_namespace_t target_namespace_cli = PF_NS_UNKNOWN;

    int max_concurrent_processes=8;
    int consumer_threads_per_set=1;
    int search_threads_per_consumer=4;
    int set_consumer_threads_per_set=0;
    int set_search_threads_per_consumer=0;
    uint16_t min_prediction = 1;

    char *sample_barcodes_filename = NULL;
    int sample_max_hamming = 1;
    int sample_max_N = 0;
    int sample_constant_offset_cli = -2;  /* sentinel */
    int sample_offset_relative_cli = 0;
    feature_arrays *sample_barcodes = NULL;

    /* EmptyDrops control */
    int skip_emptydrops = 0;
    int emptydrops_expected_cells = 0;
    int emptydrops_failure_fatal = 0;
    int emptydrops_use_fdr = 0;
    int legacy_cb_rescue = 0;
    int skip_qc_outputs = 0;

    static struct option long_options[] = {
        {"whitelist", required_argument, 0, 'w'},
        {"featurelist", required_argument, 0, 'f'},
        {"maxHammingDistance", required_argument, 0, 'm'},
        {"feature_constant_offset", required_argument, 0, 'o'},
        {"barcode_constant_offset", required_argument, 0, 'B'},
        {"threads", required_argument, 0, 't'},
        {"search_threads", required_argument, 0, 'S'},
        {"stringency", required_argument, 0, 's'},
        {"min_counts", required_argument, 0, 'i'},
        {"umi_length",required_argument , 0, 'u'},
        {"barcode_length",required_argument , 0, 'b'},
        {"directory", required_argument, 0, 'd'},
        {"debug", no_argument, 0, 'v'},
        {"as_named", no_argument, 0, 'a'},
        {"reverse_complement_whitelist", no_argument, 0, 'r'},
        {"keep_existing", no_argument, 0, 'k'},
        {"read_buffer_lines", required_argument, 0, 'R'},
        {"average_read_length", required_argument, 0, 'L'},
        {"min_posterior", required_argument, 0, 'M'},
        {"consumer_threads_per_set", required_argument, 0, 'c'},
        {"barcode_fastqs", required_argument, 0, 0},
        {"forward_fastqs", required_argument, 0, 1},
        {"reverse_fastqs", required_argument, 0, 2},
        {"max_barcode_mismatches", required_argument, 0, 3},
        {"legacy-cb-rescue", no_argument, 0, 27},
        {"legacy_cb_rescue", no_argument, 0, 27}, /* alias */
        {"feature_n", required_argument, 0, 4},
        {"barcode_n", required_argument, 0, 5},
        {"barcode_fastq_pattern", required_argument, 0, 6},
        {"forward_fastq_pattern", required_argument, 0, 7},
        {"reverse_fastq_pattern", required_argument, 0, 8},
        {"max_reads", required_argument, 0, 9},
        {"limit_search", required_argument, 0, 10},
        {"feature_limited_mode", required_argument, 0, 35},
        {"feature_limited_fallback", required_argument, 0, 11},
        {"feature_prehash_max_hamming", required_argument, 0, 33},
        {"feature_prehash_max_entries", required_argument, 0, 34},
        {"feature_prehash_memory_budget", required_argument, 0, 39},
        {"use_feature_offset_array", no_argument, 0, 22},
        {"force_individual_offsets", no_argument, 0, 22},  /* alias */
        {"strict-offset-check", no_argument, 0, 26},
        {"use_feature_anchor_search", no_argument, 0, 23},
        {"require_feature_anchor_match", no_argument, 0, 24},
        {"feature_mode_bootstrap_reads", required_argument, 0, 25},
        {"autodetect_chemistry", required_argument, 0, 30},
        {"autodetect_chemistry_reads", required_argument, 0, 31},
        {"autodetect_chemistry_min_hits", required_argument, 0, 32},
        {"filtered_barcodes", required_argument, 0, 12},
        {"min_prediction", required_argument, 0, 15},
        {"min_heatmap", required_argument, 0, 16},
        {"skip_qc_outputs", no_argument, 0, 40},
        {"translate_NXT", no_argument, 0, 17},
        {"skip_empty_drops", no_argument, 0, 18},
        {"skip_emptydrops", no_argument, 0, 18},  /* alias */
        {"emptydrops_expected_cells", required_argument, 0, 19},
        {"emptydrops_failure_fatal", no_argument, 0, 20},
        {"emptydrops_use_fdr", no_argument, 0, 21},
        {"emptydrops-use-fdr", no_argument, 0, 21},  /* alias */
        {"allow_union_whitelist", no_argument, 0, 36},
        {"allow-union-whitelist", no_argument, 0, 36},  /* alias */
        {"source_namespace", required_argument, 0, 37},
        {"source-namespace", required_argument, 0, 37},
        {"target_namespace", required_argument, 0, 38},
        {"target-namespace", required_argument, 0, 38},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    if (argc == 1) { print_usage(argv[0]); return 0; }

    while ((c = getopt_long(argc, argv, "w:b:f:m:s:S:i:t:T:o:d:u:c:vakrDB:R:L:M:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'w': strcpy(whitelist_filename, optarg); break;
            case 'b': barcode_length=atoi(optarg); barcode_code_length=(barcode_length+3)/4; break;
            case 'f':
                features=read_features_file(optarg);
                number_of_features=features->number_of_features;
                maximum_feature_length=features->max_length;
                feature_code_length=(maximum_feature_length+3)/4;
                feature_offsets = features->feature_offsets;
                feature_offsets_count = features->number_of_features;
                feature_anchors = features->feature_anchors;
                feature_anchor_lengths = features->feature_anchor_lengths;
                feature_suffix_anchors = features->feature_suffix_anchors;
                feature_suffix_anchor_lengths = features->feature_suffix_anchor_lengths;
                feature_anchor_count = features->number_of_features;
                if (feature_mode_bootstrap_reads > 0) {
                    feature_mode_offsets = malloc(sizeof(int) * features->number_of_features);
                    feature_mode_hist = calloc((size_t)features->number_of_features * (size_t)feature_mode_max_offset, sizeof(unsigned int));
                    if (!feature_mode_offsets || !feature_mode_hist) {
                        fprintf(stderr, "Failed to allocate feature mode arrays\n");
                        exit(EXIT_FAILURE);
                    }
                    for (int i = 0; i < features->number_of_features; i++) {
                        feature_mode_offsets[i] = -1;
                    }
                }
                break;
            case 'm': maxHammingDistance=atoi(optarg); break;
            case 's': stringency=(uint16_t)atoi(optarg); break;
            case 'S': set_search_threads_per_consumer=atoi(optarg); break;
            case 'i': min_counts=(uint16_t)atoi(optarg); break;
            case 'd': strcpy(directory, optarg); if (directory[strlen(directory)-1] != '/'){ strcat(directory, "/"); } break;
            case 'o': feature_constant_offset=atoi(optarg); feature_constant_offset_explicit=1; break;
            case 't': max_concurrent_processes=atoi(optarg); break;
            case 'u': umi_length=atoi(optarg); umi_code_length=(umi_length+3)/4; break;
            case 'c': set_consumer_threads_per_set=atoi(optarg); break;
            case 'v': debug = 1; break;
            case 'a': sample_flag=0; break;
            case 'k': keep_existing=1; break;
            case 'r': reverse_complement_whitelist=1; break;
            case 'B': barcode_constant_offset=atoi(optarg); break;
            case 'R': read_buffer_lines=atoi(optarg); break;
            case 'L': average_read_length=atoi(optarg); break;
            case 'M': min_posterior=(double)atof(optarg); break;
            case 0: barcodeFastqFilesString = strdup(optarg); break;
            case 1: forwardFastqFilesString = strdup(optarg); break;
            case 2: reverseFastqFilesString = strdup(optarg); break;
            case 3: max_barcode_mismatches=atoi(optarg); break;    
            case 27: legacy_cb_rescue = 1; break;
            case 4: max_feature_n=atoi(optarg); break;
            case 5: max_barcode_n=atoi(optarg); break;
            case 6: strcpy(barcode_pattern, optarg); break;
            case 7: strcpy(forward_pattern, optarg); break;
            case 8: strcpy(reverse_pattern, optarg); break;
            case 9: max_reads=atoll(optarg); break;
            case 10: limit_search = atoi(optarg); break;
            case 35:
                if (!parse_feature_limited_mode(optarg, &feature_limited_fallback_mode_cli)) {
                    fprintf(stderr, "Error: --feature_limited_mode must be one of: in_window_full, in_window_simple, full, simple, 0, 1\n");
                    return 1;
                }
                break;
            case 11:
                fprintf(stderr,
                    "WARNING: --feature_limited_fallback is a deprecated alias for --feature_limited_mode.\n"
                    "         Behavior is in-window only; --limit_search is a hard bound (no out-of-window rescue).\n"
                    "         Please migrate to --feature_limited_mode <in_window_full|in_window_simple>.\n");
                if (!parse_feature_limited_mode(optarg, &feature_limited_fallback_mode_cli)) {
                    fprintf(stderr, "Error: --feature_limited_fallback must be one of: full, simple, 0, 1\n");
                    return 1;
                }
                break;
            case 33:
                feature_prehash_max_hamming = atoi(optarg);
                if (feature_prehash_max_hamming < 0) feature_prehash_max_hamming = 0;
                if (feature_prehash_max_hamming > 2) feature_prehash_max_hamming = 2;
                break;
            case 34:
                feature_prehash_max_entries = strtoull(optarg, NULL, 10);
                break;
            case 39:
                feature_prehash_memory_budget = strtoull(optarg, NULL, 10);
                break;
            case 22: use_feature_offset_array_cli = 1; break;
            case 26: strict_offset_check_cli = 1; break;
            case 23: use_feature_anchor_search_cli = 1; break;
            case 24: require_feature_anchor_match_cli = 1; break;
            case 25: feature_mode_bootstrap_reads_cli = atoi(optarg); break;
            case 30: autodetect_chemistry_cli = atoi(optarg); break;
            case 31: autodetect_chemistry_reads_cli = atoi(optarg); break;
            case 32: autodetect_chemistry_min_hits_cli = atoi(optarg); break;
            case 12: filtered_barcodes_filename = strdup(optarg); break;
            case 15: min_prediction = atoi(optarg); break;
            case 16: min_heatmap = atoi(optarg); break;
            case 40: skip_qc_outputs = 1; break;
            case 17: translate_NXT = 1; fprintf(stderr, "translate_NXT enabled: complementing positions 8 and 9 at output/filter time.\n"); break;
            case 18: skip_emptydrops = 1; break;
            case 19: emptydrops_expected_cells = atoi(optarg); break;
            case 20: emptydrops_failure_fatal = 1; break;
            case 21: emptydrops_use_fdr = 1; break;
            case 36: allow_union_whitelist_cli = 1; break;
            case 37:
                source_namespace_cli = pf_namespace_from_string(optarg);
                if (source_namespace_cli != PF_NS_NXT && source_namespace_cli != PF_NS_TRU) {
                    fprintf(stderr, "ERROR: --source_namespace must be NXT or TRU (got '%s')\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 38:
                target_namespace_cli = pf_namespace_from_string(optarg);
                if (target_namespace_cli != PF_NS_NXT && target_namespace_cli != PF_NS_TRU) {
                    fprintf(stderr, "ERROR: --target_namespace must be NXT or TRU (got '%s')\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    if (use_feature_offset_array_cli) {
        use_feature_offset_array = 1;
    }
    if (use_feature_anchor_search_cli) {
        use_feature_anchor_search = 1;
    }
    if (require_feature_anchor_match_cli) {
        require_feature_anchor_match = 1;
        use_feature_anchor_search = 1;
    }
    if (feature_mode_bootstrap_reads_cli > 0) {
        feature_mode_bootstrap_reads = feature_mode_bootstrap_reads_cli;
        feature_mode_reads_seen = 0;
        feature_mode_bootstrap_done = 0;
    }
    if (feature_mode_bootstrap_reads > 0 && features && !feature_mode_offsets) {
        feature_mode_offsets = malloc(sizeof(int) * features->number_of_features);
        feature_mode_hist = calloc((size_t)features->number_of_features * (size_t)feature_mode_max_offset, sizeof(unsigned int));
        if (!feature_mode_offsets || !feature_mode_hist) {
            fprintf(stderr, "Failed to allocate feature mode arrays\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < features->number_of_features; i++) {
            feature_mode_offsets[i] = -1;
        }
    }
    feature_limited_fallback_mode = feature_limited_fallback_mode_cli;
    if (feature_prehash_memory_budget == 0) {
        feature_prehash_memory_budget = prehash_detect_memory_budget();
    }
    fprintf(stderr, "Feature prehash policy: max_hamming=%d max_entries=%llu memory_budget=%lluGB\n",
            feature_prehash_max_hamming,
            (unsigned long long)feature_prehash_max_entries,
            feature_prehash_memory_budget ? (unsigned long long)(feature_prehash_memory_budget / (1024ULL*1024ULL*1024ULL)) : 0ULL);
    
    /* Feature offset preflight detection */
    if (use_feature_offset_array_cli && feature_constant_offset_explicit) {
        fprintf(stderr, "Error: Cannot specify both --force_individual_offsets and --feature_constant_offset.\n");
        fprintf(stderr, "       Use --force_individual_offsets for per-feature offsets from pattern column,\n");
        fprintf(stderr, "       or --feature_constant_offset N for a single global offset.\n");
        return 1;
    }
    
    if (!use_feature_offset_array_cli && !feature_constant_offset_explicit && features && feature_offsets) {
        /* Auto-detect: scan feature_offsets for heterogeneity */
        int offset_counts[256] = {0};  /* count occurrences of each offset 0-255 */
        int max_offset_seen = -1;
        int valid_offsets = 0;
        
        for (int i = 0; i < feature_offsets_count; i++) {
            int off = feature_offsets[i];
            if (off >= 0 && off < 256) {
                offset_counts[off]++;
                valid_offsets++;
                if (off > max_offset_seen) max_offset_seen = off;
            }
        }
        
        if (valid_offsets > 0) {
            /* Find dominant offset */
            int dominant_offset = 0;
            int dominant_count = 0;
            int second_count = 0;
            
            for (int i = 0; i <= max_offset_seen; i++) {
                if (offset_counts[i] > dominant_count) {
                    second_count = dominant_count;
                    dominant_count = offset_counts[i];
                    dominant_offset = i;
                } else if (offset_counts[i] > second_count) {
                    second_count = offset_counts[i];
                }
            }
            
            /* Check for heterogeneity: second offset > 5% of dominant */
            double heterogeneity_threshold = 0.05;
            if (second_count > 0 && (double)second_count / (double)dominant_count > heterogeneity_threshold) {
                fprintf(stderr, "\n");
                fprintf(stderr, "WARNING: Multiple feature offsets detected in pattern column.\n");
                fprintf(stderr, "         Dominant offset: %d (used by %d features)\n", dominant_offset, dominant_count);
                fprintf(stderr, "         Other offsets detected:\n");
                for (int i = 0; i <= max_offset_seen; i++) {
                    if (i != dominant_offset && offset_counts[i] > 0) {
                        double pct = 100.0 * offset_counts[i] / dominant_count;
                        fprintf(stderr, "           offset %d: %d features (%.1f%%)\n", i, offset_counts[i], pct);
                    }
                }
                fprintf(stderr, "\n");
                if (strict_offset_check_cli) {
                    fprintf(stderr, "ERROR: --strict-offset-check is set. To proceed, choose one of:\n");
                    fprintf(stderr, "  1. --force_individual_offsets   Use per-feature offsets (slower for large feature sets)\n");
                    fprintf(stderr, "  2. --feature_constant_offset %d  Use dominant offset globally (faster)\n", dominant_offset);
                    fprintf(stderr, "  3. Remove --strict-offset-check to use dominant offset with warning\n");
                    fprintf(stderr, "\n");
                    return 1;
                }
                fprintf(stderr, "         Proceeding with dominant offset %d.\n", dominant_offset);
                fprintf(stderr, "         Use --force_individual_offsets for per-feature offsets, or\n");
                fprintf(stderr, "         --strict-offset-check to make this an error.\n\n");
            }
            
            /* Use dominant offset as global */
            feature_constant_offset = dominant_offset;
            fprintf(stderr, "[offset-detect] Auto-detected global offset: %d (from %d features with pattern)\n", 
                    dominant_offset, valid_offsets);
        } else {
            /* No valid offsets from pattern column - default to 0 */
            feature_constant_offset = 0;
            fprintf(stderr, "[offset-detect] No pattern offsets found, using default offset: 0\n");
        }
    } else if (feature_constant_offset == -1) {
        /* No features loaded yet or no offsets, default to 0 */
        feature_constant_offset = 0;
    }
    
    khash_t(strptr) *filtered_barcodes_hash = NULL;
    if (filtered_barcodes_filename) {
        int filtered_barcodes_found = 0;
        if (!pf_file_exists(filtered_barcodes_filename)) {
            printf("filtered_barcodes_filename: %s does not exist  \n", filtered_barcodes_filename);
            char full_path[2048];
            if (strlen(directory) > 0) {
                snprintf(full_path, sizeof(full_path), "%s%s", directory, filtered_barcodes_filename);
                printf("Will check for filtered barcodes file at %s\n", full_path);
                if (pf_file_exists(full_path)) {
                    free(filtered_barcodes_filename);
                    filtered_barcodes_filename = strdup(full_path);
                    filtered_barcodes_found = 1;
                } else {
                    fprintf(stderr, "Warning: Filtered barcodes file not found at %s or %s\n", filtered_barcodes_filename, full_path);
                }
            } else {
                fprintf(stderr, "Warning: Filtered barcodes file not found at %s\n", filtered_barcodes_filename);
            }
        }
        else{
            filtered_barcodes_found = 1;
            printf("Will use filtered barcodes file: %s\n", filtered_barcodes_filename);
        }
        
        if (filtered_barcodes_found) {
            filtered_barcodes_hash = kh_init(strptr);
            read_barcodes_into_hash(filtered_barcodes_filename, filtered_barcodes_hash);

            /* Ingress namespace normalization */
            if ((source_namespace_cli == PF_NS_NXT || source_namespace_cli == PF_NS_TRU) &&
                (target_namespace_cli == PF_NS_NXT || target_namespace_cli == PF_NS_TRU) &&
                source_namespace_cli != target_namespace_cli &&
                !allow_union_whitelist_cli) {
                int rc = pf_normalize_hash_namespace(filtered_barcodes_hash);
                if (rc < 0) {
                    fprintf(stderr,
                        "ERROR: Failed to normalize filtered barcodes from %s to %s "
                        "(memory allocation error). Aborting.\n",
                        pf_namespace_to_string(source_namespace_cli),
                        pf_namespace_to_string(target_namespace_cli));
                    return EXIT_FAILURE;
                }
                fprintf(stderr,
                    "NOTICE: Normalized %d filtered barcodes from %s to %s namespace.\n",
                    rc, pf_namespace_to_string(source_namespace_cli),
                    pf_namespace_to_string(target_namespace_cli));
            }

            if (allow_union_whitelist_cli) {
                int added = expand_hash_union_namespace(filtered_barcodes_hash);
                if (added < 0) {
                    fprintf(stderr,
                        "ERROR: --allow_union_whitelist expansion failed "
                        "(memory allocation error). Aborting.\n");
                    return EXIT_FAILURE;
                }
                if (added > 0) {
                    fprintf(stderr,
                        "WARNING: --allow_union_whitelist active. Expanded filtered barcodes "
                        "with %d NXT/TRU translations (total %u).\n"
                        "         This is legacy compat for union whitelists. "
                        "Migrate to namespace-split files for new workflows.\n",
                        added, kh_size(filtered_barcodes_hash));
                }
            }

            if (!allow_union_whitelist_cli) {
                int src_known = (source_namespace_cli == PF_NS_NXT || source_namespace_cli == PF_NS_TRU);
                int tgt_known = (target_namespace_cli == PF_NS_NXT || target_namespace_cli == PF_NS_TRU);
                if (!src_known || !tgt_known) {
                    fprintf(stderr,
                        "ERROR: Filtered barcodes loaded with incomplete namespace "
                        "metadata (source=%s, target=%s) and without "
                        "--allow_union_whitelist.  With exact-only matching, "
                        "a namespace mismatch will silently drop barcodes.\n"
                        "  Set both --source_namespace and --target_namespace, "
                        "or use --allow_union_whitelist.\n",
                        pf_namespace_to_string(source_namespace_cli),
                        pf_namespace_to_string(target_namespace_cli));
                    return EXIT_FAILURE;
                }
            }
        }
    }
    int positional_arg_count = argc - optind;
    if (positional_arg_count > 0 && (barcodeFastqFilesString != NULL || forwardFastqFilesString != NULL || reverseFastqFilesString != NULL)) {
        fprintf(stderr, "Error: Cannot specify both positional arguments and --barcode_fastqs\n");
        return 1;
    }
    if (strcmp(barcode_pattern, forward_pattern) == 0 || strcmp(barcode_pattern, reverse_pattern) == 0 || strcmp(forward_pattern, reverse_pattern) == 0){
        fprintf(stderr, "Error: Barcode, forward, and reverse patterns must be different\n");
        return 1;
    }
    


    fastq_files_collection fastq_files;
    memset(&fastq_files, 0, sizeof(fastq_files));
    if (positional_arg_count && pf_is_directory(argv[optind])){
        organize_fastq_files_by_directory(positional_arg_count, argc, argv, optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern);
    }
    else{
        organize_fastq_files_by_type(positional_arg_count, argc, argv,optind, barcodeFastqFilesString, forwardFastqFilesString, reverseFastqFilesString, &fastq_files, barcode_pattern, forward_pattern, reverse_pattern,sample_flag);
    }
    whitelist_hash = kh_init(u32ptr);
    read_whiteList(whitelist_filename, whitelist_hash, reverse_complement_whitelist);
    if (sample_barcodes_filename) {
        sample_barcodes = read_features_file(sample_barcodes_filename);
        if (!sample_barcodes) { fprintf(stderr, "Failed to load sample barcodes file %s\n", sample_barcodes_filename); exit(1);}    
        if ((sample_constant_offset_cli == -2) && sample_offset_relative_cli == 0) {
            fprintf(stderr, "Error: must specify --sample_offset or --sample_offset_rel when --sample_barcodes given\n");
            exit(1);
        }
        if (sample_constant_offset_cli != -2 && sample_offset_relative_cli != 0) {
            fprintf(stderr, "Error: specify only ONE of --sample_offset or --sample_offset_rel\n");
            exit(1);
        }
    }
    int demux_nsamples = (sample_barcodes) ? sample_barcodes->number_of_features : 1;
    initialize_unit_sizes();
    const int nSamples=fastq_files.nsamples;
    if (set_search_threads_per_consumer){
        search_threads_per_consumer=set_search_threads_per_consumer;
    }
    if (set_consumer_threads_per_set){
        consumer_threads_per_set=set_consumer_threads_per_set;
    } 
    fprintf(stderr, "Using %d consumer threads and %d search threads per consumer. Max concurrent processes %d\n", consumer_threads_per_set, search_threads_per_consumer, max_concurrent_processes);
    int concurrent_processes=0;
    atomic_int *thread_counter=mmap(NULL, sizeof(atomic_int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (thread_counter == MAP_FAILED){
        perror("Failed to allocate memory for thread counter");
        exit(EXIT_FAILURE);
    }
    atomic_init(thread_counter, 0);
    int threads_per_set = consumer_threads_per_set * (1 + search_threads_per_consumer);
    fprintf(stderr,"threads per set %d\n", threads_per_set);
    for (int index=0; index<nSamples; index++){
        const int i=fastq_files.sorted_index[index];
    
        if (concurrent_processes >= max_concurrent_processes){
            wait(NULL);
            concurrent_processes--;
        }
        pid_t pid=fork();
        if (pid< 0){
            perror("Failed to fork");
            exit(EXIT_FAILURE); 
        }
        else if (pid== 0){
            atomic_fetch_add(thread_counter, threads_per_set);
            
            char sample_directory[FILENAME_LENGTH];
            if (sample_flag){
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[i]);
                strcat(sample_directory, "/");
            }
            else{
                strcpy(sample_directory, directory);
                strcat(sample_directory, fastq_files.sample_names[0]);
                strcat(sample_directory, "/");
            }
            fprintf(stderr, "Processing sample directory %s\n", sample_directory);
            if (existing_output_skip(keep_existing, sample_directory)) exit(0);
            sample_args args;
            memset(&args, 0, sizeof(args));  // Zero-initialize all fields
            // The following are now initialized in process_files_in_sample
            // memory_pool_collection *pools=initialize_memory_pool_collection();
            // statistics stats;
            // data_structures hashes;
            // initialize_data_structures(&hashes);
            // initialize_statistics(&stats);
            args.sample_index = i;
            args.directory = sample_directory;
            if (filtered_barcodes_filename)
              args.filtered_barcodes_name = filtered_barcodes_filename;
            else
              args.filtered_barcodes_name = NULL;
            args.fastq_files = &fastq_files;
            args.features = features;
            args.maxHammingDistance = maxHammingDistance;
            args.nThreads = search_threads_per_consumer;
            args.pools = NULL;
            args.stats = NULL;
            args.hashes = NULL;
            args.stringency = stringency;
            args.min_counts = min_counts;
            args.barcode_constant_offset = barcode_constant_offset;
            args.feature_constant_offset = feature_constant_offset;
            args.read_buffer_lines = read_buffer_lines;
            args.average_read_length = average_read_length;
            args.min_posterior = min_posterior;
            args.legacy_cb_rescue = legacy_cb_rescue;
            args.consumer_threads_per_set = consumer_threads_per_set;
            args.filtered_barcodes_hash = filtered_barcodes_hash;
            args.min_prediction = min_prediction;
            args.min_heatmap = min_heatmap;
            args.demux_nsamples = demux_nsamples;
            args.sample_barcodes = sample_barcodes;
            args.sample_max_hamming = sample_max_hamming;
            args.sample_max_N = sample_max_N;
            args.sample_constant_offset = (sample_constant_offset_cli!=-2)? sample_constant_offset_cli : -1;
            args.sample_offset_relative = sample_offset_relative_cli;
            args.sample_hashes = NULL;
            args.sample_stats = NULL;
            args.sample_pools = NULL;
            
            /* EmptyDrops control from CLI flags */
            int child_error = 0;
            args.skip_emptydrops = skip_emptydrops;
            args.emptydrops_failure_fatal = emptydrops_failure_fatal;
            args.expected_cells = emptydrops_expected_cells;
            args.emptydrops_use_fdr = emptydrops_use_fdr;
            args.skip_qc_outputs = skip_qc_outputs;
            args.error_out = &child_error;
            struct chem_detect_state chem_detect_buf;
            if (autodetect_chemistry_cli) {
                if (autodetect_chemistry_reads_cli < 1) {
                    fprintf(stderr, "ERROR: --autodetect_chemistry_reads must be >= 1 (got %d)\n",
                            autodetect_chemistry_reads_cli);
                    exit(EXIT_FAILURE);
                }
                if (autodetect_chemistry_min_hits_cli < 1) {
                    fprintf(stderr, "ERROR: --autodetect_chemistry_min_hits must be >= 1 (got %d)\n",
                            autodetect_chemistry_min_hits_cli);
                    exit(EXIT_FAILURE);
                }
                memset(&chem_detect_buf, 0, sizeof(chem_detect_buf));
                chem_detect_buf.max_reads = autodetect_chemistry_reads_cli;
                chem_detect_buf.min_hits = autodetect_chemistry_min_hits_cli;
                args.chem_detect = &chem_detect_buf;
            } else {
                args.chem_detect = NULL;
            }
            
            process_files_in_sample(&args);
            // cleanup_sample is handled within process_files_in_sample
            atomic_fetch_add(thread_counter, -threads_per_set);
            exit(child_error ? 1 : 0);
        }
        concurrent_processes++;

    }
    
    /* Wait for all children and check exit status */
    int any_failed = 0;
    while (concurrent_processes > 0) {
        int status;
        pid_t child_pid = waitpid(-1, &status, 0);
        if (child_pid > 0) {
            if (WIFEXITED(status)) {
                int exit_code = WEXITSTATUS(status);
                if (exit_code != 0) {
                    fprintf(stderr, "Child process %d exited with error code %d\n", child_pid, exit_code);
                    any_failed = 1;
                }
            } else if (WIFSIGNALED(status)) {
                fprintf(stderr, "Child process %d killed by signal %d\n", child_pid, WTERMSIG(status));
                any_failed = 1;
            }
        }
        concurrent_processes--;
    }
    
    kh_destroy(u32ptr, whitelist_hash);
    free(whitelist);
    if (barcodeFastqFilesString) free(barcodeFastqFilesString);
    if (forwardFastqFilesString) free(forwardFastqFilesString);
    if (reverseFastqFilesString) free(reverseFastqFilesString);
    if (filtered_barcodes_filename) free(filtered_barcodes_filename);
    if (filtered_barcodes_hash) free_strptr_hash(filtered_barcodes_hash);
    free_feature_arrays(features);
    /* Hash is now owned by feature_arrays; global alias is stale after free.
     * Null out global to prevent accidental use-after-free. */
    feature_code_hash.h64 = NULL;
    free_fastq_files_collection(&fastq_files);
    
    if (any_failed) {
        fprintf(stderr, "One or more samples failed\n");
        return 1;
    }
    return 0;
}
