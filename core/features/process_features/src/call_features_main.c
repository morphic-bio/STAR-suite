/**
 * @file call_features_main.c
 * @brief CLI tool for feature calling (dominant mode and CR9-compatible GMM mode)
 * 
 * Usage: call_features [options] <mex_dir> <output_dir>
 */

#include "../include/call_features.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

static void print_usage(const char *prog) {
    fprintf(stderr, "\nUsage: %s [options] <mex_dir> <output_dir>\n\n", prog);
    fprintf(stderr, "Assigns features (perturbations) to cells based on UMI counts.\n\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  mex_dir           Directory containing MEX files (matrix.mtx, barcodes.txt, features.txt)\n");
    fprintf(stderr, "  output_dir        Output directory for feature calls\n\n");
    fprintf(stderr, "Mode selection:\n");
    fprintf(stderr, "      --gmm                  Use CR9-compatible GMM calling (default: dominant mode)\n\n");
    fprintf(stderr, "Dominant mode options:\n");
    fprintf(stderr, "  -m, --min_counts <int>     Minimum UMI count for a feature to be considered (default: 2)\n");
    fprintf(stderr, "  -f, --fraction <float>     Required fraction of total for dominant call (default: 0.8)\n");
    fprintf(stderr, "  -g, --margin <int>         Required margin over second-best feature (default: 1)\n\n");
    fprintf(stderr, "GMM mode options:\n");
    fprintf(stderr, "  -u, --min_umi <int>        Minimum UMI threshold for GMM calls (default: 3)\n\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "  -h, --help                 Show this help\n\n");
    fprintf(stderr, "Dominant mode output files:\n");
    fprintf(stderr, "  feature_calls.csv          Per-cell feature assignments\n");
    fprintf(stderr, "  feature_calls_summary.txt  Summary statistics\n\n");
    fprintf(stderr, "GMM mode output files (CR9-compatible):\n");
    fprintf(stderr, "  protospacer_calls_per_cell.csv    Per-cell feature calls\n");
    fprintf(stderr, "  protospacer_calls_summary.csv     Summary statistics\n");
    fprintf(stderr, "  protospacer_umi_thresholds.csv    UMI thresholds per feature\n");
    fprintf(stderr, "  protospacer_umi_thresholds.json   UMI thresholds (JSON format)\n\n");
}

int main(int argc, char *argv[]) {
    /* Parse options */
    cf_config *config = cf_config_create();
    cf_gmm_config *gmm_config = cf_gmm_config_create();
    
    if (!config || !gmm_config) {
        fprintf(stderr, "Failed to create config\n");
        return 1;
    }
    
    int use_gmm = 0;
    
    static struct option long_options[] = {
        {"min_counts", required_argument, 0, 'm'},
        {"fraction", required_argument, 0, 'f'},
        {"margin", required_argument, 0, 'g'},
        {"gmm", no_argument, 0, 'G'},
        {"min_umi", required_argument, 0, 'u'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, "m:f:g:u:Gh", long_options, NULL)) != -1) {
        switch (c) {
            case 'm':
                config->min_deduped_counts = atoi(optarg);
                break;
            case 'f':
                config->dominance_fraction = atof(optarg);
                break;
            case 'g':
                config->dominance_margin = atoi(optarg);
                break;
            case 'G':
                use_gmm = 1;
                break;
            case 'u':
                gmm_config->min_umi_threshold = atoi(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                cf_config_destroy(config);
                cf_gmm_config_destroy(gmm_config);
                return 0;
            default:
                print_usage(argv[0]);
                cf_config_destroy(config);
                cf_gmm_config_destroy(gmm_config);
                return 1;
        }
    }
    
    /* Check positional arguments */
    if (argc - optind < 2) {
        print_usage(argv[0]);
        cf_config_destroy(config);
        cf_gmm_config_destroy(gmm_config);
        return 1;
    }
    
    const char *mex_dir = argv[optind];
    const char *output_dir = argv[optind + 1];
    
    int result;
    
    if (use_gmm) {
        /* GMM mode (CR9-compatible) */
        printf("=== call_features (GMM mode) ===\n");
        printf("MEX directory:    %s\n", mex_dir);
        printf("Output directory: %s\n", output_dir);
        printf("Config:\n");
        printf("  min_umi:        %d\n\n", gmm_config->min_umi_threshold);
        
        result = cf_process_mex_dir_gmm(mex_dir, output_dir, gmm_config);
    } else {
        /* Dominant mode */
        printf("=== call_features (dominant mode) ===\n");
        printf("MEX directory:    %s\n", mex_dir);
        printf("Output directory: %s\n", output_dir);
        printf("Config:\n");
        printf("  min_counts:     %d\n", config->min_deduped_counts);
        printf("  fraction:       %.2f\n", config->dominance_fraction);
        printf("  margin:         %d\n\n", config->dominance_margin);
        
        result = cf_process_mex_dir(mex_dir, output_dir, config);
    }
    
    cf_config_destroy(config);
    cf_gmm_config_destroy(gmm_config);
    return (result == 0) ? 0 : 1;
}
