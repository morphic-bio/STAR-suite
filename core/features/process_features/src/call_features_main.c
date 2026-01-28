/**
 * @file call_features_main.c
 * @brief CLI tool for feature calling (dominant, GMM, and NB-EM modes)
 * 
 * Usage: call_features [options] <mex_dir> <output_dir>
 */

#include "../include/call_features.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

typedef enum {
    MODE_DOMINANT = 0,
    MODE_GMM = 1,
    MODE_NBEM = 2
} call_mode;

static void print_usage(const char *prog) {
    fprintf(stderr, "\nUsage: %s [options] <mex_dir> <output_dir>\n\n", prog);
    fprintf(stderr, "Assigns features (perturbations) to cells based on UMI counts.\n\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  mex_dir           Directory containing MEX files (matrix.mtx, barcodes.txt, features.txt)\n");
    fprintf(stderr, "  output_dir        Output directory for feature calls\n\n");
    fprintf(stderr, "Mode selection (mutually exclusive):\n");
    fprintf(stderr, "      --gmm                  Use CR9-compatible GMM calling\n");
    fprintf(stderr, "      --nb-em                Use NB-EM calling (SCEPTRE-style)\n");
    fprintf(stderr, "      (default: dominant mode)\n\n");
    fprintf(stderr, "Dominant mode options:\n");
    fprintf(stderr, "  -m, --min_counts <int>     Minimum UMI count for a feature to be considered (default: 2)\n");
    fprintf(stderr, "  -f, --fraction <float>     Required fraction of total for dominant call (default: 0.8)\n");
    fprintf(stderr, "  -g, --margin <int>         Required margin over second-best feature (default: 1)\n\n");
    fprintf(stderr, "GMM mode options:\n");
    fprintf(stderr, "  -u, --min_umi <int>        Minimum UMI threshold for GMM calls (default: 3)\n\n");
    fprintf(stderr, "NB-EM mode options:\n");
    fprintf(stderr, "      --moi <auto|low|high>  MOI mode (default: auto)\n");
    fprintf(stderr, "      --moi-min-umi <int>    Min UMI for guide presence in MOI calc (default: 1)\n");
    fprintf(stderr, "      --moi-pmulti <float>   Threshold for high MOI classification (default: 0.05)\n");
    fprintf(stderr, "      --prob-threshold <f>   Posterior threshold for positive call (default: 0.8)\n");
    fprintf(stderr, "      --backup-threshold <i> UMI fallback threshold when EM fails (default: 5)\n");
    fprintf(stderr, "      --max-iter <int>       Max EM iterations (default: 100)\n");
    fprintf(stderr, "      --tol <float>          Convergence tolerance (default: 1e-4)\n");
    fprintf(stderr, "      --global-phi <float>   Fixed dispersion (default: auto via MoM)\n");
    fprintf(stderr, "      --covariates <tsv>     Optional GEX covariate TSV\n\n");
    fprintf(stderr, "      --poisson-em           Use Poisson mixture instead of NB in EM\n");
    fprintf(stderr, "      --sceptre-parity       SCEPTRE-like: Poisson EM, per-guide calls, no fallback\n");
    fprintf(stderr, "      --nbem-debug           Write per-feature debug CSV\n");
    fprintf(stderr, "      --nbem-debug-max <int> Limit debug rows (default: all)\n");
    fprintf(stderr, "      --nbem-debug-csv <csv> Debug CSV path (default: output_dir/nbem_debug.csv)\n\n");
    fprintf(stderr, "      --nbem-debug-feature <name>  Log EM iterations for this feature\n");
    fprintf(stderr, "      --nbem-debug-iter-csv <csv>  Per-iteration debug CSV (default: none)\n\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "  -h, --help                 Show this help\n\n");
    fprintf(stderr, "Output files:\n");
    fprintf(stderr, "  Dominant mode:\n");
    fprintf(stderr, "    feature_calls.csv, feature_calls_summary.txt\n");
    fprintf(stderr, "  GMM mode:\n");
    fprintf(stderr, "    protospacer_calls_per_cell.csv, protospacer_calls_summary.csv,\n");
    fprintf(stderr, "    protospacer_umi_thresholds.csv, protospacer_umi_thresholds.json\n");
    fprintf(stderr, "  NB-EM mode:\n");
    fprintf(stderr, "    protospacer_calls_per_cell.csv, protospacer_calls_summary.csv,\n");
    fprintf(stderr, "    nbem_feature_params.csv, nbem_cell_posteriors.mtx, nbem_summary.json\n\n");
}

int main(int argc, char *argv[]) {
    /* Parse options */
    cf_config *config = cf_config_create();
    cf_gmm_config *gmm_config = cf_gmm_config_create();
    cf_nbem_config *nbem_config = cf_nbem_config_create();
    
    if (!config || !gmm_config || !nbem_config) {
        fprintf(stderr, "Failed to create config\n");
        return 1;
    }
    
    call_mode mode = MODE_DOMINANT;
    
    /* Long option codes for options without short equivalents */
    enum {
        OPT_GMM = 256,
        OPT_NBEM,
        OPT_MOI,
        OPT_MOI_MIN_UMI,
        OPT_MOI_PMULTI,
        OPT_PROB_THRESHOLD,
        OPT_BACKUP_THRESHOLD,
        OPT_MAX_ITER,
        OPT_TOL,
        OPT_GLOBAL_PHI,
        OPT_COVARIATES,
        OPT_POISSON_EM,
        OPT_SCEPTRE_PARITY,
        OPT_NBEM_DEBUG,
        OPT_NBEM_DEBUG_MAX,
        OPT_NBEM_DEBUG_CSV,
        OPT_NBEM_DEBUG_FEATURE,
        OPT_NBEM_DEBUG_ITER_CSV
    };
    
    static struct option long_options[] = {
        {"min_counts", required_argument, 0, 'm'},
        {"fraction", required_argument, 0, 'f'},
        {"margin", required_argument, 0, 'g'},
        {"gmm", no_argument, 0, OPT_GMM},
        {"nb-em", no_argument, 0, OPT_NBEM},
        {"min_umi", required_argument, 0, 'u'},
        {"moi", required_argument, 0, OPT_MOI},
        {"moi-min-umi", required_argument, 0, OPT_MOI_MIN_UMI},
        {"moi-pmulti", required_argument, 0, OPT_MOI_PMULTI},
        {"prob-threshold", required_argument, 0, OPT_PROB_THRESHOLD},
        {"backup-threshold", required_argument, 0, OPT_BACKUP_THRESHOLD},
        {"max-iter", required_argument, 0, OPT_MAX_ITER},
        {"tol", required_argument, 0, OPT_TOL},
        {"global-phi", required_argument, 0, OPT_GLOBAL_PHI},
        {"covariates", required_argument, 0, OPT_COVARIATES},
        {"poisson-em", no_argument, 0, OPT_POISSON_EM},
        {"sceptre-parity", no_argument, 0, OPT_SCEPTRE_PARITY},
        {"nbem-debug", no_argument, 0, OPT_NBEM_DEBUG},
        {"nbem-debug-max", required_argument, 0, OPT_NBEM_DEBUG_MAX},
        {"nbem-debug-csv", required_argument, 0, OPT_NBEM_DEBUG_CSV},
        {"nbem-debug-feature", required_argument, 0, OPT_NBEM_DEBUG_FEATURE},
        {"nbem-debug-iter-csv", required_argument, 0, OPT_NBEM_DEBUG_ITER_CSV},
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
            case OPT_GMM:
                if (mode == MODE_NBEM) {
                    fprintf(stderr, "Error: --gmm and --nb-em are mutually exclusive\n");
                    return 1;
                }
                mode = MODE_GMM;
                break;
            case OPT_NBEM:
                if (mode == MODE_GMM) {
                    fprintf(stderr, "Error: --gmm and --nb-em are mutually exclusive\n");
                    return 1;
                }
                mode = MODE_NBEM;
                break;
            case 'u':
                gmm_config->min_umi_threshold = atoi(optarg);
                break;
            case OPT_MOI:
                if (strcmp(optarg, "auto") == 0) {
                    nbem_config->moi_mode = CF_MOI_AUTO;
                } else if (strcmp(optarg, "low") == 0) {
                    nbem_config->moi_mode = CF_MOI_LOW;
                } else if (strcmp(optarg, "high") == 0) {
                    nbem_config->moi_mode = CF_MOI_HIGH;
                } else {
                    fprintf(stderr, "Error: --moi must be auto, low, or high\n");
                    return 1;
                }
                break;
            case OPT_MOI_MIN_UMI:
                nbem_config->moi_min_umi = atoi(optarg);
                break;
            case OPT_MOI_PMULTI:
                nbem_config->moi_pmulti_threshold = atof(optarg);
                break;
            case OPT_PROB_THRESHOLD:
                nbem_config->prob_threshold = atof(optarg);
                break;
            case OPT_BACKUP_THRESHOLD:
                nbem_config->backup_threshold = atoi(optarg);
                break;
            case OPT_MAX_ITER:
                nbem_config->max_iter = atoi(optarg);
                break;
            case OPT_TOL:
                nbem_config->tol = atof(optarg);
                break;
            case OPT_GLOBAL_PHI:
                nbem_config->global_phi = atof(optarg);
                break;
            case OPT_COVARIATES:
                nbem_config->covariate_tsv = optarg;
                break;
            case OPT_POISSON_EM:
                nbem_config->use_poisson = 1;
                break;
            case OPT_SCEPTRE_PARITY:
                nbem_config->use_poisson = 1;
                nbem_config->sceptre_parity = 1;
                break;
            case OPT_NBEM_DEBUG:
                nbem_config->debug = 1;
                break;
            case OPT_NBEM_DEBUG_MAX:
                nbem_config->debug_max_features = atoi(optarg);
                break;
            case OPT_NBEM_DEBUG_CSV:
                nbem_config->debug_csv = optarg;
                break;
            case OPT_NBEM_DEBUG_FEATURE:
                nbem_config->debug_feature = optarg;
                break;
            case OPT_NBEM_DEBUG_ITER_CSV:
                nbem_config->debug_iter_csv = optarg;
                break;
            case 'h':
                print_usage(argv[0]);
                cf_config_destroy(config);
                cf_gmm_config_destroy(gmm_config);
                cf_nbem_config_destroy(nbem_config);
                return 0;
            default:
                print_usage(argv[0]);
                cf_config_destroy(config);
                cf_gmm_config_destroy(gmm_config);
                cf_nbem_config_destroy(nbem_config);
                return 1;
        }
    }
    
    /* Check positional arguments */
    if (argc - optind < 2) {
        print_usage(argv[0]);
        cf_config_destroy(config);
        cf_gmm_config_destroy(gmm_config);
        cf_nbem_config_destroy(nbem_config);
        return 1;
    }
    
    const char *mex_dir = argv[optind];
    const char *output_dir = argv[optind + 1];
    
    int result;
    
    switch (mode) {
        case MODE_GMM:
            /* GMM mode (CR9-compatible) */
            printf("=== call_features (GMM mode) ===\n");
            printf("MEX directory:    %s\n", mex_dir);
            printf("Output directory: %s\n", output_dir);
            printf("Config:\n");
            printf("  min_umi:        %d\n\n", gmm_config->min_umi_threshold);
            result = cf_process_mex_dir_gmm(mex_dir, output_dir, gmm_config);
            break;
            
        case MODE_NBEM:
            /* NB-EM mode (SCEPTRE-style) */
            printf("=== call_features (NB-EM mode) ===\n");
            printf("MEX directory:    %s\n", mex_dir);
            printf("Output directory: %s\n", output_dir);
            printf("Config:\n");
            printf("  moi_mode:       %s\n", 
                   nbem_config->moi_mode == CF_MOI_AUTO ? "auto" :
                   nbem_config->moi_mode == CF_MOI_LOW ? "low" : "high");
            printf("  prob_threshold: %.2f\n", nbem_config->prob_threshold);
            printf("  backup_thresh:  %d\n", nbem_config->backup_threshold);
            printf("  max_iter:       %d\n", nbem_config->max_iter);
            printf("  tol:            %.1e\n", nbem_config->tol);
            if (nbem_config->global_phi > 0) {
                printf("  global_phi:     %.4f\n", nbem_config->global_phi);
            } else {
                printf("  global_phi:     auto (MoM)\n");
            }
            if (nbem_config->covariate_tsv) {
                printf("  covariates:     %s\n", nbem_config->covariate_tsv);
            }
            printf("\n");
            result = cf_process_mex_dir_nbem(mex_dir, output_dir, nbem_config);
            break;
            
        default:
            /* Dominant mode */
            printf("=== call_features (dominant mode) ===\n");
            printf("MEX directory:    %s\n", mex_dir);
            printf("Output directory: %s\n", output_dir);
            printf("Config:\n");
            printf("  min_counts:     %d\n", config->min_deduped_counts);
            printf("  fraction:       %.2f\n", config->dominance_fraction);
            printf("  margin:         %d\n\n", config->dominance_margin);
            result = cf_process_mex_dir(mex_dir, output_dir, config);
            break;
    }
    
    cf_config_destroy(config);
    cf_gmm_config_destroy(gmm_config);
    cf_nbem_config_destroy(nbem_config);
    return (result == 0) ? 0 : 1;
}
