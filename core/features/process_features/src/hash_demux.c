#include "../include/hash_demux.h"
#include "../include/mex_10x.h"
#include "../include/utils.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

static int normalize_token(const char *in, char *out, size_t out_len) {
    size_t j = 0;
    if (!in || !out || out_len == 0) return 0;
    for (size_t i = 0; in[i] != '\0' && j + 1 < out_len; ++i) {
        unsigned char c = (unsigned char)in[i];
        if (isalnum(c)) {
            out[j++] = (char)tolower(c);
        }
    }
    out[j] = '\0';
    return (int)j;
}

int pf_is_hash_like_feature_type(const char *feature_type) {
    char norm[128];
    if (!feature_type) return 0;
    normalize_token(feature_type, norm, sizeof(norm));
    return strcmp(norm, "multiplexingcapture") == 0
        || strcmp(norm, "hto") == 0
        || strcmp(norm, "cmo") == 0
        || strcmp(norm, "cellplexcmo") == 0;
}

static int row_matches_selector(feature_arrays *features, int idx,
                                const char *selector) {
    if (!selector || !selector[0]) return 0;
    const char *id = features->feature_ids ? features->feature_ids[idx] : features->feature_names[idx];
    const char *name = features->feature_names[idx];
    const char *ftype = (features->feature_types && features->feature_types[idx])
        ? features->feature_types[idx] : "";

    if (strncmp(selector, "feature_type:", 13) == 0) {
        char want[128];
        char have[128];
        normalize_token(selector + 13, want, sizeof(want));
        normalize_token(ftype, have, sizeof(have));
        return want[0] != '\0' && strcmp(want, have) == 0;
    }
    if (strncmp(selector, "id_prefix:", 10) == 0) {
        const char *prefix = selector + 10;
        return id && strncmp(id, prefix, strlen(prefix)) == 0;
    }
    if (strncmp(selector, "name_regex:", 11) == 0) {
        /* Minimal fallback: treat value as literal prefix on name or id. */
        const char *prefix = selector + 11;
        if (name && strncmp(name, prefix, strlen(prefix)) == 0) return 1;
        if (id && strncmp(id, prefix, strlen(prefix)) == 0) return 1;
        return 0;
    }
    return 0;
}

int pf_build_hash_feature_mask(feature_arrays *features,
                               const char *selector,
                               const char *library_feature_type,
                               unsigned char *out_mask,
                               int *n_hash_out,
                               int *n_protein_out) {
    if (!features || !out_mask) return -1;
    int n_hash = 0;
    int n_protein = 0;
    const int n = features->number_of_features;
    const int library_is_hash = pf_is_hash_like_feature_type(library_feature_type);

    for (int i = 0; i < n; ++i) {
        int is_hash = 0;
        if (selector && selector[0]) {
            is_hash = row_matches_selector(features, i, selector);
        } else if (library_is_hash) {
            is_hash = 1;
        } else {
            const char *ftype = (features->feature_types && features->feature_types[i])
                ? features->feature_types[i] : "";
            if (pf_is_hash_like_feature_type(ftype)) {
                is_hash = 1;
            } else {
                const char *id = features->feature_ids ? features->feature_ids[i] : features->feature_names[i];
                if (id && strncmp(id, "hashtag", 7) == 0) {
                    is_hash = 1;
                }
            }
        }
        out_mask[i] = (unsigned char)is_hash;
        if (is_hash) n_hash++;
        else n_protein++;
    }
    if (n_hash_out) *n_hash_out = n_hash;
    if (n_protein_out) *n_protein_out = n_protein;
    return 0;
}

static void invert_mask(const unsigned char *hash_mask, unsigned char *protein_mask, int n) {
    for (int i = 0; i < n; ++i) {
        protein_mask[i] = hash_mask[i] ? 0 : 1;
    }
}

typedef struct {
    int feature_idx;
    uint32_t count;
} hash_feat_count_t;

static int compare_feat_count_desc(const void *a, const void *b) {
    const hash_feat_count_t *aa = (const hash_feat_count_t *)a;
    const hash_feat_count_t *bb = (const hash_feat_count_t *)b;
    if (aa->count != bb->count) {
        return (aa->count < bb->count) ? 1 : -1;
    }
    return aa->feature_idx - bb->feature_idx;
}

static const char *feature_label(feature_arrays *features, int idx) {
    if (features->feature_ids && features->feature_ids[idx] && features->feature_ids[idx][0]) {
        return features->feature_ids[idx];
    }
    return features->feature_names[idx];
}

static int write_hash_demux_outputs(const pf_hash_mex_config *config,
                                    const unsigned char *hash_mask,
                                    int n_hash,
                                    pf_hash_demux_stats *stats_out) {
    if (!config || !config->base.counts || !config->base.features) return -1;
    const char *method = (config->hash_demux_method && config->hash_demux_method[0])
        ? config->hash_demux_method : PF_HASH_DEMUX_METHOD_RATIO;
    if (strcmp(method, PF_HASH_DEMUX_METHOD_RATIO) != 0) {
        fprintf(stderr, "Error: unsupported hash demux method '%s'\n", method);
        return -1;
    }

    char out_dir[4096];
    snprintf(out_dir, sizeof(out_dir), "%s", config->base.mex_output_dir);

    char assign_dir[4096];
    snprintf(assign_dir, sizeof(assign_dir), "%s", config->base.assign_output_dir);

    char barcodes_in[4096];
    snprintf(barcodes_in, sizeof(barcodes_in), "%s/barcodes.txt", assign_dir);
    FILE *bc_in = fopen(barcodes_in, "r");
    if (!bc_in) return -1;

    char assignments_path[4096];
    char singlet_path[4096];
    char doublet_path[4096];
    char negative_path[4096];
    snprintf(assignments_path, sizeof(assignments_path), "%s/hash_demux_assignments.tsv", out_dir);
    snprintf(singlet_path, sizeof(singlet_path), "%s/singlet_barcodes.tsv", out_dir);
    snprintf(doublet_path, sizeof(doublet_path), "%s/doublet_barcodes.tsv", out_dir);
    snprintf(negative_path, sizeof(negative_path), "%s/negative_barcodes.tsv", out_dir);

    FILE *assign_fp = fopen(assignments_path, "w");
    FILE *singlet_fp = fopen(singlet_path, "w");
    FILE *doublet_fp = fopen(doublet_path, "w");
    FILE *negative_fp = fopen(negative_path, "w");
    if (!assign_fp || !singlet_fp || !doublet_fp || !negative_fp) {
        if (assign_fp) fclose(assign_fp);
        if (singlet_fp) fclose(singlet_fp);
        if (doublet_fp) fclose(doublet_fp);
        if (negative_fp) fclose(negative_fp);
        fclose(bc_in);
        return -1;
    }

    fprintf(assign_fp,
            "barcode\thash_assignment\thash_classification\thash_total_umis\t"
            "hash_top_feature\thash_top_count\thash_second_feature\thash_second_count\t"
            "hash_top_ratio\n");

    pf_hash_demux_stats stats;
    memset(&stats, 0, sizeof(stats));
    stats.n_hash_features = n_hash;

    feature_arrays *features = config->base.features;
    pf_counts_result *counts = config->base.counts;
    const int n_features = features->number_of_features;

    char line[4096];
    while (fgets(line, sizeof(line), bc_in)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;

        unsigned char code[8];
        memset(code, 0, sizeof(code));
        string2code(line, barcode_length, code);
        uint32_t bkey = *(uint32_t *)code;

        hash_feat_count_t ranked[256];
        int n_ranked = 0;
        uint32_t total = 0;
        khint_t kb = kh_get(u32ptr, counts->barcode_to_deduped_hash, bkey);
        if (kb != kh_end(counts->barcode_to_deduped_hash)) {
            khash_t(u32u32) *feat_hash = (khash_t(u32u32) *)kh_val(counts->barcode_to_deduped_hash, kb);
            for (khint_t kf = kh_begin(feat_hash); kf != kh_end(feat_hash); ++kf) {
                if (!kh_exist(feat_hash, kf)) continue;
                uint32_t f_idx = kh_key(feat_hash, kf);
                if (f_idx == 0 || f_idx > (uint32_t)n_features) continue;
                if (!hash_mask[f_idx - 1]) continue;
                uint32_t val = kh_val(feat_hash, kf);
                if (val == 0) continue;
                if (n_ranked < (int)(sizeof(ranked) / sizeof(ranked[0]))) {
                    ranked[n_ranked].feature_idx = (int)f_idx - 1;
                    ranked[n_ranked].count = val;
                    n_ranked++;
                }
                total += val;
            }
        }
        qsort(ranked, (size_t)n_ranked, sizeof(ranked[0]), compare_feat_count_desc);

        const char *top_feat = "";
        const char *second_feat = "";
        uint32_t top_count = 0;
        uint32_t second_count = 0;
        if (n_ranked > 0) {
            top_feat = feature_label(features, ranked[0].feature_idx);
            top_count = ranked[0].count;
        }
        if (n_ranked > 1) {
            second_feat = feature_label(features, ranked[1].feature_idx);
            second_count = ranked[1].count;
        }

        const char *classification = "negative";
        char assignment_buf[512];
        assignment_buf[0] = '\0';
        const char *assignment = assignment_buf;
        double ratio = 0.0;
        if ((int)total < config->hash_min_total || (int)top_count < config->hash_min_top) {
            classification = "negative";
            stats.n_negative++;
            fprintf(negative_fp, "%s\n", line);
        } else if (second_count == 0 && top_count > 0) {
            classification = "singlet";
            snprintf(assignment_buf, sizeof(assignment_buf), "%s", top_feat);
            ratio = INFINITY;
            stats.n_singlet++;
            fprintf(singlet_fp, "%s\n", line);
        } else if (second_count > 0 && ((double)top_count / (double)second_count) >= config->hash_min_ratio) {
            classification = "singlet";
            snprintf(assignment_buf, sizeof(assignment_buf), "%s", top_feat);
            ratio = (double)top_count / (double)second_count;
            stats.n_singlet++;
            fprintf(singlet_fp, "%s\n", line);
        } else {
            classification = "doublet";
            if (n_ranked >= 2) {
                snprintf(assignment_buf, sizeof(assignment_buf), "%s|%s", top_feat, second_feat);
            } else {
                snprintf(assignment_buf, sizeof(assignment_buf), "%s", top_feat);
            }
            if (second_count > 0) {
                ratio = (double)top_count / (double)second_count;
            }
            stats.n_doublet++;
            fprintf(doublet_fp, "%s\n", line);
        }

        if (isinf(ratio)) {
            fprintf(assign_fp, "%s\t%s\t%s\t%u\t%s\t%u\t%s\t%u\tinf\n",
                    line, assignment, classification, total, top_feat, top_count,
                    second_feat, second_count);
        } else {
            fprintf(assign_fp, "%s\t%s\t%s\t%u\t%s\t%u\t%s\t%u\t%.6f\n",
                    line, assignment, classification, total, top_feat, top_count,
                    second_feat, second_count, ratio);
        }
    }

    fclose(bc_in);
    fclose(assign_fp);
    fclose(singlet_fp);
    fclose(doublet_fp);
    fclose(negative_fp);

    char summary_path[4096];
    snprintf(summary_path, sizeof(summary_path), "%s/hash_demux_summary.json", out_dir);
    FILE *summary_fp = fopen(summary_path, "w");
    if (summary_fp) {
        fprintf(summary_fp, "{\n");
        fprintf(summary_fp, "  \"method\": ");
        pf_fprint_json_string(summary_fp, method);
        fprintf(summary_fp, ",\n  \"hash_demux_enabled\": true,\n");
        if (config->hash_feature_selector && config->hash_feature_selector[0]) {
            fprintf(summary_fp, "  \"hash_feature_selector\": ");
            pf_fprint_json_string(summary_fp, config->hash_feature_selector);
            fprintf(summary_fp, ",\n");
        }
        fprintf(summary_fp, "  \"hash_min_total\": %d,\n", config->hash_min_total);
        fprintf(summary_fp, "  \"hash_min_top\": %d,\n", config->hash_min_top);
        fprintf(summary_fp, "  \"hash_min_ratio\": %.3f,\n", config->hash_min_ratio);
        fprintf(summary_fp, "  \"n_hash_features\": %d,\n", n_hash);
        fprintf(summary_fp, "  \"n_singlet\": %d,\n", stats.n_singlet);
        fprintf(summary_fp, "  \"n_doublet\": %d,\n", stats.n_doublet);
        fprintf(summary_fp, "  \"n_negative\": %d,\n", stats.n_negative);
        fprintf(summary_fp, "  \"hash_features\": [\n");
        int first = 1;
        for (int i = 0; i < n_features; ++i) {
            if (!hash_mask[i]) continue;
            if (!first) fprintf(summary_fp, ",\n");
            first = 0;
            fprintf(summary_fp, "    ");
            pf_fprint_json_string(summary_fp, feature_label(features, i));
        }
        fprintf(summary_fp, "\n  ]");
        if (config->base.command_line && config->base.command_line[0]) {
            fprintf(summary_fp, ",\n  \"command\": ");
            pf_fprint_json_string(summary_fp, config->base.command_line);
        }
        fprintf(summary_fp, "\n}\n");
        fclose(summary_fp);
    }

    char command_path[4096];
    snprintf(command_path, sizeof(command_path), "%s/hash_demux_command.txt", out_dir);
    FILE *cmd_fp = fopen(command_path, "w");
    if (cmd_fp) {
        fprintf(cmd_fp, "hash_demux=yes\n");
        fprintf(cmd_fp, "hash_demux_method=%s\n", method);
        if (config->hash_feature_selector && config->hash_feature_selector[0]) {
            fprintf(cmd_fp, "hash_feature_selector=%s\n", config->hash_feature_selector);
        }
        fprintf(cmd_fp, "hash_min_total=%d\n", config->hash_min_total);
        fprintf(cmd_fp, "hash_min_top=%d\n", config->hash_min_top);
        fprintf(cmd_fp, "hash_min_ratio=%.3f\n", config->hash_min_ratio);
        if (config->base.command_line) fprintf(cmd_fp, "%s\n", config->base.command_line);
        fclose(cmd_fp);
    }

    if (stats_out) *stats_out = stats;
    return 0;
}

static int copy_feature_ref_snapshot(const pf_adt_mex_config *config, const char *output_dir) {
    if (config->features && config->features->source_csv_path) {
        char snapshot[4096];
        snprintf(snapshot, sizeof(snapshot), "%s/feature_reference.csv", output_dir);
        if (pf_copy_file(config->features->source_csv_path, snapshot) != 0) {
            fprintf(stderr, "Warning: failed to copy feature reference snapshot to %s\n", snapshot);
        }
        return 0;
    }
    if (adt_feature_ref_path[0]) {
        char snapshot[4096];
        snprintf(snapshot, sizeof(snapshot), "%s/feature_reference.csv", output_dir);
        if (pf_copy_file(adt_feature_ref_path, snapshot) != 0) {
            fprintf(stderr, "Warning: failed to copy feature reference snapshot to %s\n", snapshot);
        }
    }
    return 0;
}

static int hash_outputs_required(const pf_hash_mex_config *config) {
    if (!config) return 0;
    if (config->hash_demux_mode == PF_HASH_DEMUX_YES) return 1;
    if (config->hash_feature_selector && config->hash_feature_selector[0]) return 1;
    if (pf_is_hash_like_feature_type(config->library_feature_type)) return 1;
    return 0;
}

int pf_write_adt_mex_outputs(const pf_hash_mex_config *config) {
    if (!config || !config->base.assign_output_dir || !config->base.mex_output_dir || !config->base.features) {
        return -1;
    }
    feature_arrays *features = config->base.features;
    const int n = features->number_of_features;
    unsigned char *hash_mask = calloc((size_t)n, sizeof(unsigned char));
    unsigned char *protein_mask = calloc((size_t)n, sizeof(unsigned char));
    if (!hash_mask || !protein_mask) {
        free(hash_mask);
        free(protein_mask);
        return -1;
    }

    int n_hash = 0;
    int n_protein = 0;
    if (pf_build_hash_feature_mask(features,
                                   config->hash_feature_selector,
                                   config->library_feature_type,
                                   hash_mask,
                                   &n_hash,
                                   &n_protein) != 0) {
        free(hash_mask);
        free(protein_mask);
        return -1;
    }
    invert_mask(hash_mask, protein_mask, n);

    if (n_hash == 0 && hash_outputs_required(config)) {
        fprintf(stderr,
                "Error: hash feature output required but no hash rows matched"
                " (hash_demux_mode=%d",
                config->hash_demux_mode);
        if (config->hash_feature_selector && config->hash_feature_selector[0]) {
            fprintf(stderr, ", selector='%s'", config->hash_feature_selector);
        }
        if (config->library_feature_type && config->library_feature_type[0]) {
            fprintf(stderr, ", library_feature_type='%s'", config->library_feature_type);
        }
        fprintf(stderr, ")\n");
        free(hash_mask);
        free(protein_mask);
        return -1;
    }

    if (mkdir_p(config->base.mex_output_dir) != 0) {
        free(hash_mask);
        free(protein_mask);
        return -1;
    }

    const int mixed = (n_hash > 0 && n_protein > 0);
    const int hash_only = (n_hash > 0 && n_protein == 0);
    int rc = 0;

    if (hash_only) {
        char hash_dir[4096];
        snprintf(hash_dir, sizeof(hash_dir), "%s/hash", config->base.mex_output_dir);
        if (mkdir_p(hash_dir) != 0) {
            rc = -1;
        } else {
            copy_feature_ref_snapshot(&config->base, hash_dir);
            if (pf_write_mex_10x_from_counts_subset(config->base.assign_output_dir,
                                                    hash_dir,
                                                    features,
                                                    config->base.counts,
                                                    hash_mask,
                                                    PF_HASH_FEATURE_TYPE_DEFAULT,
                                                    1) != 0) {
                fprintf(stderr, "Error: failed to write hash MEX files\n");
                rc = -1;
            }
        }
    } else if (mixed) {
        char hash_dir[4096];
        char protein_dir[4096];
        snprintf(hash_dir, sizeof(hash_dir), "%s/hash", config->base.mex_output_dir);
        snprintf(protein_dir, sizeof(protein_dir), "%s/protein", config->base.mex_output_dir);
        if (mkdir_p(hash_dir) != 0 || mkdir_p(protein_dir) != 0) {
            rc = -1;
        } else {
            copy_feature_ref_snapshot(&config->base, protein_dir);
            if (pf_write_mex_10x_from_counts_subset(config->base.assign_output_dir,
                                                    hash_dir,
                                                    features,
                                                    config->base.counts,
                                                    hash_mask,
                                                    PF_HASH_FEATURE_TYPE_DEFAULT,
                                                    1) != 0) {
                fprintf(stderr, "Error: failed to write hash MEX files\n");
                rc = -1;
            }
            if (pf_write_mex_10x_from_counts_subset(config->base.assign_output_dir,
                                                    protein_dir,
                                                    features,
                                                    config->base.counts,
                                                    protein_mask,
                                                    PF_ADT_FEATURE_TYPE_DEFAULT,
                                                    1) != 0) {
                fprintf(stderr, "Error: failed to write protein MEX files\n");
                rc = -1;
            }
        }
    } else {
        /* Protein-only legacy layout at mex_output_dir root. */
        if (pf_write_adt_protein_outputs(&config->base) != 0) {
            rc = -1;
        }
    }

    int run_demux = 0;
    if (config->hash_demux_mode == PF_HASH_DEMUX_YES) {
        run_demux = (n_hash > 0);
    } else if (config->hash_demux_mode == PF_HASH_DEMUX_AUTO) {
        run_demux = (n_hash > 0);
    }
    if (run_demux && rc == 0) {
        pf_hash_demux_stats demux_stats;
        if (write_hash_demux_outputs(config, hash_mask, n_hash, &demux_stats) != 0) {
            fprintf(stderr, "Error: hash demux output failed\n");
            rc = -1;
        }
    }

    free(hash_mask);
    free(protein_mask);
    return rc;
}
