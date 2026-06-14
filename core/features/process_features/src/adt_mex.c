#include "../include/adt_mex.h"
#include "../include/mex_10x.h"
#include "../include/utils.h"
#include "../include/globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int count_lines_in_file(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    int lines = 0;
    char buf[4096];
    while (fgets(buf, sizeof(buf), fp)) {
        if (buf[0] != '\0' && buf[0] != '\n') lines++;
    }
    fclose(fp);
    return lines;
}

static int write_provenance(const pf_adt_mex_config *config, int n_features, int n_barcodes, long long total_umis) {
    char path[4096];
    snprintf(path, sizeof(path), "%s/protein_quant_summary.json", config->mex_output_dir);

    FILE *fp = fopen(path, "w");
    if (!fp) return -1;

    fprintf(fp, "{\n");
    fprintf(fp, "  \"mode\": \"adt_mex\",\n");
    fprintf(fp, "  \"feature_type\": \"%s\",\n", PF_ADT_FEATURE_TYPE_DEFAULT);
    if (config->features && config->features->source_csv_path) {
        fprintf(fp, "  \"feature_ref_path\": \"%s\",\n", config->features->source_csv_path);
    } else if (adt_feature_ref_path[0]) {
        fprintf(fp, "  \"feature_ref_path\": \"%s\",\n", adt_feature_ref_path);
    }
    if (config->features && config->features->source_csv_fingerprint) {
        fprintf(fp, "  \"feature_ref_fingerprint\": \"%s\",\n", config->features->source_csv_fingerprint);
    }
    fprintf(fp, "  \"read_layout\": \"10x_feature_barcode\",\n");
    fprintf(fp, "  \"barcode_read\": \"R1\",\n");
    fprintf(fp, "  \"feature_read\": \"R2\",\n");
    fprintf(fp, "  \"barcode_length\": %d,\n", config->barcode_length);
    fprintf(fp, "  \"umi_length\": %d,\n", config->umi_length);
    fprintf(fp, "  \"barcode_offset\": %d,\n", config->barcode_offset);
    fprintf(fp, "  \"feature_offset\": %d,\n", config->feature_offset);
    fprintf(fp, "  \"max_hamming_distance\": %d,\n", config->max_hamming_distance);
    fprintf(fp, "  \"umi_dedup_stringency\": %d,\n", config->stringency);
    fprintf(fp, "  \"num_features\": %d,\n", n_features);
    fprintf(fp, "  \"num_barcodes\": %d,\n", n_barcodes);
    fprintf(fp, "  \"total_assigned_umis\": %lld,\n", total_umis);
    fprintf(fp, "  \"per_feature_assigned_umis\": {\n");
    if (config->features && config->counts && config->counts->total_deduped_counts) {
        for (int i = 0; i < config->features->number_of_features; ++i) {
            const char *label = config->features->feature_names[i];
            long long umi = (long long)config->counts->total_deduped_counts[i];
            fprintf(fp, "    \"%s\": %lld%s\n", label, umi, (i + 1 < config->features->number_of_features) ? "," : "");
        }
    }
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"multiomics_manifest\": {\n");
    fprintf(fp, "    \"protein.mex_dir\": \"%s\",\n", config->mex_output_dir);
    if (config->features && config->features->source_csv_path) {
        char snapshot[4096];
        snprintf(snapshot, sizeof(snapshot), "%s/feature_reference.csv", config->mex_output_dir);
        fprintf(fp, "    \"protein.feature_ref\": \"%s\",\n", snapshot);
    }
    fprintf(fp, "    \"protein.normalization\": \"clr\"\n");
    fprintf(fp, "  }");
    if (config->command_line && config->command_line[0]) {
        fprintf(fp, ",\n  \"command\": ");
        fputc('"', fp);
        for (const char *p = config->command_line; *p; ++p) {
            if (*p == '"' || *p == '\\') fputc('\\', fp);
            fputc(*p, fp);
        }
        fputc('"', fp);
    }
    fprintf(fp, "\n}\n");
    fclose(fp);
    return 0;
}

static int write_command_manifest(const pf_adt_mex_config *config) {
    char path[4096];
    snprintf(path, sizeof(path), "%s/protein_quant_command.txt", config->mex_output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "mode=adt_mex\n");
    fprintf(fp, "feature_type=%s\n", PF_ADT_FEATURE_TYPE_DEFAULT);
    if (config->command_line) fprintf(fp, "%s\n", config->command_line);
    fclose(fp);
    return 0;
}

int pf_write_adt_protein_outputs(const pf_adt_mex_config *config) {
    if (!config || !config->assign_output_dir || !config->mex_output_dir || !config->features) {
        return -1;
    }

    if (mkdir_p(config->mex_output_dir) != 0) {
        return -1;
    }

  if (config->features->source_csv_path) {
        char snapshot[4096];
        snprintf(snapshot, sizeof(snapshot), "%s/feature_reference.csv", config->mex_output_dir);
        if (pf_copy_file(config->features->source_csv_path, snapshot) != 0) {
            fprintf(stderr, "Warning: failed to copy feature reference snapshot to %s\n", snapshot);
        }
    } else if (adt_feature_ref_path[0]) {
        char snapshot[4096];
        snprintf(snapshot, sizeof(snapshot), "%s/feature_reference.csv", config->mex_output_dir);
        if (pf_copy_file(adt_feature_ref_path, snapshot) != 0) {
            fprintf(stderr, "Warning: failed to copy feature reference snapshot to %s\n", snapshot);
        }
    }

    if (pf_write_mex_10x_from_features(config->assign_output_dir,
                                       config->mex_output_dir,
                                       config->features,
                                       PF_ADT_FEATURE_TYPE_DEFAULT,
                                       1) != 0) {
        fprintf(stderr, "Error: failed to write ADT protein MEX files\n");
        return -1;
    }

    char barcodes_path[4096];
    snprintf(barcodes_path, sizeof(barcodes_path), "%s/barcodes.txt", config->assign_output_dir);
    int n_barcodes = count_lines_in_file(barcodes_path);
    if (n_barcodes < 0) n_barcodes = 0;

    long long total_umis = 0;
    if (config->counts && config->counts->total_deduped_counts) {
        for (int i = 0; i < config->features->number_of_features; ++i) {
            total_umis += config->counts->total_deduped_counts[i];
        }
    }

    if (write_provenance(config, config->features->number_of_features, n_barcodes, total_umis) != 0) {
        fprintf(stderr, "Warning: failed to write protein_quant_summary.json\n");
    }
    write_command_manifest(config);
    return 0;
}
