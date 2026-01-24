/**
 * @file heatmap.c
 * @brief Plotly-based heatmap generation (HTML + JSON)
 * 
 * Replaces the previous cairo-based PNG heatmaps with Plotly HTML+JSON outputs.
 * This removes the libcairo dependency while providing interactive heatmaps
 * and programmatic JSON data files.
 * 
 * Outputs:
 *   - Feature_types_heatmap.html + Feature_types_heatmap.json (richness/coexpression)
 *   - Feature_counts_heatmap.html + Feature_counts_heatmap.json (deduped counts)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/common.h"
#include "../include/globals.h"

/* Plasma colorscale for Plotly (subset of key points) */
static const char *PLASMA_COLORSCALE = 
    "[[0.0, 'rgb(13, 8, 135)'], "
    "[0.1, 'rgb(70, 3, 159)'], "
    "[0.2, 'rgb(114, 1, 168)'], "
    "[0.3, 'rgb(156, 23, 158)'], "
    "[0.4, 'rgb(189, 55, 134)'], "
    "[0.5, 'rgb(216, 87, 107)'], "
    "[0.6, 'rgb(237, 121, 83)'], "
    "[0.7, 'rgb(251, 159, 58)'], "
    "[0.8, 'rgb(253, 202, 38)'], "
    "[0.9, 'rgb(240, 248, 33)'], "
    "[1.0, 'rgb(240, 249, 33)']]";

/* -------------------------------------------------------------------------- *
 *  JSON helper functions                                                      *
 * -------------------------------------------------------------------------- */

static void json_escape_string(FILE *fp, const char *str) {
    fputc('"', fp);
    for (const char *p = str; *p; p++) {
        switch (*p) {
            case '"':  fprintf(fp, "\\\""); break;
            case '\\': fprintf(fp, "\\\\"); break;
            case '\n': fprintf(fp, "\\n"); break;
            case '\r': fprintf(fp, "\\r"); break;
            case '\t': fprintf(fp, "\\t"); break;
            default:   fputc(*p, fp); break;
        }
    }
    fputc('"', fp);
}

/* -------------------------------------------------------------------------- *
 *  Write JSON data file                                                       *
 * -------------------------------------------------------------------------- */

static int write_heatmap_json(const char *json_file,
                              const char *title,
                              feature_arrays *features,
                              int num_cols,
                              int *filter_mask,
                              int *column_sums,
                              int **matrix_data,
                              int is_deduped) {
    FILE *fp = fopen(json_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create JSON file %s\n", json_file);
        return -1;
    }
    
    fprintf(fp, "{\n");
    fprintf(fp, "  \"title\": ");
    json_escape_string(fp, title);
    fprintf(fp, ",\n");
    
    /* Feature labels (row labels) */
    fprintf(fp, "  \"feature_labels\": [");
    int first = 1;
    for (int i = 0; i < features->number_of_features; i++) {
        if (!filter_mask[i]) continue;
        if (!first) fprintf(fp, ", ");
        json_escape_string(fp, features->feature_names[i]);
        first = 0;
    }
    fprintf(fp, "],\n");
    
    /* Column labels (x-axis: 1, 2, 3, ...) */
    fprintf(fp, "  \"column_labels\": [");
    for (int j = 0; j < num_cols; j++) {
        if (j > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", j + 1);
    }
    fprintf(fp, "],\n");
    
    /* Column sums (for bar graph) */
    fprintf(fp, "  \"column_sums\": [");
    for (int j = 0; j < num_cols; j++) {
        if (j > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", column_sums[j]);
    }
    fprintf(fp, "],\n");
    
    /* Matrix data (2D array) */
    fprintf(fp, "  \"matrix\": [\n");
    first = 1;
    for (int i = 0; i < features->number_of_features; i++) {
        if (!filter_mask[i]) continue;
        if (!first) fprintf(fp, ",\n");
        fprintf(fp, "    [");
        for (int j = 0; j < num_cols; j++) {
            if (j > 0) fprintf(fp, ", ");
            /* matrix_data is 1-indexed with row i+1, col j+1 */
            int val = matrix_data[i + 1][j + 1];
            fprintf(fp, "%d", val);
        }
        fprintf(fp, "]");
        first = 0;
    }
    fprintf(fp, "\n  ],\n");
    
    /* Metadata */
    fprintf(fp, "  \"type\": \"%s\",\n", is_deduped ? "deduped_counts" : "coexpression");
    fprintf(fp, "  \"num_features\": %d,\n", features->number_of_features);
    
    int filtered_count = 0;
    for (int i = 0; i < features->number_of_features; i++) {
        if (filter_mask[i]) filtered_count++;
    }
    fprintf(fp, "  \"num_filtered_features\": %d,\n", filtered_count);
    fprintf(fp, "  \"num_columns\": %d\n", num_cols);
    fprintf(fp, "}\n");
    
    fclose(fp);
    return 0;
}

/* -------------------------------------------------------------------------- *
 *  Write HTML file with embedded Plotly                                       *
 * -------------------------------------------------------------------------- */

static int write_heatmap_html(const char *html_file,
                              const char *title,
                              const char *json_file_basename,
                              feature_arrays *features,
                              int num_cols,
                              int *filter_mask,
                              int *column_sums,
                              int max_column_sum,
                              int **matrix_data,
                              int max_value) {
    FILE *fp = fopen(html_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create HTML file %s\n", html_file);
        return -1;
    }
    
    /* HTML header */
    fprintf(fp, "<!DOCTYPE html>\n<html>\n<head>\n");
    fprintf(fp, "  <title>%s</title>\n", title);
    fprintf(fp, "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n");
    fprintf(fp, "  <style>\n");
    fprintf(fp, "    body { font-family: Arial, sans-serif; margin: 20px; }\n");
    fprintf(fp, "    .container { display: flex; flex-direction: column; }\n");
    fprintf(fp, "    #bargraph { width: 100%%; height: 200px; }\n");
    fprintf(fp, "    #heatmap { width: 100%%; height: 600px; }\n");
    fprintf(fp, "  </style>\n");
    fprintf(fp, "</head>\n<body>\n");
    fprintf(fp, "  <h2>%s</h2>\n", title);
    fprintf(fp, "  <p>Data file: <a href=\"%s\">%s</a></p>\n", json_file_basename, json_file_basename);
    fprintf(fp, "  <div class=\"container\">\n");
    fprintf(fp, "    <div id=\"bargraph\"></div>\n");
    fprintf(fp, "    <div id=\"heatmap\"></div>\n");
    fprintf(fp, "  </div>\n");
    fprintf(fp, "  <script>\n");
    
    /* Generate JavaScript data arrays */
    
    /* Column labels (x values) */
    fprintf(fp, "var x_labels = [");
    for (int j = 0; j < num_cols; j++) {
        if (j > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", j + 1);
    }
    fprintf(fp, "];\n");
    
    /* Feature labels (y values) */
    fprintf(fp, "var y_labels = [");
    int first = 1;
    for (int i = 0; i < features->number_of_features; i++) {
        if (!filter_mask[i]) continue;
        if (!first) fprintf(fp, ", ");
        fprintf(fp, "'");
        /* Escape single quotes in feature names */
        for (const char *p = features->feature_names[i]; *p; p++) {
            if (*p == '\'') fprintf(fp, "\\'");
            else fputc(*p, fp);
        }
        fprintf(fp, "'");
        first = 0;
    }
    fprintf(fp, "];\n");
    
    /* Column sums for bar graph */
    fprintf(fp, "var col_sums = [");
    for (int j = 0; j < num_cols; j++) {
        if (j > 0) fprintf(fp, ", ");
        fprintf(fp, "%d", column_sums[j]);
    }
    fprintf(fp, "];\n");
    
    /* Matrix data (z values) - rows are features, columns are counts/richness */
    fprintf(fp, "var z_data = [\n");
    first = 1;
    for (int i = 0; i < features->number_of_features; i++) {
        if (!filter_mask[i]) continue;
        if (!first) fprintf(fp, ",\n");
        fprintf(fp, "  [");
        for (int j = 0; j < num_cols; j++) {
            if (j > 0) fprintf(fp, ", ");
            int val = matrix_data[i + 1][j + 1];
            fprintf(fp, "%d", val);
        }
        fprintf(fp, "]");
        first = 0;
    }
    fprintf(fp, "\n];\n");
    
    /* Bar graph trace */
    fprintf(fp, "\n// Bar graph of column sums\n");
    fprintf(fp, "var bar_trace = {\n");
    fprintf(fp, "  x: x_labels,\n");
    fprintf(fp, "  y: col_sums,\n");
    fprintf(fp, "  type: 'bar',\n");
    fprintf(fp, "  marker: { color: 'rgb(51, 102, 204)' },\n");
    fprintf(fp, "  name: 'Column Totals'\n");
    fprintf(fp, "};\n");
    
    fprintf(fp, "var bar_layout = {\n");
    fprintf(fp, "  title: 'Column Totals',\n");
    fprintf(fp, "  xaxis: { title: 'Column' },\n");
    fprintf(fp, "  yaxis: { title: 'Total Count', range: [0, %d] },\n", max_column_sum > 0 ? max_column_sum : 1);
    fprintf(fp, "  margin: { t: 40, b: 40, l: 60, r: 20 },\n");
    fprintf(fp, "  showlegend: false\n");
    fprintf(fp, "};\n");
    
    fprintf(fp, "Plotly.newPlot('bargraph', [bar_trace], bar_layout);\n");
    
    /* Heatmap trace */
    fprintf(fp, "\n// Heatmap\n");
    fprintf(fp, "var heatmap_trace = {\n");
    fprintf(fp, "  x: x_labels,\n");
    fprintf(fp, "  y: y_labels,\n");
    fprintf(fp, "  z: z_data,\n");
    fprintf(fp, "  type: 'heatmap',\n");
    fprintf(fp, "  colorscale: %s,\n", PLASMA_COLORSCALE);
    fprintf(fp, "  zmin: 0,\n");
    fprintf(fp, "  zmax: %d,\n", max_value > 0 ? max_value : 1);
    fprintf(fp, "  colorbar: {\n");
    fprintf(fp, "    title: 'Count',\n");
    fprintf(fp, "    thickness: 20\n");
    fprintf(fp, "  },\n");
    fprintf(fp, "  hovertemplate: 'Feature: %%{y}<br>Column: %%{x}<br>Count: %%{z}<extra></extra>'\n");
    fprintf(fp, "};\n");
    
    fprintf(fp, "var heatmap_layout = {\n");
    fprintf(fp, "  title: '%s',\n", title);
    fprintf(fp, "  xaxis: { title: 'Feature Count / Richness', side: 'bottom' },\n");
    fprintf(fp, "  yaxis: { title: 'Feature', autorange: 'reversed' },\n");
    fprintf(fp, "  margin: { t: 50, b: 80, l: 150, r: 80 }\n");
    fprintf(fp, "};\n");
    
    fprintf(fp, "Plotly.newPlot('heatmap', [heatmap_trace], heatmap_layout);\n");
    
    fprintf(fp, "  </script>\n");
    fprintf(fp, "</body>\n</html>\n");
    
    fclose(fp);
    return 0;
}

/* -------------------------------------------------------------------------- *
 *  Public API: generate_heatmap (coexpression / richness heatmap)             *
 * -------------------------------------------------------------------------- */

void generate_heatmap(const char *directory,
                      feature_arrays *features,
                      int **coexpression_histograms)
{
    char html_file[1024];
    char json_file[1024];
    const char *sep = (directory[strlen(directory) - 1] == '/') ? "" : "/";
    
    snprintf(html_file, sizeof(html_file), "%s%sFeature_types_heatmap.html", directory, sep);
    snprintf(json_file, sizeof(json_file), "%s%sFeature_types_heatmap.json", directory, sep);
    
    const int num_rows = features->number_of_features;
    
    /* Find max_name_length (for future use) */
    int max_name_length = 0;
    for (int i = 0; i < num_rows; i++) {
        int len = strlen(features->feature_names[i]);
        if (len > max_name_length) max_name_length = len;
    }
    
    /* Create filter mask and calculate num_cols */
    int *filter_mask = malloc(num_rows * sizeof(int));
    if (!filter_mask) return;
    
    int filtered_rows = 0;
    int num_cols = 0;
    
    /* coexpression_histograms is 1-indexed with column 0 for totals */
    for (int i = 0; i < num_rows; i++) {
        filter_mask[i] = (coexpression_histograms[i + 1][0] > 0);
        if (filter_mask[i]) filtered_rows++;
        if (coexpression_histograms[i + 1][0] > num_cols) {
            num_cols = coexpression_histograms[i + 1][0];
        }
    }
    
    if (filtered_rows == 0 || num_cols == 0) {
        fprintf(stderr, "Warning: No data for coexpression heatmap\n");
        free(filter_mask);
        return;
    }
    
    /* Calculate column sums */
    int *column_sums = calloc(num_cols, sizeof(int));
    if (!column_sums) {
        free(filter_mask);
        return;
    }
    
    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; j++) {
            column_sums[j - 1] += coexpression_histograms[i + 1][j];
        }
    }
    
    /* Correct for overcounting: richness j overcounts by factor of j */
    for (int j = 1; j <= num_cols; j++) {
        if (j > 0) column_sums[j - 1] /= j;
    }
    
    int max_column_sum = 0;
    for (int j = 0; j < num_cols; j++) {
        if (column_sums[j] > max_column_sum) max_column_sum = column_sums[j];
    }
    
    /* Find max value for color scaling */
    int max_value = 0;
    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; j++) {
            int v = coexpression_histograms[i + 1][j];
            if (v > max_value) max_value = v;
        }
    }
    
    /* Write JSON data file */
    write_heatmap_json(json_file, "Feature Types (Richness) Heatmap",
                       features, num_cols, filter_mask, column_sums,
                       coexpression_histograms, 0);
    
    /* Write HTML file */
    write_heatmap_html(html_file, "Feature Types (Richness) Heatmap",
                       "Feature_types_heatmap.json",
                       features, num_cols, filter_mask, column_sums,
                       max_column_sum, coexpression_histograms, max_value);
    
    fprintf(stderr, "Feature types heatmap saved to: %s\n", html_file);
    fprintf(stderr, "Feature types data saved to: %s\n", json_file);
    
    free(column_sums);
    free(filter_mask);
}

/* -------------------------------------------------------------------------- *
 *  Public API: generate_deduped_heatmap (feature counts heatmap)              *
 * -------------------------------------------------------------------------- */

void generate_deduped_heatmap(const char *directory,
                              feature_arrays *features,
                              int **deduped_histograms,
                              int max_deduped_count,
                              int *total_deduped_counts,
                              int histogram_minimum_counts)
{
    char html_file[1024];
    char json_file[1024];
    const char *sep = (directory[strlen(directory) - 1] == '/') ? "" : "/";
    
    snprintf(html_file, sizeof(html_file), "%s%sFeature_counts_heatmap.html", directory, sep);
    snprintf(json_file, sizeof(json_file), "%s%sFeature_counts_heatmap.json", directory, sep);
    
    if (histogram_minimum_counts == -1) {
        histogram_minimum_counts = 1;
    }
    
    const int num_rows = features->number_of_features;
    
    /* Create filter mask */
    int *filter_mask = malloc(num_rows * sizeof(int));
    if (!filter_mask) return;
    
    int filtered_rows = 0;
    for (int i = 0; i < num_rows; i++) {
        filter_mask[i] = (total_deduped_counts[i] > histogram_minimum_counts);
        if (filter_mask[i]) filtered_rows++;
    }
    
    if (filtered_rows == 0) {
        fprintf(stderr, "No features passed the minimum count threshold for the deduplicated heatmap.\n");
        free(filter_mask);
        return;
    }
    
    int num_cols = max_deduped_count;
    
    /* Calculate column sums (starting from index 1) */
    long *column_sums_long = calloc(num_cols + 1, sizeof(long));
    if (!column_sums_long) {
        free(filter_mask);
        return;
    }
    
    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; j++) {
            column_sums_long[j] += deduped_histograms[i + 1][j];
        }
    }
    
    /* Convert to int array (drop unused column 0) */
    int *column_sums = malloc(num_cols * sizeof(int));
    if (!column_sums) {
        free(column_sums_long);
        free(filter_mask);
        return;
    }
    
    long max_column_sum = 0;
    for (int j = 0; j < num_cols; j++) {
        column_sums[j] = (int)column_sums_long[j + 1];
        if (column_sums_long[j + 1] > max_column_sum) {
            max_column_sum = column_sums_long[j + 1];
        }
    }
    free(column_sums_long);
    
    /* Find max value for color scaling */
    int max_value = 0;
    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; j++) {
            int v = deduped_histograms[i + 1][j];
            if (v > max_value) max_value = v;
        }
    }
    
    /* Write JSON data file */
    write_heatmap_json(json_file, "Feature Counts (Deduplicated) Heatmap",
                       features, num_cols, filter_mask, column_sums,
                       deduped_histograms, 1);
    
    /* Write HTML file */
    write_heatmap_html(html_file, "Feature Counts (Deduplicated) Heatmap",
                       "Feature_counts_heatmap.json",
                       features, num_cols, filter_mask, column_sums,
                       (int)max_column_sum, deduped_histograms, max_value);
    
    fprintf(stderr, "Feature counts heatmap saved to: %s\n", html_file);
    fprintf(stderr, "Feature counts data saved to: %s\n", json_file);
    
    free(column_sums);
    free(filter_mask);
}
