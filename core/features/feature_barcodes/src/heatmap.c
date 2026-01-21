#include <cairo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/common.h"
#include "../include/globals.h"
#include "plasma_colormap_16.h"
#include "plasma_colormap_64.h"
#include "plasma_colormap_256.h"
#include "plasma_colormap_1024.h"

#define CELL_SIZE 10  // Size of each cell in the heatmap
#define BAR_WIDTH 20  // Width of the color bar
#define BASE_PADDING 10  // Base padding for additional adjustments
#define BAR_GRAPH_HEIGHT 100   // Height of the bar graph area

const double (*select_colormap(int dynamic_range, int *colormap_size))[3] {
    if (dynamic_range < 20) {
        *colormap_size = 16;
        return plasma_colormap_16;
    } else if (dynamic_range < 100) {
        *colormap_size = 64;
        return plasma_colormap_64;
    } else if (dynamic_range < 1000) {
        *colormap_size = 256;
        return plasma_colormap_256;
    } else {
        *colormap_size = 1024;
        return plasma_colormap_1024;
    }
}


// Function to map a normalized value (0–1) to an RGB color using a colormap
void value_to_color(double normalized_value, double *r, double *g, double *b, const double colormap[][3], int colormap_size) {
    // Clamp the normalized value to [0, 1]
    if (normalized_value < 0.0) normalized_value = 0.0;
    if (normalized_value > 1.0) normalized_value = 1.0;

    // Scale the normalized value to the colormap size
    int index = (int)(normalized_value * (colormap_size - 1)); // Map to 0–(colormap_size-1)

    // Map to RGB values from the provided colormap
    *r = colormap[index][0];
    *g = colormap[index][1];
    *b = colormap[index][2];
}


// Function to normalize values for coloring
double normalize(int value, int max_value) {
    return (double)value / max_value;
}

/* -------------------------------------------------------------------------- *
 *  Shared drawing logic                                                      *
 * -------------------------------------------------------------------------- */
typedef int (*value_cb)(int row, int col, void *ctx);

/* central renderer – used by both public helpers */
static void draw_heatmap_core(const char *output_file,
                              feature_arrays *features,
                              int left_padding,
                              int right_padding,
                              int num_cols,
                              int filtered_rows,
                              int *column_sums,
                              int max_column_sum,
                              int max_value,
                              int *filter_mask,
                              value_cb value_fn,
                              void *value_ctx,
                              int clear_background)
{
    int width  = left_padding + num_cols * CELL_SIZE + BAR_WIDTH + right_padding;
    int height = BASE_PADDING + BAR_GRAPH_HEIGHT + BASE_PADDING + 20 +
                 filtered_rows * CELL_SIZE + BASE_PADDING;

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
                                                          width, height);
    cairo_t *cr = cairo_create(surface);

    if (clear_background) {          /* only needed by the deduped map    */
        cairo_set_source_rgb(cr, 1, 1, 1);
        cairo_paint(cr);
    }

    /* ---------- bar graph (column sums) ---------- */
    for (int j = 0; j < num_cols; j++) {
        double bar_h = 0.0;
        if (max_column_sum > 0)
            bar_h = (double)column_sums[j] / max_column_sum * BAR_GRAPH_HEIGHT;

        cairo_set_source_rgb(cr, 0.20, 0.40, 0.80);
        cairo_rectangle(cr, left_padding + j * CELL_SIZE,
                        BASE_PADDING + BAR_GRAPH_HEIGHT - bar_h,
                        CELL_SIZE, bar_h);
        cairo_fill(cr);
    }

    /* side labels (max / mid) */
    cairo_set_font_size(cr, 10);
    cairo_set_source_rgb(cr, 0, 0, 0);
    char lbl[32];
    snprintf(lbl, sizeof lbl, "%d", max_column_sum);
    cairo_move_to(cr, left_padding - 40, BASE_PADDING + 5);
    cairo_show_text(cr, lbl);

    snprintf(lbl, sizeof lbl, "%d", max_column_sum / 2);
    cairo_move_to(cr, left_padding - 40, BASE_PADDING + BAR_GRAPH_HEIGHT / 2 + 5);
    cairo_show_text(cr, lbl);

    /* x-axis labels every 5th column */
    for (int j = 0; j < num_cols; j++) {
        if (j == 0 || (j + 1) % 5 == 0) {
            snprintf(lbl, sizeof lbl, "%d", j + 1);
            cairo_move_to(cr, left_padding + j * CELL_SIZE + CELL_SIZE / 4,
                          BASE_PADDING + BAR_GRAPH_HEIGHT + 15);
            cairo_show_text(cr, lbl);
        }
    }

    /* pick a colour-map based on dynamic range */
    int cmap_sz;
    const double (*plasma)[3] = select_colormap(max_value, &cmap_sz);

    /* ---------- coloured matrix ---------- */
    int draw_row = 0;
    for (int i = 0; i < features->number_of_features; ++i) {
        if (!filter_mask[i]) continue;

        for (int j = 0; j < num_cols; ++j) {
            double v      = value_fn(i, j, value_ctx);
            double inten  = normalize((int)v, max_value);
            double r, g, b;
            value_to_color(inten, &r, &g, &b, plasma, cmap_sz);

            cairo_set_source_rgb(cr, r, g, b);
            cairo_rectangle(cr,
                            left_padding + j * CELL_SIZE,
                            BASE_PADDING + BAR_GRAPH_HEIGHT + 20
                                + draw_row * CELL_SIZE,
                            CELL_SIZE, CELL_SIZE);
            cairo_fill(cr);
        }

        /* row label */
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
                               CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10);
        cairo_move_to(cr, BASE_PADDING,
                      BASE_PADDING + BAR_GRAPH_HEIGHT + 20
                          + draw_row * CELL_SIZE + CELL_SIZE / 2);
        cairo_show_text(cr, features->feature_names[i]);
        draw_row++;
    }

    /* ---------- colour-bar ---------- */
    for (int i = 0; i < height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20; i++) {
        double norm = 1.0 - (double)i /
                      (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20);
        double r, g, b;
        value_to_color(norm, &r, &g, &b, plasma, cmap_sz);

        cairo_set_source_rgb(cr, r, g, b);
        cairo_rectangle(cr,
                        width - right_padding - BAR_WIDTH,
                        BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i,
                        BAR_WIDTH, 1);
        cairo_fill(cr);
    }

    cairo_set_font_size(cr, 10);
    for (int i = 0;
         i <= height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20;
         i += 10 * CELL_SIZE) {
        double norm = 1.0 - (double)i /
                      (height - BAR_GRAPH_HEIGHT - BASE_PADDING * 3 - 20);
        int v = (int)(norm * max_value);
        snprintf(lbl, sizeof lbl, "%d", v);
        cairo_move_to(cr,
                      width - right_padding + 5,
                      BASE_PADDING + BAR_GRAPH_HEIGHT + 20 + i + 4);
        cairo_show_text(cr, lbl);
    }

    cairo_surface_write_to_png(surface, output_file);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

/* ---------- small helpers that fetch the actual cell value ---------- */
static int coexpression_value(int row, int col, void *ctx)
{
    int **hist = ctx;                 /* co-expression matrix */
    return hist[row + 1][col + 1];
}

void generate_heatmap(const char *directory,
                      feature_arrays *features,
                      int **coexpression_histograms)
{
    char output_file[1024];
    if (directory[strlen(directory) - 1] == '/') {
        snprintf(output_file, sizeof(output_file), "%sFeature_types_heatmap.png", directory);
    } else {
        snprintf(output_file, sizeof(output_file), "%s/Feature_types_heatmap.png", directory);
    }
    //find max_name_length
    int max_name_length=0;
    const int num_rows = features->number_of_features;
    for (int i=0; i<num_rows; i++){
        if (strlen(features->feature_names[i])>max_name_length){
            max_name_length=strlen(features->feature_names[i]);
        }
    }

    // Calculate dynamic padding
    int left_padding = BASE_PADDING + max_name_length * 7; // 7 pixels per character as a rough estimate
    int right_padding = BASE_PADDING + 50; // Allow space for color bar labels
    int filtered_rows = 0;
    int filter_mask[num_rows];
    int num_cols = 0;
    
    //create a filter mask - necessary later for color map and calculate max column
    //remember that coexpression_histograms is 1 indexed with the 0 column saved for totals
    for (int i = 0; i < num_rows; i++) {
       filter_mask[i]=coexpression_histograms[i+1][0]>0;
         if (filter_mask[i]){
              filtered_rows++;
         }
       if (coexpression_histograms[i+1][0] > num_cols){
           num_cols=coexpression_histograms[i+1][0];
       }
    }
    // calculate column sums and find the maximum column sum
    int *column_sums = calloc(num_cols, sizeof(int));
    if (!column_sums) return;

    for (int i = 0; i < num_rows; i++) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; j++) { // j is total features
            column_sums[j-1] += coexpression_histograms[i+1][j];
        }
    }

    // Correct for overcounting for the top histogram
    for (int j = 1; j <= num_cols; j++) {
        // a richness of j overcounts by a factor of j
        if (j > 0) column_sums[j-1] /= j;
    }

    int max_column_sum = 0;
    for (int j = 0; j < num_cols; j++) {
        if (column_sums[j] > max_column_sum) {
            max_column_sum = column_sums[j];
        }
    }
    int max_value = 0;
    for (int i = 0; i < num_rows; ++i) {
        if (!filter_mask[i]) continue;
        for (int j = 1; j <= num_cols; ++j) {
            int v = coexpression_histograms[i + 1][j];
            if (v > max_value) max_value = v;
        }
    }
    // Adjust canvas dimensions to include the bar graph
//    int width = left_padding + num_cols * CELL_SIZE + BAR_WIDTH + right_padding;
//;    int height = BASE_PADDING + BAR_GRAPH_HEIGHT + BASE_PADDING + 20 + filtered_rows * CELL_SIZE + BASE_PADDING; // +20 for x-axis labels

    draw_heatmap_core(output_file, features,
                      left_padding, right_padding,
                      num_cols, filtered_rows,
                      column_sums,                /* bar-graph */
                      max_column_sum,
                      max_value,                  /* ← correct scale for colours */
                      filter_mask,
                      coexpression_value, coexpression_histograms,
                      0);
    free(column_sums);
} 


void generate_deduped_heatmap(const char *directory,
                              feature_arrays *features,
                              int **deduped_histograms,
                              int max_deduped_count,
                              int *total_deduped_counts,
                              int histogram_minimum_counts)
{
    char output_file[1024];

    if (histogram_minimum_counts == -1) {
        histogram_minimum_counts = 1;
    }

    if (directory[strlen(directory) - 1] == '/') {
        snprintf(output_file, sizeof(output_file), "%sFeature_counts_heatmap.png", directory);
    } else {
        snprintf(output_file, sizeof(output_file), "%s/Feature_counts_heatmap.png", directory);
    }

    int max_name_length = 0;
    const int num_rows = features->number_of_features;
    for (int i = 0; i < num_rows; i++) {
        if (strlen(features->feature_names[i]) > max_name_length) {
            max_name_length = strlen(features->feature_names[i]);
        }
    }

    int left_padding = BASE_PADDING + max_name_length * 7;
    int right_padding = BASE_PADDING + 50;
    
    int filtered_rows = 0;
    int filter_mask[num_rows];
    size_t num_cols = max_deduped_count;

    for (int i = 0; i < num_rows; i++) {
        filter_mask[i] = (total_deduped_counts[i] > histogram_minimum_counts);
        if (filter_mask[i]) {
            filtered_rows++;
        }
    }

    if (filtered_rows == 0) {
        fprintf(stderr, "No features passed the minimum count threshold for the deduplicated heatmap.\n");
        return;
    }

    long column_sums[num_cols + 1];
    memset(column_sums, 0, sizeof(column_sums));
    
    for (int i = 0; i < num_rows; i++) {
        if(filter_mask[i]){
            for(guint j=1; j <= num_cols; ++j){
                column_sums[j] += deduped_histograms[i+1][j];
            }
        }
    }

    /* recompute max_column_sum **excluding** the unused j==0 column */
    long max_column_sum = 0;
    for (size_t j = 1; j <= num_cols; ++j)        /* start at 1 */
        if (column_sums[j] > max_column_sum) max_column_sum = column_sums[j];

    /* find the maximum single-cell value (skip index 0) */
    int max_value = 0;
    for (int i = 0; i < num_rows; ++i) {
        if (!filter_mask[i]) continue;
        for (guint j = 1; j <= num_cols; ++j) {     /* start at 1 */
            int v = deduped_histograms[i + 1][j];
            if (v > max_value) max_value = v;
        }
    }

    /* drop the unused first column */
    int effective_cols = (int)num_cols;
    int *col_sums_int  = malloc(sizeof(int) * effective_cols);
    for (int j = 0; j < effective_cols; ++j)
        col_sums_int[j] = (int)column_sums[j + 1];

    draw_heatmap_core(output_file, features,
                      left_padding, right_padding,
                      effective_cols, filtered_rows,
                      col_sums_int,
                      max_column_sum,             /* bar-graph */
                      max_value,                  /* ← correct colour scale    */
                      filter_mask,
                      coexpression_value, deduped_histograms,
                      1);

    free(col_sums_int);
    fprintf(stderr, "Feature counts heatmap saved to: %s\n", output_file);
} 



