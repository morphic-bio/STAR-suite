#include "../include/plot_histogram.h"
#include "../include/globals.h"
// EMfit.h removed - EM functionality no longer needed
#include <math.h>
#include <string.h>
#include <stdlib.h>

// Helper function to print uint32_t array as JavaScript
static void print_js_array_uint32(FILE *fp, const char *var_name, const uint32_t *data, int len) {
    fprintf(fp, "var %s = [", var_name);
    for(int i=0; i<len; ++i) {
        fprintf(fp, "%u%s", data[i], (i==len-1?"":","));
    }
    fprintf(fp, "];\n");
}

void plot_simple_histogram(const char *directory,
                           const char *filename_leaf,
                           const char *title,
                           const char *xtitle,
                           const char *ytitle,
                           vec_u32_t *hist_data) {
    if (!hist_data || vec_u32_size(hist_data) == 0) {
        fprintf(stderr, "Warning: Empty histogram data provided for plotting '%s'\n", title);
        return;
    }

    char filename[FILENAME_LENGTH];
    sprintf(filename, "%s/%s", directory, filename_leaf);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create file %s\n", filename);
        return;
    }

    int hist_len = vec_u32_size(hist_data);
    uint32_t *data = hist_data->a;

    // --- Write HTML and Plotly JavaScript ---
    fprintf(fp, "<!DOCTYPE html><html><head><title>%s</title>", title);
    fprintf(fp, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script></head>");
    fprintf(fp, "<body><div id=\"plot\" style=\"width:100%%;height:600px;\"></div><script>\n");

    // --- Generate Data Variables ---
    print_js_array_uint32(fp, "hist_y", data, hist_len);
    fprintf(fp, "var x_axis = [");
    for(int i=0; i<hist_len; ++i) fprintf(fp, "%d%s", i, (i==hist_len-1?"":","));
    fprintf(fp, "];\n");
    
    // --- Histogram Trace ---
    fprintf(fp, "var hist_trace = { x: x_axis, y: hist_y, type: 'bar', name: 'Frequency' };\n");

    // --- Combine Data into a single array for Plotly ---
    fprintf(fp, "var data = [hist_trace];\n");

    // --- Define Layout ---
    int last_bin_with_data = 0;
    for (int i = hist_len - 1; i >= 1; i--) {
        if (data[i] >= 1) { last_bin_with_data = i; break; }
    }
    int x_axis_max = (last_bin_with_data > 0) ? ((last_bin_with_data + 4) / 5 * 5) : 5;
     if (x_axis_max == 0) x_axis_max = 5;

    fprintf(fp, "var layout = {\n");
    fprintf(fp, "    title: { text: '%s', x: 0.5 },\n", title);
    fprintf(fp, "    xaxis: { title: '%s', range: [0, %d] },\n", xtitle, x_axis_max);
    fprintf(fp, "    yaxis: { title: '%s', type: 'linear' },\n", ytitle);
    fprintf(fp, "    showlegend: true,\n");
    fprintf(fp, "    plot_bgcolor: 'rgba(240,240,240,0.95)',\n");
    fprintf(fp, "    paper_bgcolor: 'rgba(240,240,240,0.95)'\n");
    fprintf(fp, "};\n");
    
    // --- Final Plotting Call ---
    fprintf(fp, "Plotly.newPlot('plot', data, layout);\n");
    fprintf(fp, "</script></body></html>\n");

    fclose(fp);
    fprintf(stderr, "Histogram plot saved to: %s\n", filename);
}