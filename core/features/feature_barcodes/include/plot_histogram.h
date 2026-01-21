#ifndef PLOT_HISTOGRAM_H
#define PLOT_HISTOGRAM_H

#include "common.h"
#include <glib.h>
// EMfit.h removed - EM functionality no longer needed

// EM histogram plotting functions removed - no longer needed

void plot_simple_histogram(const char *directory,
                           const char *filename_leaf,
                           const char *title,
                           const char *xtitle,
                           const char *ytitle,
                           GArray *hist_data);

#endif // PLOT_HISTOGRAM_H
