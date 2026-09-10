#ifndef STAR_SCRNA_TRACE_H
#define STAR_SCRNA_TRACE_H
#include "scrna_api.h"
#include "OrdMagStage.h"
#include "SparseCountView.h"
// Internal C++ result extension. Existing public C structures and functions stay ABI compatible.
struct ScrnaTrace {
    SimpleEmptyDropsResult ordmag;
    std::vector<uint32_t> retainIndices;
};
int scrnaEmptyDropsTrace(const scrna_matrix_input*, const scrna_ed_config*,
    const uint8_t*, uint32_t bootstrapThreads, scrna_ed_result*, ScrnaTrace*, const SparseCountView* = nullptr);
#endif
