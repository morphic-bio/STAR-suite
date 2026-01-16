#ifndef H_SampleMatrixData
#define H_SampleMatrixData

#include "IncludeDefine.h"
#include <string>
#include <vector>

// Structure to hold per-sample matrix data for OrdMag/EmptyDrops pipeline
struct SampleMatrixData {
    std::vector<uint32_t> nUMIperCB;              // UMI counts per cell
    std::vector<uint32_t> countCellGeneUMI;        // Sparse matrix: [geneIdx, cellIdx, count] triples
    std::vector<uint32_t> countCellGeneUMIindex;   // Index for each cell
    std::vector<uint32_t> nGenePerCB;              // Genes per cell
    std::vector<std::string> barcodes;             // Cell barcodes
    std::vector<std::string> features;             // Gene names
    uint32_t nCells;
    uint32_t nGenes;
    uint32_t countMatStride;                       // Stride for sparse matrix (3)
};

#endif
