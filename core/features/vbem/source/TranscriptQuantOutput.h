#ifndef CODE_TranscriptQuantOutput
#define CODE_TranscriptQuantOutput

#include "em_types.h"
#include <string>

class Transcriptome;  // Forward declaration

// Write quantification results to Salmon-compatible quant.sf format
void writeQuantSF(const EMResult& result, const TranscriptState& state, const std::string& filename);

// Write gene-level quantification to quant.genes.sf format (Legacy mode: raw counts)
// Returns: 0 on success, 1 on file open error, 2 if MissingGeneID genes encountered
int writeQuantGeneSF(
    const EMResult& result,
    const TranscriptState& state,
    const Transcriptome& tr,
    const std::string& filename
);

// Write gene-level quantification with tximport-style lengthScaledTPM counts
// Returns: 0 on success, 1 on file open error, 2 if MissingGeneID genes encountered
int writeQuantGeneSFTximport(
    const EMResult& result,
    const TranscriptState& state,
    const Transcriptome& tr,
    const std::string& filename
);

#endif // CODE_TranscriptQuantOutput

