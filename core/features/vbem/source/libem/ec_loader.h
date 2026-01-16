#ifndef EC_LOADER_H
#define EC_LOADER_H

#include "em_types.h"
#include <string>

// Load equivalence classes from Salmon eq_classes.txt format
// Format:
//   Line 1: number of transcripts
//   Line 2: number of equivalence classes
//   Lines 3 to (3+n_transcripts-1): transcript names
//   Remaining lines: EC definitions as <n_txp> <idx1> <idx2> ... <count>
ECTable load_ec_file(const std::string& ec_path, std::vector<std::string>& transcript_names);

// Load transcript lengths from Salmon quant.sf format (TSV with Name, Length columns)
// If quant_path is empty, lengths will be set to 1.0 (placeholder)
void load_transcript_lengths(TranscriptState& state, const std::string& quant_path);

// Load effective lengths from simple TSV (transcript_index <tab> eff_length)
// Overrides effective lengths in state
void load_effective_lengths(TranscriptState& state, const std::string& eff_lengths_path);

#endif // EC_LOADER_H
