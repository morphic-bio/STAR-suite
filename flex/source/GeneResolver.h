#ifndef CODE_GeneResolver
#define CODE_GeneResolver

#include <vector>
#include <cstdint>

#include "FlexGdna.h"

// Probe-aware gene resolver shared utility
// Extracted from CRKeyAggregator to allow reuse in both bam_to_counts and inline hash paths
//
// IMPORTANT: Multigene fanout is removed - this resolver returns a single gene or 0.
// Callers should insert exactly one quartet key per read group after resolution.

struct CandidateView {
    bool isGenomic;                      // false = probe alignment, true = genomic alignment
    uint16_t geneIdx15;                  // resolved gene if available (0 if unset)
    std::vector<uint16_t> zgGeneIdx15;   // parsed ZG genes mapped to 15-bit indices
    int mapq = -1;                       // alignment MAPQ if available; -1 means unknown
    int asScore = 0;                     // alignment score analogue (AS) if available
    int nm = -1;                         // mismatch count analogue (NM) if available
    FlexGdnaRegion probeRegion = FlexGdnaUnknown; // exact probe-contig region; genomic candidates are unknown
};

// Resolve a single gene from a group of candidates using their STAR scores.
// Returns: resolved gene index (15-bit), or 0 if ambiguous/unresolved/dropped
//
// Logic:
// - Rank all eligible probe and genomic alignments together by AS/NM.
// - Return the gene when every alignment in the best-score tier supports the
//   same gene.
// - Return 0 when equal-best evidence supports different genes.
// - Lower-scoring conflicting evidence does not veto a better alignment.
//
// Note: Multigene fanout is removed - this resolver returns a single gene or 0.
//
// Limitations:
// - Inline path: Uses gene→probe mapping cache (may diverge from bam_to_counts which uses ZG tags)
// - Candidate construction: Inline path builds one CandidateView per read (matches read grouping)
// - Probe/genomic fidelity: Inline path approximates isGenomic via probe index lookup (not ZG tag parsing)
//
// Parameters:
//   candidates: vector of candidate views (can be empty)
// Returns: resolved geneIdx15 (0 if ambiguous/no gene/drop)
uint16_t resolveGeneFromCandidates(const std::vector<CandidateView>& candidates);

#endif
