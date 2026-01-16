#ifndef CODE_GeneResolver
#define CODE_GeneResolver

#include <vector>
#include <cstdint>

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
    bool probeCigarOk = true;            // for probe alignments: true if CIGAR passes probe QA
};

// Resolve a single gene from a group of candidates using probe-first policy
// Returns: resolved gene index (15-bit), or 0 if ambiguous/unresolved/dropped
//
// Logic:
// - If probe alignments exist:
//   - All probe genes agree (single distinct gene) → returns that gene
//   - Probe genes disagree (multiple distinct) → returns 0 (drop)
//   - Ignores genomic hits when probe exists
// - If no probe alignments:
//   - Exactly one genomic assignment → returns that gene
//   - Multiple genomic genes → falls through to fallback (don't blanket-drop)
//   - No genes → falls through to fallback
// - Fallback: unique-max ZG tie-break across all candidates
//   - Unique max → returns that gene
//   - Tie → returns 0 (drop)
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
