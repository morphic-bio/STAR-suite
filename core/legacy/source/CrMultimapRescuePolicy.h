#ifndef STAR_CR_MULTIMAP_RESCUE_POLICY_H
#define STAR_CR_MULTIMAP_RESCUE_POLICY_H

#include <cstdint>
#include <string>
#include <vector>

namespace cr_multimap_rescue {

// NA is alignment evidence without a countable GTF feature. It participates
// in alignment-score comparison and uniqueness, but can never be emitted as a
// gene assignment.
enum class Annotation : std::uint8_t {
    NA = 0,
    Intronic = 1,
    Exonic = 2
};

enum class EvidenceMode : std::uint8_t {
    Compatibility = 0,
    Decoy = 1
};

enum class Failure : std::uint8_t {
    None = 0,
    NoCountableBest = 1,
    NaBestTie = 2,
    ConflictingBestGenes = 3,
    MultiGeneBestAlignment = 4,
    MultipleCompatibilityCandidates = 5,
    IntronicFallbackOff = 6
};

struct AlignmentEvidence {
    Annotation annotation = Annotation::NA;
    std::vector<std::uint32_t> genes;
    std::int64_t score = 0;
};

struct Decision {
    bool rescued = false;
    bool intronicFallback = false;
    std::uint64_t winnerAlignIndex = 0;
    std::uint32_t geneIndex = UINT32_MAX;
    std::uint64_t exonicCount = 0;
    std::uint64_t intronicCount = 0;
    std::uint64_t naCount = 0;
    Failure failure = Failure::None;
};

// The reference formatter excludes ordinary pseudogene biotypes from the
// countable feature axis. Keep the predicate here as defense in depth for
// custom indexes. Missing/NA biotypes are non-countable evidence.
bool biotypeIsCountable(const std::string& biotype);

// Compatibility preserves the historical exon-first CR rescue. Decoy mode
// first retains only best-score alignments, then requires every tied alignment
// to support exactly the same single countable gene. NA participates in that
// uniqueness test but is never a possible winner.
Decision evaluate(const std::vector<AlignmentEvidence>& evidence,
                  bool allowIntronicFallback,
                  EvidenceMode mode);

} // namespace cr_multimap_rescue

#endif
