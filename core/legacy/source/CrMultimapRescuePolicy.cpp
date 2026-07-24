#include "CrMultimapRescuePolicy.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <unordered_set>

namespace cr_multimap_rescue {

bool biotypeIsCountable(const std::string& biotype) {
    if (biotype.empty() || biotype == "MissingGeneType" || biotype == "NA") {
        return false;
    }

    static const std::unordered_set<std::string> retainedImmunePseudogenes = {
        "IG_V_pseudogene",
        "IG_J_pseudogene",
        "IG_C_pseudogene",
        "TR_V_pseudogene",
        "TR_J_pseudogene"
    };
    if (retainedImmunePseudogenes.count(biotype) != 0) {
        return true;
    }

    std::string normalized = biotype;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    return normalized.find("pseudogene") == std::string::npos;
}

static Decision evaluateCompatibility(const std::vector<AlignmentEvidence>& evidence,
                                      bool allowIntronicFallback) {
    Decision decision;
    std::uint64_t exonicWinner = 0;
    std::uint64_t intronicWinner = 0;

    for (std::uint64_t index = 0; index < evidence.size(); ++index) {
        if (evidence[index].annotation == Annotation::Exonic) {
            ++decision.exonicCount;
            exonicWinner = index;
        } else if (evidence[index].annotation == Annotation::Intronic) {
            ++decision.intronicCount;
            intronicWinner = index;
        } else {
            ++decision.naCount;
        }
    }

    if (decision.exonicCount == 1) {
        decision.rescued = true;
        decision.winnerAlignIndex = exonicWinner;
        if (evidence[exonicWinner].genes.size() == 1) {
            decision.geneIndex = evidence[exonicWinner].genes.front();
        }
        return decision;
    }
    if (allowIntronicFallback && decision.exonicCount == 0 && decision.intronicCount == 1) {
        decision.rescued = true;
        decision.intronicFallback = true;
        decision.winnerAlignIndex = intronicWinner;
        if (evidence[intronicWinner].genes.size() == 1) {
            decision.geneIndex = evidence[intronicWinner].genes.front();
        }
        return decision;
    }

    if (!allowIntronicFallback && decision.exonicCount == 0 && decision.intronicCount == 1) {
        decision.failure = Failure::IntronicFallbackOff;
    } else if (decision.exonicCount == 0 && decision.intronicCount == 0) {
        decision.failure = Failure::NoCountableBest;
    } else {
        decision.failure = Failure::MultipleCompatibilityCandidates;
    }
    return decision;
}

static Decision evaluateAnnotatedBest(const std::vector<AlignmentEvidence>& evidence,
                                      bool allowIntronicFallback) {
    Decision decision;
    if (evidence.empty()) {
        decision.failure = Failure::NoCountableBest;
        return decision;
    }

    std::int64_t bestScore = std::numeric_limits<std::int64_t>::min();
    for (const AlignmentEvidence& item : evidence) {
        if (item.annotation != Annotation::NA && !item.genes.empty()) {
            bestScore = std::max(bestScore, item.score);
        }
    }
    if (bestScore == std::numeric_limits<std::int64_t>::min()) {
        decision.failure = Failure::NoCountableBest;
        return decision;
    }

    std::uint32_t commonGene = UINT32_MAX;
    std::uint64_t exonicWinner = evidence.size();
    std::uint64_t intronicWinner = evidence.size();
    bool hasMultiGeneBest = false;
    bool hasConflictingBestGenes = false;
    for (std::uint64_t index = 0; index < evidence.size(); ++index) {
        const AlignmentEvidence& item = evidence[index];
        if (item.score != bestScore || item.annotation == Annotation::NA
            || item.genes.empty()) {
            continue;
        }
        if (item.annotation == Annotation::Exonic) {
            ++decision.exonicCount;
            if (exonicWinner == evidence.size()) exonicWinner = index;
        } else if (item.annotation == Annotation::Intronic) {
            ++decision.intronicCount;
            if (intronicWinner == evidence.size()) intronicWinner = index;
        }

        if (item.genes.size() != 1) {
            hasMultiGeneBest = true;
        } else if (commonGene == UINT32_MAX) {
            commonGene = item.genes.front();
        } else if (commonGene != item.genes.front()) {
            hasConflictingBestGenes = true;
        }
    }

    // Diagnose the same best-score annotated evidence set identically
    // regardless of STAR's retained-alignment ordering. Annotation-free
    // alignments are deliberately absent from this decision.
    if (hasMultiGeneBest) {
        decision.failure = Failure::MultiGeneBestAlignment;
        return decision;
    }
    if (hasConflictingBestGenes) {
        decision.failure = Failure::ConflictingBestGenes;
        return decision;
    }

    if (commonGene == UINT32_MAX) {
        decision.failure = Failure::NoCountableBest;
        return decision;
    }
    if (exonicWinner != evidence.size()) {
        decision.rescued = true;
        decision.winnerAlignIndex = exonicWinner;
        decision.geneIndex = commonGene;
        return decision;
    }
    if (allowIntronicFallback && intronicWinner != evidence.size()) {
        decision.rescued = true;
        decision.intronicFallback = true;
        decision.winnerAlignIndex = intronicWinner;
        decision.geneIndex = commonGene;
        return decision;
    }

    decision.failure = Failure::IntronicFallbackOff;
    return decision;
}

Decision evaluate(const std::vector<AlignmentEvidence>& evidence,
                  bool allowIntronicFallback,
                  EvidenceMode mode) {
    return mode == EvidenceMode::AnnotatedBest
        ? evaluateAnnotatedBest(evidence, allowIntronicFallback)
        : evaluateCompatibility(evidence, allowIntronicFallback);
}

} // namespace cr_multimap_rescue
