#include "CrMultimapRescuePolicy.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

using cr_multimap_rescue::AlignmentEvidence;
using cr_multimap_rescue::Annotation;
using cr_multimap_rescue::EvidenceMode;
using cr_multimap_rescue::Failure;

int check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAIL: " << label << '\n';
        return 1;
    }
    return 0;
}

AlignmentEvidence evidence(Annotation annotation, std::int64_t score,
                           std::initializer_list<std::uint32_t> genes) {
    AlignmentEvidence result;
    result.annotation = annotation;
    result.score = score;
    result.genes.assign(genes.begin(), genes.end());
    return result;
}

int testBiotypes() {
    int failed = 0;
    failed += check(cr_multimap_rescue::biotypeIsCountable("protein_coding"),
                    "protein_coding is countable");
    failed += check(cr_multimap_rescue::biotypeIsCountable("lncRNA"),
                    "lncRNA is countable");
    failed += check(cr_multimap_rescue::biotypeIsCountable("protein_coding_LoF"),
                    "2024-A protein_coding_LoF remains countable");
    failed += check(cr_multimap_rescue::biotypeIsCountable("IG_V_pseudogene"),
                    "documented immune pseudogene exception remains countable");
    failed += check(!cr_multimap_rescue::biotypeIsCountable("processed_pseudogene"),
                    "processed pseudogene is non-countable");
    failed += check(!cr_multimap_rescue::biotypeIsCountable("unprocessed_pseudogene"),
                    "unprocessed pseudogene is non-countable");
    failed += check(!cr_multimap_rescue::biotypeIsCountable(""),
                    "missing biotype is NA/non-countable");
    failed += check(!cr_multimap_rescue::biotypeIsCountable("MissingGeneType"),
                    "formatter NA biotype is non-countable");
    return failed;
}

int testCompatibilityMode() {
    int failed = 0;
    const auto exonThenIntron = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::Intronic, 50, {22})}, true, EvidenceMode::Compatibility);
    failed += check(exonThenIntron.rescued && !exonThenIntron.intronicFallback
                        && exonThenIntron.winnerAlignIndex == 0,
                    "compatibility mode retains exon-first policy");

    const auto intronThenExon = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {22}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::Compatibility);
    failed += check(intronThenExon.rescued && !intronThenExon.intronicFallback
                        && intronThenExon.winnerAlignIndex == 1,
                    "compatibility winner is independent of alignment order");

    const auto intronicOnly = cr_multimap_rescue::evaluate({
        evidence(Annotation::NA, 50, {}),
        evidence(Annotation::Intronic, 50, {22})}, true, EvidenceMode::Compatibility);
    failed += check(intronicOnly.rescued && intronicOnly.intronicFallback
                        && intronicOnly.winnerAlignIndex == 1,
                    "compatibility mode retains intronic fallback");
    return failed;
}

int testAnnotatedBestMode() {
    int failed = 0;

    const auto ignoredNa = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::NA, 60, {})}, true, EvidenceMode::AnnotatedBest);
    failed += check(ignoredNa.rescued && ignoredNa.geneIndex == 11
                        && ignoredNa.winnerAlignIndex == 0,
                    "annotation-free alignment does not veto retained GTF evidence");

    const auto ignoredNaReordered = cr_multimap_rescue::evaluate({
        evidence(Annotation::NA, 60, {}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::AnnotatedBest);
    failed += check(ignoredNaReordered.rescued && ignoredNaReordered.geneIndex == 11
                        && ignoredNaReordered.winnerAlignIndex == 1,
                    "ignored annotation-free alignment is order independent");

    const auto mecomScoresHigher = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 49, {11}),
        evidence(Annotation::Intronic, 50, {22})}, true, EvidenceMode::AnnotatedBest);
    failed += check(mecomScoresHigher.rescued && mecomScoresHigher.intronicFallback
                        && mecomScoresHigher.geneIndex == 22
                        && mecomScoresHigher.winnerAlignIndex == 1,
                    "higher-score intronic MECOM evidence wins over RPL22 exon");

    const auto mecomScoresHigherReordered = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {22}),
        evidence(Annotation::Exonic, 49, {11})}, true, EvidenceMode::AnnotatedBest);
    failed += check(mecomScoresHigherReordered.rescued
                        && mecomScoresHigherReordered.intronicFallback
                        && mecomScoresHigherReordered.geneIndex == 22
                        && mecomScoresHigherReordered.winnerAlignIndex == 0,
                    "higher-score MECOM result is order independent");

    const auto rpl22ScoresHigher = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::Intronic, 49, {22})}, true, EvidenceMode::AnnotatedBest);
    failed += check(rpl22ScoresHigher.rescued && !rpl22ScoresHigher.intronicFallback
                        && rpl22ScoresHigher.geneIndex == 11
                        && rpl22ScoresHigher.winnerAlignIndex == 0,
                    "higher-score RPL22 exon wins over intronic MECOM evidence");

    const auto equalScoreConflict = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::Intronic, 50, {22})}, true, EvidenceMode::AnnotatedBest);
    failed += check(!equalScoreConflict.rescued
                        && equalScoreConflict.failure == Failure::ConflictingBestGenes,
                    "equal-score RPL22/MECOM evidence remains ambiguous");

    const auto equalScoreConflictReordered = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {22}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::AnnotatedBest);
    failed += check(!equalScoreConflictReordered.rescued
                        && equalScoreConflictReordered.failure == Failure::ConflictingBestGenes,
                    "equal-score conflict is order independent");

    const auto bestMultiGene = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11, 22}),
        evidence(Annotation::Exonic, 49, {11})}, true, EvidenceMode::AnnotatedBest);
    failed += check(!bestMultiGene.rescued
                        && bestMultiGene.failure == Failure::MultiGeneBestAlignment,
                    "best annotated alignment with multiple genes is rejected");

    const auto sameGene = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {11}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::AnnotatedBest);
    failed += check(sameGene.rescued && !sameGene.intronicFallback
                        && sameGene.geneIndex == 11 && sameGene.winnerAlignIndex == 1,
                    "same-gene best-score tie prefers an exonic representative");

    const auto intronicDisabled = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {11}),
        evidence(Annotation::Intronic, 50, {11})}, false, EvidenceMode::AnnotatedBest);
    failed += check(!intronicDisabled.rescued
                        && intronicDisabled.failure == Failure::IntronicFallbackOff,
                    "annotated-best mode respects disabled intronic fallback");

    const auto noAnnotatedEvidence = cr_multimap_rescue::evaluate({
        evidence(Annotation::NA, 60, {}),
        evidence(Annotation::NA, 50, {})}, true, EvidenceMode::AnnotatedBest);
    failed += check(!noAnnotatedEvidence.rescued
                        && noAnnotatedEvidence.failure == Failure::NoCountableBest,
                    "annotation-free alignments do not invent a gene");
    return failed;
}

} // namespace

int main() {
    const int failed = testBiotypes() + testCompatibilityMode()
        + testAnnotatedBestMode();
    if (failed == 0) {
        std::cout << "PASS: CR multimap rescue policy\n";
        return 0;
    }
    return 1;
}
