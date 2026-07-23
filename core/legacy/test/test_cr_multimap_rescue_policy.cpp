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

int testDecoyMode() {
    int failed = 0;

    const auto naTie = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::NA, 50, {})}, true, EvidenceMode::Decoy);
    failed += check(!naTie.rescued && naTie.failure == Failure::NaBestTie,
                    "equal-score NA alignment vetoes false feature uniqueness");

    const auto naTieReordered = cr_multimap_rescue::evaluate({
        evidence(Annotation::NA, 50, {}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::Decoy);
    failed += check(!naTieReordered.rescued
                        && naTieReordered.failure == Failure::NaBestTie,
                    "NA veto reason is independent of best-alignment order");

    const auto lowerNa = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::NA, 49, {})}, true, EvidenceMode::Decoy);
    failed += check(lowerNa.rescued && lowerNa.geneIndex == 11
                        && lowerNa.winnerAlignIndex == 0,
                    "lower-score NA remains evidence but does not erase best annotated winner");

    const auto higherNa = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 49, {11}),
        evidence(Annotation::NA, 50, {})}, true, EvidenceMode::Decoy);
    failed += check(!higherNa.rescued && higherNa.failure == Failure::NaBestTie,
                    "higher-score NA prevents lower-score annotated rescue");

    const auto rpl22Mecom = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11}),
        evidence(Annotation::Intronic, 50, {22})}, true, EvidenceMode::Decoy);
    failed += check(!rpl22Mecom.rescued
                        && rpl22Mecom.failure == Failure::ConflictingBestGenes,
                    "equal-score RPL22-exonic/MECOM-intronic evidence remains ambiguous");

    const auto mecomRpl22 = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {22}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::Decoy);
    failed += check(!mecomRpl22.rescued
                        && mecomRpl22.failure == Failure::ConflictingBestGenes,
                    "conflicting-gene reason is independent of best-alignment order");

    const auto naAndMultiGene = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11, 22}),
        evidence(Annotation::NA, 50, {})}, true, EvidenceMode::Decoy);
    failed += check(!naAndMultiGene.rescued
                        && naAndMultiGene.failure == Failure::NaBestTie,
                    "NA veto has deterministic precedence over other best-score ambiguity");

    const auto sameGene = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {11}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::Decoy);
    failed += check(sameGene.rescued && !sameGene.intronicFallback
                        && sameGene.geneIndex == 11 && sameGene.winnerAlignIndex == 1,
                    "same-gene tied alignments are feature-unique and prefer exonic representative");

    const auto perAlignmentMultiGene = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 50, {11, 22}),
        evidence(Annotation::Exonic, 49, {11})}, true, EvidenceMode::Decoy);
    failed += check(!perAlignmentMultiGene.rescued
                        && perAlignmentMultiGene.failure == Failure::MultiGeneBestAlignment,
                    "multi-gene best alignment is not rescued");

    const auto intronicDisabled = cr_multimap_rescue::evaluate({
        evidence(Annotation::Intronic, 50, {11}),
        evidence(Annotation::Intronic, 50, {11})}, false, EvidenceMode::Decoy);
    failed += check(!intronicDisabled.rescued
                        && intronicDisabled.failure == Failure::IntronicFallbackOff,
                    "decoy mode respects disabled intronic fallback");

    const auto reordered = cr_multimap_rescue::evaluate({
        evidence(Annotation::Exonic, 49, {11}),
        evidence(Annotation::NA, 48, {}),
        evidence(Annotation::Exonic, 50, {11})}, true, EvidenceMode::Decoy);
    failed += check(reordered.rescued && reordered.winnerAlignIndex == 2
                        && reordered.geneIndex == 11,
                    "score tiebreak is independent of retained alignment order");
    return failed;
}

} // namespace

int main() {
    const int failed = testBiotypes() + testCompatibilityMode() + testDecoyMode();
    if (failed == 0) {
        std::cout << "PASS: CR multimap rescue policy\n";
        return 0;
    }
    return 1;
}
