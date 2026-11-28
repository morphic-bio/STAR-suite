#include "SplitSeq.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using splitseq::Decoder;
using splitseq::Status;

int main(int argc, char** argv) {
    assert(argc == 2);
    const std::string dir = argv[1];
    auto file = [&](const std::string& name, const std::string& text) {
        const auto path = dir + "/" + name;
        std::ofstream out(path); out << text; out.close(); assert(out);
        return path;
    };
    const auto r1 = file("r1", "ACGTACGT TGCATGCA\nTGCATGCA GCTAGCTA\n");
    const auto r2 = file("r2", "GACTGACT\n");
    const auto r3 = file("r3", "CTAGCTAG\n");
    Decoder decoder({"78_8", "48_8", "10_8"}, "0_10", {r1, r2, r3});
    std::string read(86, 'A');
    read.replace(0, 10, "ACGTACGTAC"); read.replace(78, 8, "ACGTACGT");
    read.replace(48, 8, "GACTGACT"); read.replace(10, 8, "CTAGCTAG");
    auto exact = decoder.decode(read);
    assert(exact.status == Status::Accepted && exact.correctedRounds == 0);
    assert(exact.rawCell == "ACGTACGTGACTGACTCTAGCTAG");
    assert(exact.canonicalCell == "TGCATGCAGACTGACTCTAGCTAG"); // replacement once, not chain
    assert(exact.umi == "ACGTACGTAC" && exact.cell > UINT32_MAX);
    for (size_t tier = 0; tier < 3; ++tier) {
        for (unsigned pos = 0; pos < 8; ++pos) {
            for (char base : std::string("ACGT")) {
                auto variant = read;
                const unsigned offset = decoder.regions()[tier].offset + pos;
                if (base == read[offset]) continue;
                variant[offset] = base;
                auto corrected = decoder.decode(variant);
                assert(corrected.status == Status::Accepted && corrected.cell == exact.cell);
                assert(corrected.correctedRounds == (1U << tier));
            }
        }
    }
    auto variant = read; variant[78] = 'C'; variant[48] = 'A'; variant[10] = 'A';
    assert(decoder.decode(variant).correctedRounds == 7); // independent round corrections
    assert(decoder.decode(read.substr(0,85)).status == Status::ShortRead);
    // One non-ACGT base in a round is one substitution; two are invalid.
    variant = read; variant[78] = 'N';
    auto unknown = decoder.decode(variant);
    assert(unknown.status == Status::Accepted && unknown.cell == exact.cell && unknown.correctedRounds == 1);
    assert(unknown.rawCell[0] == 'N' && unknown.canonicalCell == exact.canonicalCell);
    variant[48] = 'N'; variant[17] = 'N';
    assert(decoder.decode(variant).correctedRounds == 7 && decoder.decode(variant).cell == exact.cell);
    variant = read; variant[78] = 'N'; variant[79] = 'N';
    assert(decoder.decode(variant).status == Status::InvalidBarcode);
    variant = read; variant[78] = 'N'; variant[80] = 'A'; // N plus a substitution: two away
    assert(decoder.decode(variant).status == Status::NoMatch);
    variant = read; variant.replace(48,8,"AAAAAAAA");
    assert(decoder.decode(variant).status == Status::NoMatch);
    const auto ambiguous = file("ambiguous", "CTAGCTAG\nCTAGCTAA\n");
    const auto reversed = file("reversed", "CTAGCTAA\nCTAGCTAG\n");
    Decoder tie({"78_8","48_8","10_8"}, "0_10", {r1,r2,ambiguous});
    Decoder tieReversed({"78_8","48_8","10_8"}, "0_10", {r1,r2,reversed});
    assert(tie.decode(read).status == Status::Accepted); // exact wins over H1
    variant = read; variant[17] = 'C';
    assert(tie.decode(variant).status == Status::Ambiguous);
    assert(tieReversed.decode(variant).status == Status::Ambiguous);
    auto unknownTie = read; unknownTie[17] = 'N'; // N at the base where two entries differ
    assert(tie.decode(unknownTie).status == Status::Ambiguous);
    assert(tieReversed.decode(unknownTie).status == Status::Ambiguous);
    unknownTie = read; unknownTie[10] = 'N';      // N elsewhere: one entry fits
    assert(tie.decode(unknownTie).status == Status::Accepted && tie.decode(unknownTie).correctedRounds == 4);
    // Two observed aliases of one canonical sequence are not a conflict.
    const auto alias = file("alias", "CTAGCTAG GGGGGGGG\nCTAGCTAA GGGGGGGG\n");
    Decoder aliased({"78_8","48_8","10_8"}, "0_10", {r1,r2,alias});
    const auto resolved = aliased.decode(variant);
    assert(resolved.status == Status::Accepted && resolved.correctedRounds == 4);
    assert(resolved.canonicalCell == "TGCATGCAGACTGACTGGGGGGGG");
    auto unknownAlias = read; unknownAlias[17] = 'N'; // both fillings are aliases of one sequence
    assert(aliased.decode(unknownAlias).status == Status::Accepted);
    assert(aliased.decode(unknownAlias).canonicalCell == "TGCATGCAGACTGACTGGGGGGGG");
    // A third neighbor with a different canonical sequence is a conflict.
    const auto mixedAlias = file("mixed_alias", "CTAGCTAG GGGGGGGG\nCTAGCTAA GGGGGGGG\nCTAGCTAT TTTTTTTT\n");
    Decoder conflicted({"78_8","48_8","10_8"}, "0_10", {r1,r2,mixedAlias});
    assert(conflicted.decode(variant).status == Status::Ambiguous);
    assert(conflicted.decode(read).status == Status::Accepted);
    auto fails = [](const std::vector<std::string>& pos, const std::string& umi,
                    const std::vector<std::string>& lists) {
        try { Decoder bad(pos, umi, lists); } catch (const std::exception&) { return; }
        assert(false && "invalid configuration accepted");
    };
    for (const auto& position : {"-1_8", "1_", "1_8_bad", "0_0", "0_32", "1000001_8", "0_x"})
        fails({position,"48_8"}, "0_10", {r1,r2});
    fails({"9_8","48_8"}, "0_10", {r1,r2});
    fails({"10_8","17_8"}, "0_10", {r1,r2});
    fails({"20_16","40_16"}, "0_10", {r1,r2});
    fails({"78_8","48_8"}, "0_17", {r1,r2});
    fails({"78_8","48_8"}, "0_10", {r1});
    for (const auto& bad : {"", "ACGTACGT\nACGTACGT\n", "NNNNNNNN\n", "ACGTACGT note\n", "ACGTACGT ACGTACGT extra\n"}) {
        auto path = file("bad", bad);
        fails({"78_8","48_8"}, "0_10", {path,r2});
    }
    const auto shortList = file("short", "ACG\n");
    Decoder mixed({"10_8","20_3"}, "0_10", {r2,shortList});
    std::string mixedRead = "ACGTACGTACGACTGACTAAACG";
    assert(mixed.decode(mixedRead).canonicalCell == "GACTGACTACG");

    Decoder edit2({"78_8","48_8","10_8"}, "0_10", {r1,r2,r3}, splitseq::Policy::Edit2);
    variant = read;
    variant.replace(78,8,"CGTACGTA"); // deletion + insertion at fixed coordinates
    variant.replace(48,8,"NNCTGACT");
    variant.replace(10,8,"NNAGCTAG");
    const auto six = edit2.decode(variant);
    assert(six.status == Status::Accepted && six.totalDistance == 6 && six.cell == exact.cell);
    assert(six.correctedRounds == 7 && six.umi == read.substr(0,10));
    variant[50] = 'N';
    assert(edit2.decode(variant).status == Status::NoMatch);
    Decoder tie2({"78_8","48_8","10_8"}, "0_10", {r1,r2,reversed}, splitseq::Policy::Edit2);
    variant = read;
    variant[17] = 'C';
    const auto tied = tie2.decode(variant);
    assert(tied.status == Status::Ambiguous && tied.totalDistance == 1);
    assert(tied.matches[2].alternatives.size() == 1 && tied.canonicalCell.empty());
    const auto alias2 = Decoder({"78_8","48_8","10_8"}, "0_10", {r1,r2,alias}, splitseq::Policy::Edit2);
    assert(alias2.decode(variant).status == Status::Accepted);
    // Two equally close primers for one well versus one primer for another.
    Decoder aliases({"0_3","3_3"}, "6_3",
                    {file("bayes_r1", "AAA AAA\nAAC AAA\nAAT AAT\n"), file("bayes_r2", "CGT\n")},
                    splitseq::Policy::Edit2);
    auto rawTie = aliases.decode("AAGCGTAGC");
    splitseq::ExactCounts priors;
    auto post = aliases.resolveTie(rawTie, "IIIIII", priors, 0.975);
    assert(rawTie.status == Status::Ambiguous && std::abs(post.best - 2.0/3) < 1e-12);
    auto exactPrimer = aliases.decode("AACCGTNNN");
    for (int i=0; i<78; ++i) splitseq::countExact(exactPrimer, priors);
    post = aliases.resolveTie(rawTie, "IIIIII", priors, 0.975);
    assert(rawTie.canonicalCell == "AAACGT" && std::abs(post.best - 80.0/81) < 1e-12);
    // With an equal joint prior, a low-quality mismatch can break a tie.
    Decoder qualDiff({"0_3","3_3"}, "6_3",
                     {file("q_diff", "ACA\nGCT\n"), file("q_diff2", "CAA\n")}, splitseq::Policy::Edit2);
    auto qualityTie = qualDiff.decode("ACTCAAAGC");
    auto otherTie = qualityTie;
    auto indistinguishable = qualityTie;
    auto flat = qualDiff.resolveTie(indistinguishable, "IIIIII", {}, 0.0);
    assert(std::abs(flat.best - 0.5) < 1e-12 && indistinguishable.status == Status::Ambiguous);
    qualDiff.resolveTie(qualityTie, "II!III", {}, 0.975);
    qualDiff.resolveTie(otherTie, "!IIIII", {}, 0.975);
    assert(qualityTie.canonicalCell == "ACACAA" && otherTie.canonicalCell == "GCTCAA");
    splitseq::ExactCounts full, prefix;
    auto train = [&](const Decoder& kit, const std::string& raw, unsigned total, unsigned early) {
        const auto observation = kit.decode(raw);
        for (unsigned i = 0; i < total; ++i)
            splitseq::countExact(observation, full, i < early ? &prefix : nullptr);
    };
    train(qualDiff, "ACACAAAGC", 3, 2);
    train(qualDiff, "GCTCAAAGC", 2, 1);
    auto backup = qualDiff.decode("ACTCAAAGC");
    auto rescued = qualDiff.resolveTie(backup, "DIIIII", full, 0.95, &prefix);
    assert(rescued.frequencyRescued && backup.canonicalCell == "ACACAA");
    assert(rescued.cells.at(backup.cell) < rescued.best); // weaker posterior can win the backup
    auto primary = qualDiff.decode("ACTCAAAGC");
    assert(!qualDiff.resolveTie(primary, "!IIIII", full, 0.95, &prefix).frequencyRescued);
    assert(primary.canonicalCell == "GCTCAA");
    // Sum compatible raw primers for each physical well without pseudocounts.
    full.clear(); prefix.clear();
    train(aliases, "AAACGTAGC", 2, 1);
    train(aliases, "AACCGTAGC", 3, 1);
    train(aliases, "AATCGTAGC", 4, 2);
    rawTie = aliases.decode("AAGCGTAGC");
    post = aliases.resolveTie(rawTie, "IIIIII", full, 0.95, &prefix);
    assert(rawTie.status == Status::Ambiguous && !post.frequencyRescued); // prefix is 2:2
    train(aliases, "AACCGTAGC", 1, 1);
    post = aliases.resolveTie(rawTie, "IIIIII", full, 0.95, &prefix);
    assert(post.frequencyRescued && rawTie.canonicalCell == "AAACGT");
    train(aliases, "AATCGTAGC", 8, 0); // later observations reverse the full-library leader
    rawTie = aliases.decode("AAGCGTAGC");
    post = aliases.resolveTie(rawTie, "IIIIII", full, 0.95, &prefix);
    assert(post.frequencyRescued && rawTie.canonicalCell == "AAACGT");
    splitseq::ExactCounts noPrefixSupport;
    rawTie = aliases.decode("AAGCGTAGC");
    post = aliases.resolveTie(rawTie, "IIIIII", full, 0.95, &noPrefixSupport);
    assert(!post.frequencyRescued && rawTie.status == Status::Ambiguous);
    const auto small = file("fourEdit2","AAA\n");
    Decoder four({"0_3","3_3","6_3","9_3"}, "12_3", {small,small,small,small}, splitseq::Policy::Edit2);
    assert(four.decode("AACACCAACACCACG").totalDistance == 6);
    assert(four.decode("AACACCAACACCACG").status == Status::Accepted);
    assert(four.decode("ACCACCAACACCACG").status == Status::NoMatch); // 2+2+1+2
    std::cout << "PASS split-seq decoding, independent corrections, replacement, ambiguity and validation\n";
}
