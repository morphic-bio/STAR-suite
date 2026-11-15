#include "SplitSeq.h"
#include <cassert>
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
    std::cout << "PASS split-seq decoding, independent corrections, replacement, ambiguity and validation\n";
}
