#include "PipSeq.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace pipseq;

int main(int argc, char** argv) {
    assert(argc == 2);
    const std::string dir = argv[1];
    std::vector<std::string> paths;
    const std::string tiers[] = {"ACGTACGT", "GACTGC", "CTAGCA", "TGCATGCA"};
    for (unsigned t = 0; t < 4; ++t) {
        paths.push_back(dir + "/tier" + std::to_string(t) + ".txt");
        std::ofstream out(paths.back()); out << tiers[t] << '\n';
    }
    Decoder decoder(paths);
    const std::string body = "ACGTACGTATGGACTGCGAGCTAGCATCGAGTGCATGCA";
    assert(body.size() == 39);
    for (unsigned phase = 0; phase < 4; ++phase) {
        for (unsigned bin = 0; bin < 64; ++bin) {
            std::string bi;
            bi += "ACGT"[(bin >> 4) & 3]; bi += "ACGT"[(bin >> 2) & 3]; bi += "ACGT"[bin & 3];
            const auto d = decoder.decode(std::string(phase, 'G') + body + bi + "TTTT");
            assert(d.status == DecodeStatus::Accepted && d.phase == int(phase) && d.bin == int(bin));
            assert(d.mismatches == 0 && Decoder::cellSequence(d.cell) == "ACGTACGTGACTGCCTAGCATGCATGCA");
        }
    }
    auto changed = body;
    for (unsigned start : {0u, 11u, 20u, 31u}) changed[start] = changed[start] == 'A' ? 'C' : 'A';
    const auto corrected = decoder.decode(changed + "AAA");
    assert(corrected.status == DecodeStatus::Accepted && corrected.mismatches == 4 && corrected.bin == 0);
    assert(decoder.decode(body + "NNN").status == DecodeStatus::Accepted);
    assert(decoder.decode(body + "NNN").bin == -1);
    assert(decoder.decode(body).bin == -1);
    assert(decoder.decode(body.substr(0, 30)).status == DecodeStatus::InvalidBarcode);
    auto bad = body; bad[0] = 'N';
    assert(decoder.decode(bad + "AAA").status != DecodeStatus::Accepted);
    auto linker = body; linker.replace(8, 3, "CCC");
    assert(decoder.decode(linker + "CCC").status == DecodeStatus::Accepted); // linkers are not an extra gate
    // Two tier-1 entries one substitution apart, in both list orders.
    for (const char* list : {"CCGTACGT\nACGTACGT\n", "ACGTACGT\nCCGTACGT\n"}) {
        { std::ofstream out(paths[0]); out << list; }
        const std::string firstListed = std::string(list).substr(0, 8);
        Decoder reject(paths), first(paths, MultiMatch::FirstListed);
        for (const Decoder* d : {&reject, &first}) {
            // An exact entry is never rewritten to its listed neighbor.
            const auto a = d->decode(body + "AAA");
            assert(a.status == DecodeStatus::Accepted && a.mismatches == 0);
            assert(Decoder::cellSequence(a.cell).substr(0, 8) == "ACGTACGT");
            auto other = body; other[0] = 'C';
            const auto c = d->decode(other + "AAA");
            assert(c.status == DecodeStatus::Accepted && c.mismatches == 0);
            assert(Decoder::cellSequence(c.cell).substr(0, 8) == "CCGTACGT");
            // One substitution from exactly one entry is corrected to it.
            auto single = body; single[7] = 'A';
            const auto s = d->decode(single + "AAA");
            assert(s.status == DecodeStatus::Accepted && s.mismatches == 1);
            assert(Decoder::cellSequence(s.cell).substr(0, 8) == "ACGTACGT");
        }
        // One substitution from both entries: rejected, or first listed.
        auto both = body; both[0] = 'G';
        assert(reject.decode(both + "AAA").status == DecodeStatus::Ambiguous);
        const auto f = first.decode(both + "AAA");
        assert(f.status == DecodeStatus::Accepted && f.mismatches == 1);
        assert(Decoder::cellSequence(f.cell).substr(0, 8) == firstListed);
        // A tier beyond one substitution is a non-match, not ambiguous.
        auto far = both; far.replace(31, 8, "AAAAAAAA");
        assert(reject.decode(far + "AAA").status == DecodeStatus::NoMatch);
    }
    for (unsigned t = 0; t < 4; ++t) { std::ofstream out(paths[t]); out << std::string(t == 0 || t == 3 ? 8 : 6, 'A') << '\n'; }
    assert(Decoder(paths).decode(std::string(45, 'A')).status == DecodeStatus::PhaseTie);
    { std::ofstream out(paths[0]); out << "AAAAAAAA\nAAAAAAAA\n"; }
    bool failed = false;
    try { Decoder invalid(paths); } catch (const std::runtime_error&) { failed = true; }
    assert(failed);

    Group g;
    g.add(0, 100, false); g.add(0, 100, false); g.add(0, 101, false);
    g.add(0, 100, true); g.add(1, 100, false);
    assert(g.reads == 5 && g.fragments.size() == 4 && g.occupiedBins() == 2);
    Group other; other.add(0, 100, false); other.add(63, 100, false);
    g.merge(other);
    assert(g.reads == 7 && g.fragments.size() == 5 && g.occupiedBins() == 3);
    assert(transcripts(g, 2) == 3);
    Group boundary;
    for (unsigned b = 0; b < 10; ++b) { boundary.add(b, 100, false); boundary.add(b, 200, false); }
    assert(transcripts(boundary, 1.35) == 10);
    boundary.add(10, 100, false); boundary.add(10, 200, false);
    assert(transcripts(boundary, 1.35) == 16);
    assert(transcripts(boundary, 100) == 11);
    assert(transcripts(boundary, 1) == 22);
    assert(expectedMolecules(1) == 1 && std::abs(expectedMolecules(2) - (1 + 64.0 / 63)) < 1e-12);
    Groups low; low[{0, 0}] = g;
    assert(calibrate(low, 0).source == "unused_low_occupancy");
    Groups fit; fit[{0, 0}] = boundary;
    const auto c = calibrate(fit, 0);
    assert(c.groups == 1 && c.raw > 1 && c.raw < 2 && c.effective == c.raw);
    Group onePerBin; for (unsigned b = 0; b < 5; ++b) onePerBin.add(b, 100, false);
    fit.clear(); fit[{0, 0}] = onePerBin;
    assert(calibrate(fit, 0).raw < 1 && calibrate(fit, 0).effective == 1);
    Group saturated; for (unsigned b = 0; b < 64; ++b) saturated.add(b, 100, false);
    assert(saturated.occupiedBins() == 64 && transcripts(saturated, 1) == 64);
    fit.clear(); fit[{0, 0}] = saturated;
    failed = false; try { calibrate(fit, 0); } catch (const std::runtime_error&) { failed = true; }
    assert(failed && calibrate(fit, 2).source == "supplied");
    failed = false; try { transcripts(g, 0.5); } catch (const std::runtime_error&) { failed = true; }
    assert(failed);
    std::cout << "PASS PIP decoder and estimator unit checks\n";
}
