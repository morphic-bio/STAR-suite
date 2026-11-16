#include "PipSeq.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

    // Fragment store: identical (cell, gene, BI, IMI) collapse; strand and
    // bin are part of the key; stores filled separately merge to the same.
    const auto summaries = [](const Fragments& f) {
        std::vector<GroupSummary> out;
        f.forEachGroup([&](const GroupSummary& g) { out.push_back(g); });
        return out;
    };
    Fragments g;
    g.add(7, 3, 0, 100, false); g.add(7, 3, 0, 100, false); g.add(7, 3, 0, 101, false);
    g.add(7, 3, 0, 100, true); g.add(7, 3, 1, 100, false);
    g.finalize();
    auto s = summaries(g);
    assert(s.size() == 1 && s[0].cell == 7 && s[0].gene == 3);
    assert(s[0].reads == 5 && s[0].fragments == 4 && s[0].bins == 2);
    Fragments other; other.add(7, 3, 0, 100, false); other.add(7, 3, 63, 100, false);
    g.absorb(other);
    s = summaries(g);
    assert(other.size() == 0 && s.size() == 1);
    assert(s[0].reads == 7 && s[0].fragments == 5 && s[0].bins == 3);
    assert(transcripts(s[0].bins, s[0].fragments, 2) == 3);
    // Groups come out ordered by cell, then gene, whatever the input order.
    Fragments order;
    order.add(9, 1, 0, 5, false); order.add(2, 8, 0, 5, false); order.add(2, 4, 0, 5, false); order.add(9, 0, 1, 5, true);
    order.finalize();
    s = summaries(order);
    assert(s.size() == 4 && s[0].cell == 2 && s[0].gene == 4 && s[1].gene == 8 && s[2].cell == 9 && s[2].gene == 0 && s[3].gene == 1);
    // Many reads over few fragments: compaction during add() keeps the answer.
    Fragments heavy, split1, split2;
    for (unsigned i = 0; i < 3000000; ++i) {
        const unsigned bin = i % 7; const uint64_t at = 1000 + i % 11; const bool rev = i % 2;
        heavy.add(i % 3, i % 5, bin, at, rev);
        (i % 4 ? split1 : split2).add(i % 3, i % 5, bin, at, rev);
    }
    heavy.finalize();
    split1.absorb(split2);
    const auto h = summaries(heavy), m = summaries(split1);
    assert(h.size() == 15 && m.size() == 15 && heavy.size() == split1.size());
    uint64_t reads = 0;
    for (size_t i = 0; i < h.size(); ++i) {
        assert(h[i].cell == m[i].cell && h[i].gene == m[i].gene && h[i].reads == m[i].reads);
        assert(h[i].fragments == m[i].fragments && h[i].bins == m[i].bins);
        reads += h[i].reads;
    }
    assert(reads == 3000000);
    bool rangeFailed = false;
    try { Fragments bad; bad.add(0, 0, 64, 1, false); } catch (const std::runtime_error&) { rangeFailed = true; }
    assert(rangeFailed);

    const auto filled = [](unsigned bins, unsigned startsPerBin) {
        Fragments f;
        for (unsigned b = 0; b < bins; ++b)
            for (unsigned k = 0; k < startsPerBin; ++k) f.add(0, 0, b, 100 * (k + 1), false);
        f.finalize();
        return f;
    };
    // Estimator edges: bins <= 10 give the bin count; above ten,
    // max(bins, min(fragments, floor(fragments / IPM))).
    assert(transcripts(10, 20, 1.35) == 10);
    assert(transcripts(11, 22, 1.35) == 16);
    assert(transcripts(11, 22, 100) == 11);
    assert(transcripts(11, 22, 1) == 22);
    assert(transcripts(64, 64, 1) == 64);
    assert(expectedMolecules(1) == 1 && std::abs(expectedMolecules(2) - (1 + 64.0 / 63)) < 1e-12);
    assert(calibrate(g, 0).source == "unused_low_occupancy");
    const auto c = calibrate(filled(11, 2), 0);
    assert(c.groups == 1 && c.raw > 1 && c.raw < 2 && c.effective == c.raw);
    assert(calibrate(filled(5, 1), 0).raw < 1 && calibrate(filled(5, 1), 0).effective == 1);
    const auto saturated = filled(64, 1);
    assert(summaries(saturated)[0].bins == 64);
    failed = false; try { calibrate(saturated, 0); } catch (const std::runtime_error&) { failed = true; }
    assert(failed && calibrate(saturated, 2).source == "supplied");
    failed = false; try { transcripts(3, 5, 0.5); } catch (const std::runtime_error&) { failed = true; }
    assert(failed);
    // The worked example of the published documentation ("Transcript
    // Counting", figure): one barcode and gene, 19 distinct BI/IMI pairs in
    // 12 bins from 111 reads, correction factor 1.35, recorded count 14.
    {
        struct Row { const char* bi; uint64_t imi; unsigned reads; };
        const Row rows[] = {
            {"ACC", 140633101, 12}, {"ATA", 140633107, 3}, {"AGT", 140631790, 1}, {"AGT", 140632057, 4},
            {"AGT", 140632995, 10}, {"TTA", 140631790, 7}, {"TCG", 140633098, 3}, {"TCG", 140633101, 1},
            {"TCT", 140633108, 1}, {"CAA", 140633109, 3}, {"CGT", 140633455, 9}, {"CCG", 140632099, 5},
            {"CCG", 140632125, 16}, {"CCG", 140632128, 20}, {"CCG", 140632675, 1}, {"GGA", 140633106, 8},
            {"GTC", 140632098, 4}, {"GCC", 140632057, 1}, {"GCC", 140632126, 2}};
        const auto code = [](char ch) { return unsigned(ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2 : 3); };
        Fragments documented;
        for (const Row& r : rows)
            for (unsigned k = 0; k < r.reads; ++k)
                documented.add(1, 1, 16 * code(r.bi[0]) + 4 * code(r.bi[1]) + code(r.bi[2]), r.imi, false);
        documented.finalize();
        const auto d = summaries(documented);
        assert(d.size() == 1 && d[0].reads == 111 && d[0].fragments == 19 && d[0].bins == 12);
        assert(transcripts(d[0].bins, d[0].fragments, 1.35) == 14);
        // Calibrated on this group alone: 12 bins are expected from
        // sum_{j<12} 64/(64-j) = 13.17 molecules, so IPM = 19/13.17 = 1.44.
        const auto dc = calibrate(documented, 0);
        assert(dc.groups == 1 && std::abs(expectedMolecules(12) - 13.174204) < 1e-6);
        assert(std::abs(dc.raw - 19 / expectedMolecules(12)) < 1e-12 && dc.effective == dc.raw);
        assert(transcripts(12, 19, dc.effective) == 13);
    }

    // cDNA trimming. Sequences are numeric: A, C, G, T = 0..3, N = 4.
    const auto numeric = [](const std::string& text) {
        std::string out;
        for (char ch : text) out += char(ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2 : ch == 'T' ? 3 : 4);
        return out;
    };
    const std::string kit = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"; // template-switch oligo; holds both kit adapters
    const std::string tso = numeric(kit);
    const std::string cdna = "CGTACCGTTAGCATCGGATCCATGCTAGCTAGGCTAACGTTACGGATC";
    // One adapter given whole: it may hang off the read start by all but 12 bases.
    const auto clip5 = [&](const std::string& text) {
        const std::string n = numeric(text);
        return clip5pAdapter(n.data(), n.size(), tso.data(), tso.size(), 0.1);
    };
    assert(clip5(cdna) == 0);                                    // no adapter
    assert(clip5(kit + cdna) == 30);                             // complete, at the start
    assert(clip5("GGA" + kit + cdna) == 33);                     // complete, after other bases
    assert(clip5(kit.substr(18) + cdna) == 12);                  // last 12 bases at the start
    assert(clip5(kit.substr(19) + cdna) == 0);                   // 11 bases: below the overlap
    auto oneWrong = kit; oneWrong[5] = 'T';
    assert(clip5(oneWrong + cdna) == 30);                        // 1 of 30 mismatched
    auto fourWrong = kit; for (unsigned i : {2u, 9u, 16u, 23u}) fourWrong[i] = fourWrong[i] == 'A' ? 'C' : 'A';
    assert(clip5(fourWrong + cdna) == 0);                        // 4 of 30: above floor(0.1 * 30)
    auto partWrong = kit.substr(18); partWrong[3] = partWrong[3] == 'A' ? 'C' : 'A';
    assert(clip5(partWrong + cdna) == 12);                       // 1 of 12 allowed
    partWrong[7] = partWrong[7] == 'A' ? 'C' : 'A';
    assert(clip5(partWrong + cdna) == 0);                        // 2 of 12 not allowed
    auto withN = kit; withN[10] = 'N';
    assert(clip5(withN + cdna) == 30);                           // N is a mismatch
    assert(clip5(kit.substr(0, 20)) == 0);                       // the adapter never hangs off the read end
    assert(clip5(kit) == 30 && clip5("") == 0);
    // The kit's two 12-base adapters, as a list: each must be complete.
    const std::string adapters = numeric(kitAdapters5p);
    assert(std::string(kitAdapters5p) == kit.substr(18) + "," + kit.substr(11, 12));
    const auto clipKit = [&](const std::string& text) {
        const std::string n = numeric(text);
        return clip5pAdapter(n.data(), n.size(), adapters.data(), adapters.size(), 0.1);
    };
    assert(clipKit(cdna) == 0);
    assert(clipKit(kit + cdna) == 30);                           // both found; the later end wins the tie
    assert(clipKit("GGA" + kit + cdna) == 33);
    assert(clipKit(kit.substr(18) + cdna) == 12);                // template-switch adapter alone
    assert(clipKit(kit.substr(19) + cdna) == 0);                 // one base short
    assert(clipKit(kit.substr(0, 23) + cdna) == 23);             // SMART PCR adapter alone
    auto tailWrong = kit; tailWrong[20] = tailWrong[20] == 'A' ? 'C' : 'A';
    assert(clipKit(tailWrong + cdna) == 30);                     // one mismatch in each adapter
    tailWrong = kit; tailWrong[24] = 'C'; tailWrong[27] = 'C';
    assert(clipKit(tailWrong + cdna) == 23);                     // only the SMART PCR adapter is left
    const auto clip3 = [&](const std::string& text, char base = 'A') {
        const std::string n = numeric(text);
        return clip3pHomopolymer(n.data(), n.size(), base == 'A' ? 0 : 2, 0.1);
    };
    const std::string tail = "GCTTAGGCATCG";
    assert(clip3(cdna) == 0 && clip3(cdna, 'G') == 0);
    assert(clip3(cdna + std::string(10, 'A')) == 0);             // 10 A: two non-A in the last window
    assert(clip3(cdna + std::string(11, 'A')) == 12);            // the base before 11 A is the one mismatch
    assert(clip3(cdna + std::string(12, 'A')) == 12);
    assert(clip3(cdna + std::string(30, 'A') + tail) == 30 + tail.size()); // everything from the first A window
    assert(clip3(cdna + "AAAAACAAAAAA") == 12);                  // 11 of 12 A, no clean window
    assert(clip3(cdna + "AAAAACAAAAAAAAAAAA") == 12);            // a clean window wins over the earlier 11 of 12
    assert(clip3(cdna + "AAAACAAAACAAAA") == 0);                 // two non-A in every window
    assert(clip3(std::string(12, 'A')) == 12 && clip3("AAAA") == 0);
    assert(clip3(cdna + std::string(15, 'G'), 'G') == 15 && clip3(cdna + std::string(15, 'G')) == 0);
    assert(std::string(kitAdapters3p) == "polyA,polyG");
    // Tags: the barcode and BI as read at a phase.
    assert(Decoder::rawCell("GG" + body + "ACT", 2) == "ACGTACGTGACTGCCTAGCATGCATGCA" && Decoder::rawBin("GG" + body + "ACT", 2) == "ACT");
    assert(Decoder::rawCell(body, -1) == "-" && Decoder::rawBin(body, 0) == "-" && Decoder::rawCell(body.substr(0, 38), 0) == "-");
    std::cout << "PASS PIP decoder, fragment store, estimator and trimming unit checks\n";
}
