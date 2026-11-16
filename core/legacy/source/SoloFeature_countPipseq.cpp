#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "PipSeq.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include <iomanip>
#include <limits>

// One aligned read of a PIP library: called from SoloReadFeature::record().
void SoloReadFeature::recordPipseq(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, ReadAnnotations &readAnnot) {
    const pipseq::Decoded& decoded = soloBar.pipDecoded;
    if (decoded.status != pipseq::DecodeStatus::Accepted) { ++pipReadStats[0]; return; }
    if (decoded.bin < 0) { ++pipReadStats[1]; return; }
    if (nTr == 0) { ++pipReadStats[2]; return; }
    if (nTr != 1) { ++pipReadStats[3]; return; }
    const auto& genes = readAnnot.annotFeatures[featureType].fSet;
    if (genes.size() != 1 || alignOut[0]->nExons == 0) { ++pipReadStats[4]; return; }
    const Transcript& a = *alignOut[0];
    const bool reverse = a.Str != 0;
    const unsigned block = reverse ? a.nExons - 1 : 0;
    // The aligned 5' base: clipping is not extrapolated; strand is part
    // of the key. Internal splice structure does not change this IMI.
    const uint64_t coordinate = a.exons[block][EX_G] +
        (reverse ? a.exons[block][EX_L] - 1 : 0);
    pipFragments.add(decoded.cell, *genes.begin(), decoded.bin, coordinate, reverse);
    ++pipReadStats[5];
}

void SoloFeature::countPipseq() {
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of PIPseq counting error: " + message + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    pipseq::Fragments fragments;
    std::array<uint64_t, 8> barcodeStats{};
    std::array<uint64_t, 6> readStats{};
    try {
        for (int t = 0; t < P.runThreadN; ++t) {
            auto* sr = RAchunk[t]->RA->soloRead;
            auto* rf = sr->readFeat[pSolo.featureInd[featureType]];
            fragments.absorb(rf->pipFragments);
            for (size_t i = 0; i < barcodeStats.size(); ++i) barcodeStats[i] += sr->readBar->pipStats[i];
            for (size_t i = 0; i < readStats.size(); ++i) readStats[i] += rf->pipReadStats[i];
        }
    } catch (const std::exception& e) { fail(e.what()); }
    const auto checked32 = [&](uint64_t n) -> uint32_t {
        if (n > std::numeric_limits<uint32_t>::max()) fail("count or matrix offset exceeds STAR's 32-bit matrix capacity");
        return static_cast<uint32_t>(n);
    };
    pipseq::Calibration calibration;
    try { calibration = pipseq::calibrate(fragments, pSolo.pipseqIPM); }
    catch (const std::exception& e) { fail(e.what()); }

    outputPrefix = P.outFileNamePrefix + pSolo.outFileNames[0] + SoloFeatureTypes::Names[featureType] + '/';
    outputPrefixFiltered = outputPrefix + "filtered/";
    createDirectory(outputPrefix, P.runDirPerm, "PIP output directory", P);
    ofstream& ledger = ofstrOpen(outputPrefix + "Pipseq.groups.tsv", ERROR_OUT, P);
    ledger << "barcode\tgene_id\treads\toccupied_bins\tunique_bi_imi\testimated_transcripts\n";

    pSolo.cbWLstr.clear(); pSolo.cbWLstrOut.clear();
    indCB.clear(); indCBwl.clear(); nUMIperCB.clear(); nGenePerCB.clear();
    nReadPerCB.clear(); nReadPerCBunique.clear(); nReadPerCBtotal.clear();
    countCellGeneUMI.clear(); countCellGeneUMIindex.clear();
    countMatStride = 2;
    nReadsInput = barcodeStats[0]; nReadsMapped = readStats[5];
    uint64_t previous = UINT64_MAX;
    uint64_t totalTranscripts = 0, totalFragments = 0, nGroups = 0;
    fragments.forEachGroup([&](const pipseq::GroupSummary& group) {
        if (group.gene >= Trans.nGe) fail("gene index exceeds loaded annotation");
        if (group.cell != previous || pSolo.cbWLstr.empty()) {
            previous = group.cell;
            const uint32_t i = checked32(pSolo.cbWLstr.size());
            pSolo.cbWLstr.push_back(pipseq::Decoder::cellSequence(group.cell));
            indCB.push_back(i); indCBwl.push_back(i);
            countCellGeneUMIindex.push_back(checked32(countCellGeneUMI.size()));
            nUMIperCB.push_back(0); nGenePerCB.push_back(0);
            nReadPerCB.push_back(0); nReadPerCBunique.push_back(0); nReadPerCBtotal.push_back(0);
        }
        const uint32_t count = checked32(pipseq::transcripts(group.bins, group.fragments, calibration.effective));
        countCellGeneUMI.push_back(group.gene); countCellGeneUMI.push_back(count);
        nUMIperCB.back() = checked32(uint64_t(nUMIperCB.back()) + count);
        ++nGenePerCB.back();
        nReadPerCB.back() = checked32(uint64_t(nReadPerCB.back()) + group.reads);
        nReadPerCBunique.back() = nReadPerCBtotal.back() = nReadPerCB.back();
        totalTranscripts += count; totalFragments += group.fragments; ++nGroups;
        ledger << pSolo.cbWLstr.back() << '\t' << Trans.geID[group.gene] << '\t' << group.reads
               << '\t' << group.bins << '\t' << group.fragments << '\t' << count << '\n';
    });
    ledger.close();
    if (!ledger) fail("could not write Pipseq.groups.tsv");
    countCellGeneUMIindex.push_back(checked32(countCellGeneUMI.size()));
    pSolo.cbWLsize = pSolo.cbWLstr.size();
    pSolo.cbWLstrOut = pSolo.cbWLstr;
    nCB = checked32(pSolo.cbWLsize);
    nReadPerCBmax = nReadPerCB.empty() ? 0 : *std::max_element(nReadPerCB.begin(), nReadPerCB.end());
    filteredCells.reset(nCB);

    ofstream& metrics = ofstrOpen(outputPrefix + "Pipseq.metrics.tsv", ERROR_OUT, P);
    metrics << std::setprecision(17) << "metric\tvalue\n"
        << "mode\tIllumina\ncount_units\testimated_transcripts\n"
        << "imi_policy\tstrand_aware_aligned_5prime_base\n"
        << "phase_tie_policy\treject\nbarcode_correction\texact_then_one_substitution_per_tier\n"
        << "barcode_multi_match\t" << (pSolo.pipseqMultiMatch == "FirstListed" ? "first_listed" : "reject") << "\n"
        << "cdna_trim\t" << (pSolo.pipseqTrim == "Kit" ? "kit" : "none") << "\n"
        << "calibration_population\tall_assigned_cell_gene_groups_before_cell_filtering\n"
        << "ipm_source\t" << calibration.source << "\nipm_raw\t" << calibration.raw
        << "\nipm_effective\t" << calibration.effective << "\nipm_calibration_groups\t" << calibration.groups
        << "\nipm_floor_applied\t" << (calibration.source == "estimated" && calibration.raw < 1)
        << "\nobserved_barcodes_with_counts\t" << nCB
        << "\nunique_bi_imi\t" << totalFragments << "\nestimated_transcripts\t" << totalTranscripts << '\n';
    const char* barcodeNames[] = {"reads_total", "barcode_exact", "barcode_corrected", "barcode_invalid",
        "barcode_nonmatch", "barcode_phase_tie", "bi_invalid_with_valid_barcode", "barcode_ambiguous"};
    const char* readNames[] = {"reads_rejected_barcode", "reads_rejected_bi", "reads_unmapped",
        "reads_multimapped", "reads_no_unique_gene", "reads_counted"};
    for (size_t i = 0; i < barcodeStats.size(); ++i) metrics << barcodeNames[i] << '\t' << barcodeStats[i] << '\n';
    for (size_t i = 0; i < readStats.size(); ++i) metrics << readNames[i] << '\t' << readStats[i] << '\n';
    metrics.close();
    if (!metrics) fail("could not write Pipseq.metrics.tsv");
    P.inOut->logMain << "[PIP] groups=" << nGroups << " barcodes=" << nCB
        << " IPM=" << calibration.effective << " calibration_groups=" << calibration.groups
        << " estimated_transcripts=" << totalTranscripts << '\n';
    fragments.clear();
    outputResults(false, outputPrefix + "raw/");
    if (pSolo.cellFilter.type[0] != "None") {
        if (nCB) cellFiltering();
        // The shared caller returns without writing when it retains no cells.
        // PIP promises a valid empty filtered matrix in this case as well.
        if (filteredCells.nCells == 0) outputResults(true, outputPrefixFiltered);
    }

    // The ordinary Solo summary reads counters this mode does not fill, so
    // write the same file from the PIP counters. Counts are estimated
    // transcripts; saturation is 1 - estimated transcripts / counted reads,
    // the formula of the published description (molecules over counted reads).
    {
        const string name = SoloFeatureTypes::Names[featureType];
        const double reads = barcodeStats[0] ? double(barcodeStats[0]) : 1.0;
        ofstream& summary = ofstrOpen(outputPrefix + "Summary.csv", ERROR_OUT, P);
        summary << "Number of Reads," << barcodeStats[0] << '\n'
            << "Reads With Valid Barcodes," << (barcodeStats[1] + barcodeStats[2]) / reads << '\n'
            << "Sequencing Saturation," << (readStats[5] ? 1.0 - double(totalTranscripts) / readStats[5] : 0.0) << '\n'
            << "Reads Mapped to Genome: Unique+Multiple," << double(g_statsAll.mappedReadsU + g_statsAll.mappedReadsM) / reads << '\n'
            << "Reads Mapped to Genome: Unique," << double(g_statsAll.mappedReadsU) / reads << '\n'
            << "Reads Mapped to " << name << ": Unique " << name << ',' << readStats[5] / reads << '\n'
            << "Count Units,estimated_transcripts\n"
            << "Unique BI/IMI Fragments," << totalFragments << '\n'
            << "Estimated Transcripts," << totalTranscripts << '\n'
            << "IMIs per Molecule," << calibration.effective << '\n'
            << "Barcodes With Counts," << nCB << '\n';
        if (pSolo.cellFilter.type[0] != "None") {
            summary << "Estimated Number of Cells," << filteredCells.nCells << '\n'
                << "Estimated Transcripts in Cells," << filteredCells.nUMIinCells << '\n'
                << "Median Estimated Transcripts per Cell," << filteredCells.medianUMIperCell << '\n'
                << "Median " << name << " per Cell," << filteredCells.medianGenePerCell << '\n';
        }
        summary.close();
        if (!summary) fail("could not write Summary.csv");
    }
}
