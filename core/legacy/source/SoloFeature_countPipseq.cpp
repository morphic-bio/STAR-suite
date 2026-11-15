#include "SoloFeature.h"
#include "PipSeq.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include <iomanip>
#include <limits>

void SoloFeature::countPipseq() {
    pipseq::Groups groups;
    std::array<uint64_t, 8> barcodeStats{};
    std::array<uint64_t, 6> readStats{};
    for (int t = 0; t < P.runThreadN; ++t) {
        auto* sr = RAchunk[t]->RA->soloRead;
        auto* rf = sr->readFeat[pSolo.featureInd[featureType]];
        for (auto& entry : rf->pipGroups) groups[entry.first].merge(entry.second);
        rf->pipGroups.clear();
        for (size_t i = 0; i < barcodeStats.size(); ++i) barcodeStats[i] += sr->readBar->pipStats[i];
        for (size_t i = 0; i < readStats.size(); ++i) readStats[i] += rf->pipReadStats[i];
    }
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of PIPseq counting error: " + message + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    const auto checked32 = [&](uint64_t n) -> uint32_t {
        if (n > std::numeric_limits<uint32_t>::max()) fail("count or matrix offset exceeds STAR's 32-bit matrix capacity");
        return static_cast<uint32_t>(n);
    };
    pipseq::Calibration calibration;
    try { calibration = pipseq::calibrate(groups, pSolo.pipseqIPM); }
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
    uint64_t totalTranscripts = 0, totalFragments = 0;
    for (const auto& entry : groups) {
        const uint64_t cell = entry.first.first;
        const uint32_t gene = entry.first.second;
        const pipseq::Group& group = entry.second;
        if (gene >= Trans.nGe) fail("gene index exceeds loaded annotation");
        if (cell != previous) {
            previous = cell;
            const uint32_t i = checked32(pSolo.cbWLstr.size());
            pSolo.cbWLstr.push_back(pipseq::Decoder::cellSequence(cell));
            indCB.push_back(i); indCBwl.push_back(i);
            countCellGeneUMIindex.push_back(checked32(countCellGeneUMI.size()));
            nUMIperCB.push_back(0); nGenePerCB.push_back(0);
            nReadPerCB.push_back(0); nReadPerCBunique.push_back(0); nReadPerCBtotal.push_back(0);
        }
        const uint32_t count = checked32(pipseq::transcripts(group, calibration.effective));
        countCellGeneUMI.push_back(gene); countCellGeneUMI.push_back(count);
        nUMIperCB.back() = checked32(uint64_t(nUMIperCB.back()) + count);
        ++nGenePerCB.back();
        nReadPerCB.back() = checked32(uint64_t(nReadPerCB.back()) + group.reads);
        nReadPerCBunique.back() = nReadPerCBtotal.back() = nReadPerCB.back();
        totalTranscripts += count; totalFragments += group.fragments.size();
        ledger << pSolo.cbWLstr.back() << '\t' << Trans.geID[gene] << '\t' << group.reads
               << '\t' << group.occupiedBins() << '\t' << group.fragments.size() << '\t' << count << '\n';
    }
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
    P.inOut->logMain << "[PIP] groups=" << groups.size() << " barcodes=" << nCB
        << " IPM=" << calibration.effective << " calibration_groups=" << calibration.groups
        << " estimated_transcripts=" << totalTranscripts << '\n';
    groups.clear();
    outputResults(false, outputPrefix + "raw/");
    if (pSolo.cellFilter.type[0] != "None") {
        if (nCB) cellFiltering();
        // The shared caller returns without writing when it retains no cells.
        // PIP promises a valid empty filtered matrix in this case as well.
        if (filteredCells.nCells == 0) outputResults(true, outputPrefixFiltered);
    }

    // The ordinary Solo summary reads counters this mode does not fill, so
    // write the same file from the PIP counters. Counts are estimated
    // transcripts; saturation is 1 - distinct BI/IMI fragments / counted reads.
    {
        const string name = SoloFeatureTypes::Names[featureType];
        const double reads = barcodeStats[0] ? double(barcodeStats[0]) : 1.0;
        ofstream& summary = ofstrOpen(outputPrefix + "Summary.csv", ERROR_OUT, P);
        summary << "Number of Reads," << barcodeStats[0] << '\n'
            << "Reads With Valid Barcodes," << (barcodeStats[1] + barcodeStats[2]) / reads << '\n'
            << "Sequencing Saturation," << (readStats[5] ? 1.0 - double(totalFragments) / readStats[5] : 0.0) << '\n'
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
