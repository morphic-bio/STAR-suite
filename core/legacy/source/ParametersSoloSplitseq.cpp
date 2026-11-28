#include "Parameters.h"
#include "SplitSeq.h"
#include "ErrorWarning.h"
#include <cmath>
#include <sstream>
#include <stdexcept>

void ParametersSolo::initializeSplitseq() {
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of split-seq parameter error: " + message + "\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    const bool edit2 = splitseqMode == "FixedEdit2";
    if (!edit2) {
        for (const auto* param : pP->parArray)
            if (param->inputLevel > 0 && (param->nameString == "soloSplitseqMinPosterior" ||
                                          param->nameString == "soloSplitseqFrequencyBackup"))
                fail("--" + param->nameString + " requires --soloSplitseqMode FixedEdit2");
    }
    if (splitseqMode == "None") {
        if (splitseqWhitelist != vector<string>{"-"} || splitseqCBposition != vector<string>{"-"} || splitseqUMIposition != "-")
            fail("--soloSplitseqWhitelist/CBposition/UMIposition require --soloSplitseqMode Fixed or FixedEdit2");
        return;
    }
    if (splitseqMode != "Fixed" && !edit2) fail("--soloSplitseqMode must be None, Fixed or FixedEdit2");
    if (edit2) {
        if (splitseqFrequencyBackup != "yes" && splitseqFrequencyBackup != "no")
            fail("--soloSplitseqFrequencyBackup must be yes or no");
        pP->inOut->logMain << "[SPLIT effective] soloSplitseqFrequencyBackup = " << splitseqFrequencyBackup << '\n';
        // Validate the whole token. A malformed cutoff interpreted as zero would admit calls.
        std::istringstream input(splitseqMinPosteriorInput);
        if (!(input >> splitseqMinPosterior) || !input.eof() || !std::isfinite(splitseqMinPosterior) ||
            splitseqMinPosterior < 0.0 || splitseqMinPosterior > 1.0)
            fail("--soloSplitseqMinPosterior must be a probability from 0 to 1");
        pP->inOut->logMain << "[SPLIT effective] soloSplitseqMinPosterior = " << splitseqMinPosteriorInput << '\n';
    }
    if (pipseqMode != "None") fail("split-seq and PIPseq modes are mutually exclusive");
    if (pP->runMode != "alignReads" || pP->readFilesTypeN != 1 || pP->readNends != 2 || barcodeReadIn != 0)
        fail("requires alignReads, two FASTQ reads (cDNA then barcode/UMI), and --soloBarcodeMate 0");
    if (pP->outSAMtype != vector<string>{"None"})
        fail("initial split-seq RNA implementation requires --outSAMtype None");
    if (flexModeStr != "no" || pP->twoPass.yes || pP->runRestart.type != 0 ||
        (!pP->pfMulti.pfMultiConfig.empty() && pP->pfMulti.pfMultiConfig != "-"))
        fail("Flex, two-pass, restart and pfMulti guide integration are not supported with split-seq RNA mode");

    try { splitseqDecoder = std::make_shared<splitseq::Decoder>(splitseqCBposition, splitseqUMIposition, splitseqWhitelist,
                    edit2 ? splitseq::Policy::Edit2 : splitseq::Policy::Hamming1); }
    catch (const std::exception& e) { fail(e.what()); }

    const auto parameter = [&](const string& name, const string& wanted, bool owned) {
        assayParameter("split-seq", "split-seq", "SPLIT", name, wanted, owned);
    };
    // Feed the existing sparse no-global-whitelist spool with canonical packed
    // tuples. Per-round references are enforced by the decoder, before counting.
    parameter("soloType", "CB_UMI_Simple", true);
    parameter("soloCBwhitelist", "None", true);
    parameter("soloCBtype", "Sequence", true);
    parameter("soloCBstart", "1", true);
    parameter("soloCBlen", std::to_string(splitseqDecoder->cellLength()), true);
    parameter("soloUMIstart", "1", true);
    parameter("soloUMIlen", std::to_string(splitseqDecoder->umiRegion().length), true);
    parameter("soloBarcodeReadLength", "0", true);
    parameter("soloCBposition", "-", true);
    parameter("soloUMIposition", "-", true);
    parameter("soloCBmatchWLtype", "Exact", true); // Resolved by the per-round decoder.
    parameter("soloInlineHashMode", "no", true);
    parameter("soloInlineCBCorrection", "no", true);
    parameter("soloSkipProcessing", "no", true);
    parameter("soloWriteKeysBin", "no", true);
    parameter("soloCellReadStats", "None", true);
    parameter("soloUMICorrection", "none", true); // Keep standard soloUMIdedup; no optional postprocessor.
    parameter("soloMultiMappers", "Unique", false);
    parameter("soloFeatures", "GeneFull", false);
    parameter("soloCellFilter", "EmptyDrops_CR", false);
    if (featureIn.size() != 1 || (featureIn[0] != "Gene" && featureIn[0] != "GeneFull" &&
                                  featureIn[0] != "GeneFull_Ex50pAS" && featureIn[0] != "GeneFull_ExonOverIntron"))
        fail("requires exactly one --soloFeatures Gene, GeneFull, GeneFull_Ex50pAS or GeneFull_ExonOverIntron");
    // The shared GEX resolver uses the GEX feature to select automatic intronic fallback.
    parameter("soloCrGexFeature", featureIn.front() == "Gene" ? "Gene" : "GeneFull", false);
    parameter("soloCrMultimapRescue", "yes", false);
    cbWLsize = 0;
    for (size_t i = 0; i < splitseqDecoder->regions().size(); ++i)
        pP->inOut->logMain << "[SPLIT] round " << i + 1 << " entries=" << splitseqDecoder->tierSize(i) << '\n';
    if (edit2) {
        const auto path = pP->outFileNamePrefix + "Splitseq.ties.tsv";
        try { splitseqTieWriter = std::make_shared<splitseq::TieWriter>(path, pP->outFileTmp + "Splitseq.ties.pending"); }
        catch (const std::ios_base::failure& e) { fail("cannot open tie trace " + path + ": " + e.what()); }
    }
    pP->inOut->logMain << "[SPLIT] fixed layout; "
        << (edit2 ? "minimum edit distance, <=2 per round and <=6 total; Bayesian ties from full-input exact oligo counts; "
                  : "per-round exact then unique Hamming-1 correction; ")
        << "replacement applied once; observed tuples only; standard UMI counting; one sublibrary per invocation\n";
}
