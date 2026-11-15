#include "Parameters.h"
#include "SplitSeq.h"
#include "ErrorWarning.h"
#include <stdexcept>

void ParametersSolo::initializeSplitseq() {
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of split-seq parameter error: " + message + "\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    if (splitseqMode == "None") {
        if (splitseqWhitelist != vector<string>{"-"} || splitseqCBposition != vector<string>{"-"} || splitseqUMIposition != "-")
            fail("--soloSplitseqWhitelist/CBposition/UMIposition require --soloSplitseqMode Fixed");
        return;
    }
    if (splitseqMode != "Fixed") fail("--soloSplitseqMode must be None or Fixed");
    if (pipseqMode != "None") fail("split-seq and PIPseq modes are mutually exclusive");
    if (pP->runMode != "alignReads" || pP->readFilesTypeN != 1 || pP->readNends != 2 || barcodeReadIn != 0)
        fail("requires alignReads, two FASTQ reads (cDNA then barcode/UMI), and --soloBarcodeMate 0");
    if (pP->outSAMtype != vector<string>{"None"})
        fail("initial split-seq RNA implementation requires --outSAMtype None");
    if (flexModeStr != "no" || pP->twoPass.yes || pP->runRestart.type != 0 ||
        (!pP->pfMulti.pfMultiConfig.empty() && pP->pfMulti.pfMultiConfig != "-"))
        fail("Flex, two-pass, restart and pfMulti guide integration are not supported with split-seq RNA mode");

    try { splitseqDecoder = std::make_shared<splitseq::Decoder>(splitseqCBposition, splitseqUMIposition, splitseqWhitelist); }
    catch (const std::exception& e) { fail(e.what()); }

    const auto parameter = [&](const string& name, const string& wanted, bool owned) {
        for (auto* p : pP->parArray) {
            if (p->nameString != name) continue;
            if (p->inputLevel > 0) {
                std::ostringstream out; out << *p;
                std::istringstream values(out.str());
                string token, actual;
                while (values >> token) { if (!actual.empty()) actual += ' '; actual += token; }
                if (owned && actual != wanted) fail("--" + name + " conflicts with split-seq mode (effective value: " + wanted + ")");
            } else {
                std::istringstream in(wanted);
                p->inputValues(in);
            }
            pP->inOut->logMain << "[SPLIT effective] " << name << " = " << *p << '\n';
            return;
        }
        fail("internal parameter registration missing: " + name);
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
    parameter("soloCrMultimapRescue", "no", true);
    parameter("soloCellReadStats", "None", true);
    parameter("soloUMICorrection", "none", true); // Keep standard soloUMIdedup; no optional postprocessor.
    parameter("soloMultiMappers", "Unique", false);
    parameter("soloFeatures", "GeneFull", false);
    parameter("soloCellFilter", "EmptyDrops_CR", false);
    if (featureIn != vector<string>{"Gene"} && featureIn != vector<string>{"GeneFull"})
        fail("requires exactly one --soloFeatures Gene or GeneFull");
    cbWLsize = 0;
    for (size_t i = 0; i < splitseqDecoder->regions().size(); ++i)
        pP->inOut->logMain << "[SPLIT] round " << i + 1 << " entries=" << splitseqDecoder->tierSize(i) << '\n';
    pP->inOut->logMain << "[SPLIT] fixed layout; per-round exact then unique Hamming-1 correction; "
        "replacement applied once; observed tuples only; standard UMI counting; one sublibrary per invocation\n";
}
