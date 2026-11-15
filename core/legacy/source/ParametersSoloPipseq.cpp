#include "Parameters.h"
#include "PipSeq.h"
#include "ErrorWarning.h"
#include <cmath>
#include <stdexcept>

void ParametersSolo::initializePipseq() {
    const auto fail = [&](const string& message) {
        exitWithError("EXITING because of PIPseq parameter error: " + message + "\n",
                      std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    if (pipseqMode == "None") {
        if (pipseqIPM != 0 || pipseqMultiMatch != "Reject" || (!pipseqWhitelist.empty() && pipseqWhitelist != vector<string>{"-"}))
            fail("--soloPipseqWhitelist/IPM/MultiMatch require --soloPipseqMode Illumina");
        return;
    }
    if (pipseqMode != "Illumina") fail("--soloPipseqMode must be None or Illumina");
    if (pP->runMode != "alignReads" || pP->readFilesTypeN != 1 || pP->readNends != 2 || barcodeReadIn != 0)
        fail("requires alignReads with --readFilesIn R2.fastq R1.fastq and --soloBarcodeMate 0");
    if (pP->outSAMtype != vector<string>{"None"})
        fail("initial RNA implementation requires --outSAMtype None; PIP BAM tags are not implemented");
    if (flexModeStr != "no" || pP->twoPass.yes || pP->runRestart.type != 0 ||
        (!pP->pfMulti.pfMultiConfig.empty() && pP->pfMulti.pfMultiConfig != "-"))
        fail("Flex, two-pass, restart and pfMulti guide integration are not supported with PIP RNA mode");
    if (!std::isfinite(pipseqIPM) || (pipseqIPM != 0 && pipseqIPM < 1))
        fail("--soloPipseqIPM must be 0 (automatic) or a finite value >= 1");
    if (pipseqMultiMatch != "Reject" && pipseqMultiMatch != "FirstListed")
        fail("--soloPipseqMultiMatch must be Reject or FirstListed");

    // Assay-owned settings bypass ordinary CB/UMI algorithms. An explicit
    // conflicting value is an error, not a silently ignored request.
    const auto parameter = [&](const string& name, const string& wanted, bool owned) {
        for (auto* p : pP->parArray) {
            if (p->nameString != name) continue;
            if (p->inputLevel > 0) {
                std::ostringstream out; out << *p;
                std::istringstream values(out.str());
                string token, actual;
                while (values >> token) { if (!actual.empty()) actual += ' '; actual += token; }
                if (owned && actual != wanted) fail("--" + name + " conflicts with PIP mode (effective value: " + wanted + ")");
            } else {
                std::istringstream in(wanted);
                p->inputValues(in);
            }
            pP->inOut->logMain << "[PIP effective] " << name << " = " << *p << '\n';
            return;
        }
        fail("internal parameter registration missing: " + name);
    };
    parameter("soloType", "CB_UMI_Simple", true);
    parameter("soloCBwhitelist", "None", true);
    parameter("soloCBtype", "Sequence", true);
    parameter("soloCBstart", "1", true);
    parameter("soloCBlen", "28", true);
    parameter("soloUMIstart", "40", true);
    parameter("soloUMIlen", "3", true);
    parameter("soloBarcodeReadLength", "0", true);
    parameter("soloCBmatchWLtype", "Exact", true);
    parameter("soloUMIdedup", "NoDedup", true); // Output column only; PIP estimator owns counting.
    parameter("soloUMIfiltering", "-", true);
    parameter("soloUMICorrection", "none", true);
    parameter("soloInlineHashMode", "no", true);
    parameter("soloInlineCBCorrection", "no", true);
    parameter("soloSkipProcessing", "no", true);
    parameter("soloWriteKeysBin", "no", true);
    parameter("soloCrMultimapRescue", "no", true);
    parameter("soloMultiMappers", "Unique", true);
    parameter("soloCellReadStats", "None", true);
    parameter("soloFeatures", "GeneFull", false);
    parameter("soloCellFilter", "None", false);
    if (featureIn != vector<string>{"Gene"} && featureIn != vector<string>{"GeneFull"})
        fail("requires exactly one --soloFeatures Gene or GeneFull");
    try {
        pipseqDecoder = std::make_shared<pipseq::Decoder>(pipseqWhitelist,
            pipseqMultiMatch == "FirstListed" ? pipseq::MultiMatch::FirstListed : pipseq::MultiMatch::Reject);
    }
    catch (const std::exception& e) { fail(e.what()); }
    cbWLsize = 0;
    pP->inOut->logMain << "[PIP] Illumina BI/IMI estimator; ordinary UMI correction/filtering disabled; "
        "RNA only; unique genomic alignments; exact tier match first, then one substitution ("
        << (pipseqMultiMatch == "FirstListed" ? "several candidates: first listed" : "several candidates: rejected")
        << "); phase ties rejected; "
        "IMI=strand-aware aligned 5-prime base; sample=one biological library per invocation\n";
}
