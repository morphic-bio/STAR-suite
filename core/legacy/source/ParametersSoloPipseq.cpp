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
        if (pipseqIPM != 0 || pipseqMultiMatch != "Reject" || pipseqTrim != "Kit" ||
            (!pipseqWhitelist.empty() && pipseqWhitelist != vector<string>{"-"}))
            fail("--soloPipseqWhitelist/IPM/MultiMatch/Trim require --soloPipseqMode Illumina");
        return;
    }
    if (pipseqMode != "Illumina") fail("--soloPipseqMode must be None or Illumina");
    if (pP->runMode != "alignReads" || pP->readFilesTypeN != 1 || pP->readNends != 2 || barcodeReadIn != 0)
        fail("requires alignReads with --readFilesIn R2.fastq R1.fastq and --soloBarcodeMate 0");
    if (pP->outSAMattrPresent.UB || pP->outSAMattrPresent.CY || pP->outSAMattrPresent.UY ||
        pP->outSAMattrPresent.sM || pP->outSAMattrPresent.sQ)
        fail("--outSAMattributes UB, CY, UY, sM and sQ are not defined in PIP mode; CR, CB, UR, sS, GX, GN and sF are");
    if (flexModeStr != "no" || pP->twoPass.yes || pP->runRestart.type != 0 ||
        (!pP->pfMulti.pfMultiConfig.empty() && pP->pfMulti.pfMultiConfig != "-"))
        fail("Flex, two-pass, restart and pfMulti guide integration are not supported with PIP RNA mode");
    if (!std::isfinite(pipseqIPM) || (pipseqIPM != 0 && pipseqIPM < 1))
        fail("--soloPipseqIPM must be 0 (automatic) or a finite value >= 1");
    if (pipseqMultiMatch != "Reject" && pipseqMultiMatch != "FirstListed")
        fail("--soloPipseqMultiMatch must be Reject or FirstListed");
    if (pipseqTrim != "Kit" && pipseqTrim != "None")
        fail("--soloPipseqTrim must be Kit or None");

    // Assay-owned settings bypass ordinary CB/UMI algorithms. An explicit
    // conflicting value is an error, not a silently ignored request.
    const auto parameter = [&](const string& name, const string& wanted, bool owned) {
        assayParameter("PIPseq", "PIP", "PIP", name, wanted, owned);
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
    parameter("soloCellFilter", "EmptyDrops_CR", false);
    if (featureIn.size() != 1 || (featureIn[0] != "Gene" && featureIn[0] != "GeneFull" &&
                                  featureIn[0] != "GeneFull_Ex50pAS" && featureIn[0] != "GeneFull_ExonOverIntron"))
        fail("requires exactly one --soloFeatures Gene, GeneFull, GeneFull_Ex50pAS or GeneFull_ExonOverIntron");
    if (pipseqTrim == "Kit") {
        // cDNA trimming as in the published pipeline description. The clip
        // type is the mode's; its sequences and thresholds can be replaced.
        parameter("clipAdapterType", "PIPseq", true);
        parameter("clip5pNbases", "1", false);
        parameter("clip5pAdapterSeq", pipseq::kitAdapters5p, false);
        parameter("clip5pAdapterMMp", "0.1", false);
        parameter("clip3pAdapterSeq", pipseq::kitAdapters3p, false);
        parameter("clip3pAdapterMMp", "0.1", false);
    }
    // Maximum intron length given for STAR in the published pipeline description.
    parameter("alignIntronMax", "25000", false);
    try {
        pipseqDecoder = std::make_shared<pipseq::Decoder>(pipseqWhitelist,
            pipseqMultiMatch == "FirstListed" ? pipseq::MultiMatch::FirstListed : pipseq::MultiMatch::Reject);
    }
    catch (const std::exception& e) { fail(e.what()); }
    cbWLsize = 0;
    pP->inOut->logMain << "[PIP] Illumina BI/IMI estimator; ordinary UMI correction/filtering disabled; "
        "RNA only; unique genomic alignments; exact tier match first, then one substitution ("
        << (pipseqMultiMatch == "FirstListed" ? "several candidates: first listed" : "several candidates: rejected")
        << "); phase ties rejected; cDNA trimming " << (pipseqTrim == "Kit" ? "by the mode (leading base, 5-prime adapters, poly-A/G)" : "not applied by the mode")
        << "; IMI=strand-aware aligned 5-prime base; sample=one biological library per invocation\n";
}
