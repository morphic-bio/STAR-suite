#include "SoloFeature.h"
#include "SequenceFuns.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include <cstdlib>
#include <fstream>

void SoloFeature::writeReadIdTagTable() {
    // Check STAR_EMIT_CBUB_TABLE env var (static, cached once)
    static bool tableEmitChecked = false;
    static bool tableEmitEnabled = false;
    if (!tableEmitChecked) {
        const char* envVal = std::getenv("STAR_EMIT_CBUB_TABLE");
        tableEmitEnabled = (envVal != nullptr && std::string(envVal) != "0" && std::string(envVal) != "");
        if (tableEmitEnabled) {
            P.inOut->logMain << "CB/UB table export enabled (STAR_EMIT_CBUB_TABLE)" << endl;
        }
        tableEmitChecked = true;
    }
    
    if (!tableEmitEnabled) {
        return;  // Table export disabled
    }
    
    // Hard error if packedReadInfo is empty (Flex modes may not populate it)
    // NOTE: --flex yes auto-sets --soloFlexMinimalMemory yes and --soloInlineCBCorrection yes,
    //       which leave packedReadInfo empty. This is a validation/testing feature that requires
    //       packedReadInfo, so it's incompatible with production Flex defaults.
    if (packedReadInfo.data.empty()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: STAR_EMIT_CBUB_TABLE enabled but packedReadInfo is empty\n";
        errOut << "NOTE: This occurs when:\n";
        errOut << "  - --soloFlexMinimalMemory yes (default with --flex yes)\n";
        errOut << "  - --soloInlineCBCorrection yes (default with --flex yes)\n";
        errOut << "  - --soloInlineHashMode yes (inline hash mode)\n";
        errOut << "SOLUTION: For validation/testing, override Flex defaults:\n";
        errOut << "  --flex yes --soloFlexMinimalMemory no --soloInlineCBCorrection no\n";
        errOut << "  OR disable STAR_EMIT_CBUB_TABLE (this is a testing feature, not for production Flex runs)";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Output path: alongside BAM, not in Solo.out/Gene/ (that dir may not exist in inline/flex mode)
    string tablePath = P.outFileNamePrefix + "Aligned.out.cb_ub.readId.tsv";
    
    ofstream &tableOut = ofstrOpen(tablePath, ERROR_OUT, P);
    
    // Write header
    tableOut << "readId\tcb\tub\tstatus\n";
    
    // Iterate packedReadInfo and decode matching addBAMtags logic exactly
    for (uint32_t readId = 0; readId < packedReadInfo.data.size(); ++readId) {
        uint8_t status = packedReadInfo.getStatus(readId);
        uint32_t cbIdx = packedReadInfo.getCB(readId);
        uint32_t umiPacked = packedReadInfo.getUMI(readId);
        
        string cbStr = "-";
        string ubStr = "-";
        
        // Decode matching addBAMtags logic exactly:
        // status==0: cb="-", ub="-"
        // status==1: cb=decoded, ub=decoded
        // status==2: cb=decoded, ub="-"
        if (status == 0) {
            // Both absent
            cbStr = "-";
            ubStr = "-";
        } else if (status == 1) {
            // Both present
            if (cbIdx < pSolo.cbWLstr.size()) {
                cbStr = pSolo.cbWLstr[cbIdx];
            } else {
                cbStr = "-";  // Out of range, treat as absent
            }
            // Decode UMI using length-aware decoder
            if (pSolo.umiL > 0 && pSolo.umiL <= 16) {
                ubStr = convertNuclInt64toString(umiPacked, pSolo.umiL);
            } else {
                ubStr = "-";  // Invalid UMI length
            }
        } else if (status == 2) {
            // CB present, UB absent (invalid UMI)
            if (cbIdx < pSolo.cbWLstr.size()) {
                cbStr = pSolo.cbWLstr[cbIdx];
            } else {
                cbStr = "-";  // Out of range
            }
            ubStr = "-";
        }
        
        // Write TSV row
        tableOut << readId << "\t" << cbStr << "\t" << ubStr << "\t" << static_cast<int>(status) << "\n";
    }
    
    tableOut.close();
    P.inOut->logMain << "Wrote CB/UB readId table: " << tablePath << " (" << packedReadInfo.data.size() << " reads)" << endl;
}

