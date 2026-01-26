#include "SoloFeature.h"
#include "CbUbTagInjector.h"
#include "SequenceFuns.h"
#include "BAMfunctions.h"
#include "ErrorWarning.h"
#include <cstring>

// Legacy overload: extracts iRead from trailing 8 bytes (legacy bin sorter convention with Y-bit encoding)
void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1) {
    // Extract iReadWithY from trailing 8 bytes (legacy convention: bit63=Y-flag)
    uint64_t iReadWithY = *(uint64_t*)(bam0 + size0);
    
    // Mask Y-bit and extract readId (bits[63:32])
    uint64_t iRead = (iReadWithY & ~(1ULL << 63));
    uint32_t readId = static_cast<uint32_t>(iRead >> 32);
    
    addBAMtags(bam0, size0, bam1, readId);
}

// Main implementation: explicit readId for samtools sorter (uint32 readId, no shifting needed)
//
// This function uses CbUbTagInjector for CB/UB tag injection.
// PackedReadInfo is used to look up CB/UMI values by readId.
// CB and UB tags are independent - either can be emitted without the other.
void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1, uint32_t readId) {
    // Early return if tags not requested
    if (!pSolo.samAttrYes) {
        return;
    }
    
    // Check if CB or UB are requested
    bool needCB = P.outSAMattrPresent.CB;
    bool needUB = P.outSAMattrPresent.UB;
    if (!needCB && !needUB) {
        return;
    }
    
    // Check packedReadInfo is initialized
    bool isFlexMode = (pSolo.inlineHashMode || pSolo.inlineCBCorrection);
    bool isMinimalMemory = (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode);
    
    if (packedReadInfo.data.empty()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: PackedReadInfo not initialized for CB/UB injection\n";
        if (isFlexMode || isMinimalMemory) {
            errOut << "NOTE: Flex mode detected (inlineHashMode=" << (pSolo.inlineHashMode ? "true" : "false");
            if (pSolo.inlineCBCorrection) errOut << ", inlineCBCorrection=true";
            if (isMinimalMemory) errOut << ", soloFlexMinimalMemory=true";
            errOut << ")\n";
            errOut << "SOLUTION: PackedReadInfo must be populated for tag injection. ";
            if (isMinimalMemory) {
                errOut << "Disable --soloFlexMinimalMemory or ensure PackedReadInfo is populated for tag injection.";
            } else if (pSolo.inlineCBCorrection) {
                errOut << "Inline CB correction mode requires PackedReadInfo for tag injection.";
            } else {
                errOut << "Ensure PackedReadInfo is initialized before tag injection.";
            }
        } else {
            errOut << "SOLUTION: This should not happen - please report this error";
        }
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Check readId in range
    if (readId >= packedReadInfo.data.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: readId out of range for PackedReadInfo\n";
        errOut << "readId=" << readId << " packedReadInfo.data.size()=" << packedReadInfo.data.size() << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Lookup from PackedReadInfo
    uint32_t cbIdx = packedReadInfo.getCB(readId);
    uint32_t umiPacked = packedReadInfo.getUMI(readId);
    uint8_t status = packedReadInfo.getStatus(readId);
    
    // Convert PackedReadInfo data to CbUbTagInjector format:
    // - cbIdxPlus1: 1-based (add 1 to 0-based cbIdx, or 0 if status==0)
    // - umiValid: true only if status==1
    uint32_t cbIdxPlus1 = (status != 0 && cbIdx < pSolo.cbWLstr.size()) ? (cbIdx + 1) : 0;
    bool umiValid = (status == 1);
    
    // CB bounds check (for error reporting)
    if (needCB && status != 0 && cbIdx >= pSolo.cbWLstr.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: cbIdx out of range for whitelist\n";
        errOut << "cbIdx=" << cbIdx << " pSolo.cbWLstr.size()=" << pSolo.cbWLstr.size() << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Use CbUbTagInjector for CB/UB tag injection
    uint32_t sizeOut = 0;
    bool injected = CbUbTagInjector::injectTags(
        bam0, size0, bam1, sizeOut,
        cbIdxPlus1, umiPacked, umiValid,
        pSolo.cbWLstr, pSolo.umiL,
        needCB, needUB
    );
    
    if (injected) {
        // CB/UB tags were injected - update pointers
        bam0 = bam1;
        size0 = sizeOut;
    }
}
