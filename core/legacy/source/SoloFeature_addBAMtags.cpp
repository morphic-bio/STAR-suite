#include "SoloFeature.h"
#include "SequenceFuns.h"
#include "BAMfunctions.h"
#include "ErrorWarning.h"
#include <cstring>
#include <cstdlib>

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
// Policy: CB+UB together (requireCbUbTogether)
// - When true (default): UB requires CB, CB lookup failures are fatal if UB is requested.
//   This is an intentional policy to ensure data consistency.
// - When false: Allow UB-only or CB-only injection without hard-failure.
//   CB lookup is only guarded when needCB is true, and UB decode failures are non-fatal.
//
// Policy control: Set environment variable STAR_REQUIRE_CBUB_TOGETHER=no to relax the policy.
// Default: true (strict policy)
//
// Policy exception: status==2 (invalid UMI) is an intentional exception - CB is emitted, UB is skipped.
// This matches recordReadInfo semantics where CB is valid but UMI is invalid.
void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1, uint32_t readId) {
    // Policy flag: require CB+UB together
    // Can be relaxed via environment variable STAR_REQUIRE_CBUB_TOGETHER=no
    // Default: true (strict policy)
    static const bool requireCbUbTogetherDefault = true;
    static const bool requireCbUbTogether = (std::getenv("STAR_REQUIRE_CBUB_TOGETHER") != nullptr) 
        ? (std::string(std::getenv("STAR_REQUIRE_CBUB_TOGETHER")) == "no" ? false : true)
        : requireCbUbTogetherDefault;
    
    // ZI tag emission flag: controlled by STAR_EMIT_READID_TAG env var (default: off)
    static bool ziTagEmitChecked = false;
    static bool ziTagEmitEnabled = false;
    if (!ziTagEmitChecked) {
        const char* envVal = std::getenv("STAR_EMIT_READID_TAG");
        ziTagEmitEnabled = (envVal != nullptr && std::string(envVal) != "0" && std::string(envVal) != "");
        if (ziTagEmitEnabled) {
            P.inOut->logMain << "ZI tag emission enabled (STAR_EMIT_READID_TAG)" << endl;
        }
        ziTagEmitChecked = true;
    }
    
    // Early return if tags not requested (check both CB/UB and ZI)
    if (!pSolo.samAttrYes && !ziTagEmitEnabled) {
        return;
    }
    
    // Check if CB or UB are requested
    bool needCB = P.outSAMattrPresent.CB;
    bool needUB = P.outSAMattrPresent.UB;
    // If neither CB/UB nor ZI are requested, return early
    if (!needCB && !needUB && !ziTagEmitEnabled) {
        return;
    }
    
    // Policy enforcement: if requireCbUbTogether and UB is requested, require CB
    if (requireCbUbTogether && needUB && !needCB) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: UB tag requested but CB tag not requested\n";
        errOut << "SOLUTION: Request both CB and UB tags, or set environment variable STAR_REQUIRE_CBUB_TOGETHER=no";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        return;
    }
    
    // Hard error: UB requested but UMI length invalid
    if (needUB && (pSolo.umiL == 0 || pSolo.umiL > 16)) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: UB tag requested but UMI length is invalid\n";
        errOut << "pSolo.umiL=" << static_cast<int>(pSolo.umiL) << " (must be 1-16 for UB injection)\n";
        errOut << "SOLUTION: Set --soloUMIlen to a value between 1 and 16, or do not request UB tag";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        return;
    }
    
    // Hard error: Check packedReadInfo is initialized
    // Stage 2: Flex-aware check - PackedReadInfo may be empty in Flex modes
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
    
    // readId is passed directly (no shifting needed)
    
    // Hard error: Check readId in range
    if (readId >= packedReadInfo.data.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: readId out of range for PackedReadInfo\n";
        errOut << "readId=" << readId << " packedReadInfo.data.size()=" << packedReadInfo.data.size() << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Lookup from PackedReadInfo
    // Stage 2: Flex-aware - PackedReadInfo is the primary source for both Flex and non-Flex paths
    // Future: If Flex uses alternative sources (e.g., inline hash), add lookup logic here
    uint32_t cbIdx = packedReadInfo.getCB(readId);
    uint32_t umiPacked = packedReadInfo.getUMI(readId);
    uint8_t status = packedReadInfo.getStatus(readId);
    
    // Status-aware tag emission (matches recordReadInfo semantics)
    // status==0: skip both CB and UB (missing CB) - no error unless policy requires CB
    // status==1: emit both CB and UB (valid)
    // status==2: emit CB only, skip UB (invalid UMI)
    //   POLICY EXCEPTION: status==2 is an intentional exception to CB+UB together policy.
    //   Invalid UMI still yields CB tag (matches recordReadInfo semantics where CB is valid but UMI is not).
    //   This allows partial tag injection when UMI correction fails but CB is still valid.
    // Stage 2: Flex path uses same status semantics as non-Flex
    if (status == 0) {
        // Policy: if requireCbUbTogether and UB is requested, missing CB is fatal
        // BUT: Allow ZI-only emission even when CB/UB can't be emitted (for validation)
        if (requireCbUbTogether && needUB && !ziTagEmitEnabled) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: UB tag requested but CB is missing (status==0)\n";
            errOut << "SOLUTION: This is enforced by requireCbUbTogether policy. Set environment variable STAR_REQUIRE_CBUB_TOGETHER=no to allow UB-only injection.";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            return;
        }
        // Even if CB/UB aren't emitted, still emit ZI if enabled (for validation)
        if (!ziTagEmitEnabled) {
            return;  // No tags to add
        }
        // Fall through to emit ZI tag only (skip CB/UB policy check if ZI-only)
    }
    
    // Policy check for status==2: invalid UMI but valid CB
    // EXCEPTION: status==2 is allowed to emit CB-only (intentional exception to CB+UB together)
    // This matches recordReadInfo semantics where CB is valid but UMI is invalid
    if (status == 2 && requireCbUbTogether && needUB) {
        // Note: This is an intentional exception - CB will be emitted, UB will be skipped
        // Documented above in status-aware tag emission comment
    }
    
    // Decode CB with bounds check
    // Policy: if requireCbUbTogether and UB is requested, CB lookup failures are fatal
    // Otherwise, only guard CB lookup when needCB is true
    bool cbLookupRequired = needCB || (requireCbUbTogether && needUB);
    if (status != 0 && cbLookupRequired && cbIdx >= pSolo.cbWLstr.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: cbIdx out of range for whitelist\n";
        errOut << "cbIdx=" << cbIdx << " pSolo.cbWLstr.size()=" << pSolo.cbWLstr.size() << "\n";
        if (requireCbUbTogether && needUB && !needCB) {
            errOut << "NOTE: CB lookup required because UB is requested and requireCbUbTogether=true\n";
        }
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    string cbStr;
    bool emitCB = false;
    if (needCB && status != 0 && cbIdx < pSolo.cbWLstr.size()) {
        cbStr = pSolo.cbWLstr[cbIdx];
        emitCB = true;
    }
    
    // Decode UMI using length-aware decoder (only if status == 1)
    // UMI length validation already checked above (hard error if needUB && umiL invalid)
    string ubStr;
    bool emitUB = false;
    if (status == 1 && needUB) {
        // UMI length is guaranteed valid here (1-16) due to check above
        // Use convertNuclInt64toString from SequenceFuns.h
        // Do NOT skip when umiPacked == 0 (valid all-A UMI)
        ubStr = convertNuclInt64toString(umiPacked, pSolo.umiL);
        emitUB = true;
    }
    // Note: status==2 intentionally skips UB (invalid UMI) but still emits CB (policy exception)
    
    // If no tags to add (including ZI), leave bam0/size0 unchanged
    if (!emitCB && !emitUB && !ziTagEmitEnabled) {
        return;
    }
    
    // Block_size sanity check: verify size0 matches block_size + 4
    uint32_t blockSizeFromHeader = *(uint32_t*)bam0;
    uint32_t expectedSize = blockSizeFromHeader + 4;
    if (size0 != expectedSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: BAM record size mismatch\n";
        errOut << "size0=" << size0 << " expected (block_size+4)=" << expectedSize << " block_size=" << blockSizeFromHeader << "\n";
        errOut << "SOLUTION: This indicates corrupted BAM data - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Manual raw tag append (no htslib bam_aux_append - unsafe on non-owned memory)
    // BAM format: [block_size:uint32][core:32 bytes][data:variable]
    // block_size = size of core + data (doesn't include the 4-byte size field itself)
    // data includes: QNAME, CIGAR, SEQ, QUAL, aux tags
    
    // Read current block_size (already validated above)
    uint32_t blockSize = blockSizeFromHeader;
    
    // Calculate tag sizes
    uint32_t tagSize = 0;
    if (emitCB) {
        tagSize += 3 + cbStr.size() + 1;  // CB:Z:<cbStr>\0
    }
    if (emitUB) {
        tagSize += 3 + ubStr.size() + 1;  // UB:Z:<ubStr>\0
    }
    if (ziTagEmitEnabled) {
        tagSize += 7;  // ZI:i:<readId> = ['Z']['I']['I'][4-byte uint32] = 7 bytes
    }
    
    // Calculate new sizes
    uint32_t newBlockSize = blockSize + tagSize;
    uint32_t newSize = newBlockSize + 4;  // +4 for block_size field
    
    // Hard error on size overflow
    if (newSize > BAM_ATTR_MaxSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: BAM record exceeds maximum size after tag injection\n";
        errOut << "newSize=" << newSize << " BAM_ATTR_MaxSize=" << BAM_ATTR_MaxSize << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Copy existing record to bam1
    memcpy(bam1, bam0, size0);
    
    // Update block_size
    uint32_t* blockSizePtr = reinterpret_cast<uint32_t*>(bam1);
    *blockSizePtr = newBlockSize;
    
    // Append tags at end of existing data (after block_size + core + existing data)
    // Current data ends at: bam0 + 4 + blockSize
    // New data ends at: bam1 + 4 + newBlockSize
    char* tagPtr = bam1 + 4 + blockSize;
    
    if (emitCB) {
        // Append CB tag: CB:Z:<cbStr>\0
        tagPtr[0] = 'C';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, cbStr.c_str(), cbStr.size() + 1);
        tagPtr += 3 + cbStr.size() + 1;
    }
    
    if (emitUB) {
        // Append UB tag: UB:Z:<ubStr>\0
        tagPtr[0] = 'U';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, ubStr.c_str(), ubStr.size() + 1);
        tagPtr += 3 + ubStr.size() + 1;
    }
    
    if (ziTagEmitEnabled) {
        // Append ZI tag: ZI:i:<readId> (uint32 little-endian)
        tagPtr[0] = 'Z';
        tagPtr[1] = 'I';
        tagPtr[2] = 'I';  // type 'I' = uint32
        *reinterpret_cast<uint32_t*>(tagPtr + 3) = readId;  // little-endian uint32
        tagPtr += 7;
    }
    
    // Update bam0/size0 to point to new record
    bam0 = bam1;
    size0 = newSize;
}

