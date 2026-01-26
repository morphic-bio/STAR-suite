#include "CbUbTagInjector.h"
#include "SequenceFuns.h"
#include <cstring>
#include <cstdio>

bool CbUbTagInjector::injectTags(
    const char* bamIn,
    uint32_t sizeIn,
    char* bamOut,
    uint32_t& sizeOut,
    uint32_t cbIdxPlus1,
    uint32_t umiPacked,
    bool umiValid,
    const std::vector<std::string>& cbWLstr,
    uint8_t umiL,
    bool needCB,
    bool needUB
) {
    // Default: no injection, output = input
    sizeOut = sizeIn;
    
    // Early return if neither tag is requested
    if (!needCB && !needUB) {
        return false;
    }
    
    // Validate UMI length if UB is requested
    if (needUB && (umiL == 0 || umiL > 16)) {
        // Invalid UMI length - log warning once and skip UB injection
        static bool warningLogged = false;
        if (!warningLogged) {
            fprintf(stderr, "WARNING: CbUbTagInjector: UB tag requested but umiL=%u is invalid (expected 1-16). "
                           "UB tags will not be injected. Check --soloUMIlen configuration.\n", (unsigned)umiL);
            warningLogged = true;
        }
        needUB = false;
    }
    
    // Determine what to emit based on CB/UMI validity
    // cbIdxPlus1 == 0 means no CB match
    // umiValid == false means UMI contains N or is otherwise invalid
    bool hasCB = (cbIdxPlus1 > 0 && cbIdxPlus1 <= cbWLstr.size());
    bool hasUMI = umiValid && (umiL > 0);
    
    // CB and UB are independent - emit each if requested and valid
    bool emitCB = needCB && hasCB;
    bool emitUB = needUB && hasUMI;
    
    if (!emitCB && !emitUB) {
        return false;  // Nothing to inject
    }
    
    // Decode CB string from whitelist
    std::string cbStr;
    if (emitCB) {
        cbStr = cbWLstr[cbIdxPlus1 - 1];  // Convert 1-based to 0-based index
    }
    
    // Decode UMI string
    std::string ubStr;
    if (emitUB) {
        ubStr = convertNuclInt64toString(umiPacked, umiL);
    }
    
    // Validate input BAM record structure
    // BAM format: [block_size:uint32][core:32 bytes][data:variable]
    // block_size = size of core + data (doesn't include the 4-byte size field itself)
    if (sizeIn < 4) {
        return false;  // Invalid BAM record
    }
    
    uint32_t blockSizeFromHeader = *(const uint32_t*)bamIn;
    uint32_t expectedSize = blockSizeFromHeader + 4;
    if (sizeIn != expectedSize) {
        // Size mismatch - don't modify
        return false;
    }
    
    // Calculate tag sizes
    // BAM aux tag format: [tag:2 bytes][type:1 byte][value:variable]
    // For Z (string) type: value is null-terminated string
    uint32_t tagSize = 0;
    if (emitCB) {
        tagSize += 3 + cbStr.size() + 1;  // CB:Z:<cbStr>\0
    }
    if (emitUB) {
        tagSize += 3 + ubStr.size() + 1;  // UB:Z:<ubStr>\0
    }
    
    // Calculate new sizes
    uint32_t newBlockSize = blockSizeFromHeader + tagSize;
    uint32_t newSize = newBlockSize + 4;
    
    // Check size limits
    if (newSize > BAM_ATTR_MaxSize) {
        return false;  // Would exceed max size
    }
    
    // Copy existing record to output buffer
    memcpy(bamOut, bamIn, sizeIn);
    
    // Update block_size in output
    *(uint32_t*)bamOut = newBlockSize;
    
    // Append tags at end of existing data
    char* tagPtr = bamOut + 4 + blockSizeFromHeader;
    
    if (emitCB) {
        tagPtr[0] = 'C';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, cbStr.c_str(), cbStr.size() + 1);
        tagPtr += 3 + cbStr.size() + 1;
    }
    
    if (emitUB) {
        tagPtr[0] = 'U';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, ubStr.c_str(), ubStr.size() + 1);
        tagPtr += 3 + ubStr.size() + 1;
    }
    
    sizeOut = newSize;
    return true;
}

uint32_t CbUbTagInjector::maxTagOverhead(uint32_t cbLen, uint32_t umiLen) {
    // CB tag: 'C' + 'B' + 'Z' + cbLen + '\0' = 3 + cbLen + 1
    // UB tag: 'U' + 'B' + 'Z' + umiLen + '\0' = 3 + umiLen + 1
    return (3 + cbLen + 1) + (3 + umiLen + 1);
}
