#ifndef H_CbUbTagInjector
#define H_CbUbTagInjector

#include "IncludeDefine.h"
#include <string>
#include <vector>

/**
 * CbUbTagInjector: Shared helper for injecting CB:Z: and UB:Z: tags into BAM records.
 * 
 * This provides a common implementation for both:
 * - Sorted BAM output (via SoloFeature::addBAMtags)
 * - Unsorted BAM output (via BAMoutput::unsortedOneAlign)
 * 
 * CB and UB tags are independent - either can be emitted without the other.
 */
class CbUbTagInjector {
public:
    /**
     * Inject CB/UB tags into a BAM record.
     * 
     * @param bamIn        Input BAM record buffer
     * @param sizeIn       Size of input BAM record (block_size + 4)
     * @param bamOut       Output buffer (must be at least BAM_ATTR_MaxSize bytes)
     * @param sizeOut      [out] Size of output BAM record after injection
     * @param cbIdxPlus1   1-based whitelist index (0 = no CB)
     * @param umiPacked    Packed UMI value (24-bit or umiL*2 bits)
     * @param umiValid     Whether UMI is valid (false if contains N or other issues)
     * @param cbWLstr      CB whitelist strings (for decoding cbIdxPlus1)
     * @param umiL         UMI length in bases (1-16)
     * @param needCB       Whether CB tag is requested (P.outSAMattrPresent.CB)
     * @param needUB       Whether UB tag is requested (P.outSAMattrPresent.UB)
     * 
     * @return true if tags were injected (sizeOut updated), false if no injection needed
     */
    static bool injectTags(
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
    );
    
    /**
     * Compute the maximum additional size needed for CB/UB tags.
     * Useful for pre-allocating buffers.
     * 
     * @param cbLen  Maximum CB length (typically 16)
     * @param umiLen Maximum UMI length (typically 12)
     * @return Maximum additional bytes for CB+UB tags
     */
    static uint32_t maxTagOverhead(uint32_t cbLen, uint32_t umiLen);
};

#endif // H_CbUbTagInjector
