#include "SoloFeature.h"
#include "SequenceFuns.h"
#include <algorithm>

void SoloFeature::resolveCbUb(uint32_t iread, CbUbResult& result) const
{
    // Get packed values from readInfo
    result.cbIdx = getPackedCB(iread);
    result.umiPacked = getPackedUMI(iread);
    result.status = getPackedStatus(iread);

    // Initialize strings to "-" (missing/invalid marker)
    result.cb = "-";
    result.ub = "-";

    // Resolve CB string if status is valid and index is in range
    if (result.status == 1 && result.cbIdx < pSolo.cbWLstr.size()) {
        result.cb = pSolo.cbWLstr[result.cbIdx];
    }

    // Resolve UMI string if status is valid
    if (result.status == 1) {
        result.ub.clear();
        result.ub.reserve(pSolo.umiL);
        uint32_t tmp = result.umiPacked;
        
        // Decode UMI from packed representation (LSB-first, 2 bits per base)
        for (int i = pSolo.umiL - 1; i >= 0; --i) {
            uint8_t b = tmp & 0x3;
            tmp >>= 2;
            char c = 'N';
            switch (b) {
                case 0: c = 'A'; break;
                case 1: c = 'C'; break;
                case 2: c = 'G'; break;
                case 3: c = 'T'; break;
            }
            result.ub.push_back(c);
        }
        std::reverse(result.ub.begin(), result.ub.end());
    }

    // Debug logging removed temporarily - was causing segfaults in multithreaded code
    // TODO: Re-enable with proper thread-safe logging
}

