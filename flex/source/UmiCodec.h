#ifndef CODE_UmiCodec
#define CODE_UmiCodec

// Lightweight UMI encode/decode helpers for inline flex path
// Extracted from SoloTagLedger to avoid the full ledger dependency

#include <string>
#include <cstdint>
#include <climits>

// Encode a 12bp UMI string to a 24-bit packed value (2 bits per base)
// Returns UINT32_MAX on invalid input
inline uint32_t encodeUMI12(const std::string& umi) {
    if (umi.length() != 12) {
        return UINT32_MAX;
    }
    
    uint32_t packed = 0;
    for (int i = 0; i < 12; i++) {
        char base = umi[i];
        uint32_t code = 0;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else return UINT32_MAX; // Invalid base
        
        packed |= (code << (i * 2));
    }
    return packed;
}

// Decode a 24-bit packed UMI value back to a 12bp string
inline std::string decodeUMI12(uint32_t packedUR) {
    const char bases[] = "ACGT";
    std::string umi(12, 'N');
    for (int i = 0; i < 12; i++) {
        uint32_t code = (packedUR >> (i * 2)) & 3;
        umi[i] = bases[code];
    }
    return umi;
}

#endif

