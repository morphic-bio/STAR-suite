#ifndef CODE_PackedReadInfo
#define CODE_PackedReadInfo
#include <cstdint>
#include <vector>
#include <stdexcept>

struct PackedReadInfoLayout {
    uint8_t umiBits;      // 2 * umiLength
    uint8_t cbBits;       // ceil(log2(whitelistSize+1))
    uint8_t statusBits;   // fixed 3
    uint8_t reserved;     // padding
    uint8_t umiLength;    // original length
    uint32_t whitelistSize;
};

class PackedReadInfo {
public:
    PackedReadInfoLayout L{};
    std::vector<uint64_t> data; // one 64-bit word per read
    void init(uint32_t nReads, uint32_t wlSize, uint8_t umiLength);
    void set(uint32_t readId, uint32_t cbIdx, uint32_t umiPacked, uint8_t status);
    uint32_t getCB(uint32_t readId) const;
    uint32_t getUMI(uint32_t readId) const;
    uint8_t  getStatus(uint32_t readId) const;
    void setStatus(uint32_t readId, uint8_t status);
    uint64_t pack(uint32_t cbIdx, uint32_t umiPacked, uint8_t status) const;
    void unpack(uint64_t w, uint32_t &cbIdx, uint32_t &umiPacked, uint8_t &status) const;
private:
    uint64_t maskUMI() const;
    uint64_t maskCB() const;
};

uint8_t log2ceil_u32(uint32_t n);
#endif
