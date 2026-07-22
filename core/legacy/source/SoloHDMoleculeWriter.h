#ifndef H_SoloHDMoleculeWriter
#define H_SoloHDMoleculeWriter

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class ParametersSolo;

struct SoloHDMoleculeRecord {
    uint32_t row2 = 0;
    uint32_t col2 = 0;
    uint32_t featureIdx = 0;
    uint16_t umiLength = 0;
    uint16_t flags = 0;
    uint64_t umiLow = 0;
    uint64_t umiHigh = 0;
};

class SoloHDBarcodeMap {
public:
    bool initialize(const ParametersSolo &pSolo, const std::string &coordPath, std::string &error);
    bool coordForWL(uint32_t wlIdx, uint32_t &row2, uint32_t &col2) const;
    size_t validCount() const;

    static bool parseCoordinateBarcode(const std::string &barcode, uint32_t &row2, uint32_t &col2);

private:
    std::vector<uint32_t> row2_;
    std::vector<uint32_t> col2_;
    std::vector<uint8_t> valid_;
};

class SoloHDMoleculeWriter {
public:
    static const uint16_t kHeaderSize = 64;
    static const uint16_t kRecordSize = 32;
    static const uint16_t kVersion = 1;
    static const uint32_t kUmiEncoding4Bit = 1;
    static const uint16_t kFlagPrededuplicated = 1;

    SoloHDMoleculeWriter();
    ~SoloHDMoleculeWriter();

    bool open(const std::string &path, uint64_t featureCount, uint16_t flags, std::string &error);
    bool writeRecord(const SoloHDMoleculeRecord &record, std::string &error);
    bool writeRecords(const std::vector<SoloHDMoleculeRecord> &records, std::string &error);
    bool close(std::string &error);
    uint64_t recordCount() const;

    static bool packStarUmi(uint32_t packedUmi, uint16_t umiLength, uint64_t &umiLow, uint64_t &umiHigh);

private:
    std::ofstream out_;
    uint64_t featureCount_;
    uint16_t flags_;
    uint64_t recordCount_;
    bool open_;
};

#endif
