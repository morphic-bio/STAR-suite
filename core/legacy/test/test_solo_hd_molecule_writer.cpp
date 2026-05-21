#include "ParametersSolo.h"
#include "SoloHDMoleculeWriter.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace {

uint16_t readLe16(const std::vector<unsigned char> &bytes, size_t offset)
{
    return static_cast<uint16_t>(bytes[offset]) |
           static_cast<uint16_t>(bytes[offset + 1]) << 8;
}

uint32_t readLe32(const std::vector<unsigned char> &bytes, size_t offset)
{
    return static_cast<uint32_t>(bytes[offset]) |
           static_cast<uint32_t>(bytes[offset + 1]) << 8 |
           static_cast<uint32_t>(bytes[offset + 2]) << 16 |
           static_cast<uint32_t>(bytes[offset + 3]) << 24;
}

uint64_t readLe64(const std::vector<unsigned char> &bytes, size_t offset)
{
    uint64_t value = 0;
    for (int ii = 0; ii < 8; ++ii) {
        value |= static_cast<uint64_t>(bytes[offset + ii]) << (8 * ii);
    }
    return value;
}

std::vector<unsigned char> readFile(const std::string &path)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    assert(in.good());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(in)),
                                      std::istreambuf_iterator<char>());
}

} // namespace

int main()
{
    uint32_t row = 0;
    uint32_t col = 0;
    assert(SoloHDBarcodeMap::parseCoordinateBarcode("s_002um_00123_00456-1", row, col));
    assert(row == 123);
    assert(col == 456);
    assert(SoloHDBarcodeMap::parseCoordinateBarcode("s_002um_7_8", row, col));
    assert(row == 7);
    assert(col == 8);
    assert(!SoloHDBarcodeMap::parseCoordinateBarcode("AAACCC", row, col));

    ParametersSolo *pSolo = new ParametersSolo();
    pSolo->cbWLstr.push_back("AAAA");
    pSolo->cbWLstr.push_back("CCCC");
    pSolo->cbWLstr.push_back("s_002um_00003_00004-1");
    pSolo->cbWLstrOut.push_back("s_002um_00001_00002-1");
    pSolo->cbWLstrOut.push_back("CCCC");
    pSolo->cbWLstrOut.push_back("s_002um_00003_00004-1");

    SoloHDBarcodeMap directMap;
    std::string error;
    assert(directMap.initialize(*pSolo, "", error));
    assert(directMap.validCount() == 2);
    assert(directMap.coordForWL(0, row, col));
    assert(row == 1);
    assert(col == 2);
    assert(!directMap.coordForWL(1, row, col));
    assert(directMap.coordForWL(2, row, col));
    assert(row == 3);
    assert(col == 4);

    const std::string mapPath = "/tmp/star_hd_molecule_writer_coords.tsv";
    {
        std::ofstream out(mapPath.c_str());
        out << "barcode\trow2\tcol2\n";
        out << "CCCC\t5\t6\n";
    }
    SoloHDBarcodeMap fileMap;
    error.clear();
    assert(fileMap.initialize(*pSolo, mapPath, error));
    assert(fileMap.validCount() == 3);
    assert(fileMap.coordForWL(1, row, col));
    assert(row == 5);
    assert(col == 6);

    uint64_t umiLow = 0;
    uint64_t umiHigh = 0;
    const uint32_t acgtPacked = (0u << 6) | (1u << 4) | (2u << 2) | 3u;
    assert(SoloHDMoleculeWriter::packStarUmi(acgtPacked, 4, umiLow, umiHigh));
    assert(umiLow == 0x3210u);
    assert(umiHigh == 0);

    const std::string moleculePath = "/tmp/star_hd_molecule_writer_test.sshdmol";
    std::remove(moleculePath.c_str());
    SoloHDMoleculeWriter writer;
    error.clear();
    assert(writer.open(moleculePath, 100, SoloHDMoleculeWriter::kFlagPrededuplicated, error));
    SoloHDMoleculeRecord record;
    record.row2 = 5;
    record.col2 = 6;
    record.featureIdx = 7;
    record.umiLength = 4;
    record.umiLow = umiLow;
    record.umiHigh = umiHigh;
    assert(writer.writeRecord(record, error));
    assert(writer.recordCount() == 1);
    assert(writer.close(error));

    const std::vector<unsigned char> bytes = readFile(moleculePath);
    assert(bytes.size() == SoloHDMoleculeWriter::kHeaderSize + SoloHDMoleculeWriter::kRecordSize);
    assert(std::string(reinterpret_cast<const char *>(&bytes[0]), 8) == "SSHDMOL1");
    assert(readLe16(bytes, 8) == SoloHDMoleculeWriter::kVersion);
    assert(readLe16(bytes, 10) == SoloHDMoleculeWriter::kHeaderSize);
    assert(readLe16(bytes, 12) == SoloHDMoleculeWriter::kRecordSize);
    assert(readLe16(bytes, 14) == SoloHDMoleculeWriter::kFlagPrededuplicated);
    assert(readLe32(bytes, 16) == SoloHDMoleculeWriter::kUmiEncoding4Bit);
    assert(readLe64(bytes, 24) == 1);
    assert(readLe64(bytes, 32) == 100);

    const size_t rec = SoloHDMoleculeWriter::kHeaderSize;
    assert(readLe32(bytes, rec + 0) == 5);
    assert(readLe32(bytes, rec + 4) == 6);
    assert(readLe32(bytes, rec + 8) == 7);
    assert(readLe16(bytes, rec + 12) == 4);
    assert(readLe16(bytes, rec + 14) == 0);
    assert(readLe64(bytes, rec + 16) == 0x3210u);
    assert(readLe64(bytes, rec + 24) == 0);

    std::remove(mapPath.c_str());
    std::remove(moleculePath.c_str());
    std::cout << "SoloHDMoleculeWriter tests passed\n";
    return 0;
}
