#include "SoloHDMoleculeWriter.h"

#include "ParametersSolo.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace {

void writeLe16(std::ostream &out, uint16_t value)
{
    const unsigned char bytes[2] = {
        static_cast<unsigned char>(value & 0xffu),
        static_cast<unsigned char>((value >> 8) & 0xffu),
    };
    out.write(reinterpret_cast<const char *>(bytes), 2);
}

void writeLe32(std::ostream &out, uint32_t value)
{
    const unsigned char bytes[4] = {
        static_cast<unsigned char>(value & 0xffu),
        static_cast<unsigned char>((value >> 8) & 0xffu),
        static_cast<unsigned char>((value >> 16) & 0xffu),
        static_cast<unsigned char>((value >> 24) & 0xffu),
    };
    out.write(reinterpret_cast<const char *>(bytes), 4);
}

void writeLe64(std::ostream &out, uint64_t value)
{
    unsigned char bytes[8] = {};
    for (int ii = 0; ii < 8; ++ii) {
        bytes[ii] = static_cast<unsigned char>((value >> (8 * ii)) & 0xffu);
    }
    out.write(reinterpret_cast<const char *>(bytes), 8);
}

std::vector<std::string> splitLine(std::string line)
{
    for (char &ch : line) {
        if (ch == ',') {
            ch = ' ';
        }
    }
    std::istringstream in(line);
    std::vector<std::string> fields;
    std::string field;
    while (in >> field) {
        fields.push_back(field);
    }
    return fields;
}

bool parseUint32(const std::string &s, uint32_t &value)
{
    if (s.empty()) {
        return false;
    }
    char *end = nullptr;
    errno = 0;
    unsigned long parsed = std::strtoul(s.c_str(), &end, 10);
    if (errno != 0 || end == s.c_str() || *end != '\0' || parsed > std::numeric_limits<uint32_t>::max()) {
        return false;
    }
    value = static_cast<uint32_t>(parsed);
    return true;
}

void indexBarcode(const std::string &barcode,
                  uint32_t index,
                  std::unordered_map<std::string, uint32_t> &barcodeToIndex)
{
    if (!barcode.empty() && barcodeToIndex.count(barcode) == 0) {
        barcodeToIndex[barcode] = index;
    }
}

bool assignCoord(uint32_t idx,
                 uint32_t row2,
                 uint32_t col2,
                 std::vector<uint32_t> &rows,
                 std::vector<uint32_t> &cols,
                 std::vector<uint8_t> &valid,
                 std::string &error)
{
    if (idx >= valid.size()) {
        error = "internal barcode coordinate index out of range";
        return false;
    }
    if (valid[idx] != 0 && (rows[idx] != row2 || cols[idx] != col2)) {
        std::ostringstream err;
        err << "conflicting HD coordinates for whitelist index " << idx
            << ": existing=(" << rows[idx] << "," << cols[idx] << ")"
            << " new=(" << row2 << "," << col2 << ")";
        error = err.str();
        return false;
    }
    rows[idx] = row2;
    cols[idx] = col2;
    valid[idx] = 1;
    return true;
}

} // namespace

bool SoloHDBarcodeMap::parseCoordinateBarcode(const std::string &barcode, uint32_t &row2, uint32_t &col2)
{
    const std::string prefix = "s_002um_";
    if (barcode.compare(0, prefix.size(), prefix) != 0) {
        return false;
    }

    const size_t rowStart = prefix.size();
    const size_t rowEnd = barcode.find('_', rowStart);
    if (rowEnd == std::string::npos || rowEnd == rowStart || rowEnd + 1 >= barcode.size()) {
        return false;
    }

    const size_t colStart = rowEnd + 1;
    size_t colEnd = barcode.find('-', colStart);
    if (colEnd == std::string::npos) {
        colEnd = barcode.size();
    }
    if (colEnd == colStart) {
        return false;
    }

    uint32_t parsedRow = 0;
    uint32_t parsedCol = 0;
    if (!parseUint32(barcode.substr(rowStart, rowEnd - rowStart), parsedRow) ||
        !parseUint32(barcode.substr(colStart, colEnd - colStart), parsedCol)) {
        return false;
    }
    row2 = parsedRow;
    col2 = parsedCol;
    return true;
}

bool SoloHDBarcodeMap::initialize(const ParametersSolo &pSolo, const std::string &coordPath, std::string &error)
{
    const size_t n = pSolo.cbWLstr.size();
    row2_.assign(n, 0);
    col2_.assign(n, 0);
    valid_.assign(n, 0);

    std::unordered_map<std::string, uint32_t> barcodeToIndex;
    barcodeToIndex.reserve(n * 2 + 16);
    for (uint32_t ii = 0; ii < n; ++ii) {
        indexBarcode(pSolo.cbWLstr[ii], ii, barcodeToIndex);
        if (pSolo.cbWLstrOut.size() == pSolo.cbWLstr.size()) {
            indexBarcode(pSolo.cbWLstrOut[ii], ii, barcodeToIndex);
        }

        uint32_t row = 0;
        uint32_t col = 0;
        if (pSolo.cbWLstrOut.size() == pSolo.cbWLstr.size() &&
            parseCoordinateBarcode(pSolo.cbWLstrOut[ii], row, col)) {
            if (!assignCoord(ii, row, col, row2_, col2_, valid_, error)) {
                return false;
            }
        } else if (parseCoordinateBarcode(pSolo.cbWLstr[ii], row, col)) {
            if (!assignCoord(ii, row, col, row2_, col2_, valid_, error)) {
                return false;
            }
        }
    }

    if (coordPath.empty()) {
        return validCount() > 0;
    }

    std::ifstream in(coordPath.c_str());
    if (!in.good()) {
        error = "could not open HD barcode coordinate map: " + coordPath + " (" + std::strerror(errno) + ")";
        return false;
    }

    std::string line;
    uint64_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::vector<std::string> fields = splitLine(line);
        if (fields.empty()) {
            continue;
        }
        std::string header = fields[0];
        std::transform(header.begin(), header.end(), header.begin(), ::tolower);
        if (header == "barcode" || header == "cb") {
            continue;
        }

        uint32_t row = 0;
        uint32_t col = 0;
        bool parsedCoord = false;
        if (fields.size() == 1) {
            parsedCoord = parseCoordinateBarcode(fields[0], row, col);
        } else if (parseCoordinateBarcode(fields[1], row, col)) {
            parsedCoord = true;
        } else if (fields.size() >= 3 && parseUint32(fields[1], row) && parseUint32(fields[2], col)) {
            parsedCoord = true;
        }

        if (!parsedCoord) {
            std::ostringstream err;
            err << "could not parse HD coordinate map line " << lineNo << " in " << coordPath;
            error = err.str();
            return false;
        }

        std::unordered_map<std::string, uint32_t>::const_iterator it = barcodeToIndex.find(fields[0]);
        if (it == barcodeToIndex.end()) {
            continue;
        }
        if (!assignCoord(it->second, row, col, row2_, col2_, valid_, error)) {
            return false;
        }
    }

    if (!in.eof()) {
        error = "failed while reading HD barcode coordinate map: " + coordPath;
        return false;
    }
    return validCount() > 0;
}

bool SoloHDBarcodeMap::coordForWL(uint32_t wlIdx, uint32_t &row2, uint32_t &col2) const
{
    if (wlIdx >= valid_.size() || valid_[wlIdx] == 0) {
        return false;
    }
    row2 = row2_[wlIdx];
    col2 = col2_[wlIdx];
    return true;
}

size_t SoloHDBarcodeMap::validCount() const
{
    size_t n = 0;
    for (uint8_t v : valid_) {
        if (v != 0) {
            ++n;
        }
    }
    return n;
}

SoloHDMoleculeWriter::SoloHDMoleculeWriter()
    : featureCount_(0), flags_(0), recordCount_(0), open_(false)
{
}

SoloHDMoleculeWriter::~SoloHDMoleculeWriter()
{
    std::string ignored;
    close(ignored);
}

bool SoloHDMoleculeWriter::open(const std::string &path,
                                uint64_t featureCount,
                                uint16_t flags,
                                std::string &error)
{
    if (open_) {
        error = "HD molecule writer is already open";
        return false;
    }

    out_.open(path.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
    if (!out_.good()) {
        error = "could not open HD molecule output: " + path + " (" + std::strerror(errno) + ")";
        return false;
    }

    featureCount_ = featureCount;
    flags_ = flags;
    recordCount_ = 0;
    open_ = true;

    out_.write("SSHDMOL1", 8);
    writeLe16(out_, kVersion);
    writeLe16(out_, kHeaderSize);
    writeLe16(out_, kRecordSize);
    writeLe16(out_, flags_);
    writeLe32(out_, kUmiEncoding4Bit);
    writeLe32(out_, 0);
    writeLe64(out_, 0);
    writeLe64(out_, featureCount_);
    writeLe64(out_, 0);
    writeLe64(out_, 0);
    writeLe64(out_, 0);
    if (!out_.good()) {
        error = "failed while writing HD molecule header: " + path;
        return false;
    }
    return true;
}

bool SoloHDMoleculeWriter::writeRecord(const SoloHDMoleculeRecord &record, std::string &error)
{
    if (!open_) {
        error = "HD molecule writer is not open";
        return false;
    }
    if (record.umiLength == 0 || record.umiLength > 32) {
        error = "HD molecule record has invalid UMI length";
        return false;
    }
    if (featureCount_ != 0 && record.featureIdx >= featureCount_) {
        std::ostringstream err;
        err << "HD molecule feature index " << record.featureIdx
            << " is outside feature count " << featureCount_;
        error = err.str();
        return false;
    }

    writeLe32(out_, record.row2);
    writeLe32(out_, record.col2);
    writeLe32(out_, record.featureIdx);
    writeLe16(out_, record.umiLength);
    writeLe16(out_, record.flags);
    writeLe64(out_, record.umiLow);
    writeLe64(out_, record.umiHigh);
    if (!out_.good()) {
        error = "failed while writing HD molecule record";
        return false;
    }
    ++recordCount_;
    return true;
}

bool SoloHDMoleculeWriter::writeRecords(const std::vector<SoloHDMoleculeRecord> &records, std::string &error)
{
    for (const SoloHDMoleculeRecord &record : records) {
        if (!writeRecord(record, error)) {
            return false;
        }
    }
    return true;
}

bool SoloHDMoleculeWriter::close(std::string &error)
{
    if (!open_) {
        return true;
    }

    out_.seekp(24, std::ios::beg);
    writeLe64(out_, recordCount_);
    out_.flush();
    if (!out_.good()) {
        error = "failed while finalizing HD molecule header";
        open_ = false;
        out_.close();
        return false;
    }
    out_.close();
    open_ = false;
    return true;
}

uint64_t SoloHDMoleculeWriter::recordCount() const
{
    return recordCount_;
}

bool SoloHDMoleculeWriter::packStarUmi(uint32_t packedUmi,
                                       uint16_t umiLength,
                                       uint64_t &umiLow,
                                       uint64_t &umiHigh)
{
    if (umiLength == 0 || umiLength > 16) {
        return false;
    }

    umiLow = 0;
    umiHigh = 0;
    for (uint16_t ii = 0; ii < umiLength; ++ii) {
        const uint16_t shift = static_cast<uint16_t>(2 * (umiLength - 1 - ii));
        const uint64_t code = (static_cast<uint64_t>(packedUmi) >> shift) & 0x3u;
        if (ii < 16) {
            umiLow |= code << (4 * ii);
        } else {
            umiHigh |= code << (4 * (ii - 16));
        }
    }
    return true;
}
