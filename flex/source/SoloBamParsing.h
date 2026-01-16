#ifndef CODE_SoloBamParsing
#define CODE_SoloBamParsing

#include <string>
#include <vector>
#include <cstdint>

namespace SoloBamParsing {

// Parse ZG:Z tag from BAM record
// Returns empty string if tag not found or record malformed
std::string parseZG(const char* bamRecord, uint32_t recordSize);

// Parse UR:Z tag (raw UMI) from BAM record
// Returns empty string if tag not found or record malformed
std::string parseUR(const char* bamRecord, uint32_t recordSize);

// Parse CB:Z (corrected barcode) tag; falls back to CR:Z when missing
std::string parseCB(const char* bamRecord, uint32_t recordSize);

// Parse core BAM fields plus NM and NH tags
// Extracts: QNAME, MAPQ, NM (edit distance), NH (number of hits), read length, FLAG, and read bases
void parseCoreAndNM(const char* bamRecord, uint32_t recordSize,
                    std::string &qname,
                    int &mapq,
                    int &nm,
                    int &nh,
                    int &readLen,
                    uint16_t &flag,
                    std::string &readBases);

} // namespace SoloBamParsing

#endif
