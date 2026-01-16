#include "SoloBamParsing.h"
#include "Parameters.h"
#include "Solo.h"
#include "SoloFeature.h"
#include "ProbeListIndex.h"
#include <cstring>
#include <unordered_map>

namespace {

std::string parseStringTag(const char* bamRecord, uint32_t recordSize, char tag0, char tag1) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(bamRecord);
    if (recordSize < 32) return std::string();
    
    const uint8_t* ptr = p + 4; // skip block_size
    const uint8_t* end = p + recordSize;
    if (end - ptr < 32) return std::string();
    
    ptr += 4; // refID
    ptr += 4; // pos
    uint8_t l_read_name = *ptr++;
    ptr += 1; // mapq
    ptr += 2; // bin
    uint16_t n_cigar_op = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    ptr += 2; // flag
    int32_t l_seq = *reinterpret_cast<const int32_t*>(ptr); ptr += 4;
    ptr += 4; // next_refID
    ptr += 4; // next_pos
    ptr += 4; // tlen
    ptr += l_read_name;
    ptr += static_cast<size_t>(n_cigar_op) * 4u;
    ptr += static_cast<size_t>((l_seq + 1) / 2);
    ptr += static_cast<size_t>(l_seq);
    
    while (ptr + 3 <= end) {
        char t0 = static_cast<char>(ptr[0]);
        char t1 = static_cast<char>(ptr[1]);
        char tt = static_cast<char>(ptr[2]);
        ptr += 3;
        
        if (t0 == tag0 && t1 == tag1 && tt == 'Z') {
            const char* s = reinterpret_cast<const char*>(ptr);
            size_t maxlen = static_cast<size_t>(end - ptr);
            size_t len = strnlen(s, maxlen);
            return std::string(s, len);
        }
        
        switch (tt) {
            case 'A': case 'c': case 'C': ptr += 1; break;
            case 's': case 'S': ptr += 2; break;
            case 'i': case 'I': case 'f': ptr += 4; break;
            case 'Z': case 'H': {
                while (ptr < end && *ptr != 0) ptr++;
                if (ptr < end) ptr++;
                break;
            }
            case 'B': {
                if (ptr + 5 > end) { ptr = end; break; }
                char et = static_cast<char>(*ptr++);
                uint32_t n = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
                size_t esz = (et=='c'||et=='C')?1: (et=='s'||et=='S')?2:4;
                ptr += esz * static_cast<size_t>(n);
                break;
            }
            default: ptr = end; break;
        }
    }
    
    return std::string();
}

} // namespace

namespace SoloBamParsing {

std::string parseZG(const char* bamRecord, uint32_t recordSize) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(bamRecord);
    if (recordSize < 32) return std::string();
    
    const uint8_t* ptr = p + 4; // skip block_size
    const uint8_t* end = p + recordSize;
    if (end - ptr < 32) return std::string();
    
    ptr += 4; // refID
    ptr += 4; // pos
    uint8_t l_read_name = *ptr++;
    ptr += 1; // mapq
    ptr += 2; // bin
    uint16_t n_cigar_op = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    ptr += 2; // flag
    int32_t l_seq = *reinterpret_cast<const int32_t*>(ptr); ptr += 4;
    ptr += 4; // next_refID
    ptr += 4; // next_pos
    ptr += 4; // tlen
    ptr += l_read_name;
    ptr += static_cast<size_t>(n_cigar_op) * 4u;
    ptr += static_cast<size_t>((l_seq + 1) / 2);
    ptr += static_cast<size_t>(l_seq);
    
    // Parse optional tags to find ZG:Z
    while (ptr + 3 <= end) {
        char t0 = static_cast<char>(ptr[0]);
        char t1 = static_cast<char>(ptr[1]);
        char tt = static_cast<char>(ptr[2]);
        ptr += 3;
        
        if (t0 == 'Z' && t1 == 'G' && tt == 'Z') {
            const char* s = reinterpret_cast<const char*>(ptr);
            size_t maxlen = static_cast<size_t>(end - ptr);
            size_t len = strnlen(s, maxlen);
            return std::string(s, len);
        }
        
        // Skip over this tag's value
        switch (tt) {
            case 'A': case 'c': case 'C': ptr += 1; break;
            case 's': case 'S': ptr += 2; break;
            case 'i': case 'I': case 'f': ptr += 4; break;
            case 'Z': case 'H': {
                while (ptr < end && *ptr != 0) ptr++;
                if (ptr < end) ptr++;
                break;
            }
            case 'B': {
                if (ptr + 5 > end) { ptr = end; break; }
                char et = static_cast<char>(*ptr++);
                uint32_t n = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
                size_t esz = (et=='c'||et=='C')?1: (et=='s'||et=='S')?2:4;
                ptr += esz * static_cast<size_t>(n);
                break;
            }
            default: ptr = end; break;
        }
    }
    
    return std::string();
}

std::string parseUR(const char* bamRecord, uint32_t recordSize) {
    return parseStringTag(bamRecord, recordSize, 'U', 'R');
}

std::string parseCB(const char* bamRecord, uint32_t recordSize) {
    std::string cb = parseStringTag(bamRecord, recordSize, 'C', 'B');
    if (cb.empty()) {
        cb = parseStringTag(bamRecord, recordSize, 'C', 'R'); // fallback to raw tag
    }
    return cb;
}

void parseCoreAndNM(const char* bamRecord, uint32_t recordSize,
                    std::string &qname,
                    int &mapq,
                    int &nm,
                    int &nh,
                    int &readLen,
                    uint16_t &flag,
                    std::string &readBases)
{
    const uint8_t* p = reinterpret_cast<const uint8_t*>(bamRecord);
    const uint8_t* ptr = p + 4; // skip block_size
    const uint8_t* end = p + recordSize;
    
    int32_t refID = *reinterpret_cast<const int32_t*>(ptr); ptr += 4; (void)refID;
    int32_t pos = *reinterpret_cast<const int32_t*>(ptr); ptr += 4; (void)pos;
    uint8_t l_read_name = *ptr++;
    mapq = static_cast<int>(*ptr++);
    uint16_t bin = *reinterpret_cast<const uint16_t*>(ptr); (void)bin; ptr += 2;
    uint16_t n_cigar_op = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    flag = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    int32_t l_seq = *reinterpret_cast<const int32_t*>(ptr); ptr += 4;
    readLen = l_seq;
    ptr += 4; // next_refID
    ptr += 4; // next_pos
    ptr += 4; // tlen
    
    // Extract QNAME
    if (ptr + l_read_name > end) return;
    qname.assign(reinterpret_cast<const char*>(ptr), l_read_name - 1); // exclude trailing NUL
    ptr += l_read_name;
    
    // Skip CIGAR
    ptr += static_cast<size_t>(n_cigar_op) * 4u;
    
    // Decode sequence
    const uint8_t* seqPtr = ptr;
    size_t seqBytes = static_cast<size_t>((l_seq + 1) / 2);
    auto nibbleToBase = [](uint8_t code) -> char {
        static const char lut[16] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
        char b = lut[code & 0xF];
        return (b == '=' ? 'N' : b);
    };
    
    readBases.clear();
    readBases.reserve(static_cast<size_t>(l_seq));
    for (size_t i = 0; i < seqBytes; i++) {
        uint8_t byte = seqPtr[i];
        char b1 = nibbleToBase(byte >> 4);
        readBases.push_back(b1);
        if ((i * 2 + 1) < static_cast<size_t>(l_seq)) {
            char b2 = nibbleToBase(byte & 0xF);
            readBases.push_back(b2);
        }
    }
    ptr += seqBytes;
    
    // Skip quality scores
    ptr += static_cast<size_t>(l_seq);
    
    // Parse NM:i and NH:i tags
    nm = -1;
    nh = 1;
    while (ptr + 3 <= end) {
        char t0 = static_cast<char>(ptr[0]);
        char t1 = static_cast<char>(ptr[1]);
        char tt = static_cast<char>(ptr[2]);
        ptr += 3;
        
        if (t0 == 'N' && t1 == 'M' && (tt == 'i' || tt == 'I')) {
            if (ptr + 4 <= end) nm = *reinterpret_cast<const int32_t*>(ptr);
            ptr += 4;
            continue;
        }
        if (t0 == 'N' && t1 == 'H' && (tt == 'i' || tt == 'I')) {
            if (ptr + 4 <= end) nh = *reinterpret_cast<const int32_t*>(ptr);
            return;
        }
        
        // Skip over this tag's value
        switch (tt) {
            case 'A': case 'c': case 'C': ptr += 1; break;
            case 's': case 'S': ptr += 2; break;
            case 'i': case 'I': case 'f': ptr += 4; break;
            case 'Z': case 'H': {
                while (ptr < end && *ptr != 0) ptr++;
                if (ptr < end) ptr++;
                break;
            }
            case 'B': {
                if (ptr + 5 > end) { ptr = end; break; }
                char et = static_cast<char>(*ptr++);
                uint32_t n = *reinterpret_cast<const uint32_t*>(ptr); ptr += 4;
                size_t esz = (et=='c'||et=='C')?1: (et=='s'||et=='S')?2:4;
                ptr += esz * static_cast<size_t>(n);
                break;
            }
            default: ptr = end; break;
        }
    }
}

} // namespace SoloBamParsing
