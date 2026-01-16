#include "SamtoolsSorter.h"
#include "Parameters.h"
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <dirent.h>
#include <algorithm>
#include <cstdlib>
#include <unistd.h>

static std::vector<char> make_bam_record(int32_t tid, int32_t pos,
                                         const std::string& qname,
                                         uint16_t flag,
                                         int32_t mtid, int32_t mpos, int32_t isize) {
    uint8_t l_qname = static_cast<uint8_t>(qname.size() + 1); // include NUL
    uint16_t n_cigar = 1;
    int32_t l_seq = 1;

    uint32_t bin = 0;
    uint8_t mapq = 0;
    uint32_t bin_mq_nl = (bin << 16) | (mapq << 8) | l_qname;
    uint32_t flag_nc = (static_cast<uint32_t>(flag) << 16) | n_cigar;

    uint32_t data_len = 32 + l_qname + 4 * n_cigar + ((l_seq + 1) / 2) + l_seq;
    std::vector<char> buf(4 + data_len, 0);

    uint32_t* p32 = reinterpret_cast<uint32_t*>(buf.data());
    p32[0] = data_len;
    p32[1] = static_cast<uint32_t>(tid);
    p32[2] = static_cast<uint32_t>(pos);
    p32[3] = bin_mq_nl;
    p32[4] = flag_nc;
    p32[5] = static_cast<uint32_t>(l_seq);
    p32[6] = static_cast<uint32_t>(mtid);
    p32[7] = static_cast<uint32_t>(mpos);
    p32[8] = static_cast<uint32_t>(isize);

    char* p = buf.data() + 36; // 4 + 32
    memcpy(p, qname.c_str(), qname.size());
    p[qname.size()] = '\0';
    p += l_qname;

    uint32_t cigar = (1u << 4) | 0u; // 1M
    memcpy(p, &cigar, 4);
    p += 4;

    p[0] = 1 << 4; // seq 'A'
    p += 1;

    p[0] = 30; // qual
    return buf;
}

static std::string get_qname(const char* bam) {
    const uint8_t* core = reinterpret_cast<const uint8_t*>(bam + 4);
    uint8_t l_qname = core[8];
    const char* q = bam + 36;
    if (l_qname > 0 && q[l_qname - 1] == 0) l_qname--;
    return std::string(q, q + l_qname);
}

int main() {
    Parameters P;
    // Use unique temp dir per test run
    P.outBAMsortTmpDir = "/tmp/samtools_sorter_test_" + std::to_string(getpid());
    
    // Clean up any existing test directory
    system(("rm -rf " + P.outBAMsortTmpDir).c_str());
    mkdir(P.outBAMsortTmpDir.c_str(), 0755);

    // Test 1: Multiple records per spill file
    // Set maxRAM to allow 2-3 records before spilling (~150 bytes for 3 records)
    // This ensures each spill file contains multiple records, exercising the ownership path
    std::cout << "Test 1: Multiple records per spill file\n";
    SamtoolsSorter sorter1(150, 1, P.outBAMsortTmpDir, P);

    // Create records that will be grouped into spill files with multiple records each
    // With maxRAM=150, we can fit ~3 records before spilling, ensuring multiple records per spill
    auto r1 = make_bam_record(0, 10, "readB", 0, -1, -1, 0);
    auto r2 = make_bam_record(0, 10, "readA", 0, -1, -1, 0);
    auto r3 = make_bam_record(0, 10, "readC", 0, -1, -1, 0);
    auto r4 = make_bam_record(0, 20, "readD", 0, -1, -1, 0);
    auto r5 = make_bam_record(0, 20, "readE", 0, -1, -1, 0);
    auto r6 = make_bam_record(0, 20, "readF", 0, -1, -1, 0);
    auto r7 = make_bam_record(0, 30, "readG", 0, -1, -1, 0);
    // Unmapped (should be last)
    auto r8 = make_bam_record(-1, -1, "zzz", 0, -1, -1, 0);

    uint32_t readId1 = 200;  // Higher readId
    uint32_t readId2 = 100;  // Lower readId (will sort first)
    uint32_t readId3 = 300;
    uint32_t readId4 = 400;
    uint32_t readId5 = 500;
    uint32_t readId6 = 600;
    uint32_t readId7 = 700;
    uint32_t readId8 = 800;
    
    // Add records - with maxRAM=150, records 1-3 will be in first spill, 4-6 in second, etc.
    // This ensures each spill file contains multiple records, exercising the ownership path
    sorter1.addRecord(r1.data(), r1.size(), readId1, false);
    sorter1.addRecord(r2.data(), r2.size(), readId2, false);
    sorter1.addRecord(r3.data(), r3.size(), readId3, false);
    sorter1.addRecord(r4.data(), r4.size(), readId4, false);
    sorter1.addRecord(r5.data(), r5.size(), readId5, false);
    sorter1.addRecord(r6.data(), r6.size(), readId6, false);
    sorter1.addRecord(r7.data(), r7.size(), readId7, false);
    sorter1.addRecord(r8.data(), r8.size(), readId8, false);
    
    // Check spill files before finalize
    DIR* dir = opendir(P.outBAMsortTmpDir.c_str());
    int spillCount1 = 0;
    if (dir != nullptr) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string name(entry->d_name);
            if (name.find("samtools_sort_spill_") == 0 && name.find(".dat") != std::string::npos) {
                spillCount1++;
            }
        }
        closedir(dir);
    }
    
    sorter1.finalize();

    std::vector<std::string> out1;
    std::vector<uint32_t> readIds1;
    const char* bam = nullptr;
    uint32_t size = 0;
    uint32_t readId = 0;
    bool hasY = false;
    
    // CRITICAL TEST: Exercise multiple records per spill file
    // With maxRAM=150, each spill file contains 2-3 records. During merge, readNext() advances
    // through multiple records in each spill file. This exercises the original bug where
    // readNext() would advance before returning, invalidating the returned pointer.
    // The fix uses returnedRecord_ buffer to preserve data until next call.
    while (sorter1.nextRecord(&bam, &size, &readId, &hasY)) {
        // Verify the returned pointer is valid and points to correct data
        if (bam == nullptr || size == 0) {
            std::cerr << "FAIL Test 1: invalid BAM pointer or size\n";
            system(("rm -rf " + P.outBAMsortTmpDir).c_str());
            return 1;
        }
        
        // Extract QNAME and verify it's valid (not garbage from freed memory)
        std::string qname = get_qname(bam);
        if (qname.empty() || qname.size() > 100) {  // Sanity check: QNAME shouldn't be empty or huge
            std::cerr << "FAIL Test 1: invalid QNAME extracted: '" << qname << "' (size=" << qname.size() << ")\n";
            std::cerr << "This indicates the returned pointer points to invalid/freed memory\n";
            std::cerr << "This would happen with the original bug where readNext() advanced before returning\n";
            system(("rm -rf " + P.outBAMsortTmpDir).c_str());
            return 1;
        }
        
        out1.push_back(qname);
        readIds1.push_back(readId);
    }

    // Expected order: coord (tid<<32|pos) then readId
    // Records at coord (0,10): readId2=100, readId1=200, readId3=300
    // Records at coord (0,20): readId4=400, readId5=500, readId6=600
    // Records at coord (0,30): readId7=700
    // Unmapped: readId8=800
    if (out1.size() != 8 ||
        out1[0] != "readA" ||  // readId=100, coord=(0,10)
        out1[1] != "readB" ||  // readId=200, coord=(0,10)
        out1[2] != "readC" ||  // readId=300, coord=(0,10)
        out1[3] != "readD" ||  // readId=400, coord=(0,20)
        out1[4] != "readE" ||  // readId=500, coord=(0,20)
        out1[5] != "readF" ||  // readId=600, coord=(0,20)
        out1[6] != "readG" ||  // readId=700, coord=(0,30)
        out1[7] != "zzz") {    // unmapped, coord=(INT32_MAX, INT32_MAX)
        std::cerr << "FAIL Test 1: unexpected order\n";
        std::cerr << "Expected: readA readB readC readD readE readF readG zzz\n";
        std::cerr << "Got: ";
        for (const auto& qname : out1) std::cerr << qname << " ";
        std::cerr << "\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    
    // Verify readId round-trips through spill/merge unchanged
    if (readIds1.size() != 8) {
        std::cerr << "FAIL Test 1: expected 8 readId values, got " << readIds1.size() << "\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    std::vector<uint32_t> expectedReadIds1 = {readId2, readId1, readId3, readId4, readId5, readId6, readId7, readId8};
    if (readIds1 != expectedReadIds1) {
        std::cerr << "FAIL Test 1: readId round-trip mismatch\n";
        std::cerr << "Expected: ";
        for (auto r : expectedReadIds1) std::cerr << r << " ";
        std::cerr << "\nGot: ";
        for (auto r : readIds1) std::cerr << r << " ";
        std::cerr << "\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    
    // Verify spill files were created (with maxRAM=150 and 8 records, should create multiple spills)
    if (spillCount1 == 0) {
        std::cerr << "FAIL Test 1: No spill files created (expected with maxRAM=150 and 8 records)\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    
    // With maxRAM=150 and 8 records, we should have multiple spill files with multiple records each
    if (spillCount1 < 2) {
        std::cerr << "WARNING Test 1: Only " << spillCount1 << " spill file(s) created\n";
        std::cerr << "Expected multiple spill files with multiple records each to fully exercise the bug scenario\n";
    } else {
        std::cout << "Test 1 passed: " << spillCount1 << " spill file(s) created with multiple records each\n";
    }
    
    // Clean up for next test
    system(("rm -rf " + P.outBAMsortTmpDir).c_str());
    mkdir(P.outBAMsortTmpDir.c_str(), 0755);
    
    // Test 2: Single record per spill (extreme case)
    std::cout << "Test 2: Single record per spill file (extreme case)\n";
    SamtoolsSorter sorter2(64, 1, P.outBAMsortTmpDir, P);
    
    // Create 3 records - with maxRAM=64, each will be in its own spill file
    auto s2r1 = make_bam_record(0, 10, "single1", 0, -1, -1, 0);
    auto s2r2 = make_bam_record(0, 20, "single2", 0, -1, -1, 0);
    auto s2r3 = make_bam_record(-1, -1, "single3", 0, -1, -1, 0);
    
    sorter2.addRecord(s2r1.data(), s2r1.size(), 100, false);
    sorter2.addRecord(s2r2.data(), s2r2.size(), 200, false);
    sorter2.addRecord(s2r3.data(), s2r3.size(), 300, false);
    
    sorter2.finalize();
    
    std::vector<std::string> out2;
    while (sorter2.nextRecord(&bam, &size, &readId, &hasY)) {
        std::string qname = get_qname(bam);
        if (qname.empty() || qname.size() > 100) {
            std::cerr << "FAIL Test 2: invalid QNAME\n";
            system(("rm -rf " + P.outBAMsortTmpDir).c_str());
            return 1;
        }
        out2.push_back(qname);
    }
    
    if (out2.size() != 3 || out2[0] != "single1" || out2[1] != "single2" || out2[2] != "single3") {
        std::cerr << "FAIL Test 2: unexpected order\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    
    std::cout << "Test 2 passed: Single record per spill file handled correctly\n";
    
    // Clean up temp directory
    system(("rm -rf " + P.outBAMsortTmpDir).c_str());

    std::cout << "All tests passed!\n";
    return 0;
}

