#include <iostream>
#include <cassert>
#include <cstdint>
#include "../source/CellRangerFormatter.h"
#include <fstream>
#include <cstdlib> // For system()

int main() {
    // Create a test file with known content
    const char* testFile = "/tmp/test_cksum_input.txt";
    const char* testContent = "The quick brown fox jumps over the lazy dog";
    
    // Write test file
    std::ofstream out(testFile);
    out << testContent;
    out.close();
    
    // Compute cksum
    uint32_t computedCrc;
    uint64_t computedSize;
    if (!CellRangerFormatter::computeCksumFile(testFile, computedCrc, computedSize)) {
        std::cerr << "✗ cksum test FAILED: could not compute cksum" << std::endl;
        return 1;
    }
    
    // Compute expected cksum using system cksum command for verification
    // Expected: cksum of "The quick brown fox jumps over the lazy dog\n"
    // We'll compute it using system cksum and compare
    std::string cmd = "cksum " + std::string(testFile) + " | awk '{print $1 \" \" $2}'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "✗ cksum test FAILED: could not run system cksum" << std::endl;
        return 1;
    }
    
    char buffer[128];
    std::string result = "";
    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL) {
            result += buffer;
        }
    }
    pclose(pipe);
    
    // Parse expected values from system cksum output
    uint32_t expectedCrc = 0;
    uint64_t expectedSize = 0;
    if (sscanf(result.c_str(), "%u %lu", &expectedCrc, &expectedSize) != 2) {
        std::cerr << "✗ cksum test FAILED: could not parse system cksum output" << std::endl;
        return 1;
    }
    
    std::cout << "Test file: " << testFile << std::endl;
    std::cout << "Computed cksum: " << computedCrc << " size: " << computedSize << std::endl;
    std::cout << "Expected cksum: " << expectedCrc << " size: " << expectedSize << std::endl;
    
    if (computedCrc == expectedCrc && computedSize == expectedSize) {
        std::cout << "✓ cksum test PASSED" << std::endl;
        return 0;
    } else {
        std::cerr << "✗ cksum test FAILED: mismatch" << std::endl;
        std::cerr << "  Computed: crc=" << computedCrc << " size=" << computedSize << std::endl;
        std::cerr << "  Expected: crc=" << expectedCrc << " size=" << expectedSize << std::endl;
        return 1;
    }
}

