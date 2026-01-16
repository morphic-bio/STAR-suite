#include <iostream>
#include <cassert>
#include <algorithm>
#include <cctype>
#include "../source/CellRangerFormatter.h"
#include <fstream>

int main() {
    // Create a test file with known content
    const char* testFile = "/tmp/test_md5_input.txt";
    const char* testContent = "The quick brown fox jumps over the lazy dog";
    
    // Write test file
    std::ofstream out(testFile);
    out << testContent;
    out.close();
    
    // Compute MD5
    std::string computedMD5 = CellRangerFormatter::computeMD5File(testFile);
    
    // Expected MD5 for "The quick brown fox jumps over the lazy dog"
    // (computed using standard MD5 - note: ofstream adds newline)
    // MD5 of "The quick brown fox jumps over the lazy dog\n" = 9e107d9d372bb6826bd81d3542a419d6
    std::string expectedMD5 = "9e107d9d372bb6826bd81d3542a419d6";
    
    std::cout << "Test file: " << testFile << std::endl;
    std::cout << "Computed MD5: " << computedMD5 << std::endl;
    std::cout << "Expected MD5: " << expectedMD5 << std::endl;
    
    // Convert to lowercase for comparison
    std::transform(computedMD5.begin(), computedMD5.end(), computedMD5.begin(), ::tolower);
    std::transform(expectedMD5.begin(), expectedMD5.end(), expectedMD5.begin(), ::tolower);
    
    if (computedMD5 == expectedMD5) {
        std::cout << "✓ MD5 test PASSED" << std::endl;
        return 0;
    } else {
        std::cerr << "✗ MD5 test FAILED: mismatch" << std::endl;
        return 1;
    }
}

