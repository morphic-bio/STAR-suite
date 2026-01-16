/**
 * test_mex_cb16.cpp - Unit test for MexWriter barcode truncation (cb_len parameter)
 * 
 * Tests:
 * 1. Barcodes are truncated to 16bp when cb_len=16
 * 2. Duplicate detection after truncation fails fast
 * 3. No truncation when cb_len=-1 (default)
 * 4. Short barcodes passed through unchanged
 * 5. Default cb_len preserves original behavior
 * 
 * Note: This tests MexWriter::writeMex() directly. The STAR inline flag plumbing
 * (--soloFlexKeepCBTag) is NOT exercised here since it requires a full STAR build.
 * For STAR inline integration testing, use a full STAR run with the flag.
 * 
 * Compile:
 *   g++ -std=c++11 -I../source -o test_mex_cb16 test_mex_cb16.cpp ../source/MexWriter.cpp
 * 
 * Run:
 *   ./test_mex_cb16
 */

#include "MexWriter.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

// Helper: read first line from file
std::string readFirstLine(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    if (f.is_open()) {
        std::getline(f, line);
    }
    return line;
}

// Helper: count lines in file
int countLines(const std::string& path) {
    std::ifstream f(path);
    int count = 0;
    std::string line;
    while (std::getline(f, line)) {
        count++;
    }
    return count;
}

// Helper: cleanup test files
void cleanup(const std::string& prefix) {
    std::remove((prefix + "barcodes.tsv").c_str());
    std::remove((prefix + "features.tsv").c_str());
    std::remove((prefix + "matrix.mtx").c_str());
}

int main() {
    int passed = 0;
    int failed = 0;
    
    std::string testDir = "/tmp/test_mex_cb16/";
    std::string mkdirCmd = "mkdir -p " + testDir;
    std::system(mkdirCmd.c_str());
    
    // Test data: 24bp barcodes (CB16 + TAG8)
    std::vector<std::string> barcodes24 = {
        "AAACCCAAGAAACACTAAGTAGAG",  // CB16: AAACCCAAGAAACACT, TAG8: AAGTAGAG
        "AAACCCAAGAAACGGAAAGTAGAG",  // CB16: AAACCCAAGAAACGGA, TAG8: AAGTAGAG
        "AAACCCAAGAAACTGAAAGTAGAG"   // CB16: AAACCCAAGAAACTGA, TAG8: AAGTAGAG
    };
    
    std::vector<MexWriter::Feature> features = {
        MexWriter::Feature("GENE1", "Gene1", "Gene Expression"),
        MexWriter::Feature("GENE2", "Gene2", "Gene Expression")
    };
    
    std::vector<MexWriter::Triplet> triplets = {
        {0, 0, 10},
        {1, 1, 20},
        {2, 0, 30}
    };
    
    //-------------------------------------------------------------------------
    // Test 1: Barcodes truncated to 16bp when cb_len=16
    //-------------------------------------------------------------------------
    {
        std::cout << "Test 1: Barcode truncation to 16bp... ";
        cleanup(testDir);
        
        int result = MexWriter::writeMex(testDir, barcodes24, features, triplets, 16);
        
        if (result == 0) {
            // Check barcodes.tsv
            std::ifstream bcFile(testDir + "barcodes.tsv");
            std::vector<std::string> writtenBarcodes;
            std::string line;
            while (std::getline(bcFile, line)) {
                writtenBarcodes.push_back(line);
            }
            
            bool allCorrect = true;
            for (const auto& bc : writtenBarcodes) {
                if (bc.length() != 16) {
                    allCorrect = false;
                    std::cerr << "\n  ERROR: barcode length " << bc.length() << " != 16: " << bc;
                }
            }
            
            if (allCorrect && writtenBarcodes.size() == 3) {
                std::cout << "PASSED\n";
                passed++;
            } else {
                std::cout << "FAILED (wrong barcode count or length)\n";
                failed++;
            }
        } else {
            std::cout << "FAILED (writeMex returned " << result << ")\n";
            failed++;
        }
    }
    
    //-------------------------------------------------------------------------
    // Test 2: No truncation when cb_len=-1 (default)
    //-------------------------------------------------------------------------
    {
        std::cout << "Test 2: No truncation when cb_len=-1... ";
        cleanup(testDir);
        
        int result = MexWriter::writeMex(testDir, barcodes24, features, triplets, -1);
        
        if (result == 0) {
            std::string firstBc = readFirstLine(testDir + "barcodes.tsv");
            
            if (firstBc.length() == 24) {
                std::cout << "PASSED\n";
                passed++;
            } else {
                std::cout << "FAILED (barcode length " << firstBc.length() << " != 24)\n";
                failed++;
            }
        } else {
            std::cout << "FAILED (writeMex returned " << result << ")\n";
            failed++;
        }
    }
    
    //-------------------------------------------------------------------------
    // Test 3: Default cb_len (no explicit parameter) = no truncation
    //-------------------------------------------------------------------------
    {
        std::cout << "Test 3: Default cb_len (no truncation)... ";
        cleanup(testDir);
        
        int result = MexWriter::writeMex(testDir, barcodes24, features, triplets);
        
        if (result == 0) {
            std::string firstBc = readFirstLine(testDir + "barcodes.tsv");
            
            if (firstBc.length() == 24) {
                std::cout << "PASSED\n";
                passed++;
            } else {
                std::cout << "FAILED (barcode length " << firstBc.length() << " != 24)\n";
                failed++;
            }
        } else {
            std::cout << "FAILED (writeMex returned " << result << ")\n";
            failed++;
        }
    }
    
    //-------------------------------------------------------------------------
    // Test 4: Duplicate detection after truncation (should fail)
    //-------------------------------------------------------------------------
    {
        std::cout << "Test 4: Duplicate detection after truncation... ";
        cleanup(testDir);
        
        // Two barcodes with same CB16 prefix but different tags
        std::vector<std::string> duplicateBarcodes = {
            "AAACCCAAGAAACACTAAGTAGAG",  // CB16: AAACCCAAGAAACACT
            "AAACCCAAGAAACACTDIFFERNT"   // CB16: AAACCCAAGAAACACT (same!)
        };
        
        std::vector<MexWriter::Triplet> smallTriplets = {
            {0, 0, 10},
            {1, 0, 20}
        };
        
        // Redirect stderr to suppress expected error message
        int result = MexWriter::writeMex(testDir, duplicateBarcodes, features, smallTriplets, 16);
        
        if (result == -1) {
            std::cout << "PASSED (correctly detected duplicate)\n";
            passed++;
        } else {
            std::cout << "FAILED (should have returned -1 for duplicates)\n";
            failed++;
        }
    }
    
    //-------------------------------------------------------------------------
    // Test 5: Short barcodes (< cb_len) passed through unchanged
    //-------------------------------------------------------------------------
    {
        std::cout << "Test 5: Short barcodes passed unchanged... ";
        cleanup(testDir);
        
        std::vector<std::string> shortBarcodes = {
            "AAACCCAA",         // 8bp (shorter than 16)
            "AAACCCAAGAAACACT"  // 16bp (exactly 16)
        };
        
        std::vector<MexWriter::Triplet> smallTriplets = {
            {0, 0, 10},
            {1, 0, 20}
        };
        
        int result = MexWriter::writeMex(testDir, shortBarcodes, features, smallTriplets, 16);
        
        if (result == 0) {
            std::ifstream bcFile(testDir + "barcodes.tsv");
            std::string bc1, bc2;
            std::getline(bcFile, bc1);
            std::getline(bcFile, bc2);
            
            if (bc1 == "AAACCCAA" && bc2 == "AAACCCAAGAAACACT") {
                std::cout << "PASSED\n";
                passed++;
            } else {
                std::cout << "FAILED (barcodes modified: '" << bc1 << "', '" << bc2 << "')\n";
                failed++;
            }
        } else {
            std::cout << "FAILED (writeMex returned " << result << ")\n";
            failed++;
        }
    }
    
    // Cleanup
    cleanup(testDir);
    std::string rmCmd = "rmdir " + testDir + " 2>/dev/null";
    std::system(rmCmd.c_str());
    
    // Summary
    std::cout << "\n=== Summary ===\n";
    std::cout << "Passed: " << passed << "/" << (passed + failed) << "\n";
    
    return (failed == 0) ? 0 : 1;
}

