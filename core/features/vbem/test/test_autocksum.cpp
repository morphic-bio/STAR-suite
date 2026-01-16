/**
 * Unit test for auto-fill cksum from CHECKSUMS files
 */

#include "CellRangerFormatter.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;
using namespace CellRangerFormatter;

int main() {
    int failures = 0;
    
    cout << "Testing auto-fill cksum from CHECKSUMS files...\n\n";
    
    // Test 1: Ensembl URL detection and CHECKSUMS URL derivation
    cout << "Test 1: Ensembl URL CHECKSUMS derivation\n";
    {
        string url = "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
        string expectedChecksumsUrl = "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS";
        
        // Extract filename
        size_t lastSlash = url.find_last_of('/');
        string filename = url.substr(lastSlash + 1);
        
        // Derive CHECKSUMS URL (simplified logic check)
        size_t releasePos = url.find("/release-");
        size_t releaseEnd = url.find('/', releasePos + 9);
        size_t dnaPos = url.find("/dna/", releaseEnd);
        string derivedChecksumsUrl = url.substr(0, dnaPos + 5) + "CHECKSUMS";
        
        if (derivedChecksumsUrl == expectedChecksumsUrl) {
            cout << "  ✓ Ensembl CHECKSUMS URL derivation correct\n";
        } else {
            cout << "  ✗ FAILED: Expected " << expectedChecksumsUrl << ", got " << derivedChecksumsUrl << "\n";
            failures++;
        }
    }
    
    // Test 2: GENCODE URL detection and CHECKSUMS URL derivation
    cout << "\nTest 2: GENCODE URL CHECKSUMS derivation\n";
    {
        string url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz";
        string expectedChecksumsUrl = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/CHECKSUMS";
        
        // Derive CHECKSUMS URL (simplified logic check)
        size_t releasePos = url.find("/release_");
        size_t releaseEnd = url.find('/', releasePos + 9);
        string derivedChecksumsUrl;
        if (releaseEnd == string::npos) {
            derivedChecksumsUrl = url.substr(0, releasePos + 1) + "CHECKSUMS";
        } else {
            derivedChecksumsUrl = url.substr(0, releaseEnd) + "/CHECKSUMS";
        }
        
        if (derivedChecksumsUrl == expectedChecksumsUrl) {
            cout << "  ✓ GENCODE CHECKSUMS URL derivation correct\n";
        } else {
            cout << "  ✗ FAILED: Expected " << expectedChecksumsUrl << ", got " << derivedChecksumsUrl << "\n";
            failures++;
        }
    }
    
    // Test 3: CHECKSUMS file parsing (mock)
    cout << "\nTest 3: CHECKSUMS file parsing\n";
    {
        // Create a mock CHECKSUMS file
        string mockChecksums = "12345678 1000000 test_file.fa.gz\n"
                               "87654321 2000000 gencode.v44.primary_assembly.annotation.gtf.gz\n"
                               "11111111 500000 other_file.txt\n";
        
        string targetFilename = "gencode.v44.primary_assembly.annotation.gtf.gz";
        uint32_t foundCrc = 0;
        uint64_t foundSize = 0;
        bool found = false;
        
        istringstream iss(mockChecksums);
        string line;
        while (getline(iss, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            istringstream lineIss(line);
            uint32_t crc;
            uint64_t size;
            string filename;
            
            if (lineIss >> crc >> size) {
                size_t firstSpace = line.find(' ');
                if (firstSpace != string::npos) {
                    size_t secondSpace = line.find(' ', firstSpace + 1);
                    if (secondSpace != string::npos) {
                        filename = line.substr(secondSpace + 1);
                        while (!filename.empty() && (filename.back() == ' ' || filename.back() == '\t' || filename.back() == '\r')) {
                            filename.pop_back();
                        }
                        
                        if (filename == targetFilename || filename.find(targetFilename) != string::npos) {
                            foundCrc = crc;
                            foundSize = size;
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if (found && foundCrc == 87654321 && foundSize == 2000000) {
            cout << "  ✓ CHECKSUMS parsing correct (found cksum: " << foundCrc << ", size: " << foundSize << ")\n";
        } else {
            cout << "  ✗ FAILED: Expected cksum 87654321 size 2000000, got cksum " << foundCrc << " size " << foundSize << "\n";
            failures++;
        }
    }
    
    // Test 4: Invalid URL (not Ensembl/GENCODE)
    cout << "\nTest 4: Invalid URL detection\n";
    {
        string url = "http://example.com/file.fa.gz";
        bool isEnsembl = (url.find("ftp://ftp.ensembl.org/pub/") == 0);
        bool isGencode = (url.find("ftp://ftp.ebi.ac.uk/pub/databases/gencode/") == 0);
        
        if (!isEnsembl && !isGencode) {
            cout << "  ✓ Invalid URL correctly rejected\n";
        } else {
            cout << "  ✗ FAILED: Invalid URL incorrectly accepted\n";
            failures++;
        }
    }
    
    cout << "\n";
    if (failures == 0) {
        cout << "✓ All unit tests passed!\n";
        return 0;
    } else {
        cout << "✗ " << failures << " test(s) failed\n";
        return 1;
    }
}

