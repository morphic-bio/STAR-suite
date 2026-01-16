#include "TranscriptomeFasta.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <array>

bool TranscriptomeFasta::writeTranscriptomeFasta(const Genome& genome, const GTF& gtf, 
                                                  const std::string& outPath,
                                                  Parameters& P,
                                                  uint32_t lineWidth,
                                                  bool overwrite) {
    if (!gtf.gtfYes) {
        P.inOut->logMain << "Skipping transcriptome FASTA generation: no GTF file provided\n";
        return false;
    }
    
    // Check if file exists
    ifstream testFile(outPath.c_str());
    if (testFile.good()) {
        testFile.close();
        if (!overwrite) {
            P.inOut->logMain << "Transcriptome FASTA file already exists: " << outPath << "\n";
            P.inOut->logMain << "Skipping generation (use --genomeGenerateTranscriptomeOverwrite Yes to overwrite)\n";
            return false;
        }
        P.inOut->logMain << "Overwriting existing transcriptome FASTA file: " << outPath << "\n";
    }
    
    // Atomic write: write to temporary file first
    string tmpPath = outPath + ".tmp";
    ofstream outFile(tmpPath.c_str(), ios::out | ios::binary);
    if (outFile.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open transcriptome FASTA output file: " << tmpPath << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    }
    
    // Build transcript -> exons map and compute transcript start/end positions
    // GTF_transcriptGeneSJ.cpp sorts exons by (trStart, trEnd, trId, exStart, exEnd) where:
    // - trStart = start of first exon (when sorted by start)
    // - trEnd = end of last exon (when sorted by start) - NOT max(exon ends)!
    // We must match this ordering exactly for byte-parity with transcriptInfo.tab
    map<uint64, vector<array<uint64, 2>>> transcriptExons; // transcript index -> vector of (start, end) pairs
    vector<array<uint64, 3>> transcriptStartEndId; // (trStart, trEnd, trId) for sorting
    
    for (uint64 iex = 0; iex < gtf.exonN; iex++) {
        uint64 trid = gtf.exonLoci[iex][0]; // exT = 0
        uint64 exStart = gtf.exonLoci[iex][1]; // exS = 1
        uint64 exEnd = gtf.exonLoci[iex][2]; // exE = 2
        transcriptExons[trid].push_back({exStart, exEnd});
    }
    
    // Compute transcript start/end for each transcript
    // Must match GTF_transcriptGeneSJ.cpp logic exactly:
    // - trStart = min exon start (same as start of first exon after sort)
    // - trEnd = end of the last exon when exons are sorted by start position
    uint64 transcriptsWithNoExons = 0;
    for (uint64 trid = 0; trid < gtf.transcriptID.size(); trid++) {
        if (transcriptExons.find(trid) == transcriptExons.end()) {
            transcriptsWithNoExons++;
            continue; // Skip transcripts with no exons
        }
        
        vector<array<uint64, 2>>& exons = transcriptExons[trid];
        
        // Sort exons by start position first
        sort(exons.begin(), exons.end(), [](const array<uint64, 2>& a, const array<uint64, 2>& b) {
            return a[0] < b[0];
        });
        
        // trStart = start of first exon, trEnd = end of last exon (after sorting by start)
        uint64 trStart = exons.front()[0];
        uint64 trEnd = exons.back()[1];  // End of LAST exon in start-sorted order
        transcriptStartEndId.push_back({trStart, trEnd, trid});
    }
    
    // Sort transcripts by (trStart, trEnd, trId) to match transcriptInfo.tab ordering
    // This matches the qsort with funCompareArrays<uint64,5> on (trStart, trEnd, trId, exStart, exEnd)
    // Since we're sorting transcripts (not exons), exStart/exEnd are implicit (first exon for ties)
    sort(transcriptStartEndId.begin(), transcriptStartEndId.end(), 
         [](const array<uint64, 3>& a, const array<uint64, 3>& b) {
             if (a[0] != b[0]) return a[0] < b[0]; // Sort by trStart
             if (a[1] != b[1]) return a[1] < b[1]; // Then by trEnd
             return a[2] < b[2]; // Then by trId (for deterministic ordering)
         });
    
    // Nucleotide code to character conversion
    const char nuclChar[] = "ACGT";
    
    uint64 transcriptsWritten = 0;
    uint64 totalBases = 0;
    uint64 skippedTranscripts = 0;
    
    // Process transcripts in sorted order (to match transcriptInfo.tab ordering)
    for (const auto& trInfo : transcriptStartEndId) {
        uint64 trid = trInfo[2];
        // transcriptStartEndId only contains transcripts with exons, so no need to check
        
        // Exons already sorted by start position during trEnd computation above
        vector<array<uint64, 2>>& exons = transcriptExons[trid];
        
        // Extract transcript sequence
        string transcriptSeq;
        bool skipTranscript = false;
        
        for (const auto& exon : exons) {
            uint64 exStart = exon[0];
            uint64 exEnd = exon[1];
            
            // Find chromosome for this exon
            uint chrIdx = genome.chrBin[exStart >> genome.pGe.gChrBinNbits];
            
            // Validate chromosome index
            if (chrIdx >= genome.nChrReal) {
                P.inOut->logMain << "WARNING: Exon at genome position " << exStart 
                                 << " maps to invalid chromosome index " << chrIdx 
                                 << ", skipping transcript " << gtf.transcriptID[trid] << "\n";
                skipTranscript = true;
                break;
            }
            
            // Check bounds
            uint64 chrStartPos = genome.chrStart[chrIdx];
            uint64 chrEndPos = genome.chrStart[chrIdx] + genome.chrLength[chrIdx];
            
            if (exStart < chrStartPos || exEnd >= chrEndPos) {
                P.inOut->logMain << "WARNING: Exon coordinates out of bounds for chromosome " 
                                 << genome.chrName[chrIdx] << ", skipping transcript " 
                                 << gtf.transcriptID[trid] << "\n";
                skipTranscript = true;
                break;
            }
            
            // Extract sequence from genome
            for (uint64 pos = exStart; pos <= exEnd; pos++) {
                uint8_t code = static_cast<uint8_t>(genome.G[pos]);
                if (code < 4) {
                    transcriptSeq += nuclChar[code];
                } else {
                    transcriptSeq += 'N'; // N or other
                }
            }
        }
        
        if (skipTranscript) {
            skippedTranscripts++;
            continue;
        }
        
        // Handle strand: reverse-complement if negative strand
        // transcriptStrand: 1 = '+', 2 = '-', 0 = '.'
        if (gtf.transcriptStrand[trid] == 2) {
            // Negative strand: reverse exon order and reverse-complement sequence
            // Note: exons are already in genomic order, so we reverse-complement the whole sequence
            revComplementNucleotides(transcriptSeq);
        }
        
        // Write FASTA header (transcript_id without version)
        string transcriptID = gtf.transcriptID[trid];
        // Remove version suffix if present (e.g., "ENST00000123456.1" -> "ENST00000123456")
        size_t dotPos = transcriptID.find('.');
        if (dotPos != string::npos) {
            transcriptID = transcriptID.substr(0, dotPos);
        }
        
        outFile << ">" << transcriptID << "\n";
        
        // Write sequence with line wrapping
        uint64 seqLen = transcriptSeq.length();
        for (uint64 i = 0; i < seqLen; i += lineWidth) {
            uint64 chunkLen = min((uint64)lineWidth, seqLen - i);
            outFile.write(transcriptSeq.c_str() + i, chunkLen);
            outFile << "\n";
        }
        
        transcriptsWritten++;
        totalBases += seqLen;
    }
    
    outFile.close();
    
    // Atomic rename
    if (rename(tmpPath.c_str(), outPath.c_str()) != 0) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not rename temporary transcriptome FASTA file from " 
               << tmpPath << " to " << outPath << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    }
    
    P.inOut->logMain << "Transcriptome FASTA generation completed:\n";
    P.inOut->logMain << "  Output file: " << outPath << "\n";
    P.inOut->logMain << "  Transcripts written: " << transcriptsWritten << "\n";
    P.inOut->logMain << "  Transcripts skipped (no exons): " << transcriptsWithNoExons << "\n";
    P.inOut->logMain << "  Transcripts skipped (extraction issues): " << skippedTranscripts << "\n";
    P.inOut->logMain << "  Total bases: " << totalBases << "\n";
    
    return true;
}
