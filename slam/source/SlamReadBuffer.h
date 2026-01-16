#ifndef SLAM_READ_BUFFER_H
#define SLAM_READ_BUFFER_H

#include <cstdint>
#include <vector>
#include <string>
#include <set>

// Per-position data within a buffered read
struct SlamBufferedPosition {
    uint32_t readPos;           // Position within concatenated read
    uint64_t genomicPos;        // Genomic position
    uint8_t refBase;            // Reference base (0-3)
    uint8_t readBase;           // Read base (0-3)
    uint8_t qual;               // Quality score (Phred)
    bool secondMate;            // Is this from mate 2?
    bool overlap;               // Is this position in PE overlap region?
};

// Per-read buffered data - stores everything needed to replay counting
struct SlamBufferedRead {
    // Read identification
    std::string readName;       // For debugging/SNP detection
    
    // Read geometry
    uint32_t readLength0;       // Length of mate 1
    uint32_t readLength1;       // Length of mate 2 (0 if single-end)
    bool isMinus;               // Strand (for T→C direction)
    bool oppositeStrand;        // Opposite to gene strand
    
    // Gene assignments
    std::vector<uint32_t> geneIds;    // Assigned genes
    double weight;                     // Assignment weight
    bool isIntronic;                   // Intronic classification
    uint32_t fileIndex = 0;            // Source file index (for trimScope=per-file)
    
    // Per-position data for mismatch/transition counting
    std::vector<SlamBufferedPosition> positions;
    
    // SNP detection data (if enabled)
    std::vector<uint32_t> mismatchGenomicPositions;  // For SNP buffering
    uint16_t nT;                // T count (for EM)
    uint8_t k;                  // Conversion count (for EM)
};

// Buffer for storing reads during auto-trim variance collection
class SlamReadBuffer {
public:
    SlamReadBuffer(uint64_t maxReads = 1000000);
    
    // Add a read to the buffer
    // Returns true if buffer accepted the read (not full)
    bool addRead(SlamBufferedRead&& read);
    
    // Check if buffer is full
    bool isFull() const { return reads_.size() >= maxReads_; }
    
    // Get number of buffered reads
    uint64_t size() const { return reads_.size(); }
    
    // Get max capacity
    uint64_t capacity() const { return maxReads_; }
    
    // Access buffered reads for replay
    const std::vector<SlamBufferedRead>& reads() const { return reads_; }
    
    // Clear buffer (for per-file mode)
    void clear();
    
    // Estimate memory usage
    uint64_t memoryBytes() const;
    
private:
    std::vector<SlamBufferedRead> reads_;
    uint64_t maxReads_;
};

#endif // SLAM_READ_BUFFER_H
