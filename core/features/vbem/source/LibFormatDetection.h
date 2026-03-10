#ifndef LIB_FORMAT_DETECTION_H
#define LIB_FORMAT_DETECTION_H

#include "Transcript.h"
#include "ec_builder.h"
#include <string>
#include <map>
#include <cstdint>

struct TranscriptPairGeometry {
    bool valid = false;
    bool read1_forward = true;
    bool read2_forward = false;
    bool read1_leftmost = true;
    int32_t read1_left = -1;
    int32_t read1_right = -1;
    int32_t read2_left = -1;
    int32_t read2_right = -1;
    int32_t read1_fiveprime = -1;
    int32_t read2_fiveprime = -1;
};

// Library format detector for auto-detection during EC building
class LibFormatDetector {
public:
    LibFormatDetector(int window_size);
    
    // Vote on library format from a Transcript alignment
    void vote(const Transcript* tr);
    
    // Finalize detection and return the detected format
    // Exits with error if detection fails or is ambiguous
    LibraryFormat finalizeOrFail(std::ostream& logStream);
    
    // Get window size
    int windowSize() const { return window_size_; }
    
private:
    int window_size_;
    int total_votes_;
    std::map<uint8_t, int> format_votes_;  // Map from format ID to vote count
    
    // Determine library format from a Transcript object
    LibraryFormat observeFormatFromTranscript(const Transcript* tr);
};

// Helper functions

// Convert LibraryFormat to human-readable string (e.g., "IU", "ISF", "ISR", "U")
std::string formatName(const LibraryFormat& fmt);

// Extract paired-end mate geometry from a STAR Transcript alignment.
// Returns false when a proper paired-end layout is not available.
bool deriveTranscriptPairGeometry(const Transcript* tr, TranscriptPairGeometry& geometry);

// Determine LibraryFormat from extracted mate geometry.
LibraryFormat observeFormatFromGeometry(const TranscriptPairGeometry& geometry);

// Determine LibraryFormat from mate positions and orientations (ported from Salmon)
LibraryFormat hitType(int32_t pos1, bool fwd1, int32_t pos2, bool fwd2);

// Parse user-specified library format string
// Case-insensitive; exits for unrecognized input (validation should catch first)
LibraryFormat parseLibFormat(const std::string& s);

#endif // LIB_FORMAT_DETECTION_H
