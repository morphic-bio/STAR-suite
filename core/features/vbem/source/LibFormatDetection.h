#ifndef LIB_FORMAT_DETECTION_H
#define LIB_FORMAT_DETECTION_H

#include "Transcript.h"
#include "ec_builder.h"
#include <string>
#include <map>
#include <cstdint>

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

// Determine LibraryFormat from mate positions and orientations (ported from Salmon)
LibraryFormat hitType(int32_t pos1, bool fwd1, int32_t pos2, bool fwd2);

// Parse user-specified library format string
// Case-insensitive; exits for unrecognized input (validation should catch first)
LibraryFormat parseLibFormat(const std::string& s);

#endif // LIB_FORMAT_DETECTION_H

