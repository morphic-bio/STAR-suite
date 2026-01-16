#include "LibFormatDetection.h"
#include "IncludeDefine.h"
#include <iostream>
#include <cstdlib>
#include <algorithm>

LibFormatDetector::LibFormatDetector(int window_size)
    : window_size_(window_size), total_votes_(0) {
}

void LibFormatDetector::vote(const Transcript* tr) {
    // Check if this is a proper PE pair first
    if (tr->readNmates != 2) return;
    
    // Find read1 and read2 to ensure both mates are present
    uint iExRead1 = UINT_MAX;
    uint iExRead2 = UINT_MAX;
    for (uint i = 0; i < tr->nExons; i++) {
        if (tr->exons[i][EX_iFrag] == 0 && iExRead1 == UINT_MAX) {
            iExRead1 = i;
        }
        if (tr->exons[i][EX_iFrag] == 1 && iExRead2 == UINT_MAX) {
            iExRead2 = i;
        }
    }
    
    // If either mate is missing, skip voting (do not vote IU)
    if (iExRead1 == UINT_MAX || iExRead2 == UINT_MAX) {
        return;
    }
    
    LibraryFormat fmt = observeFormatFromTranscript(tr);
    uint8_t fmt_id = fmt.typeId();
    format_votes_[fmt_id]++;
    total_votes_++;
}

LibraryFormat LibFormatDetector::finalizeOrFail(std::ostream& logStream) {
    if (total_votes_ == 0) {
        logStream << "\n"
            << "EXITING because of FATAL ERROR: Library format auto-detection failed.\n"
            << "No proper paired-end alignments found in first " 
            << window_size_ << " reads.\n"
            << "SOLUTION: Specify library type explicitly with:\n"
            << "    --quantVBLibType IU   (unstranded)\n"
            << "    --quantVBLibType ISR  (stranded, first read reverse)\n"
            << "    --quantVBLibType ISF  (stranded, first read forward)\n\n";
        exit(1);  // EXIT BEFORE RETURNING - workers never spawn
    }
    
    // Find winner and check for ambiguity
    uint8_t winner_id = 0;
    int max_votes = 0, second_votes = 0;
    
    for (const auto& kv : format_votes_) {
        if (kv.second > max_votes) {
            second_votes = max_votes;
            max_votes = kv.second;
            winner_id = kv.first;
        } else if (kv.second > second_votes) {
            second_votes = kv.second;
        }
    }
    
    if (second_votes > 0 && (double)second_votes / max_votes > 0.9) {
        logStream << "\n"
            << "EXITING because of FATAL ERROR: Library format auto-detection ambiguous.\n"
            << "Top format: " << max_votes << " votes, "
            << "second: " << second_votes << " votes (within 10%).\n"
            << "SOLUTION: Specify library type explicitly with --quantVBLibType\n\n";
        exit(1);  // EXIT BEFORE RETURNING
    }
    
    // formatName() is defined in the same module (LibFormatDetection.cpp)
    int unknown_votes = 0;
    logStream << "Library format votes: ";
    for (const auto& kv : format_votes_) {
        LibraryFormat fmt = LibraryFormat::formatFromID(kv.first);
        const std::string name = formatName(fmt);
        if (name == "UNKNOWN") {
            unknown_votes += kv.second;
        }
        logStream << name << "(" << static_cast<int>(kv.first) << ")="
                  << kv.second << " ";
    }
    logStream << "\n";

    const double winner_frac = static_cast<double>(max_votes) / static_cast<double>(total_votes_);
    const double unknown_frac = static_cast<double>(unknown_votes) / static_cast<double>(total_votes_);
    LibraryFormat winner_fmt = LibraryFormat::formatFromID(winner_id);
    if (winner_frac < 0.85) {
        logStream << "WARNING: Auto-detect winner is weak (" << max_votes << "/"
                  << total_votes_ << " = " << winner_frac << ").\n";
    }
    if (winner_fmt.orientation == ReadOrientation::AWAY ||
        winner_fmt.orientation == ReadOrientation::SAME ||
        formatName(winner_fmt) == "UNKNOWN") {
        logStream << "WARNING: Auto-detected library format is outward/same-strand ("
                  << formatName(winner_fmt) << "). Check mate order and library prep.\n";
    }
    if (unknown_frac > 0.15) {
        logStream << "WARNING: UNKNOWN-format votes are high (" << unknown_votes << "/"
                  << total_votes_ << " = " << unknown_frac << ").\n";
    }
    
    return LibraryFormat::formatFromID(winner_id);
}

LibraryFormat LibFormatDetector::observeFormatFromTranscript(const Transcript* tr) {
    // Find read1 and read2 positions using EX_iFrag
    // Read1 has EX_iFrag==0, Read2 has EX_iFrag==1
    uint iExRead1 = UINT_MAX;
    uint iExRead2 = UINT_MAX;
    
    for (uint i = 0; i < tr->nExons; i++) {
        if (tr->exons[i][EX_iFrag] == 0 && iExRead1 == UINT_MAX) {
            iExRead1 = i;
        }
        if (tr->exons[i][EX_iFrag] == 1 && iExRead2 == UINT_MAX) {
            iExRead2 = i;
        }
    }
    
    // If either mate is missing, skip voting (do not vote IU)
    if (iExRead1 == UINT_MAX || iExRead2 == UINT_MAX) {
        // Return a sentinel that vote() can check, or throw
        // For now, return IU but vote() should check for proper pairs
        return LibraryFormat::IU();
    }
    
    // Positions in TRANSCRIPTOME coordinates (EX_G for transcriptomic alignment)
    int32_t pos1 = tr->exons[iExRead1][EX_G];
    int32_t pos2 = tr->exons[iExRead2][EX_G];
    
    // Strand determination: tr->Str is the strand of the "left" mate
    // From ReadAlign_outputTranscriptSAM.cpp:
    // - Str is the strand of the leftmost mate (0=forward, 1=reverse)
    // - leftMate = Str (0=read1, 1=read2) - which mate is leftmost
    //
    // For SAM flag generation, when outputting read1 (Mate==0):
    // - read1 strand = Str, read2 strand = 1-Str (assumes opposite strands)
    //
    // However, for same-orientation pairs (MSF/MSR), both reads are on the same strand.
    // Str alone cannot distinguish this case.
    //
    // We determine strand based on which mate is leftmost and try both hypotheses:
    // 1. Opposite strands (most common: inward/outward)
    // 2. Same strand (same-orientation: MSF/MSR)
    //
    // For now, use opposite strands as Str assumes, but note that same-orientation
    // detection will be limited. The hitType function checks `if (end1Fwd != end2Fwd)`
    // first, so if we always set opposite strands, it will never detect same-orientation.
    //
    // However, since Str doesn't encode same-orientation, we must use opposite strands
    // as the default. Same-orientation pairs are rare in typical RNA-seq libraries.
    
    bool read1_leftmost = (pos1 <= pos2);
    bool fwd1, fwd2;
    
    if (read1_leftmost) {
        // Read1 is leftmost: read1 strand = Str
        fwd1 = (tr->Str == 0);
        // Assume opposite strands (as Str does)
        fwd2 = !fwd1;
    } else {
        // Read2 is leftmost: read2 strand = Str  
        fwd2 = (tr->Str == 0);
        // Assume opposite strands
        fwd1 = !fwd2;
    }
    
    // Note: This will correctly detect inward/outward formats (ISF/ISR/OSF/OSR)
    // but may miss same-orientation formats (MSF/MSR) since Str doesn't encode them.
    // This is a limitation of STAR's Transcript structure for same-orientation detection.
    
    return hitType(pos1, fwd1, pos2, fwd2);
}

std::string formatName(const LibraryFormat& fmt) {
    // Compare against known formats to determine name
    if (fmt == LibraryFormat::IU()) {
        return "IU";
    } else if (fmt == LibraryFormat::ISF()) {
        return "ISF";
    } else if (fmt == LibraryFormat::ISR()) {
        return "ISR";
    } else if (fmt == LibraryFormat::U()) {
        return "U";
    } else {
        return "UNKNOWN";
    }
}

// Determine library format from paired-end read positions and orientations
// Ported from Salmon's hitType() function (ec_filter_cli.cpp lines 79-108)
LibraryFormat hitType(int32_t end1Start, bool end1Fwd, int32_t end2Start, bool end2Fwd) {
    // If reads come from opposite strands
    if (end1Fwd != end2Fwd) {
        // Read 1 from forward strand
        if (end1Fwd) {
            if (end1Start <= end2Start) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);  // ISF
            } else {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);     // OSF
            }
        }
        // Read 2 from forward strand
        if (end2Fwd) {
            if (end2Start <= end1Start) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);  // ISR
            } else {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);    // OSR
            }
        }
    } else {
        // Reads from same strand
        if (end1Fwd) {
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);  // MSF
        } else {
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);  // MSR
        }
    }
    // Default fallback
    return LibraryFormat::IU();
}

LibraryFormat parseLibFormat(const std::string& s) {
    std::string upper = s;
    for (auto& c : upper) c = std::toupper(c);
    
    if (upper == "IU") return LibraryFormat::IU();
    if (upper == "ISF") return LibraryFormat::ISF();
    if (upper == "ISR") return LibraryFormat::ISR();
    if (upper == "U") return LibraryFormat::U();
    
    // "A" should use auto-detect path, not parseLibFormat
    // Any other value is a bug (validation should have caught it)
    std::cerr << "FATAL ERROR: Invalid library format '" << s 
              << "'. This should have been caught during parameter validation.\n";
    exit(1);
}
