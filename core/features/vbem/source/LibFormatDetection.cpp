#include "LibFormatDetection.h"
#include "IncludeDefine.h"
#include <array>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>

LibFormatDetector::LibFormatDetector(int window_size)
    : window_size_(window_size), total_votes_(0) {
}

void LibFormatDetector::vote(const Transcript* tr, uint32_t read1_len, uint32_t read2_len) {
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
    
    LibraryFormat fmt = observeFormatFromTranscript(tr, read1_len, read2_len);
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
    
    const auto voteCount = [&](const LibraryFormat& fmt) -> int {
        auto it = format_votes_.find(fmt.typeId());
        return it == format_votes_.end() ? 0 : it->second;
    };

    const int inward_forward = voteCount(LibraryFormat::ISF());
    const int inward_reverse = voteCount(LibraryFormat::ISR());
    const int inward_unstranded = voteCount(LibraryFormat::IU());
    const int inward_total = inward_forward + inward_reverse + inward_unstranded;

    const int outward_forward = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA));
    const int outward_reverse = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS));
    const int outward_unstranded = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U));
    const int outward_total = outward_forward + outward_reverse + outward_unstranded;

    const int same_forward = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S));
    const int same_reverse = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A));
    const int same_unstranded = voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U));
    const int same_total = same_forward + same_reverse + same_unstranded;

    // Find exact-format winner for diagnostics. The selected auto type below
    // follows Salmon's aggregate orientation/strandedness thresholds.
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
    
    LibraryFormat detected_fmt = LibraryFormat::IU();
    const int orientation_total = inward_total + outward_total + same_total;
    const int stranded_forward =
        inward_forward + outward_forward +
        voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S));
    const int stranded_reverse =
        inward_reverse + outward_reverse +
        voteCount(LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A));
    const int stranded_total = stranded_forward + stranded_reverse;

    ReadOrientation detected_orientation = ReadOrientation::TOWARD;
    bool detected_same = false;
    if (orientation_total > 0) {
        if (inward_total >= outward_total && inward_total >= same_total) {
            detected_orientation = ReadOrientation::TOWARD;
        } else if (outward_total >= inward_total && outward_total >= same_total) {
            detected_orientation = ReadOrientation::AWAY;
        } else {
            detected_orientation = ReadOrientation::SAME;
            detected_same = true;
        }
    }

    ReadStrandedness detected_strand = ReadStrandedness::U;
    if (stranded_total > 0) {
        const double ratio_fw =
            static_cast<double>(stranded_forward) / static_cast<double>(stranded_total);
        if (ratio_fw < 0.3) {
            detected_strand = detected_same ? ReadStrandedness::A : ReadStrandedness::AS;
        } else if (ratio_fw < 0.7) {
            detected_strand = ReadStrandedness::U;
        } else {
            detected_strand = detected_same ? ReadStrandedness::S : ReadStrandedness::SA;
        }
    }
    detected_fmt = LibraryFormat(ReadType::PAIRED_END, detected_orientation, detected_strand);

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

    if (!(detected_fmt == winner_fmt)) {
        logStream << "Auto-detect collapsed exact-format votes to "
                  << formatName(detected_fmt)
                  << " using Salmon aggregate orientation/strandedness thresholds.\n";
    }

    return detected_fmt;
}

LibraryFormat LibFormatDetector::observeFormatFromTranscript(const Transcript* tr,
                                                             uint32_t read1_len,
                                                             uint32_t read2_len) {
    TranscriptPairGeometry geometry;
    if (!deriveTranscriptPairGeometry(tr, geometry)) {
        return LibraryFormat::IU();
    }
    return observeFormatFromGeometry(geometry, read1_len, read2_len);
}

std::string formatName(const LibraryFormat& fmt) {
    // Compare against known formats to determine name
    if (fmt == LibraryFormat::IU()) {
        return "IU";
    } else if (fmt == LibraryFormat::ISF()) {
        return "ISF";
    } else if (fmt == LibraryFormat::ISR()) {
        return "ISR";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA)) {
        return "OSF";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS)) {
        return "OSR";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)) {
        return "OU";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)) {
        return "MSF";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)) {
        return "MSR";
    } else if (fmt == LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)) {
        return "MU";
    } else if (fmt == LibraryFormat::U()) {
        return "U";
    } else {
        return "UNKNOWN";
    }
}

bool deriveTranscriptPairGeometry(const Transcript* tr, TranscriptPairGeometry& geometry) {
    geometry = TranscriptPairGeometry{};

    if (tr == nullptr || tr->readNmates != 2 || tr->nExons == 0) {
        return false;
    }

    std::array<bool, 2> seen = {false, false};
    std::array<int32_t, 2> left = {
        std::numeric_limits<int32_t>::max(),
        std::numeric_limits<int32_t>::max()
    };
    std::array<int32_t, 2> right = {
        std::numeric_limits<int32_t>::min(),
        std::numeric_limits<int32_t>::min()
    };

    for (uint i = 0; i < tr->nExons; ++i) {
        const uint mate = tr->exons[i][EX_iFrag];
        if (mate > 1) {
            continue;
        }
        const int32_t exon_left = static_cast<int32_t>(tr->exons[i][EX_G]);
        const int32_t exon_right =
            static_cast<int32_t>(tr->exons[i][EX_G] + tr->exons[i][EX_L] - 1);
        seen[mate] = true;
        left[mate] = std::min(left[mate], exon_left);
        right[mate] = std::max(right[mate], exon_right);
    }

    if (!seen[0] || !seen[1]) {
        return false;
    }

    geometry.read1_left = left[0];
    geometry.read1_right = right[0];
    geometry.read2_left = left[1];
    geometry.read2_right = right[1];
    geometry.read1_leftmost = geometry.read1_left <= geometry.read2_left;

    // Match STAR's transcriptome BAM emission: alignBAM marks a segment as
    // forward when its internal mate label equals tr->Str.
    geometry.read1_forward = (0 == tr->Str);
    geometry.read2_forward = (1 == tr->Str);

    geometry.read1_fiveprime = geometry.read1_forward
        ? geometry.read1_left
        : geometry.read1_right;
    geometry.read2_fiveprime = geometry.read2_forward
        ? geometry.read2_left
        : geometry.read2_right;
    geometry.valid = true;
    return true;
}

LibraryFormat observeFormatFromGeometry(const TranscriptPairGeometry& geometry,
                                        uint32_t read1_len,
                                        uint32_t read2_len) {
    if (!geometry.valid) {
        return LibraryFormat::IU();
    }
    (void)read1_len;
    (void)read2_len;
    return hitType(geometry.read1_left, geometry.read1_forward,
                   geometry.read2_left, geometry.read2_forward);
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

LibraryFormat hitType(int32_t end1Start, bool end1Fwd, uint32_t len1,
                      int32_t end2Start, bool end2Fwd, uint32_t len2,
                      bool canDovetail) {
    if (end1Fwd != end2Fwd) {
        if (end1Fwd) {
            int32_t stretch = canDovetail ? static_cast<int32_t>(len2) : 0;
            if (end1Start <= end2Start + stretch) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
            }
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
        }
        if (end2Fwd) {
            int32_t stretch = canDovetail ? static_cast<int32_t>(len1) : 0;
            if (end2Start <= end1Start + stretch) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
            }
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
        }
    } else {
        if (end1Fwd) {
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
        }
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
    }
    return LibraryFormat::IU();
}

LibraryFormat parseLibFormat(const std::string& s) {
    std::string upper = s;
    for (auto& c : upper) c = std::toupper(c);
    
    if (upper == "IU") return LibraryFormat::IU();
    if (upper == "ISF") return LibraryFormat::ISF();
    if (upper == "ISR") return LibraryFormat::ISR();
    if (upper == "OU") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U);
    if (upper == "OSF") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
    if (upper == "OSR") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
    if (upper == "MU") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U);
    if (upper == "MSF") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
    if (upper == "MSR") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
    if (upper == "U") return LibraryFormat::U();
    
    // "A" should use auto-detect path, not parseLibFormat
    // Any other value is a bug (validation should have caught it)
    std::cerr << "FATAL ERROR: Invalid library format '" << s 
              << "'. This should have been caught during parameter validation.\n";
    exit(1);
}
