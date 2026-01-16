#include "ReadAlign.h"
#include "SlamQuant.h"
#include "SlamCompat.h"
#include "SlamReadBuffer.h"
#include <sstream>

namespace {
struct GenomicInterval {
    uint64_t start;
    uint64_t end;
};

inline void addInterval(std::vector<GenomicInterval>& intervals, uint64_t start, uint64_t end) {
    if (start <= end) {
        intervals.push_back({start, end});
    }
}

inline bool containsPos(const std::vector<GenomicInterval>& intervals, size_t& idx, uint64_t pos) {
    while (idx < intervals.size() && pos > intervals[idx].end) {
        ++idx;
    }
    return idx < intervals.size() && pos >= intervals[idx].start && pos <= intervals[idx].end;
}

inline bool isOppositeStrand(const Transcriptome& tr, const std::set<uint32_t>& geneIds, uint8_t readStr) {
    bool anySense = false;
    bool anyOpposite = false;
    for (uint32_t geneId : geneIds) {
        if (geneId >= tr.geStr.size()) {
            continue;
        }
        uint8_t geneStr = tr.geStr[geneId];
        if (geneStr == 0) {
            continue;
        }
        bool sense = (geneStr == static_cast<uint8_t>(readStr + 1));
        anySense = anySense || sense;
        anyOpposite = anyOpposite || !sense;
        if (anySense && anyOpposite) {
            break;
        }
    }
    if (anySense) {
        return false;
    }
    return anyOpposite;
}

inline std::string buildReadLoc(const Transcript& trOut, const Genome& genOut) {
    if (trOut.Chr >= genOut.chrName.size()) {
        return "";
    }
    std::ostringstream oss;
    oss << genOut.chrName[trOut.Chr];
    oss << (trOut.Str == 1 ? "-" : "+");
    oss << ":";
    uint64 chrStart = (trOut.Chr < genOut.chrStart.size()) ? genOut.chrStart[trOut.Chr] : 0;
    for (uint iex = 0; iex < trOut.nExons; ++iex) {
        uint64 gStart = trOut.exons[iex][EX_G];
        if (gStart >= chrStart) {
            gStart -= chrStart;
        }
        uint64 len = trOut.exons[iex][EX_L];
        uint64 gEnd = gStart + len;
        if (iex > 0) {
            oss << "|";
        }
        oss << gStart << "-" << gEnd;
    }
    return oss.str();
}

} // namespace

bool ReadAlign::slamCollect(const Transcript& trOut, const std::set<uint32_t>& geneIds, double weight, bool isIntronic) {
    if (slamQuant == nullptr || weight <= 0.0) {
        return false;
    }
    bool debugEnabled = slamQuant->debugEnabled();
    bool debugReadMatch = false;
    if (debugEnabled) {
        debugReadMatch = slamQuant->debugReadMatch(readName);
    }
    if (geneIds.empty()) {
        if (!isIntronic) {
            slamQuant->diagnostics().readsZeroGenes++;
        }
        if (debugEnabled && debugReadMatch) {
            std::string name(readName ? readName : "");
            size_t end = name.find_first_of(" \t");
            if (end != std::string::npos) {
                name = name.substr(0, end);
            }
            if (!name.empty() && name[0] == '@') {
                name.erase(0, 1);
            }
            SlamDebugReadRecord rec;
            rec.readName = name;
            rec.readLoc = buildReadLoc(trOut, genOut);
            rec.status = SlamDebugDropReason::NoGenes;
            rec.readLength = static_cast<uint32_t>(readLength[0] + readLength[1]);
            slamQuant->debugLogRead(rec);
        }
        return false;
    }

    // Record read for variance analysis (only during detection pass)
    bool varianceCollecting = false;
    if (slamQuant && slamQuant->varianceAnalysisEnabled() && P.quant.slam.autoTrimDetectionPass) {
        varianceCollecting = slamQuant->recordVarianceRead();
    }
    
    char* R = Read1[trOut.roStr == 0 ? 0 : 2];
    bool isMinus = (trOut.Str == 1);
    bool oppositeStrand = isOppositeStrand(*chunkTr, geneIds, trOut.Str);
    
    bool debugGeneMatch = false;
    if (debugEnabled && slamQuant->debugGenesEnabled()) {
        for (uint32_t geneId : geneIds) {
            if (slamQuant->debugGeneEnabled(geneId)) {
                debugGeneMatch = true;
                break;
            }
        }
    }
    bool debugThisRead = debugReadMatch || debugGeneMatch;
    std::string readNameStr;
    std::string readLocStr;
    uint32_t readLen = static_cast<uint32_t>(readLength[0] + readLength[1]);
    if (debugEnabled && debugThisRead) {
        readNameStr = readName ? readName : "";
        size_t end = readNameStr.find_first_of(" \t");
        if (end != std::string::npos) {
            readNameStr = readNameStr.substr(0, end);
        }
        if (!readNameStr.empty() && readNameStr[0] == '@') {
            readNameStr.erase(0, 1);
        }
        readLocStr = buildReadLoc(trOut, genOut);
    }
    if (P.quant.slam.strandness == 1 && oppositeStrand) {
        slamQuant->diagnostics().readsDroppedStrandness++;
        if (debugEnabled) {
            for (uint32_t geneId : geneIds) {
                if (slamQuant->debugGeneEnabled(geneId)) {
                    slamQuant->debugCountDrop(geneId, SlamDebugDropReason::Strandness);
                }
                if (debugReadMatch || slamQuant->debugGeneEnabled(geneId)) {
                    SlamDebugReadRecord rec;
                    rec.readName = readNameStr;
                    rec.readLoc = readLocStr;
                    rec.geneId = geneId;
                    rec.intronic = isIntronic;
                    rec.oppositeStrand = oppositeStrand;
                    rec.weight = weight;
                    rec.readLength = readLen;
                    rec.status = SlamDebugDropReason::Strandness;
                    slamQuant->debugLogRead(rec);
                }
            }
        }
        return false;
    }
    if (P.quant.slam.strandness == 2 && !oppositeStrand) {
        slamQuant->diagnostics().readsDroppedStrandness++;
        if (debugEnabled) {
            for (uint32_t geneId : geneIds) {
                if (slamQuant->debugGeneEnabled(geneId)) {
                    slamQuant->debugCountDrop(geneId, SlamDebugDropReason::Strandness);
                }
                if (debugReadMatch || slamQuant->debugGeneEnabled(geneId)) {
                    SlamDebugReadRecord rec;
                    rec.readName = readNameStr;
                    rec.readLoc = readLocStr;
                    rec.geneId = geneId;
                    rec.intronic = isIntronic;
                    rec.oppositeStrand = oppositeStrand;
                    rec.weight = weight;
                    rec.readLength = readLen;
                    rec.status = SlamDebugDropReason::Strandness;
                    slamQuant->debugLogRead(rec);
                }
            }
        }
        return false;
    }
    SlamMismatchCategory category = isIntronic ? SlamMismatchCategory::Intronic : SlamMismatchCategory::Exonic;
    SlamMismatchCategory senseCategory = isIntronic ? SlamMismatchCategory::IntronicSense : SlamMismatchCategory::ExonicSense;

    // Optional dump record for external re-quant
    SlamBufferedRead dumpRead;
    const bool dumpEnabled = slamQuant->dumpEnabled() && !slamQuant->dumpBufferFull();
    if (dumpEnabled) {
        std::string name(readName ? readName : "");
        size_t end = name.find_first_of(" \t");
        if (end != std::string::npos) {
            name = name.substr(0, end);
        }
        if (!name.empty() && name[0] == '@') {
            name.erase(0, 1);
        }
        dumpRead.readName = name;
        dumpRead.readLength0 = static_cast<uint32_t>(readLength[0]);
        dumpRead.readLength1 = static_cast<uint32_t>(readLength[1]);
        dumpRead.isMinus = isMinus;
        dumpRead.oppositeStrand = oppositeStrand;
        dumpRead.geneIds.assign(geneIds.begin(), geneIds.end());
        dumpRead.weight = weight;
        dumpRead.isIntronic = isIntronic;
        dumpRead.fileIndex = static_cast<uint32_t>(P.quant.slam.currentFileIndex);
    }
    bool snpDetect = slamQuant->snpDetectEnabled() && !isIntronic;
    std::vector<uint32_t> mismatchPositions;
    std::vector<uint32_t> debugConvReadPos;
    std::vector<uint64_t> debugConvGenPos;
    bool capturePositions = debugThisRead && !isIntronic;
    if (capturePositions) {
        debugConvReadPos.reserve(8);
        debugConvGenPos.reserve(8);
    }

    uint16_t nT = 0;
    uint16_t k = 0;

    if (slamSnpMask != nullptr) {
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 len = trOut.exons[iex][EX_L];
            for (uint64 ii = 0; ii < len; ++ii) {
                uint64 gpos = gStart + ii;
                if (slamSnpMask->contains(gpos)) {
                    slamQuant->diagnostics().readsDroppedSnpMask++;
                    if (debugEnabled) {
                        for (uint32_t geneId : geneIds) {
                            if (slamQuant->debugGeneEnabled(geneId)) {
                                slamQuant->debugCountDrop(geneId, SlamDebugDropReason::SnpMask);
                            }
                            if (debugReadMatch || slamQuant->debugGeneEnabled(geneId)) {
                                SlamDebugReadRecord rec;
                                rec.readName = readNameStr;
                                rec.readLoc = readLocStr;
                                rec.geneId = geneId;
                                rec.intronic = isIntronic;
                                rec.oppositeStrand = oppositeStrand;
                                rec.weight = weight;
                                rec.readLength = readLen;
                                rec.status = SlamDebugDropReason::SnpMask;
                                slamQuant->debugLogRead(rec);
                            }
                        }
                    }
                    return false; // discard read if it overlaps SNP mask
                }
            }
        }
    }

    std::vector<GenomicInterval> mate1Intervals;
    std::vector<GenomicInterval> mate2Intervals;
    bool hasMate = (readLength[1] > 0);
    if (hasMate) {
        mate1Intervals.reserve(trOut.nExons);
        mate2Intervals.reserve(trOut.nExons);
        uint32_t split = readLength[0];
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 rStart = trOut.exons[iex][EX_R];
            uint64 len = trOut.exons[iex][EX_L];
            uint64 rEnd = rStart + len;
            if (rEnd <= split) {
                addInterval(mate1Intervals, gStart, gStart + len - 1);
            } else if (rStart > split) {
                addInterval(mate2Intervals, gStart, gStart + len - 1);
            } else {
                if (rStart < split) {
                    uint64 len1 = split - rStart;
                    addInterval(mate1Intervals, gStart, gStart + len1 - 1);
                    uint64 len2 = len - len1;
                    if (len2 > 0) {
                        addInterval(mate2Intervals, gStart + len1, gStart + len1 + len2 - 1);
                    }
                }
            }
        }
    }
    size_t mate1Idx = 0;
    size_t mate2Idx = 0;

    for (uint iex = 0; iex < trOut.nExons; ++iex) {
        uint64 gStart = trOut.exons[iex][EX_G];
        uint64 rStart = trOut.exons[iex][EX_R];
        uint64 len = trOut.exons[iex][EX_L];

        for (uint64 ii = 0; ii < len; ++ii) {
            uint64 gpos = gStart + ii;
            uint8_t r1 = static_cast<uint8_t>(R[rStart + ii]);
            uint8_t g1 = static_cast<uint8_t>(genOut.G[gpos]);
            if (r1 > 3 || g1 > 3) {
                continue;
            }
            if (snpDetect) {
                // For SLAM-seq SNP detection, only count T→C conversions
                // STAR base encoding: A=0, C=1, G=2, T=3
                // T→C: genomic T (3) → read C (1)
                // A→G: genomic A (0) → read G (2) (opposite strand equivalent)
                bool isTtoC = (g1 == 3 && r1 == 1);
                bool isAtoG = (g1 == 0 && r1 == 2);
                bool isConv = (isTtoC || isAtoG);
                bool isAnyMis = (g1 != r1);
                if (slamQuant->snpSiteDebugEnabled()) {
                    slamQuant->debugSnpSiteObserve(gpos, isAnyMis, isConv, weight, trOut.primaryFlag, trOut.mapq);
                }
                slamQuant->recordSnpObservation(gpos, isAnyMis, isConv);
            }
            uint64 rposRaw = rStart + ii;
            if (readLength[0] > 0 && rposRaw == readLength[0]) {
                continue; // skip spacer between mates
            }
            bool secondMate = (readLength[1] > 0 && rposRaw > readLength[0]);
            uint32_t readPos = static_cast<uint32_t>(secondMate ? rposRaw - 1 : rposRaw);
            bool overlap = false;
            if (hasMate) {
                if (secondMate) {
                    overlap = containsPos(mate1Intervals, mate1Idx, gpos);
                } else {
                    overlap = containsPos(mate2Intervals, mate2Idx, gpos);
                }
            }
            
            // Get quality score for this position (needed for variance and buffering)
            uint8_t qual = 30; // Default quality if not available
            if (Qual0 && !secondMate && readPos < readLength[0] && Qual0[0] && readPos < strlen(Qual0[0])) {
                qual = static_cast<uint8_t>(Qual0[0][readPos] - 33); // Convert ASCII to Phred
            } else if (Qual0 && secondMate && readLength[1] > 0 && Qual0[1]) {
                uint32_t mateLocalPos = readPos - static_cast<uint32_t>(readLength[0]);
                if (mateLocalPos < strlen(Qual0[1])) {
                    qual = static_cast<uint8_t>(Qual0[1][mateLocalPos] - 33);
                }
            }
            
            // Record variance stats for auto-trim (before trim filtering, only during detection pass)
            // Only collect if varianceCollecting is true (read was recorded and under maxReads limit)
            if (varianceCollecting) {
                // Determine if this is a T base and if it's a T→C conversion
                bool isT = false;
                bool isTc = false;
                if (!isIntronic) {
                    if (!isMinus) {
                        isT = (g1 == 3); // T
                        isTc = (g1 == 3 && r1 == 1); // T→C
                    } else {
                        isT = (g1 == 0); // A (complement of T)
                        isTc = (g1 == 0 && r1 == 2); // A→G (complement of T→C)
                    }
                }
                
                slamQuant->recordVariancePosition(readPos, qual, isT, isTc);
            }

            // Buffer for external re-quant dump (before trim filtering)
            if (dumpEnabled) {
                SlamBufferedPosition bp;
                bp.readPos = readPos;
                bp.genomicPos = gpos;
                bp.refBase = g1;
                bp.readBase = r1;
                bp.qual = qual;
                bp.secondMate = secondMate;
                bp.overlap = overlap;
                dumpRead.positions.push_back(bp);
            }
            
            // Compat position filtering (overlap and trim)
            if (slamCompat) {
                // Check overlap skip first (only when compatIgnoreOverlap is enabled)
                if (slamCompat->cfg().ignoreOverlap && overlap) {
                    slamQuant->diagnostics().compatPositionsSkippedOverlap++;
                    continue;
                }
                // Convert concatenated readPos to mate-local coordinates
                uint32_t mateLocalPos;
                uint32_t mateLen;
                if (secondMate) {
                    mateLocalPos = readPos - static_cast<uint32_t>(readLength[0]);
                    mateLen = static_cast<uint32_t>(readLength[1]);
                } else {
                    mateLocalPos = readPos;
                    mateLen = static_cast<uint32_t>(readLength[0]);
                }
                // Check trim guards (applies even in overlap regions if ignoreOverlap is off)
                if (!slamCompat->compatShouldCountPos(mateLocalPos, mateLen)) {
                    slamQuant->diagnostics().compatPositionsSkippedTrim++;
                    continue;
                }
            }
            
            bool skipMismatch = overlap && secondMate;
            if (!skipMismatch) {
                slamQuant->addTransitionBase(category, readPos, secondMate, overlap, oppositeStrand, g1, r1, weight);
                if (!oppositeStrand) {
                    slamQuant->addTransitionBase(senseCategory, readPos, secondMate, overlap, false, g1, r1, weight);
                }
            }
            
            if (!isIntronic) {
                if (!isMinus) {
                    if (g1 == 3) { // T
                        ++nT;
                        if (r1 == 1) { // C
                            ++k;
                            if (capturePositions) {
                                debugConvReadPos.push_back(readPos);
                                debugConvGenPos.push_back(gpos);
                            }
                            if (snpDetect) {
                                mismatchPositions.push_back(static_cast<uint32_t>(gpos));
                            }
                        }
                    }
                } else {
                    if (g1 == 0) { // A
                        ++nT;
                        if (r1 == 2) { // G
                            ++k;
                            if (capturePositions) {
                                debugConvReadPos.push_back(readPos);
                                debugConvGenPos.push_back(gpos);
                            }
                            if (snpDetect) {
                                mismatchPositions.push_back(static_cast<uint32_t>(gpos));
                            }
                        }
                    }
                }
            }
        }
    }

    // Commit dump record (after position buffering)
    if (dumpEnabled) {
        slamQuant->bufferDumpRead(std::move(dumpRead));
    }

    uint8_t k8 = static_cast<uint8_t>(k > 255 ? 255 : k);
    if (debugEnabled && slamQuant->debugGenesEnabled()) {
        for (uint32_t geneId : geneIds) {
            slamQuant->debugAddAssignment(geneId, weight, isIntronic, oppositeStrand, nT, k8);
        }
    }
    std::string convReadPosStr;
    std::string convGenPosStr;
    if (debugEnabled && debugThisRead && capturePositions) {
        auto joinU32 = [](const std::vector<uint32_t>& vec) -> std::string {
            std::ostringstream oss;
            for (size_t i = 0; i < vec.size(); ++i) {
                if (i > 0) {
                    oss << ",";
                }
                oss << vec[i];
            }
            return oss.str();
        };
        auto joinU64 = [](const std::vector<uint64_t>& vec) -> std::string {
            std::ostringstream oss;
            for (size_t i = 0; i < vec.size(); ++i) {
                if (i > 0) {
                    oss << ",";
                }
                oss << vec[i];
            }
            return oss.str();
        };
        convReadPosStr = joinU32(debugConvReadPos);
        convGenPosStr = joinU64(debugConvGenPos);
    }
    bool snpBuffered = (snpDetect && !mismatchPositions.empty());
    if (debugEnabled && debugThisRead) {
        for (uint32_t geneId : geneIds) {
            if (debugReadMatch || slamQuant->debugGeneEnabled(geneId)) {
                SlamDebugReadRecord rec;
                rec.readName = readNameStr;
                rec.readLoc = readLocStr;
                rec.geneId = geneId;
                rec.intronic = isIntronic;
                rec.oppositeStrand = oppositeStrand;
                rec.weight = weight;
                rec.nT = nT;
                rec.k = k8;
                rec.readLength = readLen;
                rec.status = SlamDebugDropReason::None;
                rec.snpBuffered = snpBuffered;
                rec.convReadPos = convReadPosStr;
                rec.convGenPos = convGenPosStr;
                slamQuant->debugLogRead(rec);
            }
        }
    }

    if (!isIntronic) {
        if (snpDetect && !mismatchPositions.empty()) {
            for (uint32_t geneId : geneIds) {
                slamQuant->bufferSnpRead(geneId, nT, mismatchPositions, weight);
            }
        } else {
            for (uint32_t geneId : geneIds) {
                slamQuant->addRead(geneId, nT, k8, weight);
            }
        }
    }
    slamQuant->diagnostics().readsProcessed++;
    return true;
}
