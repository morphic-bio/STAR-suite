#include "ReadAlign.h"
#include "SlamQuant.h"
#include "SlamCompat.h"
#include "SlamReadBuffer.h"
#include <algorithm>
#include <cstring>
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

inline bool containsPosAny(const std::vector<GenomicInterval>& intervals, uint64_t pos) {
    for (const GenomicInterval& interval : intervals) {
        if (pos >= interval.start && pos <= interval.end) {
            return true;
        }
    }
    return false;
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

struct SlamConsensusObservation {
    uint64_t genomicPos = 0;
    uint32_t readPos = 0;
    uint32_t mateLocalPos = 0;
    uint32_t mateLen = 0;
    uint8_t refBase = 0;
    uint8_t readBase = 0;
    uint8_t qual = 0;
    bool secondMate = false;
    bool overlap = false;
    bool trimPass = true;
    bool isT = false;
    bool isTc = false;
    bool isSnpConv = false;
};

inline bool chooseConsensusObservation(const std::vector<SlamConsensusObservation>& observations,
                                       size_t begin, size_t end, bool requireTrimPass,
                                       SlamConsensusObservation& chosen) {
    if (begin >= end) {
        return false;
    }

    int agreedBase = -1;
    for (size_t i = begin; i < end; ++i) {
        int base = static_cast<int>(observations[i].readBase);
        if (agreedBase < 0) {
            agreedBase = base;
        } else if (agreedBase != base) {
            return false;
        }
    }

    const SlamConsensusObservation* best = nullptr;
    for (size_t i = begin; i < end; ++i) {
        const SlamConsensusObservation& obs = observations[i];
        if (requireTrimPass && !obs.trimPass) {
            continue;
        }
        if (best == nullptr || obs.qual > best->qual ||
            (obs.qual == best->qual && best->secondMate && !obs.secondMate)) {
            best = &obs;
        }
    }
    if (best == nullptr) {
        return false;
    }
    chosen = *best;
    return true;
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

    std::vector<GenomicInterval> mateIntervals[2];
    bool hasMate = (readLength[1] > 0);
    if (hasMate) {
        mateIntervals[0].reserve(trOut.nExons);
        mateIntervals[1].reserve(trOut.nExons);
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint mate = trOut.exons[iex][EX_iFrag];
            if (mate > 1) {
                continue;
            }
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 len = trOut.exons[iex][EX_L];
            addInterval(mateIntervals[mate], gStart, gStart + len - 1);
        }
        for (int mate = 0; mate < 2; ++mate) {
            std::sort(mateIntervals[mate].begin(), mateIntervals[mate].end(),
                      [](const GenomicInterval& a, const GenomicInterval& b) {
                          if (a.start != b.start) {
                              return a.start < b.start;
                          }
                          return a.end < b.end;
                      });
        }
    } else if (readLength[0] > 0) {
        mateIntervals[0].reserve(trOut.nExons);
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 len = trOut.exons[iex][EX_L];
            addInterval(mateIntervals[0], gStart, gStart + len - 1);
        }
    }

    auto mateInfoForPosition = [&](uint exonIndex, uint64_t rposRaw,
                                   uint32_t& readPos, uint32_t& mateLocalPos,
                                   uint32_t& mateLen, bool& secondMate) -> bool {
        uint mate = hasMate ? trOut.exons[exonIndex][EX_iFrag] : 0;
        if (mate > 1) {
            return false;
        }
        const uint leftMate = hasMate ? trOut.Str : 0;
        const uint32_t leftMateLen = static_cast<uint32_t>(readLength[leftMate]);
        if (hasMate && rposRaw == leftMateLen) {
            return false;
        }
        uint64_t mateOffset = 0;
        if (hasMate && mate != leftMate) {
            mateOffset = static_cast<uint64_t>(leftMateLen) + 1;
        }
        if (rposRaw < mateOffset) {
            return false;
        }
        uint64_t orientedMatePos = rposRaw - mateOffset;
        mateLen = static_cast<uint32_t>(readLength[mate]);
        if (orientedMatePos >= mateLen || mateLen == 0) {
            return false;
        }

        bool mateReversed = false;
        if (hasMate) {
            mateReversed = (mate == 0) ? (trOut.Str == 1) : (trOut.Str == 0);
        } else {
            mateReversed = (trOut.Str == 1);
        }
        mateLocalPos = mateReversed ?
            (mateLen - 1 - static_cast<uint32_t>(orientedMatePos)) :
            static_cast<uint32_t>(orientedMatePos);
        secondMate = (mate == 1);
        readPos = secondMate ?
            static_cast<uint32_t>(readLength[0]) + mateLocalPos :
            mateLocalPos;
        return true;
    };

    auto qualityForPosition = [&](bool secondMate, uint32_t mateLocalPos,
                                  uint32_t mateLen) -> uint8_t {
        uint8_t qual = 30;
        const uint mate = secondMate ? 1 : 0;
        if (Qual0 && Qual0[mate] && mateLocalPos < mateLen && mateLocalPos < strlen(Qual0[mate])) {
            qual = static_cast<uint8_t>(Qual0[mate][mateLocalPos] - 33);
        }
        return qual;
    };

    auto overlapForPosition = [&](bool secondMate, uint64_t gpos) -> bool {
        if (!hasMate) {
            return false;
        }
        const int mate = secondMate ? 1 : 0;
        return containsPosAny(mateIntervals[1 - mate], gpos);
    };

    std::vector<SlamConsensusObservation> observations;
    observations.reserve(static_cast<size_t>(readLength[0] + readLength[1]));

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
            uint64 rposRaw = rStart + ii;
            uint32_t readPos = 0;
            uint32_t mateLocalPos = 0;
            uint32_t mateLen = 0;
            bool secondMate = false;
            if (!mateInfoForPosition(iex, rposRaw, readPos, mateLocalPos, mateLen, secondMate)) {
                continue;
            }
            bool overlap = overlapForPosition(secondMate, gpos);
            
            // Get quality score for this position (needed for variance and buffering)
            uint8_t qual = qualityForPosition(secondMate, mateLocalPos, mateLen);

            bool isT = false;
            bool isTc = false;
            if (!isIntronic) {
                if (!isMinus) {
                    isT = (g1 == 3);
                    isTc = (g1 == 3 && r1 == 1);
                } else {
                    isT = (g1 == 0);
                    isTc = (g1 == 0 && r1 == 2);
                }
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

            bool trimPass = true;
            if (slamCompat) {
                if (slamCompat->cfg().ignoreOverlap && overlap) {
                    trimPass = false;
                } else if (!slamCompat->compatShouldCountPos(mateLocalPos, mateLen,
                                                            secondMate ? 1u : 0u)) {
                    trimPass = false;
                }
            }

            SlamConsensusObservation obs;
            obs.genomicPos = gpos;
            obs.readPos = readPos;
            obs.mateLocalPos = mateLocalPos;
            obs.mateLen = mateLen;
            obs.refBase = g1;
            obs.readBase = r1;
            obs.qual = qual;
            obs.secondMate = secondMate;
            obs.overlap = overlap;
            obs.trimPass = trimPass;
            obs.isT = isT;
            obs.isTc = isTc;
            obs.isSnpConv = (g1 == 3 && r1 == 1) || (g1 == 0 && r1 == 2);
            observations.push_back(obs);
        }
    }

    std::sort(observations.begin(), observations.end(),
              [](const SlamConsensusObservation& a, const SlamConsensusObservation& b) {
                  if (a.genomicPos != b.genomicPos) {
                      return a.genomicPos < b.genomicPos;
                  }
                  if (a.secondMate != b.secondMate) {
                      return !a.secondMate;
                  }
                  return a.qual > b.qual;
              });

    for (size_t begin = 0; begin < observations.size();) {
        size_t end = begin + 1;
        while (end < observations.size() &&
               observations[end].genomicPos == observations[begin].genomicPos) {
            ++end;
        }

        SlamConsensusObservation untrimmedObs;
        if (chooseConsensusObservation(observations, begin, end, false, untrimmedObs)) {
            if (varianceCollecting) {
                if (hasMate && readLength[1] > 0 && slamQuant->varianceSeparateMates()) {
                    uint8_t mateVarIndex = untrimmedObs.secondMate ? 1 : 0;
                    slamQuant->recordVariancePosition(untrimmedObs.mateLocalPos, mateVarIndex,
                                                      untrimmedObs.qual, untrimmedObs.isT,
                                                      untrimmedObs.isTc);
                } else {
                    slamQuant->recordVariancePosition(untrimmedObs.readPos, untrimmedObs.qual,
                                                      untrimmedObs.isT, untrimmedObs.isTc);
                }
            }
            if (snpDetect) {
                bool isAnyMis = (untrimmedObs.refBase != untrimmedObs.readBase);
                if (slamQuant->snpSiteDebugEnabled()) {
                    slamQuant->debugSnpSiteObserve(untrimmedObs.genomicPos, isAnyMis,
                                                   untrimmedObs.isSnpConv, weight,
                                                   trOut.primaryFlag, trOut.mapq);
                }
                slamQuant->recordSnpObservation(untrimmedObs.genomicPos, isAnyMis,
                                                untrimmedObs.isSnpConv);
            }
        }

        SlamConsensusObservation obs;
        if (chooseConsensusObservation(observations, begin, end, true, obs)) {
            slamQuant->addTransitionBase(category, obs.readPos, obs.secondMate, obs.overlap,
                                         oppositeStrand, obs.refBase, obs.readBase, weight);
            if (!oppositeStrand) {
                slamQuant->addTransitionBase(senseCategory, obs.readPos, obs.secondMate,
                                             obs.overlap, false, obs.refBase, obs.readBase, weight);
            }

            if (!isIntronic && obs.isT) {
                ++nT;
                if (obs.isTc) {
                    ++k;
                    if (capturePositions) {
                        debugConvReadPos.push_back(obs.readPos);
                        debugConvGenPos.push_back(obs.genomicPos);
                    }
                    if (snpDetect) {
                        mismatchPositions.push_back(static_cast<uint32_t>(obs.genomicPos));
                    }
                }
            }
        } else if (slamCompat) {
            bool anyTrimmed = false;
            bool anyOverlapSkipped = false;
            for (size_t i = begin; i < end; ++i) {
                if (!observations[i].trimPass) {
                    anyTrimmed = true;
                    anyOverlapSkipped = anyOverlapSkipped ||
                        (slamCompat->cfg().ignoreOverlap && observations[i].overlap);
                }
            }
            if (anyOverlapSkipped) {
                slamQuant->diagnostics().compatPositionsSkippedOverlap++;
            } else if (anyTrimmed) {
                slamQuant->diagnostics().compatPositionsSkippedTrim++;
            }
        }

        begin = end;
    }

    // Commit dump record (after position buffering)
    if (dumpEnabled) {
        slamQuant->bufferDumpRead(std::move(dumpRead));
    }

    if (debugEnabled && slamQuant->debugGenesEnabled()) {
        for (uint32_t geneId : geneIds) {
            slamQuant->debugAddAssignment(geneId, weight, isIntronic, oppositeStrand, nT, k);
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
                rec.k = k;
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
                slamQuant->addRead(geneId, nT, k, weight);
            }
        }
    }
    slamQuant->diagnostics().readsProcessed++;
    return true;
}
