#include "SlamCompat.h"
#include "Transcriptome.h"
#include "Transcript.h"
#include "IncludeDefine.h"
#include <algorithm>

SlamCompat::SlamCompat(const Transcriptome& tr, const SlamCompatConfig& cfg)
    : cfg_(cfg)
{
    buildCaches(tr);
}

SlamCompat::SlamCompat(const SlamCompatConfig& cfg,
                       std::vector<std::vector<uint32_t>>&& geneToTranscripts,
                       std::vector<std::vector<std::pair<uint64_t, uint64_t>>>&& transcriptIntrons)
    : cfg_(cfg)
    , geneToTranscripts_(std::move(geneToTranscripts))
    , transcriptIntrons_(std::move(transcriptIntrons))
{
}

void SlamCompat::buildCaches(const Transcriptome& tr) {
    // Build geneToTranscripts_ from trGene
    geneToTranscripts_.clear();
    geneToTranscripts_.resize(tr.nGe);
    
    for (uint32_t trId = 0; trId < tr.nTr; ++trId) {
        uint32_t geneId = tr.trGene[trId];
        if (geneId < tr.nGe) {
            geneToTranscripts_[geneId].push_back(trId);
        }
    }
    
    // Build transcriptIntrons_ from exSE gaps
    // NOTE: exSE stores positions relative to transcript start (trS)
    // We convert to genomic coordinates here for comparison with alignment positions
    transcriptIntrons_.clear();
    transcriptIntrons_.resize(tr.nTr);
    
    for (uint32_t trId = 0; trId < tr.nTr; ++trId) {
        uint16_t exN = tr.trExN[trId];
        uint32_t exI = tr.trExI[trId];
        
        if (exN < 2) {
            // No introns for single-exon transcripts
            continue;
        }
        
        // Get transcript start (genomic coordinate)
        uint64_t trStart = tr.trS[trId];
        
        // Extract intron intervals (gaps between exons) in GENOMIC coordinates
        for (uint16_t i = 0; i < exN - 1; ++i) {
            // exSE stores relative positions; add trStart for genomic
            uint64_t exonEndGenomic = trStart + tr.exSE[2 * (exI + i) + 1];
            uint64_t nextExonStartGenomic = trStart + tr.exSE[2 * (exI + i + 1)];
            
            if (nextExonStartGenomic > exonEndGenomic + 1) {
                // Intron: [exonEnd+1, nextExonStart-1] in genomic coordinates
                transcriptIntrons_[trId].emplace_back(exonEndGenomic + 1, nextExonStartGenomic - 1);
            }
        }
    }
}

double SlamCompat::compatOverlapWeight(double baseWeight, size_t readLevelGeneCount) const {
    if (!cfg_.overlapWeight || readLevelGeneCount <= 1) return baseWeight;
    return baseWeight / static_cast<double>(readLevelGeneCount);
}

bool SlamCompat::compatShouldCountPos(uint32_t mateLocalPos, uint32_t mateLen) const {
    // Guard against underflow: if trim >= mateLen, skip all positions
    uint32_t trim5 = static_cast<uint32_t>(cfg_.trim5p);
    uint32_t trim3 = static_cast<uint32_t>(cfg_.trim3p);
    if (trim5 + trim3 >= mateLen) return false;
    
    // 5' trim guard (from mate start)
    if (mateLocalPos < trim5) return false;
    // 3' trim guard (from mate end) - safe now due to guard above
    if (mateLocalPos >= mateLen - trim3) return false;
    return true;
}

bool SlamCompat::compatIsIntronic(const Transcript& aln,
                                  const std::set<uint32_t>& geneCandidates,
                                  std::set<uint32_t>& outIntronicGenes) const {
    if (!cfg_.intronic) return false;
    
    // Gate: must be single-part (unspliced) read
    if (aln.nExons != 1) return false;
    
    uint64_t alnStart = aln.exons[0][EX_G];
    uint64_t alnEnd = alnStart + aln.exons[0][EX_L] - 1;
    
    for (uint32_t geneId : geneCandidates) {
        if (geneId >= geneToTranscripts_.size()) continue;
        
        // Count transcripts for this gene where alignment intersects an intron
        size_t intronHitCount = countTranscriptIntronHits(geneId, alnStart, alnEnd);
        if (intronHitCount > 1) {
            outIntronicGenes.insert(geneId);
        }
    }
    return !outIntronicGenes.empty();
}

size_t SlamCompat::countTranscriptIntronHits(uint32_t geneId, uint64_t alnStart, uint64_t alnEnd) const {
    if (geneId >= geneToTranscripts_.size()) return 0;
    
    size_t hitCount = 0;
    for (uint32_t trId : geneToTranscripts_[geneId]) {
        if (trId >= transcriptIntrons_.size()) continue;
        
        // Check if alignment intersects any intron of this transcript
        for (const auto& intron : transcriptIntrons_[trId]) {
            // Check overlap: [alnStart, alnEnd] intersects [intron.first, intron.second]
            if (alnStart <= intron.second && alnEnd >= intron.first) {
                hitCount++;
                break;  // Count transcript once if it has any intersecting intron
            }
        }
    }
    return hitCount;
}

uint32_t SlamCompat::computeExonOverlap(const Transcript& aln, uint64_t trStart,
                                        uint16_t trExN, const uint32_t* trExSE) const {
    uint32_t totalOverlap = 0;
    
    // For each alignment block
    for (uint iex = 0; iex < aln.nExons; ++iex) {
        uint64_t alnStart = aln.exons[iex][EX_G];
        uint64_t alnEnd = alnStart + aln.exons[iex][EX_L] - 1;
        
        // Find overlap with transcript exons
        // trExSE is relative to trStart (genomic coordinates)
        for (uint16_t trEx = 0; trEx < trExN; ++trEx) {
            uint64_t trExStart = trStart + trExSE[2 * trEx];
            uint64_t trExEnd = trStart + trExSE[2 * trEx + 1];
            
            // Compute overlap (both in genomic coordinates)
            uint64_t overlapStart = std::max(alnStart, trExStart);
            uint64_t overlapEnd = std::min(alnEnd, trExEnd);
            
            if (overlapStart <= overlapEnd) {
                totalOverlap += static_cast<uint32_t>(overlapEnd - overlapStart + 1);
            }
        }
    }
    
    return totalOverlap;
}
