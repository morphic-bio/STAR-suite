#include "ReadAlign.h"
#include "IncludeDefine.h"
#include "SequenceFuns.h"
#include "SoloReadFeature.h"
#include "SoloReadFeature_record_shared.h"
#include "SoloFeatureTypes.h"
#include "SampleDetector.h"
#include "serviceFuns.cpp"

#include <cstring>

bool ReadAlign::flexHashCacheValidateSyntheticPair(const char* r2seq, uint32_t lenR2, const char* r1seq, uint32_t lenR1,
                                                   uint16_t expectedGeneIdx15) {
    const bool savedSynth = hashCacheSynthProbe_;
    hashCacheSynthProbe_ = true;

    if (r2seq == nullptr || r1seq == nullptr || lenR2 == 0 || lenR1 == 0) {
        hashCacheSynthProbe_ = savedSynth;
        return false;
    }
    if (P.readNmates != 2) {
        hashCacheSynthProbe_ = savedSynth;
        return false;
    }

    std::strncpy(readName, "@hashCacheSynth", DEF_readNameLengthMax - 1);
    readName[DEF_readNameLengthMax - 1] = '\0';
    readNameMates[0] = readName;

    readLength[0] = lenR2;
    readLength[1] = lenR1;
    readLengthOriginal[0] = lenR2;
    readLengthOriginal[1] = lenR1;

    for (uint32_t i = 0; i < lenR2; ++i) {
        Read0[0][i] = r2seq[i];
        Qual0[0][i] = 'A';
    }
    Read0[0][lenR2] = '\0';
    for (uint32_t i = 0; i < lenR1; ++i) {
        Read0[1][i] = r1seq[i];
        Qual0[1][i] = 'A';
    }
    Read0[1][lenR1] = '\0';

    convertNucleotidesToNumbers(Read0[0], Read1[0], lenR2);
    convertNucleotidesToNumbers(Read0[1], Read1[1], lenR1);

    Lread = lenR2 + lenR1 + 1;
    readLengthPairOriginal = readLengthOriginal[0] + readLengthOriginal[1] + 1;
    if (Lread > DEF_readSeqLengthMax) {
        hashCacheSynthProbe_ = savedSynth;
        return false;
    }

    Read1[0][readLength[0]] = MARK_FRAG_SPACER_BASE;
    complementSeqNumbers(Read1[1], Read1[0] + readLength[0] + 1, readLength[1]);
    for (uint ii = 0; ii < readLength[1] / 2; ii++) {
        std::swap(Read1[0][Lread - ii - 1], Read1[0][ii + readLength[0] + 1]);
    }

    complementSeqNumbers(Read1[0], Read1[1], Lread);
    for (uint ii = 0; ii < Lread; ii++) {
        Read1[2][Lread - ii - 1] = Read1[1][ii];
    }

    outFilterMismatchNmaxTotal = std::min(P.outFilterMismatchNmax,
                                          (uint)(P.outFilterMismatchNoverReadLmax * (readLength[0] + readLength[1])));

    readAnnot.reset();

    if (P.pGe.gType == 101) {
        mapOneReadSpliceGraph();
    } else {
        mapOneRead();
    }

    peOverlapMergeMap();
    multMapSelect();
    mappedFilter();
    transformGenome();

    if (!peOv.yes) {
        chimericDetection();
    }

    if (P.pCh.out.bam && chimRecord) {
        hashCacheSynthProbe_ = savedSynth;
        return false;
    }

    waspMap();
    outputAlignments();

    hashCacheSynthProbe_ = savedSynth;

    if (unmapType >= 0) {
        return false;
    }
    if (nTr == 0) {
        return false;
    }

    if (soloRead == nullptr || soloRead->readBar == nullptr || soloRead->readFeat == nullptr) {
        return false;
    }
    SoloReadBarcode &soloBar = *soloRead->readBar;
    if (soloBar.cbMatch < 0) {
        return false;
    }

    // Match outputReadCB_flex tag policy (sample whitelist / require-match drops tagIdx==0)
    uint8_t tagIdx = 0;
    if (soloBar.detectedSampleToken != 0xFF) {
        uint16_t sampleIdx = SampleDetector::sampleIndexForToken(soloBar.detectedSampleToken);
        if (sampleIdx > 0) {
            tagIdx = static_cast<uint8_t>(sampleIdx & 0x1F);
        }
    }
    const bool dropUnmatchedTag =
        ((!P.pSolo.sampleWhitelistPath.empty()) || P.pSolo.sampleRequireMatch) && (tagIdx == 0);
    if (dropUnmatchedTag) {
        return false;
    }

    // Match hash-screen path: ReadAlign_oneRead uses SoloFeatureTypes::Gene only (not GeneFull / other gene-like features).
    if (!P.pSolo.featureYes[SoloFeatureTypes::Gene] || P.pSolo.featureInd[SoloFeatureTypes::Gene] < 0) {
        return false;
    }
    SoloReadFeature *geneSoloFeat = soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
    const int32_t geneFeatureType = SoloFeatureTypes::Gene;

    const auto &readGe = readAnnot.annotFeatures[geneFeatureType].fSet;
    if (readGe.size() == 0) {
        return false;
    }
    if (readGe.size() > 1 && !P.pSolo.multiMap.yes.multi) {
        return false;
    }

    ReadSoloFeatures reFe;
    reFe.alignOut = trMult;
    reFe.indAnnotTr = 0;

    FlexGeneInlineResolveResult res = flexResolveGeneIdx15_inlineResolver(
        geneSoloFeat, soloBar, reFe, readAnnot, geneFeatureType, (uint64_t)-1);

    return res.geneIdx15 != 0 && res.geneIdx15 == expectedGeneIdx15;
}
