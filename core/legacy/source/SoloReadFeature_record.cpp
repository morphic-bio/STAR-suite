#include "SoloReadFeature.h"
#include "SoloReadFeature_record_shared.h"
#include <cstdlib>

namespace {
bool nonFlexHashBridgeRecordMode(const ParametersSolo &pSolo)
{
    return pSolo.inlineHashMode && !pSolo.flexMode && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
}
}

// Dispatcher: routes to base (upstream-compatible) or Flex implementation
void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (pSolo.pipseqDecoder) {
        const pipseq::Decoded& decoded = soloBar.pipDecoded;
        if (decoded.status != pipseq::DecodeStatus::Accepted) { ++pipReadStats[0]; return; }
        if (decoded.bin < 0) { ++pipReadStats[1]; return; }
        if (nTr == 0) { ++pipReadStats[2]; return; }
        if (nTr != 1) { ++pipReadStats[3]; return; }
        const auto& genes = readAnnot.annotFeatures[featureType].fSet;
        if (genes.size() != 1 || alignOut[0]->nExons == 0) { ++pipReadStats[4]; return; }
        const Transcript& a = *alignOut[0];
        const bool reverse = a.Str != 0;
        const unsigned block = reverse ? a.nExons - 1 : 0;
        // The aligned 5' base: clipping is not extrapolated; strand is part
        // of the key. Internal splice structure does not change this IMI.
        const uint64_t coordinate = a.exons[block][EX_G] +
            (reverse ? a.exons[block][EX_L] - 1 : 0);
        pipGroups[{decoded.cell, *genes.begin()}].add(decoded.bin, coordinate, reverse);
        ++pipReadStats[5];
        return;
    }
    if (nonFlexHashBridgeRecordMode(pSolo)) {
        record_base(this, soloBar, nTr, alignOut, iRead, readAnnot);
    } else if (pSolo.flexMode) {
        record_flex(this, soloBar, nTr, alignOut, iRead, readAnnot);
    } else {
        record_base(this, soloBar, nTr, alignOut, iRead, readAnnot);
    }
}
