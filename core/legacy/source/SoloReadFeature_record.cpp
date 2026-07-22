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
    if (nonFlexHashBridgeRecordMode(pSolo)) {
        record_base(this, soloBar, nTr, alignOut, iRead, readAnnot);
    } else if (pSolo.flexMode) {
        record_flex(this, soloBar, nTr, alignOut, iRead, readAnnot);
    } else {
        record_base(this, soloBar, nTr, alignOut, iRead, readAnnot);
    }
}
