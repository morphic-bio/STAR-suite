#include "SoloReadFeature.h"
#include "SoloReadFeature_record_shared.h"

// Dispatcher: routes to base (upstream-compatible) or Flex implementation
void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (pSolo.flexMode || pSolo.inlineHashMode) {
        record_flex(this, soloBar, nTr, alignOut, iRead, readAnnot);
    } else {
        record_base(this, soloBar, nTr, alignOut, iRead, readAnnot);
    }
}
