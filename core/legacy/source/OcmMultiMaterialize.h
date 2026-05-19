#ifndef OCM_MULTI_MATERIALIZE_H
#define OCM_MULTI_MATERIALIZE_H

#include "Parameters.h"
#include "IncludeDefine.h"

namespace OcmMultiMaterialize {

/** Classify corrected cell barcode by OCM overhang (bases 8-9, 1-based). Returns OB1-OB4 or empty. */
string classifyBarcodeTag(const string& barcode);

} // namespace OcmMultiMaterialize

/** Materialize OCM per-sample and Cell Ranger multi-compatible MEX outputs. */
int runOcmMultiMaterialize(Parameters& P);

#endif // OCM_MULTI_MATERIALIZE_H
