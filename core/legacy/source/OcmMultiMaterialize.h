#ifndef OCM_MULTI_MATERIALIZE_H
#define OCM_MULTI_MATERIALIZE_H

#include "Parameters.h"
#include "IncludeDefine.h"

namespace OcmMultiMaterialize {

/** Classify corrected cell barcode by OCM overhang (bases 8-9, 1-based). Returns OB1-OB4 or empty. */
string classifyBarcodeTag(const string& barcode);

/** True when --ocmMultiBarcodeMode flex is active. */
bool isFlexBarcodeMode(const Parameters& P);

/** Fixed Flex-compatible TAG8 sequence for an OCM id (OB1-OB4); empty if unknown. */
string tag8ForOcmId(const string& ocmId);

/** Fixed Flex-compatible TAG8 sequence implied by CB16 bases 8-9; empty if unknown. */
string tag8ForBarcode(const string& barcode);

/** Append the OCM TAG8 suffix used by FlexFilter to a raw CB16 barcode; returns empty if unclassified. */
string appendFlexTag8(const string& barcode);

/** Strip an OCM TAG8 suffix from a CB16+TAG8 barcode, preserving any GEM suffix. */
string stripFlexTag8(const string& barcode);

} // namespace OcmMultiMaterialize

/** Materialize OCM per-sample and Cell Ranger multi-compatible MEX outputs. */
int runOcmMultiMaterialize(Parameters& P);

#endif // OCM_MULTI_MATERIALIZE_H
