#ifndef H_SoloReadFeature_record_shared
#define H_SoloReadFeature_record_shared

#include "IncludeDefine.h"
#include "SoloReadBarcode.h"
#include "ReadAnnotations.h"
#include "SoloReadFeatureStats.h"

// Forward declarations
class Transcript;
class SoloReadFeature;
class SoloReadFlagClass;

// Shared data structure used by both base and Flex implementations
class ReadSoloFeatures {
public:
    uint32 gene;
    vector<uint32> geneMult;
    vector<array<uint64,2>> sj;
    bool sjAnnot;
    uint32 indAnnotTr; //index of the annotated transcript
    Transcript **alignOut;
};

// Base implementation function declarations
uint32 outputReadCB_base(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
                         const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                         SoloReadFeature *soloReadFeat = nullptr);

// Flex implementation function declarations  
uint32 outputReadCB_flex(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
                         const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                         SoloReadFeature *soloReadFeat = nullptr);

// Base record implementation (upstream-compatible)
void record_base(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);

// Flex record implementation
void record_flex(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);

#endif
