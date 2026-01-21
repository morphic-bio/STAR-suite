#include "SoloRead.h"
#include "SoloReadFeature.h"

SoloRead::SoloRead(Parameters &Pin, int32 iChunkIn) :  iChunk(iChunkIn), P(Pin), pSolo(P.pSolo)
{
    readBar = new SoloReadBarcode(P);
    readFeat = nullptr;
    
    if (pSolo.type==0)
        return;
    if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
        return;
    
    readFeat = new SoloReadFeature*[pSolo.nFeatures];

    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii] = new SoloReadFeature(pSolo.features[ii], P, iChunk);
};

SoloRead::~SoloRead() {
    if (readFeat != nullptr) {
        for (uint32 ii = 0; ii < pSolo.nFeatures; ++ii) {
            delete readFeat[ii];
        }
        delete[] readFeat;
        readFeat = nullptr;
    }
    delete readBar;
    readBar = nullptr;
}

void SoloRead::readFlagReset()
{
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        readFeat[ii]->readFlag.flag = 0;
};
