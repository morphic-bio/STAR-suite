#include "soloInputFeatureUMI.h"
#include "SoloBinarySpool.h"
#include "SoloReadFeature.h"
#include "binarySearch2.h"

bool soloInputFeatureUMI(fstream *strIn, int32 featureType, bool readInfoYes, bool binarySpool, array<vector<uint64>,2> &sjAll, uint64 &iread, 
                            int32 &cbmatch, uint32 &feature, uint64 &umi, vector<uint32> &featVecU32, SoloReadFlagClass &readFlagCounts)
{
    if (binarySpool) {
        return SoloBinarySpool::readRecordHeader(*strIn, readInfoYes, umi, iread, readFlagCounts.flag, feature, cbmatch);
    }

    if (!(*strIn >> umi)) //end of file
        return false;

    if (readInfoYes) {
        *strIn >> iread >> readFlagCounts.flag;
    };

    switch (featureType) {
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron :
            *strIn >> feature;
            break;

        case SoloFeatureTypes::SJ :
            uint64 sj[2];
            *strIn >> sj[0] >> sj[1];
            feature=(uint32) binarySearch2(sj[0],sj[1],sjAll[0].data(),sjAll[1].data(),sjAll[0].size());
            break;

        case SoloFeatureTypes::Transcript3p :
            feature=0;
            uint32 ntr, in1;
            *strIn >> ntr;
            featVecU32.resize(2*ntr);
            for (uint32 ii=0; ii<2*ntr; ii++) {
                *strIn >> in1;
                featVecU32[ii]=in1;
            };
            break;
        };

    *strIn >> cbmatch;

    return true;
};

bool soloInputFeatureUMI(const SoloBinarySpool::MemoryBuffer &buffer, size_t &offset, int32 featureType, bool readInfoYes, bool binarySpool, array<vector<uint64>,2> &sjAll, uint64 &iread,
                            int32 &cbmatch, uint32 &feature, uint64 &umi, vector<uint32> &featVecU32, SoloReadFlagClass &readFlagCounts)
{
    if (binarySpool) {
        return SoloBinarySpool::readRecordHeader(buffer, offset, readInfoYes, umi, iread, readFlagCounts.flag, feature, cbmatch);
    }

    (void)buffer;
    (void)offset;
    (void)featureType;
    (void)readInfoYes;
    (void)sjAll;
    (void)iread;
    (void)cbmatch;
    (void)feature;
    (void)umi;
    (void)featVecU32;
    (void)readFlagCounts;
    return false;
};
