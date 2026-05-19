#ifndef VELOCYTO_MEX_WRITER_H
#define VELOCYTO_MEX_WRITER_H

#include "IncludeDefine.h"
#include "MexWriter.h"
#include "Parameters.h"

namespace VelocytoMexWriter {

struct VelocytoMexData {
    vector<string> features;
    vector<string> featureNames;
    vector<string> featureTypes;
    vector<string> barcodes;
    vector<MexWriter::Triplet> spliced;
    vector<MexWriter::Triplet> unspliced;
    vector<MexWriter::Triplet> ambiguous;
};

struct VelocytoRunMex {
    VelocytoMexData raw;
    VelocytoMexData filtered;
};

bool soloVelocytoRawReady(const string& soloOut);
VelocytoRunMex loadSoloVelocytoMex(const string& soloOut);

VelocytoMexData readVelocytoMex(const string& mexDir);
VelocytoMexData subsetVelocytoColumns(const VelocytoMexData& data,
                                      const vector<uint32_t>& colIndices);
int writeVelocytoGzDir(const string& outputDir, const VelocytoMexData& data);

int materializeRunVelocytoMex(const string& runDir,
                              const string& outputRoot,
                              ofstream& logStream);
int runVelocytoMexMaterialize(Parameters& P);

} // namespace VelocytoMexWriter

#endif // VELOCYTO_MEX_WRITER_H
