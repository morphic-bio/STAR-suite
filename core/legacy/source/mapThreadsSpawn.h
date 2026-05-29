#ifndef CODE_mapThreadsSpawn
#define CODE_mapThreadsSpawn
#include "Parameters.h"
#include "ReadAlignChunk.h"
#include <string>
void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk);
bool flexPipelineActivationGuard(Parameters &P, std::string *reason = nullptr, bool logMessages = true);
bool flexNoGenomeCountOnlyActivationGuard(Parameters &P, std::string *reason = nullptr);
void runFlexNoGenomeCountOnly(Parameters &P);

#endif
