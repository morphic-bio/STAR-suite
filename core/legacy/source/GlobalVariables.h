#ifndef GLOBAL_VARIABLES_DEF
#define GLOBAL_VARIABLES_DEF

#include "ThreadControl.h"
#include <atomic>
#include <string>
#include <vector>

// Forward declaration
class SamtoolsSorter;

extern Stats g_statsAll;
extern ThreadControl g_threadChunks;
extern std::atomic<uint64_t> g_bamRecordIndex;
extern SamtoolsSorter* g_samtoolsSorter;
extern SamtoolsSorter* g_unsortedTagBuffer;  // For unsorted BAM CB/UB tag injection (noSort mode)
extern const std::vector<std::string>* g_bamHeaderChrNames;
extern const std::vector<uint>* g_bamHeaderChrLengths;

#endif
