#include "GlobalVariables.h"
#include "SamtoolsSorter.h"
Stats g_statsAll;//global mapping statistics
ThreadControl g_threadChunks;
std::atomic<uint64_t> g_bamRecordIndex{0};
SamtoolsSorter* g_samtoolsSorter = nullptr;

