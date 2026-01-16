#ifndef THREAD_CONTROL_DEF
#define THREAD_CONTROL_DEF

#include "ReadAlignChunk.h"
#include <pthread.h>
#include <atomic>

#define MAX_chunkOutBAMposition 100000

class ThreadControl {
public:
    bool threadBool;

    pthread_t *threadArray;
    pthread_mutex_t mutexInRead, mutexOutSAM, mutexOutBAM1, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout;
    pthread_mutex_t mutexStats, mutexLogMain, mutexBAMsortBins, mutexError;
    pthread_mutex_t mutexOutYFastq[MAX_N_MATES], mutexOutNoYFastq[MAX_N_MATES];  // Y/noY FASTQ output mutexes per mate
    
    // Auto-trim synchronization
    pthread_barrier_t autoTrimBarrier;      // Barrier for file boundary synchronization
    pthread_mutex_t mutexAutoTrim;          // Mutex for auto-trim computation
    std::atomic<bool> autoTrimPending;      // Flag indicating trim computation needed
    std::atomic<uint32_t> autoTrimFileIndex; // File index that triggered trim
    std::atomic<int> autoTrimThreadsDone;   // Threads that completed current file
    bool autoTrimBarrierInitialized;        // Whether barrier was initialized

    uint chunkInN,chunkOutN;

    ThreadControl();

    static void* threadRAprocessChunks(void *RAchunk) {
        ( (ReadAlignChunk*) RAchunk )->processChunks();
        pthread_exit(0);
        return NULL;
    };
};

#endif

