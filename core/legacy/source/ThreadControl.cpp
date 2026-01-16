#include "ThreadControl.h"

ThreadControl::ThreadControl() {
    chunkInN=0;
    chunkOutN=0;
    // Initialize auto-trim synchronization state
    autoTrimPending.store(false);
    autoTrimFileIndex.store(0);
    autoTrimThreadsDone.store(0);
    autoTrimBarrierInitialized = false;
//     chunkOutBAMposition=new uint [MAX_chunkOutBAMposition];
};