#include "Parameters.h"
#include "ErrorWarning.h"
#include "input/CbqInputModule.h"
#include "input/FastxInputModule.h"
#include "input/BgzfStarAdapter.h"
#include "input/BgzfPipeGroup.h"
#include "input/FastxMateReaders.h"
#include <fstream>
#include <sys/stat.h>
#include <cerrno>
#include <csignal>
#include <sys/wait.h>

namespace {
// Terminate and reap readFilesCommand helper children to avoid lingering
// processes at shutdown.
void terminateReadCommandChildren(Parameters& P) {
    for (uint imate=0; imate<MAX_N_MATES; imate++) {
        pid_t pid = P.readFilesCommandPID[imate];
        if (pid <= 0) {
            continue;
        }

        int status = 0;
        pid_t wpid = waitpid(pid, &status, WNOHANG);
        if (wpid == pid) {
            P.readFilesCommandPID[imate] = 0;
            continue; // already exited and reaped
        }
        if (wpid == -1 && errno == ECHILD) {
            P.readFilesCommandPID[imate] = 0;
            continue; // not our child anymore
        }

        if (kill(pid, SIGKILL) == -1 && errno != ESRCH) {
            P.readFilesCommandPID[imate] = 0;
            continue;
        }

        while (waitpid(pid, &status, 0) == -1 && errno == EINTR) {
        }
        P.readFilesCommandPID[imate] = 0;
    }
}
}

void Parameters::closeReadsFiles() {
    if (fastxMateReaders) {
        // Stop the mate readers while their streams are still open. Helper
        // children are ended first, so a reader waiting on a FIFO sees EOF.
        terminateReadCommandChildren(*this);
        fastxMateReaders->stopAndJoin();
        if (fastxMateReaders->started()) {
            inOut->logMain << fastxMateReaders->summary();
        }
        fastxMateReaders.reset();
    }

    bgzfCoreInputAdapter.reset();
    bgzfCoreExhausted = false;
    bgzfCoreLaneIndex = 0;

    if (fastxInputActive && fastxInputModule) {
        fastxInputModule->close();
        fastxInputPendingRecordValid = false;
        fastxInputExhausted = false;
        fastxInputPendingRecord.reset();
        fastxInputLastLoggedLane = -1;
    }

    if (cbqInputActive && cbqInputModule) {
        cbqInputModule->close();
        cbqInputExhausted = false;
        cbqInputLastLoggedLane = -1;
        cbqInputPendingBatch.reset();
        cbqInputPendingBatchOffset = 0;
    }

    // Close all potential read streams (not just readFilesIn.size()).
    for (uint imate=0; imate<MAX_N_MATES; imate++) {
        if (inOut->readIn[imate].is_open()) {
            inOut->readIn[imate].close();
        }
    }

    if (bgzfPipes) {
        bgzfPipes->join();
        inOut->logMain << bgzfPipes->summary();
        const string error = bgzfPipes->error();
        bgzfPipes.reset();
        if (!error.empty()) exitWithError("BGZF input failed: " + error + "\n", std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
    }

    terminateReadCommandChildren(*this);
};
