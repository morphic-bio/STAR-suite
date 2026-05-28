#include "Parameters.h"
#include "ErrorWarning.h"
#include "input/CbqInputModule.h"
#include "input/FastxInputModule.h"
#include <fstream>
#include <sys/stat.h>
#include <cerrno>
#include <csignal>
#include <sys/wait.h>

void Parameters::closeReadsFiles() {
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
    }

    // Close all potential read streams (not just readFilesIn.size()).
    for (uint imate=0; imate<MAX_N_MATES; imate++) {
        if (inOut->readIn[imate].is_open()) {
            inOut->readIn[imate].close();
        }
    }

    // Terminate and reap readFilesCommand helper children to avoid lingering
    // processes at shutdown.
    for (uint imate=0; imate<MAX_N_MATES; imate++) {
        pid_t pid = readFilesCommandPID[imate];
        if (pid <= 0) {
            continue;
        }

        int status = 0;
        pid_t wpid = waitpid(pid, &status, WNOHANG);
        if (wpid == pid) {
            readFilesCommandPID[imate] = 0;
            continue; // already exited and reaped
        }
        if (wpid == -1 && errno == ECHILD) {
            readFilesCommandPID[imate] = 0;
            continue; // not our child anymore
        }

        if (kill(pid, SIGKILL) == -1 && errno != ESRCH) {
            readFilesCommandPID[imate] = 0;
            continue;
        }

        while (waitpid(pid, &status, 0) == -1 && errno == EINTR) {
        }
        readFilesCommandPID[imate] = 0;
    }
};
