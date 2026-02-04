/*
    functions that handle errors and warnings
*/
#include "ErrorWarning.h"
#include "TimeFunctions.h"
#include "GlobalVariables.h"
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

namespace {

std::string deriveBaseDir(const Parameters &P) {
    if (!P.outFileNamePrefixAutoRoot.empty()) {
        return P.outFileNamePrefixAutoRoot;
    }
    if (!P.outFileNamePrefix.empty()) {
        size_t pos = P.outFileNamePrefix.rfind('/');
        if (pos != std::string::npos) {
            return P.outFileNamePrefix.substr(0, pos + 1);
        }
    }
    return std::string("./");
}

std::string respawnScriptPath(const Parameters &P, int nextRetry) {
    if (!P.batchRespawnScript.empty() && P.batchRespawnScript != "-") {
        return P.batchRespawnScript;
    }
    std::ostringstream path;
    path << deriveBaseDir(P) << "batch_respawn_retry" << nextRetry << ".sh";
    return path.str();
}

}

void exitWithError(string messageOut, ostream &streamOut1, ostream &streamOut2, int errorInt, Parameters &P)
{
    if (P.runThreadN>1) 
        pthread_mutex_lock(&g_threadChunks.mutexError);
    time_t timeCurrent;
    time( &timeCurrent);
    if (streamOut1.good()) {
        streamOut1 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };
    if (streamOut2.good()) {
        streamOut2 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };

    // Batch-mode respawn (optional)
    if (P.batchMode && P.batchOnError == "respawn") {
        if (P.batchRetryCount < P.batchMaxRetries && P.batchResumeHasList && !P.batchResumeFastqListR1.empty()) {
            const int nextRetry = P.batchRetryCount + 1;
            std::ostringstream cmd;
            cmd << P.commandLineFull
                << " --readFilesIn \"" << P.batchResumeFastqListR1 << "\"";
            if (P.batchPaired && !P.batchResumeFastqListR2.empty()) {
                cmd << " \"" << P.batchResumeFastqListR2 << "\"";
            }
            cmd << " --batchRetryCount " << nextRetry;

            const std::string scriptPath = respawnScriptPath(P, nextRetry);
            std::ofstream script(scriptPath.c_str());
            if (script.good()) {
                script << "#!/usr/bin/env bash\n";
                script << "set -euo pipefail\n";
                script << "while kill -0 " << getpid() << " 2>/dev/null; do sleep 2; done\n";
                script << "exec " << cmd.str() << "\n";
                script.close();
                ::chmod(scriptPath.c_str(), 0755);

                if (streamOut2.good()) {
                    streamOut2 << "BATCH RESPAWN: retry " << nextRetry
                               << " of " << P.batchMaxRetries
                               << " using script " << scriptPath << "\n";
                    streamOut2 << "BATCH RESPAWN CMD: " << cmd.str() << "\n";
                }
                std::ostringstream sh;
                sh << "/bin/bash " << scriptPath << " &";
                std::system(sh.str().c_str());
            } else if (streamOut2.good()) {
                streamOut2 << "BATCH RESPAWN: failed to write script " << scriptPath << "\n";
            }
        } else if (streamOut2.good()) {
            streamOut2 << "BATCH RESPAWN: not attempted (retryCount=" << P.batchRetryCount
                       << " maxRetries=" << P.batchMaxRetries << ")\n";
        }
    }

    delete P.inOut; //to close files
//     if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexError);
    exit(errorInt);
};

void warningMessage(string messageOut, ostream &streamOut1, ostream &streamOut2, Parameters &P)
{
    if (P.runThreadN>1) 
        pthread_mutex_lock(&g_threadChunks.mutexError);
    if (streamOut1.good()) {
        streamOut1 << "!!!!! WARNING: " << messageOut <<endl;
    };
    if (streamOut2.good()) {
        streamOut2 << "!!!!! WARNING: " << messageOut <<endl;
    };
    if (P.runThreadN>1) 
        pthread_mutex_unlock(&g_threadChunks.mutexError);    
};
