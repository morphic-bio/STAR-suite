#include "Parameters.h"
#include "ErrorWarning.h"
#include "input/CbqInputModule.h"
#include "input/FastxInputModule.h"
#include <fstream>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cctype>
#include <algorithm>
#include <limits>
#include <zlib.h>

namespace {
bool endsWithCaseInsensitiveLocal(const string& value, const string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    const size_t offset = value.size() - suffix.size();
    for (size_t ii = 0; ii < suffix.size(); ++ii) {
        const unsigned char c1 = static_cast<unsigned char>(value[offset + ii]);
        const unsigned char c2 = static_cast<unsigned char>(suffix[ii]);
        if (std::tolower(c1) != std::tolower(c2)) {
            return false;
        }
    }
    return true;
}

bool isFastqPathLocal(const string& value) {
    return endsWithCaseInsensitiveLocal(value, ".fastq") ||
           endsWithCaseInsensitiveLocal(value, ".fq") ||
           endsWithCaseInsensitiveLocal(value, ".fastq.gz") ||
           endsWithCaseInsensitiveLocal(value, ".fq.gz");
}

bool writeAll(int fd, const char* data, size_t size) {
    while (size > 0) {
        const ssize_t written = ::write(fd, data, size);
        if (written < 0) {
            if (errno == EINTR) {
                continue;
            }
            return false;
        }
        data += written;
        size -= static_cast<size_t>(written);
    }
    return true;
}

bool streamGzipFileToFd(const string& gzipPath, int outFd, string& errorOut) {
    gzFile gz = gzopen(gzipPath.c_str(), "rb");
    if (gz == nullptr) {
        errorOut = "could not open gzip file: " + gzipPath;
        return false;
    }

    char buffer[1 << 20];
    int nRead = 0;
    while ((nRead = gzread(gz, buffer, sizeof(buffer))) > 0) {
        if (!writeAll(outFd, buffer, static_cast<size_t>(nRead))) {
            errorOut = "write() failed while streaming gzip input: " + gzipPath + ": " + strerror(errno);
            gzclose(gz);
            return false;
        }
    }

    if (nRead < 0) {
        int zerr = Z_OK;
        const char* zmsg = gzerror(gz, &zerr);
        errorOut = "gzread failed for " + gzipPath + ": " + (zmsg != nullptr ? string(zmsg) : string("unknown zlib error"));
        gzclose(gz);
        return false;
    }

    if (gzclose(gz) != Z_OK) {
        errorOut = "gzclose failed for " + gzipPath;
        return false;
    }

    return true;
}

[[noreturn]] void streamGzipMateToFifo(const vector<string>& mateFiles, const string& fifoPath) {
    const int outFd = open(fifoPath.c_str(), O_WRONLY);
    if (outFd < 0) {
        cerr << "EXITING: internal gzip helper could not open FIFO for write: " << fifoPath
             << " (" << strerror(errno) << ")\n";
        _exit(1);
    }

    for (uint32 ifile = 0; ifile < mateFiles.size(); ++ifile) {
        const string marker = "FILE " + to_string(ifile) + "\n";
        if (!writeAll(outFd, marker.c_str(), marker.size())) {
            cerr << "EXITING: internal gzip helper failed writing FILE marker for " << mateFiles[ifile]
                 << " (" << strerror(errno) << ")\n";
            close(outFd);
            _exit(1);
        }

        string errMsg;
        if (!streamGzipFileToFd(mateFiles[ifile], outFd, errMsg)) {
            cerr << "EXITING: internal gzip helper failed: " << errMsg << "\n";
            close(outFd);
            _exit(1);
        }
    }

    close(outFd);
    _exit(0);
}

string lowerCopyLocal(string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

bool cbqCoreRangeGateReject(Parameters& P, string& reason) {
    auto reject = [&](const string& message) {
        reason = message;
        return true;
    };

    const string mode = lowerCopyLocal(P.readFilesCbqRangeMode);
    if (mode == "off") {
        return reject("disabled by --readFilesCbqRangeMode off");
    }
    if (mode != "auto" && mode != "range") {
        return reject("unrecognized --readFilesCbqRangeMode=" + P.readFilesCbqRangeMode);
    }
    P.readFilesCbqRangeMode = mode;

    if (P.readFilesTypeN != 20 || !P.cbqInputActive) {
        return reject("input is not active Binseq/CBQ");
    }
    if (P.runThreadN <= 1) {
        return reject("runThreadN <= 1");
    }
    if (P.outSAMtype.empty()) {
        return reject("outSAMtype is empty");
    }
    const bool noAlignmentOutput = (P.outSAMtype[0] == "None" &&
                                    !P.outSAMbool && !P.outBAMcoord && !P.outBAMunsorted);
    const bool sortedBamOnly = (P.outSAMtype[0] == "BAM" &&
                                P.outBAMcoord && !P.outBAMunsorted && !P.outSAMbool);
    if (!noAlignmentOutput && !sortedBamOnly) {
        return reject("CBQ range mode currently supports --outSAMtype None or BAM SortedByCoordinate without SAM/Unsorted side output");
    }
    if (P.emitYNoYyes || P.emitYNoYFastqyes || P.emitYNoYCbqyes) {
        return reject("Y/noY sidecar emission is order-dependent");
    }
    if (P.outSAMorder == "PairedKeepInputOrder") {
        return reject("--outSAMorder PairedKeepInputOrder is order-dependent");
    }
    if (P.batchMode || P.batchModeInt != 0 || P.quant.slam.batchMode || P.quant.slam.batchModeInt != 0) {
        return reject("batch mode is not supported by initial CBQ range mode");
    }
    if (P.quant.slam.autoTrimDetectionPass || P.quant.slam.perFileProcessing || P.quant.slam.skipToFileIndex > 0) {
        return reject("SLAM per-file or auto-trim pass is not supported by initial CBQ range mode");
    }
    if (P.twoPass.yes || P.sjdbInsert.yes) {
        return reject("two-pass or run-time SJ insertion is not supported by initial CBQ range mode");
    }
    if (P.outFilterBySJoutStage != 0 || P.outFilterType == "BySJout") {
        return reject("two-stage SJ filtering is not supported by initial CBQ range mode");
    }
    if (P.readFilesNames.empty() || P.readFilesNames.front().empty()) {
        return reject("no CBQ input lanes");
    }
    reason.clear();
    return false;
}

bool prepareCbqCoreRangeTasks(Parameters& P,
                              const star::input::InputSourcePlan& cbqInputPlan,
                              string& reason) {
    P.cbqRangeTasks.clear();
    P.cbqRangeActive = false;
    P.cbqRangeFallbackReason.clear();
    P.cbqRangeTotalRecords = 0;
    if (!P.cbqRangeNextTask) {
        P.cbqRangeNextTask.reset(new std::atomic<uint32_t>(0));
    }
    if (!P.cbqRangeNextChunk) {
        P.cbqRangeNextChunk.reset(new std::atomic<uint32_t>(0));
    }
    P.cbqRangeNextTask->store(0);
    P.cbqRangeNextChunk->store(0);

    if (cbqCoreRangeGateReject(P, reason)) {
        return false;
    }

    const uint32 laneCount = static_cast<uint32>(cbqInputPlan.mate_files.front().size());
    vector<uint64> laneCounts(laneCount, 0);
    vector<uint64> laneStarts(laneCount + 1, 0);
    for (uint32 lane = 0; lane < laneCount; ++lane) {
        star::input::CbqInputModule laneReader;
        string inputError;
        if (!laneReader.configure(cbqInputPlan, &inputError) ||
            !laneReader.open_range(lane, 0, std::numeric_limits<uint64>::max(), &inputError)) {
            reason = inputError.empty()
                ? ("could not open indexed CBQ range for lane " + to_string(lane))
                : inputError;
            return false;
        }
        laneCounts[lane] = laneReader.current_lane_record_count();
        laneStarts[lane + 1] = laneStarts[lane] + laneCounts[lane];
        laneReader.close();
    }

    const uint64 inputRecords = laneStarts.back();
    const uint64 totalRecords = P.readMapNumber == static_cast<uint>(-1)
        ? inputRecords
        : std::min<uint64>(inputRecords, static_cast<uint64>(P.readMapNumber));
    if (totalRecords == 0) {
        reason = "CBQ input contains no records";
        return false;
    }
    P.cbqRangeTotalRecords = totalRecords;

    const uint64 targetRangeCount = std::min<uint64>(static_cast<uint64>(P.runThreadN), totalRecords);
    const uint64 targetRecordsPerRange = (totalRecords + targetRangeCount - 1) / targetRangeCount;
    uint64 taskOrdinal = 0;
    uint64 globalStart = 0;
    while (globalStart < totalRecords) {
        const uint64 globalEnd = std::min<uint64>(totalRecords, globalStart + targetRecordsPerRange);
        uint64 cursor = globalStart;
        while (cursor < globalEnd) {
            const auto laneIt = std::upper_bound(laneStarts.begin(), laneStarts.end(), cursor);
            if (laneIt == laneStarts.begin()) {
                reason = "internal CBQ range planning error: invalid lane start";
                return false;
            }
            uint32 lane = static_cast<uint32>(std::distance(laneStarts.begin(), laneIt) - 1);
            if (lane >= laneCount) {
                reason = "internal CBQ range planning error: lane exceeds lane count";
                return false;
            }
            const uint64 laneEndGlobal = laneStarts[lane] + laneCounts[lane];
            const uint64 take = std::min<uint64>(globalEnd, laneEndGlobal) - cursor;
            if (take == 0) {
                reason = "internal CBQ range planning error: empty task";
                return false;
            }

            Parameters::CbqRangeTask task;
            task.laneIndex = lane;
            task.firstRecord = cursor - laneStarts[lane];
            task.recordCount = take;
            task.globalFirst = cursor;
            task.taskOrdinal = taskOrdinal++;
            P.cbqRangeTasks.push_back(task);
            cursor += take;
        }
        globalStart = globalEnd;
    }

    if (P.cbqRangeTasks.empty()) {
        reason = "CBQ range planner produced no tasks";
        return false;
    }
    P.cbqRangeActive = true;
    reason.clear();
    return true;
}

void fatalCbqRangeMode(Parameters& P, const string& reason) {
    ostringstream errOut;
    errOut << "EXITING because of fatal input ERROR: --readFilesCbqRangeMode range could not be activated.\n";
    errOut << reason << "\n";
    errOut << "SOLUTION: use indexed CBQ inputs with supported order-independent settings, or set --readFilesCbqRangeMode auto/off.\n";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
}
}

void Parameters::openReadsFiles() 
{
    // Reset FIFO list to avoid stale entries when reopening (e.g. SLAM auto-trim detection pass)
    readFilesInTmp.clear();
    // Check number of mates BEFORE opening files
    // Use readNends if available (set during readFilesInit), otherwise count readFilesIn
    uint readFilesNmates = readNends;
    if (readFilesNmates == 0) {
        // Fallback: count from readFilesIn or readFilesNames
        readFilesNmates = (readFilesCommandString=="") ? readFilesIn.size() : readFilesNames.size();
    }
    if (readFilesNmates > MAX_N_MATES) {
        ostringstream errOut;
        errOut << "EXITING: because of fatal INPUT error: number of read mates (" << readFilesNmates
               << ") exceeds MAX_N_MATES=" << MAX_N_MATES << "\n";
        errOut << "SOLUTION: --readFilesIn expects 1-3 read mates (R1/R2/optional barcode read). "
               << "Multiple lanes must be comma-separated per mate, e.g. "
               << "--readFilesIn R2_L1.fastq,R2_L2.fastq R1_L1.fastq,R1_L2.fastq.\n"
               << "Do not include index reads (I1/I2) in --readFilesIn.\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (readFilesTypeN == 1 && fastxInputActive) {
        for (uint imate = 0; imate < MAX_N_MATES; imate++) {
            readFilesCommandPID[imate] = 0;
            if (inOut->readIn[imate].is_open()) {
                inOut->readIn[imate].close();
            }
        }

        if (emitYNoYFastqyes) {
            for (uint32 imate = 0; imate < readFilesNames.size(); ++imate) {
                for (const auto& fastxName : readFilesNames[imate]) {
                    if (!isFastqPathLocal(fastxName)) {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL INPUT ERROR: --emitYNoYFastq currently requires FASTQ input files.\n";
                        errOut << "Offending --readFilesIn entry: " << fastxName << "\n";
                        errOut << "SOLUTION: provide .fastq/.fq files, optionally gzip-compressed as .fastq.gz/.fq.gz, or disable --emitYNoYFastq.\n";
                        // TODO: Y-removal/Y-noY FASTQ emission is needed for other input formats.
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    }
                }
            }
        }

        vector<string> fastxReadGroups;
        if (!readFilesNames.empty() && outSAMattrRG.size() == readFilesNames.front().size()) {
            fastxReadGroups = outSAMattrRG;
        }
        star::input::InputSourcePlan fastxInputPlan =
            star::input::make_fastx_input_source_plan(
                readFilesNames,
                fastxReadGroups,
                readFilesCommandString,
                readFilesPrefixFinal,
                readFilesUseInternalGzip);

        string inputContractError;
        fastxInputModule.reset(new star::input::FastxInputModule());
        if (!fastxInputModule->configure(fastxInputPlan, &inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal input ERROR: invalid Fastx input source plan at open\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        if (!fastxInputModule->open(&inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal input ERROR: could not open Fastx input module\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
        }
        readFilesIndex = 0;
        fastxInputPendingRecordValid = false;
        fastxInputExhausted = false;
        fastxInputPendingRecord.reset();
        fastxInputLastLoggedLane = -1;
        return;
    }

    if (readFilesTypeN == 20 && cbqInputActive) {
        for (uint imate = 0; imate < MAX_N_MATES; imate++) {
            readFilesCommandPID[imate] = 0;
            if (inOut->readIn[imate].is_open()) {
                inOut->readIn[imate].close();
            }
        }

        vector<string> cbqReadGroups;
        if (!readFilesNames.empty() && outSAMattrRG.size() == readFilesNames.front().size()) {
            cbqReadGroups = outSAMattrRG;
        }
        star::input::InputSourcePlan cbqInputPlan =
            star::input::make_cbq_input_source_plan(
                readFilesNames,
                cbqReadGroups,
                readNends);

        string inputContractError;
        cbqInputModule.reset(new star::input::CbqInputModule());
        if (!cbqInputModule->configure(cbqInputPlan, &inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal input ERROR: invalid Binseq/CBQ input source plan at open\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        string cbqRangeReason;
        const bool cbqRangePrepared = prepareCbqCoreRangeTasks(*this, cbqInputPlan, cbqRangeReason);
        const string cbqRangeMode = lowerCopyLocal(readFilesCbqRangeMode);
        if (cbqRangePrepared) {
            inOut->logMain << "CBQ indexed range reader: active with "
                           << cbqRangeTasks.size() << " tasks across "
                           << readFilesN << " lanes and "
                           << cbqRangeTotalRecords << " records\n";
            readFilesIndex = 0;
            cbqInputExhausted = false;
            cbqInputLastLoggedLane = -1;
            cbqInputPendingBatch.reset();
            cbqInputPendingBatchOffset = 0;
            return;
        }
        cbqRangeFallbackReason = cbqRangeReason;
        if (cbqRangeMode == "range" || (cbqRangeMode != "auto" && cbqRangeMode != "off")) {
            fatalCbqRangeMode(*this, cbqRangeReason);
        }
        if (cbqRangeMode == "auto") {
            inOut->logMain << "CBQ indexed range reader: not active ("
                           << (cbqRangeReason.empty() ? "gate rejected command" : cbqRangeReason)
                           << "); using shared CBQ reader\n";
        }
        if (!cbqInputModule->open(&inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal input ERROR: could not open Binseq/CBQ input module\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
        }
        readFilesIndex = 0;
        cbqInputExhausted = false;
        cbqInputLastLoggedLane = -1;
        cbqInputPendingBatch.reset();
        cbqInputPendingBatchOffset = 0;
        return;
    }

    if (readFilesCommandString=="") {//read from file
        for (uint ii=0;ii<readFilesIn.size();ii++) {//open readIn files
            readFilesCommandPID[ii]=0;//no command process IDs
            if ( inOut->readIn[ii].is_open() ) inOut->readIn[ii].close();

            string rfName=readFilesPrefixFinal + readFilesIn.at(ii);

            inOut->readIn[ii].open(rfName.c_str()); //try to open the Sequences file right away, exit if failed
            if (inOut->readIn[ii].fail()) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal input ERROR: could not open readFilesIn=" << rfName <<"\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    } else {//create fifo files, execute pre-processing command

         vector<string> readsCommandFileName;
         const bool batchModeRequested = (batchMode || batchModeInt != 0 || quant.slam.batchModeInt != 0);
             const bool filterBrokenPipe = (batchModeRequested &&
                                        (readFilesCommandString.find("zcat") != std::string::npos ||
                                         readFilesCommandString.find("gzcat") != std::string::npos));
         uint imate;
         for (imate=0;imate<readFilesNames.size();imate++) {//open readIn files
            ostringstream sysCom;
            sysCom << outFileTmp <<"tmp.fifo.read"<<imate+1;
            readFilesInTmp.push_back(sysCom.str());
            remove(readFilesInTmp.at(imate).c_str());
            if (mkfifo(readFilesInTmp.at(imate).c_str(), S_IRUSR | S_IWUSR ) != 0) {
                exitWithError("Exiting because of *FATAL ERROR*: could not create FIFO file " + readFilesInTmp.at(imate) + "\n"
                            + "SOLUTION: check the if run directory supports FIFO files.\n"
                            + "If run partition does not support FIFO (e.g. Windows partitions FAT, NTFS), "
                            + "re-run on a Linux partition, or point --outTmpDir to a Linux partition.\n"
                            , std::cerr, inOut->logMain, EXIT_CODE_FIFO, *this);
            };

            inOut->logMain << "\n   Input read files for mate "<< imate+1 <<" :\n";

            // Log and validate input files before spawning helper path.
            for (uint32 ifile = 0; ifile < readFilesN; ifile++) {
                if (system(("ls -lL " + readFilesNames[imate][ifile] + " > " + outFileTmp + "/readFilesIn.info 2>&1").c_str()) != 0) {
                    warningMessage(" Could not ls " + readFilesNames[imate][ifile], std::cerr, inOut->logMain, *this);
                }

                ifstream readFilesIn_info((outFileTmp + "/readFilesIn.info").c_str());
                inOut->logMain << readFilesIn_info.rdbuf();

                // Try to open files early and fail with actionable message.
                ifstream rftry(readFilesNames[imate][ifile].c_str());
                if (!rftry.good()) {
                    exitWithError("EXITING: because of fatal INPUT file error: could not open read file: " +
                                   readFilesNames[imate][ifile] +
                                   "\nSOLUTION: check that this file exists and has read permision.\n",
                                   std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
                rftry.close();
            }

            if (readFilesUseInternalGzip) {
                inOut->logMain << "NOTE: mate " << (imate + 1)
                               << " using internal gzip FIFO helper: " << readFilesInTmp.at(imate) << "\n";
                readFilesCommandPID[imate] = 0;
                ostringstream errOut;
                const pid_t pid = fork();
                switch (pid) {
                    case -1:
                        errOut << "EXITING: because of fatal EXECUTION error: Failed forking internal gzip helper\n";
                        errOut << errno << ": " << strerror(errno) << "\n";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                        break;

                    case 0:
                        // Child: stream .gz FASTQ to FIFO using in-process zlib (no external zcat).
                        streamGzipMateToFifo(readFilesNames.at(imate), readFilesInTmp.at(imate));
                        break;

                    default:
                        readFilesCommandPID[imate] = pid;
                }

                inOut->readIn[imate].open(readFilesInTmp.at(imate).c_str());
                if (inOut->readIn[imate].fail()) {
                    ostringstream errOpen;
                    errOpen << "EXITING because of fatal input ERROR: could not open FIFO stream " << readFilesInTmp.at(imate) << "\n";
                    exitWithError(errOpen.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
                continue;
            }

            readsCommandFileName.push_back(outFileTmp+"/readsCommand_read" + to_string(imate+1));
            fstream readsCommandFile( readsCommandFileName.at(imate).c_str(), ios::out);
            readsCommandFile.close();
            readsCommandFile.open( readsCommandFileName.at(imate).c_str(), ios::in | ios::out);
            //first line in the commands file
            if (filterBrokenPipe) {
                if (sysShell == "-" || sysShell.find("bash") != std::string::npos) {
                    readsCommandFile << "#!/bin/bash\n";
                } else {
                    inOut->logMain << "WARNING: batch mode zcat filter requested, but sysShell is not bash ("
                                   << sysShell << "); skipping Broken pipe filtering.\n";
                }
            } else if (sysShell!="-") {//executed via specified shell
                readsCommandFile << "#!" <<sysShell <<"\n";
            };
            readsCommandFile << "exec > \""<<readFilesInTmp.at(imate)<<"\"\n" ; // redirect stdout to temp fifo files
            for (uint32 ifile=0; ifile<readFilesN; ifile++) {
                readsCommandFile <<"echo FILE "<< ifile << "\n";
                if (filterBrokenPipe && (sysShell == "-" || sysShell.find("bash") != std::string::npos)) {
                    readsCommandFile << readFilesCommandString <<"   "<< ("\""+readFilesNames[imate][ifile]+"\"")
                                     << " 2> >(grep -v \"Broken pipe\" >&2)\n";
                } else {
                    readsCommandFile << readFilesCommandString <<"   "<< ("\""+readFilesNames[imate][ifile]+"\"") <<"\n";
                }
            };

            readsCommandFile.flush();
            readsCommandFile.seekg(0,ios::beg);
            inOut->logMain <<"\n   readsCommandsFile:\n"<<readsCommandFile.rdbuf()<<endl;
            readsCommandFile.close();
            if (filterBrokenPipe) {
                inOut->logMain << "NOTE: batch mode detected zcat; wrapping readFilesCommand with bash to filter benign "
                               << "\"Broken pipe\" warnings.\n";
            }

            chmod(readsCommandFileName.at(imate).c_str(),S_IXUSR | S_IRUSR | S_IWUSR);

            readFilesCommandPID[imate]=0;

            ostringstream errOut;
            pid_t PID=vfork();
            switch (PID) {
                case -1:
                    errOut << "EXITING: because of fatal EXECUTION error: Failed vforking readFilesCommand\n";
                    errOut << errno << ": " << strerror(errno) << "\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    break;

                case 0:
                    //this is the child
                    execlp(readsCommandFileName.at(imate).c_str(), readsCommandFileName.at(imate).c_str(), (char*) NULL);
                    exit(0);

                default:
                    //this is the father, record PID of the children
                    readFilesCommandPID[imate]=PID;
            };

            inOut->readIn[imate].open(readFilesInTmp.at(imate).c_str());
        };

    };
    readFilesIndex=0;

    if (readFilesTypeN==10) {//SAM file - skip header lines
        readSAMheader(readFilesCommandString, readFilesNames.at(0));
    };
 
};
