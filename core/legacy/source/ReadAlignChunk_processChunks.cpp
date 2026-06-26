#include "ReadAlignChunk.h"
#include "ThreadControl.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "GlobalVariables.h"
#include "FlexDebugCounters.h"
#include "IncludeDefine.h"
#include "input/CbqInputModule.h"
#include "input/FastxInputModule.h"
#include "input/CbqStarAdapter.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <pthread.h>
#include <sstream>

inline uint64 fastqReadOneLine(ifstream &streamIn, char *arrIn);
inline void removeStringEndControl(string &str);

namespace {

bool fastxChunkHasBytes(const ReadAlignChunk& chunk) {
    for (uint imate = 0; imate < chunk.P.readNends; ++imate) {
        if (chunk.chunkInSizeBytesTotal[imate] > 0) {
            return true;
        }
    }
    return false;
}

bool fastxChunkUnderTarget(const ReadAlignChunk& chunk) {
    for (uint imate = 0; imate < chunk.P.readNends; ++imate) {
        if (chunk.chunkInSizeBytesTotal[imate] >= chunk.P.chunkInSizeBytes) {
            return false;
        }
    }
    return true;
}

uint64 decimalDigitCount(uint64 value) {
    uint64 digits = 1;
    while (value >= 10) {
        value /= 10;
        ++digits;
    }
    return digits;
}

string fastxHeaderLine(const Parameters& P,
                       const star::input::InputRecord& record,
                       const star::input::InputMateRecord& mate,
                       uint64 ordinal) {
    const char headerPrefix = mate.has_quality ? '@' : '>';
    string header;
    if (P.outSAMreadIDnumber) {
        header = string(1, headerPrefix) + to_string(ordinal);
    } else {
        header = string(1, headerPrefix) + record.read_name;
    }
    header += ' ' + to_string(ordinal) + ' ' + record.read_filter + ' ' + to_string(record.lane_index);
    return header;
}

uint64 fastxMateChunkBytes(const Parameters& P,
                           const star::input::InputRecord& record,
                           const star::input::InputMateRecord& mate,
                           uint64 ordinal) {
    uint64 nBytes = fastxHeaderLine(P, record, mate, ordinal).size() + 1;
    nBytes += mate.sequence.size() + 1;
    if (mate.has_quality) {
        nBytes += 2; // "+\n"
        nBytes += mate.quality.size() + 1;
    }
    return nBytes;
}

uint64 cbqFastqHeaderLineBytes(const Parameters& P,
                               const star::input::CbqReadView& record,
                               uint64 ordinal) {
    const uint64 readNameCoreBytes = P.outSAMreadIDnumber
        ? decimalDigitCount(ordinal)
        : static_cast<uint64>(record.read_name.size);

    uint64 bytes = 1; // '@' or '>'
    bytes += readNameCoreBytes;
    bytes += 1; // space before ordinal
    bytes += decimalDigitCount(ordinal);
    bytes += 1; // space before read filter
    bytes += 1; // read filter character
    bytes += 1; // space before lane index
    bytes += decimalDigitCount(record.lane_index);
    return bytes;
}

uint64 cbqMateFastqChunkBytes(const Parameters& P,
                              const star::input::CbqReadView& record,
                              uint32 imate,
                              uint64 ordinal) {
    const star::input::CbqSegmentView& mate = record.segments[imate];
    const uint64 sequenceLength =
        static_cast<uint64>(star::input::cbq_segment_sequence_length(mate));

    uint64 bytes = cbqFastqHeaderLineBytes(P, record, ordinal) + 1; // header newline
    bytes += sequenceLength + 1; // sequence newline
    if (mate.has_quality) {
        bytes += 2; // "+\n"
        bytes += static_cast<uint64>(mate.quality.size) + 1; // quality newline
    }
    return bytes;
}

const char* inputChunkTracePath() {
    static const char* path = std::getenv("STAR_INPUT_CHUNK_TRACE");
    return (path != nullptr && path[0] != '\0') ? path : nullptr;
}

const char* inputChunkTraceSource(const Parameters& P) {
    if (P.readFilesTypeN == 20 && P.cbqInputActive) {
        return "cbq";
    }
    if (P.readFilesTypeN == 1 && P.fastxInputActive) {
        return "fastx";
    }
    return "legacy";
}

void writeInputChunkTrace(const ReadAlignChunk& chunk,
                          uint64 chunkIndex,
                          uint64 chunkReadStart,
                          uint64 chunkReadN,
                          uint64 chunkWorkBytes,
                          bool noReadsLeft,
                          int readFilesIndexOverride = -1) {
    const char* path = inputChunkTracePath();
    if (path == nullptr) {
        return;
    }

    static pthread_mutex_t traceMutex = PTHREAD_MUTEX_INITIALIZER;
    static bool headerWritten = false;

    pthread_mutex_lock(&traceMutex);
    std::ofstream out(path, std::ios::app);
    if (out.good()) {
        if (!headerWritten) {
            out << "chunk_index\tthread\tsource\tread_start\tread_end\tread_count"
                << "\tmate1_bytes\tmate2_bytes\twork_bytes\tread_files_index\tno_reads_left\n";
            headerWritten = true;
        }
        const uint64 readStart = chunkReadN == 0 ? chunkReadStart : chunkReadStart + 1;
        const uint64 readEnd = chunkReadStart + chunkReadN;
        const uint64 mate1Bytes = chunk.chunkInSizeBytesTotal[0];
        const uint64 mate2Bytes = chunk.P.readNends > 1 ? chunk.chunkInSizeBytesTotal[1] : 0;
        const int readFilesIndex = readFilesIndexOverride >= 0
            ? readFilesIndexOverride
            : chunk.P.readFilesIndex;
        out << chunkIndex << '\t'
            << chunk.iThread << '\t'
            << inputChunkTraceSource(chunk.P) << '\t'
            << readStart << '\t'
            << readEnd << '\t'
            << chunkReadN << '\t'
            << mate1Bytes << '\t'
            << mate2Bytes << '\t'
            << chunkWorkBytes << '\t'
            << readFilesIndex << '\t'
            << (noReadsLeft ? 1 : 0) << '\n';
    }
    pthread_mutex_unlock(&traceMutex);
}

bool cbqRecordFits(ReadAlignChunk& chunk,
                   const star::input::CbqReadView& record,
                   uint64 ordinal,
                   std::array<uint64, MAX_N_MATES>* mateBytes,
                   string& errorOut) {
    if (record.segment_count < chunk.P.readNends || record.segments == nullptr) {
        errorOut = "CBQ input module returned fewer segments than STAR expects";
        return false;
    }

    mateBytes->fill(0);
    for (uint imate = 0; imate < chunk.P.readNends; ++imate) {
        const uint64 bytes = cbqMateFastqChunkBytes(chunk.P, record, imate, ordinal);
        (*mateBytes)[imate] = bytes;
        if (bytes + 2 >= chunk.P.chunkInSizeBytesArray) {
            ostringstream err;
            err << "CBQ record ";
            if (record.read_name.data != nullptr && record.read_name.size > 0) {
                err.write(record.read_name.data, static_cast<std::streamsize>(record.read_name.size));
            } else {
                err << ordinal;
            }
            err << " mate " << (imate + 1)
                << " requires " << bytes << " FASTQ-equivalent bytes, larger than STAR input chunk buffer "
                << chunk.P.chunkInSizeBytesArray;
            errorOut = err.str();
            return false;
        }
        if (chunk.chunkInSizeBytesTotal[imate] + bytes + 2 >= chunk.P.chunkInSizeBytesArray) {
            return false;
        }
    }

    return true;
}

bool fastxRecordFits(ReadAlignChunk& chunk,
                     const star::input::InputRecord& record,
                     uint64 ordinal,
                     string& errorOut) {
    if (record.mates.size() < chunk.P.readNends) {
        errorOut = "Fastx input module returned fewer mates than STAR expects";
        return false;
    }

    for (uint imate = 0; imate < chunk.P.readNends; ++imate) {
        const uint64 mateBytes = fastxMateChunkBytes(chunk.P, record, record.mates[imate], ordinal);
        if (mateBytes + 2 >= chunk.P.chunkInSizeBytesArray) {
            ostringstream err;
            err << "Fastx record " << record.read_name << " mate " << (imate + 1)
                << " requires " << mateBytes << " bytes, larger than STAR input chunk buffer "
                << chunk.P.chunkInSizeBytesArray;
            errorOut = err.str();
            return false;
        }
        if (chunk.chunkInSizeBytesTotal[imate] + mateBytes + 2 >= chunk.P.chunkInSizeBytesArray) {
            return false;
        }
    }

    return true;
}

void fastxAppendLine(char* buffer, uint64& offset, const string& line) {
    if (!line.empty()) {
        std::memcpy(buffer + offset, line.data(), line.size());
        offset += line.size();
    }
    buffer[offset++] = '\n';
}

void fastxAppendRecord(ReadAlignChunk& chunk, const star::input::InputRecord& record) {
    for (uint imate = 0; imate < chunk.P.readNends; ++imate) {
        const star::input::InputMateRecord& mate = record.mates[imate];
        uint64& offset = chunk.chunkInSizeBytesTotal[imate];
        fastxAppendLine(chunk.chunkIn[imate], offset, fastxHeaderLine(chunk.P, record, mate, chunk.P.iReadAll));
        fastxAppendLine(chunk.chunkIn[imate], offset, mate.sequence);
        if (mate.has_quality) {
            chunk.chunkIn[imate][offset++] = '+';
            chunk.chunkIn[imate][offset++] = '\n';
            fastxAppendLine(chunk.chunkIn[imate], offset, mate.quality);
        }
    }
}

void fastxLogLaneStart(Parameters& P, uint32 laneIndex) {
    if (P.fastxInputLastLoggedLane == static_cast<int>(laneIndex)) {
        return;
    }
    P.fastxInputLastLoggedLane = static_cast<int>(laneIndex);

    if (P.readFilesN <= 1 && P.readFilesCommandString.empty()) {
        return;
    }

    pthread_mutex_lock(&g_threadChunks.mutexLogMain);
    P.inOut->logMain << "Starting to map file # " << laneIndex << "\n";
    for (uint imate = 0; imate < P.readFilesNames.size(); imate++) {
        if (laneIndex < P.readFilesNames.at(imate).size()) {
            P.inOut->logMain << "mate " << imate + 1 << ":   " << P.readFilesNames.at(imate).at(laneIndex) << "\n";
        }
    }
    P.inOut->logMain << flush;
    pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
}

bool fastxStopBeforeLane(Parameters& P, uint32 laneIndex, int detectionStartFileIndex) {
    if (P.quant.slam.autoTrimDetectionPass && P.quant.slam.trimScope == "first" &&
        static_cast<int>(laneIndex) > detectionStartFileIndex) {
        P.inOut->logMain << "SLAM auto-trim detection: stopping at file boundary "
                         << "(file " << detectionStartFileIndex << " -> " << laneIndex << ")\n";
        return true;
    }

    if (P.quant.slam.perFileProcessing &&
        static_cast<int>(laneIndex) > P.quant.slam.currentFileIndex) {
        P.inOut->logMain << "SLAM per-file processing: stopping at file boundary "
                         << "(completed file " << P.quant.slam.currentFileIndex << ")\n";
        return true;
    }

    return false;
}

void fastxSkipToTargetLane(ReadAlignChunk& chunk) {
    Parameters& P = chunk.P;
    bool foundTarget = false;

    while (P.quant.slam.skipToFileIndex > 0 && P.readFilesIndex < P.quant.slam.skipToFileIndex) {
        star::input::InputRecord record;
        string inputError;
        const star::input::InputStatus status = P.fastxInputModule->next_record(&record, &inputError);
        if (status == star::input::InputStatus::Error) {
            ostringstream errOut;
            errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Fastx input module while skipping to file "
                   << P.quant.slam.skipToFileIndex << "\n"
                   << inputError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        if (status == star::input::InputStatus::End) {
            P.fastxInputExhausted = true;
            break;
        }

        P.readFilesIndex = static_cast<int>(record.lane_index);
        if (P.readFilesIndex >= P.quant.slam.skipToFileIndex) {
            P.fastxInputPendingRecord.reset(new star::input::InputRecord(record));
            P.fastxInputPendingRecordValid = true;
            foundTarget = true;
            P.inOut->logMain << "SLAM per-file: reached file " << P.readFilesIndex << "\n";
            break;
        }
    }

    if (!foundTarget && P.readFilesIndex < P.quant.slam.skipToFileIndex) {
        P.inOut->logMain << "WARNING: SLAM per-file: reached end of input before finding file "
                         << P.quant.slam.skipToFileIndex << ". Current file index: "
                         << P.readFilesIndex << "\n";
    }

    P.quant.slam.skipToFileIndex = -1;
}

[[noreturn]] void exitCbqRangeError(ReadAlignChunk& chunk,
                                    const string& message,
                                    const string& detail = string()) {
    ostringstream errOut;
    errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Binseq/CBQ indexed range input\n"
           << message << "\n";
    if (!detail.empty()) {
        errOut << detail << "\n";
    }
    exitWithError(errOut.str(), std::cerr, chunk.P.inOut->logMain, EXIT_CODE_INPUT_FILES, chunk.P);
}

void mapPreparedCbqRangeChunk(ReadAlignChunk& chunk,
                              uint64 chunkReadStart,
                              uint64 chunkReadN,
                              uint64 chunkWorkBytes,
                              uint32 laneIndex) {
    Parameters& P = chunk.P;
    chunk.noReadsLeft = false;
    chunk.iChunkIn = P.cbqRangeNextChunk
        ? P.cbqRangeNextChunk->fetch_add(1)
        : 0;

    if (P.runThreadN > 1) {
        pthread_mutex_lock(&g_threadChunks.mutexInRead);
    }
    P.iReadAll += static_cast<uint>(chunkReadN);
    if (P.runThreadN > 1) {
        pthread_mutex_unlock(&g_threadChunks.mutexInRead);
    }

    writeInputChunkTrace(chunk,
                         chunk.iChunkIn,
                         chunkReadStart,
                         chunkReadN,
                         chunkWorkBytes,
                         false,
                         static_cast<int>(laneIndex));

    const bool permitEnabled = g_threadChunks.mapPermitEnabled();
    uint64_t waitNs = 0;
    if (permitEnabled) {
        waitNs = g_threadChunks.mapPermitAcquire();
    }
    const auto workStart = std::chrono::steady_clock::now();
    const uint64_t readCountBefore = chunk.RA->iRead;
    chunk.mapCbqChunk();
    const auto workEnd = std::chrono::steady_clock::now();
    const uint64_t pipelineMapChunkNs = static_cast<uint64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(workEnd - workStart).count());

    if (P.runThreadN > 1) {
        pthread_mutex_lock(&g_threadChunks.mutexStats);
    }
    g_statsAll.pipelineChunkReadBytes += chunkWorkBytes;
    g_statsAll.pipelineChunksProcessed += 1;
    g_statsAll.pipelineMapChunkNs += pipelineMapChunkNs;
    if (P.runThreadN > 1) {
        pthread_mutex_unlock(&g_threadChunks.mutexStats);
    }

    if (permitEnabled) {
        const uint64_t workNs = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(workEnd - workStart).count());
        const uint64_t readsProcessed = chunk.RA->iRead >= readCountBefore
            ? (chunk.RA->iRead - readCountBefore)
            : chunkReadN;
        g_threadChunks.mapPermitRelease(waitNs, readsProcessed, chunkWorkBytes, workNs);
    }
}

void processCbqRangeTask(ReadAlignChunk& chunk, const Parameters::CbqRangeTask& task) {
    Parameters& P = chunk.P;
    star::input::CbqInputModule rangeReader;
    string inputError;
    if (!rangeReader.configure(P.cbqInputModule->plan(), &inputError) ||
        !rangeReader.open_range(task.laneIndex, task.firstRecord, task.recordCount, &inputError)) {
        exitCbqRangeError(chunk,
                          "could not open a worker-local CBQ range reader",
                          inputError);
    }

    std::shared_ptr<star::input::CbqReadBatchView> pendingBatch;
    uint32 pendingBatchOffset = 0;
    uint64 taskOffset = 0;
    bool rangeExhausted = false;

    while (taskOffset < task.recordCount) {
        const auto chunkReadStartTime = std::chrono::steady_clock::now();
        const uint64 chunkReadStart = task.globalFirst + taskOffset;
        uint64 chunkWorkBytes = 0;

        chunk.chunkInSizeBytesTotal = {0, 0};
        chunk.cbqStarChunk.clear();
        chunk.cbqChunkReadN = 0;

        while (fastxChunkUnderTarget(chunk) &&
               taskOffset < task.recordCount &&
               !rangeExhausted) {
            if (!pendingBatch ||
                pendingBatchOffset >= pendingBatch->record_count) {
                pendingBatch.reset();
                pendingBatchOffset = 0;

                star::input::CbqReadBatchView batch;
                inputError.clear();
                const star::input::InputStatus status = rangeReader.next_batch(&batch, &inputError);
                if (status == star::input::InputStatus::Error) {
                    exitCbqRangeError(chunk,
                                      "worker-local CBQ range reader failed",
                                      inputError);
                }
                if (status == star::input::InputStatus::End) {
                    rangeExhausted = true;
                    break;
                }
                pendingBatch.reset(new star::input::CbqReadBatchView(batch));
            }

            star::input::CbqReadBatchView& batch = *pendingBatch;
            const uint32 batchStart = pendingBatchOffset;
            uint32 recordsToCopy = 0;
            while (pendingBatchOffset + recordsToCopy < batch.record_count &&
                   taskOffset + recordsToCopy < task.recordCount &&
                   fastxChunkUnderTarget(chunk)) {
                const star::input::CbqReadView& record =
                    batch.records[pendingBatchOffset + recordsToCopy];
                const uint64 ordinal = task.globalFirst + taskOffset + recordsToCopy + 1;
                std::array<uint64, MAX_N_MATES> mateBytes{};
                string fitError;
                if (!cbqRecordFits(chunk, record, ordinal, &mateBytes, fitError)) {
                    if (!fitError.empty() && chunk.cbqChunkReadN == 0 && recordsToCopy == 0) {
                        exitCbqRangeError(chunk,
                                          "CBQ record does not fit in the STAR input chunk accounting",
                                          fitError);
                    }
                    break;
                }

                ++recordsToCopy;
                for (uint imate = 0; imate < P.readNends; ++imate) {
                    chunk.chunkInSizeBytesTotal[imate] += mateBytes[imate];
                    chunkWorkBytes += mateBytes[imate];
                }
            }

            if (recordsToCopy == 0) {
                break;
            }

            star::input::CbqReadBatchView subBatch = batch;
            subBatch.records = batch.records + batchStart;
            subBatch.record_count = recordsToCopy;
            const size_t firstAppendedRead = chunk.cbqStarChunk.reads.size();
            string adapterError;
            if (!star::input::append_cbq_batch_to_star_chunk(subBatch,
                                                             recordsToCopy,
                                                             P.readNends,
                                                             &chunk.cbqStarChunk,
                                                             &adapterError)) {
                exitCbqRangeError(chunk,
                                  "CBQ range adapter append failed",
                                  adapterError);
            }
            for (uint32 irecord = 0; irecord < recordsToCopy; ++irecord) {
                star::input::CbqStarChunkRead& read =
                    chunk.cbqStarChunk.reads[firstAppendedRead + irecord];
                read.read_ordinal = task.globalFirst + taskOffset + irecord + 1;
                read.lane_index = task.laneIndex;
            }

            chunk.cbqChunkReadN += recordsToCopy;
            taskOffset += recordsToCopy;
            pendingBatchOffset += recordsToCopy;
            if (pendingBatchOffset >= batch.record_count) {
                pendingBatch.reset();
                pendingBatchOffset = 0;
            }
        }

        if (chunk.cbqChunkReadN == 0) {
            if (rangeExhausted && taskOffset < task.recordCount) {
                ostringstream detail;
                detail << "range ended at " << taskOffset
                       << " records, expected " << task.recordCount;
                exitCbqRangeError(chunk,
                                  "CBQ range reader ended before the planned task was complete",
                                  detail.str());
            }
            exitCbqRangeError(chunk,
                              "CBQ range chunk accounting produced an empty chunk before task completion");
        }

        const auto chunkReadEndTime = std::chrono::steady_clock::now();
        const uint64_t chunkReadNs = static_cast<uint64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(chunkReadEndTime - chunkReadStartTime).count());

        if (P.runThreadN > 1) {
            pthread_mutex_lock(&g_threadChunks.mutexStats);
        }
        g_statsAll.pipelineChunkReadNs += chunkReadNs;
        if (P.runThreadN > 1) {
            pthread_mutex_unlock(&g_threadChunks.mutexStats);
        }

        mapPreparedCbqRangeChunk(chunk,
                                 chunkReadStart,
                                 chunk.cbqChunkReadN,
                                 chunkWorkBytes,
                                 task.laneIndex);
    }
}

void flushCbqRangeResidualSam(ReadAlignChunk& chunk) {
    Parameters& P = chunk.P;
    if (!P.outSAMbool ||
        P.outSAMorder == "PairedKeepInputOrder" ||
        chunk.chunkOutBAMtotal == 0) {
        return;
    }

    if (P.runThreadN > 1) {
        pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    }
    P.inOut->outSAM->write(chunk.chunkOutBAM, chunk.chunkOutBAMtotal);
    P.inOut->outSAM->clear();
    if (P.runThreadN > 1) {
        pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
    }
    chunk.RA->outSAMstream->seekp(0, ios::beg);
    chunk.chunkOutBAMtotal = 0;
}

void processCbqRangeChunks(ReadAlignChunk& chunk) {
    Parameters& P = chunk.P;
    if (!P.cbqRangeNextTask) {
        exitCbqRangeError(chunk, "CBQ range task cursor is not initialized");
    }

    for (;;) {
        const uint32 taskIndex = P.cbqRangeNextTask->fetch_add(1);
        if (taskIndex >= P.cbqRangeTasks.size()) {
            break;
        }
        processCbqRangeTask(chunk, P.cbqRangeTasks[taskIndex]);
    }
    flushCbqRangeResidualSam(chunk);
    chunk.noReadsLeft = true;
}

} // namespace


void ReadAlignChunk::processChunks() {//read-map-write chunks
    // Register this thread's counters for debug statistics collection
    flexCountersRegisterThread();

    if (processRuns > 0 && P.emitYNoYFastqyes) {
        reopenYNoYFastqOutputsForReuse();
    }
    processRuns++;

    noReadsLeft=false; //true if there no more reads left in the file
    bool newFile=false; //new file marker in the input stream
    
    // Track file index for auto-trim detection (file boundary detection)
    int detectionStartFileIndex = P.readFilesIndex;
    
    // Per-file processing: skip to the target file if needed
    // This is used when we've rewound and need to skip to a specific file for mapping
    // NOTE: This requires FILE markers in the read stream (generated by readFilesCommand or
    // when multiple files are provided via comma-separated --readFilesIn). Without FILE markers,
    // trimScope=per-file cannot reliably detect file boundaries.
    if (P.quant.slam.skipToFileIndex > 0) {
        if (P.runThreadN > 1) pthread_mutex_lock(&g_threadChunks.mutexInRead);
        
        // Re-check condition after acquiring mutex (another thread may have already skipped)
        if (P.quant.slam.skipToFileIndex > 0 && P.readFilesIndex < P.quant.slam.skipToFileIndex) {
            P.inOut->logMain << "SLAM per-file: skipping to file " << P.quant.slam.skipToFileIndex 
                             << " (currently at file " << P.readFilesIndex << ")\n";

            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                ostringstream errOut;
                errOut << ERROR_OUT << " EXITING because SLAM per-file skipping is not yet implemented for Binseq/CBQ input.\n";
                errOut << "SOLUTION: use Fastx input for SLAM per-file processing or disable per-file SLAM mode.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            } else if (P.readFilesTypeN == 1 && P.fastxInputActive) {
                fastxSkipToTargetLane(*this);
            } else {
        
            // Skip reads until we reach the target file
            // For paired-end: advance ALL mate streams in sync
            bool foundFileMarker = false;
            uint64 linesScanned = 0;
            const uint64 maxLinesToScan = 100000000;  // Safety limit: 100M lines
            
            while (P.readFilesIndex < P.quant.slam.skipToFileIndex && P.inOut->readIn[0].good()) {
                string line;
                std::getline(P.inOut->readIn[0], line);
                linesScanned++;
                
                // Check for FILE marker in mate 0
                if (line.length() >= 4 && line.substr(0, 4) == "FILE") {
                    foundFileMarker = true;
                    std::istringstream iss(line);
                    string marker;
                    int fileIdx;
                    iss >> marker >> fileIdx;
                    P.readFilesIndex = fileIdx;
                    
                    // Skip the FILE marker line for all other mates (paired-end support)
                    for (uint imate = 1; imate < P.readNends; imate++) {
                        string mateLine;
                        // Each mate stream should have its own FILE marker - skip it
                        if (P.inOut->readIn[imate].good()) {
                            std::getline(P.inOut->readIn[imate], mateLine);
                            // Verify mate has matching FILE marker (optional sanity check)
                            if (mateLine.length() >= 4 && mateLine.substr(0, 4) == "FILE") {
                                // Good - mate stream is in sync
                            }
                        }
                    }
                    
                    if (P.readFilesIndex >= P.quant.slam.skipToFileIndex) {
                        P.inOut->logMain << "SLAM per-file: reached file " << P.readFilesIndex << "\n";
                        break;
                    }
                } else {
                    // Not a FILE marker - skip corresponding lines in other mates
                    for (uint imate = 1; imate < P.readNends; imate++) {
                        if (P.inOut->readIn[imate].good()) {
                            P.inOut->readIn[imate].ignore(numeric_limits<streamsize>::max(), '\n');
                        }
                    }
                }
                
                // Safety check: warn if no FILE marker found after scanning many lines
                if (!foundFileMarker && linesScanned >= maxLinesToScan) {
                    P.inOut->logMain << "WARNING: SLAM per-file: scanned " << linesScanned 
                                     << " lines without finding FILE marker. "
                                     << "trimScope=per-file requires FILE markers in read stream "
                                     << "(use comma-separated --readFilesIn or readFilesCommand).\n";
                    break;
                }
            }
            
            // Warn if we reached end of input without finding target file
            if (P.readFilesIndex < P.quant.slam.skipToFileIndex && !P.inOut->readIn[0].good()) {
                P.inOut->logMain << "WARNING: SLAM per-file: reached end of input before finding file "
                                 << P.quant.slam.skipToFileIndex << ". Current file index: " 
                                 << P.readFilesIndex << "\n";
            }
        
            // Reset skipToFileIndex after skipping
            P.quant.slam.skipToFileIndex = -1;
            }
        }  // end re-check condition
        
        if (P.runThreadN > 1) pthread_mutex_unlock(&g_threadChunks.mutexInRead);
    }
    
    if (P.readFilesTypeN == 20 && P.cbqInputActive && P.cbqRangeActive) {
        processCbqRangeChunks(*this);
    } else {
    while (!noReadsLeft) {//continue until the input EOF
            uint64_t chunkReadN = 0;
            uint64_t chunkWorkBytes = 0;
            uint64_t pipelineMutexWaitNs = 0;
            uint64_t pipelineChunkReadNs = 0;
            uint64_t pipelineChunksProcessed = 0;
            //////////////read a chunk from input files and store in memory
        if (P.outFilterBySJoutStage<2) {//read chunks from input file

            const auto mutexWaitStart = std::chrono::steady_clock::now();
        if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexInRead);
            const auto mutexAcquired = std::chrono::steady_clock::now();

            chunkInSizeBytesTotal={0,0};
            cbqStarChunk.clear();
            cbqChunkReadN = 0;
            const uint64_t chunkReadStart = P.iReadAll;
            
            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                while (fastxChunkUnderTarget(*this) &&
                       !P.cbqInputExhausted && P.iReadAll != P.readMapNumber) {
                    if (!P.cbqInputPendingBatch ||
                        P.cbqInputPendingBatchOffset >= P.cbqInputPendingBatch->record_count) {
                        P.cbqInputPendingBatch.reset();
                        P.cbqInputPendingBatchOffset = 0;

                        star::input::CbqReadBatchView batch;
                        string inputError;
                        const star::input::InputStatus status = P.cbqInputModule->next_batch(&batch, &inputError);
                        if (status == star::input::InputStatus::Error) {
                            ostringstream errOut;
                            errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Binseq/CBQ input module\n"
                                   << inputError << "\n";
                            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                        }
                        if (status == star::input::InputStatus::End) {
                            P.cbqInputExhausted = true;
                            break;
                        }
                        P.cbqInputPendingBatch.reset(new star::input::CbqReadBatchView(batch));
                    }

                    star::input::CbqReadBatchView& batch = *P.cbqInputPendingBatch;
                    P.readFilesIndex = static_cast<int>(batch.lane_index);
                    if (fastxStopBeforeLane(P, batch.lane_index, detectionStartFileIndex)) {
                        noReadsLeft = true;
                        break;
                    }

                    fastxLogLaneStart(P, batch.lane_index);
                    const uint32 batchStart = P.cbqInputPendingBatchOffset;
                    uint32 recordsToCopy = 0;
                    while (P.cbqInputPendingBatchOffset + recordsToCopy < batch.record_count &&
                           P.iReadAll != P.readMapNumber &&
                           fastxChunkUnderTarget(*this)) {
                        const star::input::CbqReadView& record =
                            batch.records[P.cbqInputPendingBatchOffset + recordsToCopy];
                        std::array<uint64, MAX_N_MATES> mateBytes{};
                        string fitError;
                        if (!cbqRecordFits(*this, record, P.iReadAll + 1, &mateBytes, fitError)) {
                            if (!fitError.empty() && cbqChunkReadN == 0 && recordsToCopy == 0) {
                                ostringstream errOut;
                                errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Binseq/CBQ STAR chunk accounting\n"
                                       << fitError << "\n";
                                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                            }
                            break;
                        }

                        ++recordsToCopy;
                        ++P.iReadAll;
                        for (uint imate = 0; imate < P.readNends; ++imate) {
                            chunkInSizeBytesTotal[imate] += mateBytes[imate];
                            chunkWorkBytes += mateBytes[imate];
                        }
                    }

                    if (recordsToCopy == 0) {
                        break;
                    }

                    star::input::CbqReadBatchView subBatch = batch;
                    subBatch.records = batch.records + batchStart;
                    subBatch.record_count = recordsToCopy;
                    string adapterError;
                    if (!star::input::append_cbq_batch_to_star_chunk(subBatch,
                                                                     recordsToCopy,
                                                                     P.readNends,
                                                                     &cbqStarChunk,
                                                                     &adapterError)) {
                        ostringstream errOut;
                        errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Binseq/CBQ STAR adapter\n"
                               << adapterError << "\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    }
                    cbqChunkReadN += recordsToCopy;
                    P.cbqInputPendingBatchOffset += recordsToCopy;
                    if (P.cbqInputPendingBatchOffset >= batch.record_count) {
                        P.cbqInputPendingBatch.reset();
                        P.cbqInputPendingBatchOffset = 0;
                    }
                }
            } else if (P.readFilesTypeN == 1 && P.fastxInputActive) {
                while (fastxChunkUnderTarget(*this)) {
                    if (P.iReadAll == P.readMapNumber) {
                        break;
                    }
                    if (P.fastxInputExhausted && !P.fastxInputPendingRecordValid) {
                        break;
                    }

                    star::input::InputRecord record;
                    if (P.fastxInputPendingRecordValid) {
                        record = *P.fastxInputPendingRecord;
                        P.fastxInputPendingRecordValid = false;
                        P.fastxInputPendingRecord.reset();
                    } else {
                        string inputError;
                        const star::input::InputStatus status = P.fastxInputModule->next_record(&record, &inputError);
                        if (status == star::input::InputStatus::Error) {
                            ostringstream errOut;
                            errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Fastx input module\n"
                                   << inputError << "\n";
                            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                        }
                        if (status == star::input::InputStatus::End) {
                            P.fastxInputExhausted = true;
                            break;
                        }
                    }

                    P.readFilesIndex = static_cast<int>(record.lane_index);
                    if (fastxStopBeforeLane(P, record.lane_index, detectionStartFileIndex)) {
                        noReadsLeft = true;
                        break;
                    }

                    string fitError;
                    if (!fastxRecordFits(*this, record, static_cast<uint64>(P.iReadAll) + 1, fitError)) {
                        if (fastxChunkHasBytes(*this) && fitError.empty()) {
                            P.fastxInputPendingRecord.reset(new star::input::InputRecord(record));
                            P.fastxInputPendingRecordValid = true;
                            break;
                        }

                        ostringstream errOut;
                        errOut << ERROR_OUT << " EXITING because of FATAL ERROR in Fastx input buffering\n";
                        if (!fitError.empty()) {
                            errOut << fitError << "\n";
                        } else {
                            errOut << "Fastx record does not fit in the STAR input chunk buffer\n";
                        }
                        errOut << "SOLUTION: increase --limitIObufferSize or reduce maximum read/header length.\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    }

                    fastxLogLaneStart(P, record.lane_index);
                    P.iReadAll++;
                    P.readFilesIndex = static_cast<int>(record.lane_index);
                    fastxAppendRecord(*this, record);
                }
            } else {
            while (chunkInSizeBytesTotal[0] < P.chunkInSizeBytes && chunkInSizeBytesTotal[1] < P.chunkInSizeBytes && P.inOut->readIn[0].good() && P.inOut->readIn[1].good()) {
                char nextChar=P.inOut->readIn[0].peek();
                if (P.iReadAll==P.readMapNumber) {//do not read any more reads
                    break;
                    
                ///////////////////////////////////////////////////////////////////////////////////// SAM                        
                } else if (P.readFilesTypeN==10 && P.inOut->readIn[0].good() && P.outFilterBySJoutStage!=2) {//SAM input && not eof && not 2nd stage


                    if (nextChar=='@') {//with SAM input linest that start with @ are headers
                        P.inOut->readIn[0].ignore(DEF_readNameSeqLengthMax,'\n'); //read line and skip it
                        continue;
                    };

                    string str1;
                    P.inOut->readIn[0] >> str1;
                    if (str1=="FILE") {
                        newFile=true;
                    } else {
                        P.iReadAll++; //increment read number

                        uint64 flag1; 
                        P.inOut->readIn[0] >> flag1;
                        uint imate1=0;
                        for (uint imate=0;imate<P.readNmates;imate++) {//not readNends: this is SAM input
                            if (imate>0) {
                                string str2;
                                uint64 flag2;
                                P.inOut->readIn[0] >> str2; //for imate=0 str1 was already read
                                P.inOut->readIn[0] >> flag2; //read name and flag
                                
                                if ( str1 != str2 ) {
                                    ostringstream errOut;
                                    errOut << ERROR_OUT <<" EXITING because of FATAL ERROR in input BAM file: the consecutive lines in paired-end BAM have different read IDs:\n"
                                           << str1 <<"   vs   "<< str2 << '\n'
                                           << "\n SOLUTION: fix BAM file formatting. Paired-end reads should be always consecutive lines, with exactly 2 lines per paired-end read" ;
                                    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                                };
                                
                                if (! ( ((flag1 & 0x40) && (flag2 & 0x80)) || ((flag2 & 0x40) && (flag1 & 0x80)) ) ) {
                                    ostringstream errOut;
                                    errOut << ERROR_OUT <<" EXITING because of FATAL ERROR in input BAM file: the consecutive lines in paired-end BAM have wrong mate FLAG bits:\n"
                                           << str1 <<"   "<< flag1 <<"   vs   "<< str2 <<"   "<< flag2 << '\n'
                                           << "\n SOLUTION: fix BAM file formatting. Paired-end reads should be always consecutive lines, with exactly 2 lines per paired-end read."
                                           << " Mate1 should have 0x40 bit set in the FLAG, Mate2 should have 0x80 bit set in the FLAG";
                                    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                                };
                                
                                str1 = str2;   //used below for both mates
                                flag1 = flag2; //used below for both mates
                            };
                            char passFilterIllumina=(flag1 & 0x800 ? 'Y' : 'N');

                            if (imate==1) {//2nd line is always opposite of the 1st one
                                imate1=1-imate1;
                            } else if (P.readNmates==2 && (flag1 & 0x80)) {//not readNends: this is SAM input
                                imate1=1;
                            } else {
                                imate1=0;
                            };

                            //read ID or number
                            if (P.outSAMreadID=="Number") {
                                chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "@%llu", P.iReadAll);
                            } else {
                                chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "@%s", str1.c_str());
                            };

                            //iReadAll, passFilterIllumina, passFilterIllumina
                            chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], " %llu %c %i", P.iReadAll, passFilterIllumina, P.readFilesIndex);

                            string dummy;
                            for (int ii=3; ii<=9; ii++)
                                P.inOut->readIn[0] >> dummy; //skip fields until sequence

                            string seq1,qual1;
                            P.inOut->readIn[0]  >> seq1 >> qual1;
                            if (flag1 & 0x10) {//sequence reverse-coomplemented
                                revComplementNucleotides(seq1);
                                reverse(qual1.begin(),qual1.end());
                            };
                            
                            string attrs;
                            getline(P.inOut->readIn[0], attrs); //rest of the SAM line: str1 is now all SAM attributes - it's added to the read ID line (1st "fastq" line)
                            chunkInSizeBytesTotal[imate1] += sprintf(chunkIn[imate1] + chunkInSizeBytesTotal[imate1], "%s\n%s\n+\n%s\n", attrs.c_str(), seq1.c_str(), qual1.c_str());
                        };
                    };
                    
                ///////////////////////////////////////////////////////////////////////////////////// FASTQ    
                } else if (nextChar=='@') {//fastq, not multi-line
                    P.iReadAll++; //increment read number
                    if (P.outFilterBySJoutStage!=2) {//not the 2nd stage of the 2-stage mapping, read ID from the 1st read
                        string readID;
                        P.inOut->readIn[0] >> readID;
                        removeStringEndControl(readID);
                        if (P.outSAMreadIDnumber) {
                            readID="@"+to_string(P.iReadAll);
                        };
                        //read the second field of the read name line
                        char passFilterIllumina='N';
                        if (P.inOut->readIn[0].peek()!='\n') {//2nd field exists
                            string field2;
                            P.inOut->readIn[0] >> field2;
                            if (field2.length()>=3 && field2[1]==':' && field2[2]=='Y' && field2[3]==':' )
                                passFilterIllumina='Y';
                        };
                        
                        //add extra information to readID line
                        readID += ' '+ to_string(P.iReadAll)+' '+passFilterIllumina+' '+to_string(P.readFilesIndex);

                        //ignore the rest of the read name for both mates
                        for (uint imate=0; imate<P.readNends; imate++)
                            P.inOut->readIn[imate].ignore(DEF_readNameSeqLengthMax,'\n');

                        //copy the same readID to both mates
                        for (uint imate=0; imate<P.readNends; imate++) {
                            chunkInSizeBytesTotal[imate] += 1 + readID.copy(chunkIn[imate] + chunkInSizeBytesTotal[imate], readID.size(),0);
                            chunkIn[imate][chunkInSizeBytesTotal[imate]-1]='\n';
                        };
                    };
                    //copy 3 (4 for stage 2) lines: sequence, dummy, quality
                    for (uint imate=0; imate<P.readNends; imate++) {
                        // read 1st line for 2nd stage only
                        if (P.outFilterBySJoutStage == 2)
                            chunkInSizeBytesTotal[imate] += fastqReadOneLine(P.inOut->readIn[imate], chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                        //sequence
                        chunkInSizeBytesTotal[imate] += fastqReadOneLine(P.inOut->readIn[imate], chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                        //skip 3rd line, record '+'
                        P.inOut->readIn[imate].ignore(DEF_readNameSeqLengthMax, '\n');
                        chunkIn[imate][chunkInSizeBytesTotal[imate]] = '+';
                        chunkIn[imate][chunkInSizeBytesTotal[imate]+1] = '\n';
                        chunkInSizeBytesTotal[imate] += 2;
                        //quality
                        uint64 lenIn = fastqReadOneLine(P.inOut->readIn[imate], chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                        chunkInSizeBytesTotal[imate] += lenIn;
                    };
                } else if (nextChar=='>') {//fasta, can be multiline, which is converted to single line
                    P.iReadAll++; //increment read number
                    for (uint imate=0; imate<P.readNends; imate++) {
                        if (P.outFilterBySJoutStage!=2) {//not the 2nd stage of the 2-stage mapping

                            if (P.outSAMreadID=="Number") {
                                chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], ">%llu", P.iReadAll);
                            } else {
                                P.inOut->readIn[imate] >> (chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                                chunkInSizeBytesTotal[imate] += strlen(chunkIn[imate] + chunkInSizeBytesTotal[imate]);
                            };

                            P.inOut->readIn[imate].ignore(DEF_readNameSeqLengthMax,'\n');

                            chunkInSizeBytesTotal[imate] += sprintf(chunkIn[imate] + chunkInSizeBytesTotal[imate], " %llu %c %i \n", P.iReadAll, 'N', P.readFilesIndex);
                        };
                        
                        //read multi-line fasta
                        nextChar=P.inOut->readIn[imate].peek();
                        while (nextChar!='@' && nextChar!='>' && nextChar!=' ' && nextChar!='\n' && P.inOut->readIn[imate].good()) {
                            P.inOut->readIn[imate].getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readSeqLengthMax + 1 );
                            if (P.inOut->readIn[imate].gcount()<2) 
                                break; //no more input
                                
                            chunkInSizeBytesTotal[imate] += P.inOut->readIn[imate].gcount()-1; //-1 because \n was counted, bu wee need to remove it
                            if ( int(chunkIn[imate][chunkInSizeBytesTotal[imate]-1]) < 33 ) {//remove control char at the end if present
                                chunkInSizeBytesTotal[imate]--;
                            };
                            
                            nextChar=P.inOut->readIn[imate].peek();
                        };
                        chunkIn[imate][chunkInSizeBytesTotal[imate]]='\n';
                        chunkInSizeBytesTotal[imate] ++;
                    };
                } else if (nextChar==' ' || nextChar=='\n' || !P.inOut->readIn[0].good()) {//end of stream
                    P.inOut->logMain << "Thread #" <<iThread <<" end of input stream, nextChar="<<int(nextChar) <<endl;
                    break;
                } else {
                    string word1;
                    P.inOut->readIn[0] >> word1;
                    if (word1=="FILE") {//new file marker
                        newFile=true;
                    } else {//error
                        ostringstream errOut;
                        string str1;
                        std::getline(P.inOut->readIn[0], str1);
                        errOut << ERROR_OUT <<" EXITING because of FATAL ERROR in input reads: wrong read ID line format: the read ID lines should start with @ or > \n";
                        errOut << "Offending line for read # " << P.iReadAll+1 << "\n" << word1 <<" "<< str1 << "\n";
                        errOut << "SOLUTION: verify and correct the input read files\n";
                        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                    };
                };

                if (newFile) {
                        P.inOut->readIn[0] >> P.readFilesIndex;
                        
                        // Check for per-file stop BEFORE processing the new file
                        // This ensures we don't process any reads from the next file
                        bool shouldStop = false;
                        
                        // Auto-trim detection: stop at file boundary for trimScope=first
                        if (P.quant.slam.autoTrimDetectionPass && P.quant.slam.trimScope == "first") {
                            if (P.readFilesIndex > detectionStartFileIndex) {
                                P.inOut->logMain << "SLAM auto-trim detection: stopping at file boundary "
                                                 << "(file " << detectionStartFileIndex << " -> " << P.readFilesIndex << ")\n";
                                shouldStop = true;
                            }
                        }
                        
                        // Per-file processing: stop at file boundary during both detection and mapping
                        if (P.quant.slam.perFileProcessing) {
                            if (P.readFilesIndex > P.quant.slam.currentFileIndex) {
                                P.inOut->logMain << "SLAM per-file processing: stopping at file boundary "
                                                 << "(completed file " << P.quant.slam.currentFileIndex << ")\n";
                                shouldStop = true;
                            }
                        }
                        
                        if (shouldStop) {
                            noReadsLeft = true;
                            newFile = false;
                            break;
                        }
                        
                        // Only log "Starting to map file" if we're actually going to process it
                        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
                        P.inOut->logMain << "Starting to map file # " << P.readFilesIndex<<"\n";
                        for (uint imate=0; imate<P.readFilesNames.size(); imate++) {
                            P.inOut->logMain << "mate " <<imate+1 <<":   "<<P.readFilesNames.at(imate).at(P.readFilesIndex) <<"\n";
                            P.inOut->readIn[imate].ignore(numeric_limits<streamsize>::max(),'\n');
                        };
                        P.inOut->logMain<<flush;
                        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
                        newFile=false;
                };
            };
            }
            if (P.readFilesTypeN != 20 && P.readNends == 2U) {
                const bool m0Empty = (chunkInSizeBytesTotal[0] == 0ULL);
                const bool m1Empty = (chunkInSizeBytesTotal[1] == 0ULL);
                if (m0Empty != m1Empty) {
                    ostringstream errOut;
                    errOut << ERROR_OUT
                           << " EXITING because of FATAL INPUT ERROR: paired mates have unequal "
                              "FASTQ buffering for this chunk (truncated pair or mismatched mate file lists).\n"
                           << "chunkInSizeBytesTotal mate1(bytes)=" << chunkInSizeBytesTotal[0]
                           << ", mate2(bytes)=" << chunkInSizeBytesTotal[1] << "\n";
                    if (m0Empty && !m1Empty) {
                        errOut << "Mate 1 produced no bytes for this chunk while mate 2 still has buffered data "
                                  "(mate 1 stream likely exhausted first).\n";
                    } else {
                        errOut << "Mate 2 produced no bytes for this chunk while mate 1 still has buffered data "
                                  "(mate 2 stream likely exhausted first).\n";
                    }
                    errOut << "SOLUTION: verify R1 and R2 have the same record counts and matching order "
                              "in --readFilesIn.\n";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                }
            }
            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                if (cbqChunkReadN == 0) {
                    noReadsLeft=true;
                    iChunkIn=g_threadChunks.chunkInN;
                    g_threadChunks.chunkInN++;
                } else {
                    noReadsLeft=false;
                    iChunkIn=g_threadChunks.chunkInN;
                    g_threadChunks.chunkInN++;
                }
            } else if (chunkInSizeBytesTotal[0]==0) {
                noReadsLeft=true; //true if there no more reads left in the file
                iChunkIn=g_threadChunks.chunkInN;//to keep things consistent
                g_threadChunks.chunkInN++;
            } else {
                noReadsLeft=false;
                iChunkIn=g_threadChunks.chunkInN;
                g_threadChunks.chunkInN++;
            };

            if (!(P.readFilesTypeN == 20 && P.cbqInputActive)) {
                for (uint imate=0; imate<P.readNends; imate++)
                    chunkIn[imate][chunkInSizeBytesTotal[imate]]='\n';//extra empty line at the end of the chunks
            }
            chunkReadN = P.iReadAll - chunkReadStart;
            if (!(P.readFilesTypeN == 20 && P.cbqInputActive)) {
                for (uint imate=0; imate<P.readNends; imate++) {
                    chunkWorkBytes += chunkInSizeBytesTotal[imate];
                }
            }
            writeInputChunkTrace(*this, iChunkIn, chunkReadStart, chunkReadN, chunkWorkBytes, noReadsLeft);

            if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexInRead);

            {
                const auto chunkReadEnd = std::chrono::steady_clock::now();
                pipelineMutexWaitNs = static_cast<uint64>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(mutexAcquired - mutexWaitStart).count());
                pipelineChunkReadNs = static_cast<uint64>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(chunkReadEnd - mutexAcquired).count());
                pipelineChunksProcessed = 1;
            }

        } else {//read from one file per thread
            noReadsLeft=true;
            for (uint imate=0; imate<P.readNends; imate++) {
                RA->chunkOutFilterBySJoutFiles[imate].flush();
                RA->chunkOutFilterBySJoutFiles[imate].seekg(0,ios::beg);
                RA->readInStream[imate]=& RA->chunkOutFilterBySJoutFiles[imate];
            };
        };

        const bool permitEnabled = g_threadChunks.mapPermitEnabled();
        uint64_t waitNs = 0;
        if (permitEnabled) {
            // Pool size is configured once at startup by mapThreadsSpawn.cpp
            // (and may be re-targeted by a deliberately-ticking controller
            // such as PfPermitController). Re-targeting from every map-worker
            // acquire is policy in the hot path and silently re-clamped the
            // pool to runThreadN, undoing the wider chromapAtac budget.
            // See multiomic-atac-scrna plans/2026-04-27-atac-permits-controller-followups.md
            // Step 3.
            waitNs = g_threadChunks.mapPermitAcquire();
        }
        const auto workStart = std::chrono::steady_clock::now();
        const uint64_t readCountBefore = RA->iRead;
        if (P.readFilesTypeN == 20 && P.cbqInputActive) {
            mapCbqChunk();
        } else {
            mapChunk();
        }
        const auto workEnd = std::chrono::steady_clock::now();
        const uint64_t pipelineMapChunkNs = static_cast<uint64>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(workEnd - workStart).count());
        if (pipelineChunksProcessed > 0) {
            if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexStats);
            g_statsAll.pipelineMutexWaitNs += pipelineMutexWaitNs;
            g_statsAll.pipelineChunkReadNs += pipelineChunkReadNs;
            g_statsAll.pipelineChunkReadBytes += chunkWorkBytes;
            g_statsAll.pipelineChunksProcessed += pipelineChunksProcessed;
            g_statsAll.pipelineMapChunkNs += pipelineMapChunkNs;
            if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexStats);
        }
        if (permitEnabled) {
            const uint64_t workNs = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(workEnd - workStart).count());
            const uint64_t readsProcessed = RA->iRead >= readCountBefore ? (RA->iRead - readCountBefore) : chunkReadN;
            g_threadChunks.mapPermitRelease(waitNs, readsProcessed, chunkWorkBytes, workNs);
        }

        if (iThread==0 && P.runThreadN>1 && P.outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
            chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
        };

    };//cycle over input chunks
    }

    // Close Y/noY FASTQ gzip streams even if no reads were processed in this thread
    if (P.emitYNoYFastqyes && P.emitYNoYFastqCompression == "gz") {
        const uint32 yFastqEmitCount = P.yFastqEmitReadCount();
        for (uint32 imate = 0; imate < yFastqEmitCount; imate++) {
            if (RA->chunkOutYFastqGz[imate] != nullptr) {
                gzclose(RA->chunkOutYFastqGz[imate]);
                RA->chunkOutYFastqGz[imate] = nullptr;
            }
            if (RA->chunkOutNoYFastqGz[imate] != nullptr) {
                gzclose(RA->chunkOutNoYFastqGz[imate]);
                RA->chunkOutNoYFastqGz[imate] = nullptr;
            }
        }
    }

    // Skip output operations during detection-only passes.
    if (P.outFilterBySJoutStage!=1 && RA->iRead>0 &&
        !P.quant.slam.autoTrimDetectionPass &&
        !P.quant.transcriptVB.inDetectionMode) {//not the first stage of the 2-stage mapping, and not detection-only
        if (P.outBAMunsorted && chunkOutBAMunsorted!=NULL) chunkOutBAMunsorted->unsortedFlush();
        if (P.outBAMcoord) chunkOutBAMcoord->coordFlush();
        if (chunkOutBAMquant!=NULL) chunkOutBAMquant->unsortedFlush();

        //the thread is finished mapping reads, concatenate the temp files into output files
        if (P.pCh.segmentMin>0) {
            chunkFstreamCat (RA->chunkOutChimSAM, P.inOut->outChimSAM, P.runThreadN>1, g_threadChunks.mutexOutChimSAM);
            chunkFstreamCat (*RA->chunkOutChimJunction, P.inOut->outChimJunction, P.runThreadN>1, g_threadChunks.mutexOutChimJunction);
        };
        if (P.outReadsUnmapped=="Fastx" ) {
            if (P.runThreadN>1)
                pthread_mutex_lock(&g_threadChunks.mutexOutUnmappedFastx);

            for (uint ii=0;ii<P.readNends;ii++) {
                chunkFstreamCat (RA->chunkOutUnmappedReadsStream[ii],P.inOut->outUnmappedReadsStream[ii], false, g_threadChunks.mutexOutUnmappedFastx);
            };

            if (P.runThreadN>1)
                pthread_mutex_unlock(&g_threadChunks.mutexOutUnmappedFastx);
        };
        
        // Concatenate Y/noY FASTQ thread outputs
        if (P.emitYNoYFastqyes && P.emitYNoYFastqCompression != "gz") {
            const uint32 yFastqEmitCount = P.yFastqEmitReadCount();
            for (uint32 imate = 0; imate < yFastqEmitCount; imate++) {
                chunkFstreamCat(RA->chunkOutYFastqStream[imate], P.inOut->outYFastqStream[imate], P.runThreadN > 1, g_threadChunks.mutexOutYFastq[imate]);
                chunkFstreamCat(RA->chunkOutNoYFastqStream[imate], P.inOut->outNoYFastqStream[imate], P.runThreadN > 1, g_threadChunks.mutexOutNoYFastq[imate]);
            }
        }
    };
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexLogMain);
    P.inOut->logMain << "Completed: thread #" <<iThread <<endl;
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
};

inline uint64 fastqReadOneLine(ifstream &streamIn, char *arrIn)
{
    uint64 lenIn;
    streamIn.getline(arrIn, DEF_readNameSeqLengthMax+1 );
    lenIn = streamIn.gcount(); //=seqLength+1: includes \0 but not \n. We will replace \0 with \n
    
    if ( int(arrIn[lenIn-2]) < 33 ) {//remove control char at the end if present
        --lenIn;
    };
    
    arrIn[lenIn-1]='\n'; //replace \0 with \n
    return lenIn; //lenIn contains \n at the end
};

inline void removeStringEndControl(string &str)
{//removes control character (including space) from the end of the string
    if (int(str.back())<33)
        str.pop_back();
};
