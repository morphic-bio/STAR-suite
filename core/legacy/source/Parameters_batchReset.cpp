// Parameters_batchReset.cpp
// Reset functions for batch mode (--batchMode / --slamBatchMode)

#include "Parameters.h"
#include "ErrorWarning.h"
#include "InOutStreams.h"
#include "TimeFunctions.h"
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>

// Helper: create directory recursively
static bool ensureDirExists(const std::string& path, mode_t mode) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) {
        return S_ISDIR(st.st_mode);
    }
    // Try to create
    if (mkdir(path.c_str(), mode) == 0 || errno == EEXIST) {
        return true;
    }
    // Try creating parent first
    size_t lastSlash = path.rfind('/');
    if (lastSlash != std::string::npos && lastSlash > 0) {
        if (ensureDirExists(path.substr(0, lastSlash), mode)) {
            return mkdir(path.c_str(), mode) == 0 || errno == EEXIST;
        }
    }
    return false;
}

// Helper: basename/path/FASTQ tagging (mirrors Parameters.cpp)
static std::string pathBasenameLocal(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) {
        return path;
    }
    return path.substr(pos + 1);
}

static std::string outputDirFromPrefixLocal(const std::string& prefix) {
    size_t pos = prefix.find_last_of("/\\");
    if (pos == std::string::npos) {
        return "";
    }
    return prefix.substr(0, pos + 1);
}

static bool hasReadTokenLocal(const std::string& name) {
    return name.rfind("_R1") != std::string::npos ||
           name.rfind("_R2") != std::string::npos ||
           name.rfind("_r1") != std::string::npos ||
           name.rfind("_r2") != std::string::npos;
}

static std::string insertTagBeforeReadTokenLocal(const std::string& name, const std::string& tag) {
    size_t pos = std::string::npos;
    auto updatePos = [&](size_t candidate) {
        if (candidate != std::string::npos) {
            pos = (pos == std::string::npos) ? candidate : std::max(pos, candidate);
        }
    };
    updatePos(name.rfind("_R1"));
    updatePos(name.rfind("_R2"));
    updatePos(name.rfind("_r1"));
    updatePos(name.rfind("_r2"));
    if (pos == std::string::npos) {
        return std::string();
    }
    return name.substr(0, pos) + tag + name.substr(pos);
}

static bool endsWithLocal(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string adjustCompressionExtLocal(const std::string& name, const std::string& compression) {
    if (compression == "gz") {
        if (!endsWithLocal(name, ".gz")) {
            return name + ".gz";
        }
        return name;
    }
    if (endsWithLocal(name, ".gz")) {
        return name.substr(0, name.size() - 3);
    }
    return name;
}

// Extract sample name from a FASTQ path (exported for use in STAR.cpp)
std::string extractSampleNameFromFastq(const std::string& path) {
    // Find the filename
    size_t lastSlash = path.rfind('/');
    std::string filename = (lastSlash == std::string::npos) ? path : path.substr(lastSlash + 1);
    
    // Remove common suffixes in order
    std::vector<std::string> suffixes = {
        ".fastq.gz", ".fq.gz", ".fastq", ".fq",
        "_R1_001", "_R2_001", "_R1", "_R2",
        "_1", "_2"
    };
    
    for (const auto& suffix : suffixes) {
        if (filename.size() > suffix.size()) {
            size_t pos = filename.rfind(suffix);
            if (pos != std::string::npos && pos + suffix.size() == filename.size()) {
                filename = filename.substr(0, pos);
            }
        }
    }
    
    // Also try to remove _S##_ pattern (Illumina sample index)
    size_t sPos = filename.rfind("_S");
    if (sPos != std::string::npos && sPos > 0) {
        // Check if followed by digits
        size_t checkPos = sPos + 2;
        while (checkPos < filename.size() && std::isdigit(filename[checkPos])) {
            checkPos++;
        }
        if (checkPos > sPos + 2 && (checkPos == filename.size() || filename[checkPos] == '_')) {
            filename = filename.substr(0, sPos);
        }
    }
    
    return filename;
}

void Parameters::resetForBatchSample(int sampleIndex, const std::string& sampleName) {
    // Update batch tracking state
    quant.slam.batchCurrentIndex = sampleIndex;
    if (slamDumpGlobalCount) {
        slamDumpGlobalCount->store(0);
    }

    const std::string batchLabel = quant.slam.yes ? "SLAM Batch Mode" : "Batch Mode";
    
    // Log the sample transition
    inOut->logMain << "\n" << timeMonthDayTime() 
                   << " ===== " << batchLabel << ": Processing sample " << (sampleIndex + 1) 
                   << "/" << quant.slam.batchTotalCount 
                   << " (" << sampleName << ") =====\n" << std::flush;
    *inOut->logStdOut << timeMonthDayTime() 
                      << " ===== " << batchLabel << ": " << sampleName 
                      << " (" << (sampleIndex + 1) << "/" << quant.slam.batchTotalCount << ") =====\n" << std::flush;
    
    // Reset read file index and reopen files for new sample
    readFilesIndex = -1;
    iReadAll = 0;

    // Clear per-sample FIFO list for readFilesCommand path to avoid stale tmp paths
    readFilesInTmp.clear();
    for (uint imate = 0; imate < MAX_N_MATES; ++imate) {
        readFilesCommandPID[imate] = 0;
    }

    // Read files will be reopened by the caller after output paths/tmp dir are updated
    
    if (quant.slam.yes) {
        // Reset SLAM-specific state
        quant.slam.autoTrimComputed = false;
        quant.slam.autoTrim5p = 0;
        quant.slam.autoTrim3p = 0;
        quant.slam.autoTrimFileIndex = 0;
        quant.slam.autoTrimReplayDone = false;
        quant.slam.autoTrimDetectionPass = false;
        quant.slam.currentFileIndex = 0;
        quant.slam.varianceStddevTcRate.clear();
        quant.slam.snpErrEst = 0.0;
        quant.slam.snpErrFallbackReason.clear();

        // For samples after blank, use the blank's error rate
        if (sampleIndex > 0 && quant.slam.batchBlankProcessed && quant.slam.errorRateFromBlank) {
            // Use the captured blank error rate
            quant.slam.errorRate = quant.slam.batchBlankErrorRate;
            quant.slam.snpErrUsed = quant.slam.batchBlankErrorRate;
            quant.slam.errorRateFromBlank = false;  // Don't re-compute, use the value
            inOut->logMain << "  Using blank error rate: " << quant.slam.batchBlankErrorRate << "\n" << std::flush;
        }
    }
    
    // Update iReadAll counter (will be reset by the alignment loop)
    iReadAll = 0;
}

void Parameters::reconfigureOutputPathsForSample(const std::string& sampleName) {
    // This reconfigures output paths for the new sample
    // Uses outFileNamePrefixAuto if set, otherwise appends to existing prefix
    
    std::string newPrefix;
    std::string alignDir;
    std::string countsDir;
    std::string qcDir;
    std::string yDir;
    
    if (outFileNamePrefixAuto && !outFileNamePrefixAutoRoot.empty()) {
        // Auto mode: use structured directory layout
        outFileNamePrefixAutoSample = sampleName;
        alignDir = outFileNamePrefixAutoRoot + "alignments/" + sampleName + "/";
        countsDir = outFileNamePrefixAutoRoot + "counts/" + sampleName + "/";
        qcDir = outFileNamePrefixAutoRoot + "qc/" + sampleName + "/";
        yDir = outFileNamePrefixAutoRoot + "y_separated/" + sampleName + "/";
        newPrefix = alignDir + sampleName + "_";
    } else {
        // Manual mode: use sample name in prefix
        // Use the ORIGINAL output prefix stored at batch start (not the current one)
        std::string baseDir;
        if (!quant.slam.batchOriginalPrefix.empty()) {
            baseDir = quant.slam.batchOriginalPrefix;
        } else {
            // Fallback: find base directory from current outFileNamePrefix
            size_t lastSlash = outFileNamePrefix.rfind('/');
            baseDir = (lastSlash != std::string::npos) ? 
                outFileNamePrefix.substr(0, lastSlash + 1) : "./";
        }
        alignDir = baseDir + sampleName + "/";
        countsDir = baseDir + "counts/" + sampleName + "/";
        qcDir = baseDir + "qc/" + sampleName + "/";
        yDir = baseDir + "y_separated/" + sampleName + "/";
        newPrefix = alignDir + sampleName + "_";
    }
    
    // Create directories
    if (!ensureDirExists(alignDir, runDirPerm)) {
        std::ostringstream errOut;
        errOut << "EXITING: cannot create directory " << alignDir << ": " << strerror(errno) << "\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_RUNTIME, *this);
    }
    if (quant.slam.yes || quant.geCount.yes || quant.transcriptVB.yes) {
        ensureDirExists(countsDir, runDirPerm);
    }
    if (quant.slam.yes) {
        ensureDirExists(qcDir, runDirPerm);
    }
    if (emitNoYBAMyes || emitYNoYFastqyes) {
        ensureDirExists(yDir, runDirPerm);
    }
    
    // Update the main output prefix
    outFileNamePrefix = newPrefix;
    
    // Update temp directory for this sample
    outFileTmp = outFileNamePrefix + "_STARtmp/";
    ensureDirExists(outFileTmp, runDirPerm);

    // Update per-sample quant outputs for batch mode
    if (quant.geCount.yes) {
        if (outFileNamePrefixAuto) {
            quant.geCount.outFile = countsDir + sampleName + ".ReadsPerGene.out.tab";
        } else {
            quant.geCount.outFile = outFileNamePrefix + "ReadsPerGene.out.tab";
        }
    }
    if (quant.transcriptVB.yes) {
        if (outFileNamePrefixAuto) {
            quant.transcriptVB.outFile = countsDir + sampleName + ".quant.sf";
            quant.transcriptVB.outFileGene = countsDir + sampleName + ".quant.genes.sf";
            quant.transcriptVB.outFileGeneTximport = countsDir + sampleName + ".quant.genes.tximport.sf";
        } else {
            quant.transcriptVB.outFile = outFileNamePrefix + "quant.sf";
            quant.transcriptVB.outFileGene = outFileNamePrefix + "quant.genes.sf";
            quant.transcriptVB.outFileGeneTximport = outFileNamePrefix + "quant.genes.tximport.sf";
        }
    }

    // Update transcriptome BAM output for batch mode so each sample gets its own writer.
    if (quant.trSAM.bamYes) {
        if (outStd == "BAM_Quant") {
            outQuantBAMfileName = "-";
        } else {
            outQuantBAMfileName = newPrefix + "Aligned.toTranscriptome.out.bam";
        }
    }
    
    // Update BAM output paths
    if (outBAMunsorted) {
        outBAMfileUnsortedName = newPrefix + "Aligned.out.bam";
    }
    if (outBAMcoord) {
        outBAMfileCoordName = newPrefix + "Aligned.sortedByCoord.out.bam";
    }
    
    // Update SLAM output paths
    if (quant.slam.yes) {
        quant.slam.outFile = countsDir + sampleName + ".SlamQuant.out";
        if (quant.slam.grandSlamOut) {
            quant.slam.grandSlamOutFile = countsDir + sampleName + ".SlamQuant.grandslam.tsv";
        }
        quant.slam.slamQcJson = qcDir + sampleName + ".slam_qc.json";
        quant.slam.slamQcHtml = qcDir + sampleName + ".slam_qc.html";
        if (quant.slam.dumpBinary != "-" && quant.slam.dumpBinary != "None" && !quant.slam.dumpBinary.empty()) {
            quant.slam.dumpBinary = alignDir + sampleName + "_slam_dump.bin";
        }
        if (quant.slam.dumpWeights != "-" && quant.slam.dumpWeights != "None" && !quant.slam.dumpWeights.empty()) {
            quant.slam.dumpWeights = alignDir + sampleName + "_slam_weights.bin";
        }
    }
    
    // Update Y/noY BAM paths
    if (emitNoYBAMyes) {
        outBAMfileNoYName = yDir + sampleName + "_noY.bam";
        outBAMfileYName = yDir + sampleName + "_Y.bam";
    }

    // Update Y/noY FASTQ paths (per-sample in batch mode)
    if (emitYNoYFastqyes) {
        const bool hasYPrefix = !YFastqOutputPrefix.empty() && YFastqOutputPrefix != "-";
        const bool hasNoYPrefix = !noYFastqOutputPrefix.empty() && noYFastqOutputPrefix != "-";
        const std::string ext = (emitYNoYFastqCompression == "gz") ? ".fastq.gz" : ".fastq";
        const std::string outputDir = yDir.empty() ? outputDirFromPrefixLocal(outFileNamePrefix) : yDir;
        bool useMateFallback = false;

        if (!hasYPrefix || !hasNoYPrefix) {
            for (uint imate = 0; imate < readNmates; imate++) {
                if (readFilesNames.size() <= imate || readFilesNames[imate].empty()) {
                    useMateFallback = true;
                    break;
                }
                std::string base = pathBasenameLocal(readFilesNames[imate][0]);
                if (!hasReadTokenLocal(base)) {
                    useMateFallback = true;
                    break;
                }
            }
        }

        for (uint imate = 0; imate < readNmates; imate++) {
            if (hasYPrefix) {
                outYFastqFile[imate] = YFastqOutputPrefix + "mate" + std::to_string(imate + 1) + ext;
            } else if (!useMateFallback && readFilesNames.size() > imate && !readFilesNames[imate].empty()) {
                std::string base = pathBasenameLocal(readFilesNames[imate][0]);
                std::string tagged = insertTagBeforeReadTokenLocal(base, "_Y");
                if (!tagged.empty()) {
                    outYFastqFile[imate] = outputDir + adjustCompressionExtLocal(tagged, emitYNoYFastqCompression);
                } else {
                    outYFastqFile[imate] = outFileNamePrefix + "Y_reads.mate" + std::to_string(imate + 1) + ext;
                }
            } else {
                outYFastqFile[imate] = outFileNamePrefix + "Y_reads.mate" + std::to_string(imate + 1) + ext;
            }

            if (hasNoYPrefix) {
                outNoYFastqFile[imate] = noYFastqOutputPrefix + "mate" + std::to_string(imate + 1) + ext;
            } else if (!useMateFallback && readFilesNames.size() > imate && !readFilesNames[imate].empty()) {
                std::string base = pathBasenameLocal(readFilesNames[imate][0]);
                std::string tagged = insertTagBeforeReadTokenLocal(base, "_noY");
                if (!tagged.empty()) {
                    outNoYFastqFile[imate] = outputDir + adjustCompressionExtLocal(tagged, emitYNoYFastqCompression);
                } else {
                    outNoYFastqFile[imate] = outFileNamePrefix + "noY_reads.mate" + std::to_string(imate + 1) + ext;
                }
            } else {
                outNoYFastqFile[imate] = outFileNamePrefix + "noY_reads.mate" + std::to_string(imate + 1) + ext;
            }
        }

        // Reopen uncompressed Y/noY FASTQ streams per sample
        if (emitYNoYFastqCompression == "none") {
            for (uint imate = 0; imate < readNmates; imate++) {
                if (inOut->outYFastqStream[imate].is_open()) {
                    inOut->outYFastqStream[imate].close();
                }
                if (inOut->outNoYFastqStream[imate].is_open()) {
                    inOut->outNoYFastqStream[imate].close();
                }
                inOut->outYFastqStream[imate].open(outYFastqFile[imate].c_str());
                inOut->outNoYFastqStream[imate].open(outNoYFastqFile[imate].c_str());
                if (!inOut->outYFastqStream[imate].is_open() || !inOut->outNoYFastqStream[imate].is_open()) {
                    std::ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR: could not create Y/noY FASTQ output files\n";
                    errOut << "Solution: check that you have permission to write and disk space\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
                }
            }
        }
    }
    
    // Update log output path
    outLogFileName = newPrefix + "Log.out";

    // Re-open log files for this sample
    if (inOut->logMain.is_open()) {
        inOut->logMain.flush();
        inOut->logMain.close();
    }
    inOut->logMain.open(outLogFileName.c_str());
    if (inOut->logMain.fail()) {
        std::ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create log file: "
               << outLogFileName << "\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    }

    if (inOut->logProgress.is_open()) {
        inOut->logProgress.flush();
        inOut->logProgress.close();
    }
    inOut->logProgress.open((outFileNamePrefix + "Log.progress.out").c_str());

    if (quant.trSAM.bamYes) {
        if (inOut->outQuantBAMfile != nullptr) {
            bgzf_close(inOut->outQuantBAMfile);
            inOut->outQuantBAMfile = nullptr;
        }
        inOut->outQuantBAMfile = bgzf_open(outQuantBAMfileName.c_str(),
                                           ("w" + to_string((long long) quant.trSAM.bamCompression)).c_str());
        if (inOut->outQuantBAMfile == nullptr) {
            std::ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT ERROR: could not create transcriptome BAM file: "
                   << outQuantBAMfileName << "\n"
                   << "SOLUTION: check the path and permissions\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    }

    // Reopen BAM outputs for this sample (batch mode uses per-sample output paths)
    if (outBAMunsorted) {
        if (inOut->outBAMfileUnsorted != nullptr) {
            bgzf_close(inOut->outBAMfileUnsorted);
            inOut->outBAMfileUnsorted = nullptr;
        }
        if (inOut->outBAMfileY != nullptr) {
            bgzf_close(inOut->outBAMfileY);
            inOut->outBAMfileY = nullptr;
        }
        if (inOut->outBAMfileNoY != nullptr) {
            bgzf_close(inOut->outBAMfileNoY);
            inOut->outBAMfileNoY = nullptr;
        }

        if (pSolo.samAttrYes && !pSolo.skipProcessing) {
            // Unsorted BAM will be produced after Solo (tag injection)
            inOut->outBAMfileUnsorted = nullptr;
            inOut->outBAMfileY = nullptr;
            inOut->outBAMfileNoY = nullptr;
        } else {
            inOut->outBAMfileUnsorted = bgzf_open(outBAMfileUnsortedName.c_str(),
                                                  ("w" + to_string((long long)outBAMcompression)).c_str());
            if (inOut->outBAMfileUnsorted == nullptr) {
                std::ostringstream errOut;
                errOut << "EXITING because of fatal OUTPUT ERROR: could not create unsorted BAM file: "
                       << outBAMfileUnsortedName << "\n"
                       << "SOLUTION: check the path and permissions\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }

            if (emitNoYBAMyes) {
                inOut->outBAMfileY = bgzf_open(outBAMfileYName.c_str(),
                                               ("w" + to_string((long long)outBAMcompression)).c_str());
                if (inOut->outBAMfileY == nullptr) {
                    std::ostringstream errOut;
                    errOut << "EXITING because of fatal OUTPUT ERROR: could not create Y BAM file: "
                           << outBAMfileYName << "\n"
                           << "SOLUTION: check the path and permissions\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }

                inOut->outBAMfileNoY = bgzf_open(outBAMfileNoYName.c_str(),
                                                 ("w" + to_string((long long)outBAMcompression)).c_str());
                if (inOut->outBAMfileNoY == nullptr) {
                    std::ostringstream errOut;
                    errOut << "EXITING because of fatal OUTPUT ERROR: could not create noY BAM file: "
                           << outBAMfileNoYName << "\n"
                           << "SOLUTION: check the path and permissions\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }
        }
    }

    inOut->logMain << "  Output directory: " << alignDir << "\n";
    inOut->logMain << "  Output prefix: " << newPrefix << "\n" << std::flush;
}
