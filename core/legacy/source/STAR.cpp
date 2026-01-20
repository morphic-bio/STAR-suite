#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <limits>

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "Chain.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "SlamQuant.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "sysRemoveDir.h"
#include "BAMfunctions.h"
#include "bamSortByCoordinate.h"
#include "SamtoolsSorter.h"
#include "Transcriptome.h"
#include "signalFromBAM.h"
#include "mapThreadsSpawn.h"
#include "SjdbClass.h"
#include "sjdbInsertJunctions.h"
#include "Variation.h"
#include "Solo.h"
#include "samHeaders.h"
#include "systemFunctions.h"
#include "ProbeListIndex.h"
#include "TranscriptQuantEC.h"
#include "LibFormatDetection.h"
#include "TrimQcOutput.h"
#include "vb_engine.h"
#include "em_engine.h"
#include "ec_loader.h"
#include "gc_bias.h"
#include "fld_accumulator.h"
#include "TranscriptQuantOutput.h"
#include "SlamQuant.h"
#include "SlamDump.h"
#include "SlamVarianceAnalysis.h"
#include "SlamQcOutput.h"
#include "SnpMaskBuild.h"
#include "CrMultiProcess.h"
// Note: effective_length.h not included due to Transcriptome class name conflict
// Use wrapper function instead
#include "effective_length_wrapper.h"
#include "InlineCBCorrection.h"
#include "alignment_model.h"  // For Transcriptome and AlignmentModel
#include <memory>
#include <unordered_map>
#include <cstdlib>
#include <cerrno>
#include <cstring>

#include "twoPassRunPass1.h"

#include "htslib/htslib/sam.h"

namespace {

bool ensureDirectoryTree(const std::string &path, mode_t mode, std::string &failedPath, int &failedErrno) {
    if (path.empty()) {
        return true;
    }
    std::string normalized = path;
    while (!normalized.empty() && normalized.back() == '/') {
        normalized.pop_back();
    }
    if (normalized.empty()) {
        return true;
    }

    size_t start = 0;
    std::string current;
    if (!normalized.empty() && normalized[0] == '/') {
        current = "/";
        start = 1;
    }

    while (start <= normalized.size()) {
        size_t slash = normalized.find('/', start);
        size_t len = (slash == std::string::npos) ? normalized.size() - start : slash - start;
        std::string token = normalized.substr(start, len);
        if (!token.empty()) {
            if (current.empty()) {
                current = token;
            } else if (current == "/") {
                current += token;
            } else {
                current += "/" + token;
            }
            if (mkdir(current.c_str(), mode) != 0 && errno != EEXIST) {
                failedPath = current;
                failedErrno = errno;
                return false;
            }
        }
        if (slash == std::string::npos) {
            break;
        }
        start = slash + 1;
    }
    return true;
}

bool buildVbGeneWeights(const std::vector<const SlamReadBuffer*>& buffers,
                        const std::vector<double>& genePosterior,
                        uint64_t maxReads,
                        std::vector<double>* out,
                        std::string* err) {
    if (out == nullptr) {
        if (err) *err = "Null output for VB gene weights";
        return false;
    }
    std::vector<const SlamBufferedRead*> reads;
    reads.reserve(maxReads > 0 ? static_cast<size_t>(maxReads) : 1024);
    for (const auto* buffer : buffers) {
        if (buffer == nullptr) continue;
        for (const auto& r : buffer->reads()) {
            if (maxReads > 0 && reads.size() >= maxReads) break;
            reads.push_back(&r);
        }
        if (maxReads > 0 && reads.size() >= maxReads) break;
    }
    out->assign(reads.size(), 0.0);
    if (reads.empty()) {
        return true;
    }
    std::vector<double> scores(reads.size(), 0.0);
    std::unordered_map<std::string, double> sumByRead;
    std::unordered_map<std::string, size_t> countByRead;
    sumByRead.reserve(reads.size());
    countByRead.reserve(reads.size());

    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = *reads[i];
        double score = 0.0;
        for (uint32_t gid : r.geneIds) {
            if (gid < genePosterior.size()) {
                score += genePosterior[gid];
            }
        }
        scores[i] = score;
        std::string key = std::to_string(r.fileIndex) + "\t" + r.readName;
        sumByRead[key] += score;
        countByRead[key] += 1;
    }
    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = *reads[i];
        std::string key = std::to_string(r.fileIndex) + "\t" + r.readName;
        double denom = sumByRead[key];
        size_t count = countByRead[key];
        if (denom > 0.0) {
            (*out)[i] = scores[i] / denom;
        } else if (count > 0) {
            (*out)[i] = 1.0 / static_cast<double>(count);
        }
    }
    return true;
}

} // namespace
#include "parametersDefault.xxd"

void usage(int usageType)
{
    cout << "Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq\n";
    cout << "Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022\n\n";
    cout << "STAR version=" << STAR_VERSION << "\n";
    cout << "STAR compilation time,server,dir=" << COMPILATION_TIME_PLACE << "\n";
    cout << "For more details see:\n";
    cout << "<https://github.com/alexdobin/STAR>\n";
    cout << "<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>\n";

    if (usageType == 0)
    { // brief
        cout << "\nTo list all parameters, run STAR --help\n";
    }
    else if (usageType == 1)
    { // full
        cout.write(reinterpret_cast<char *>(parametersDefault),
                   parametersDefault_len);
    };
    exit(0);
};

int main(int argInN, char *argIn[])
{
    // If no argument is given, or the first argument is either '-h' or '--help', run usage()
    if (argInN == 1)
    {
        usage(0);
    }
    else if (argInN == 2 && (strcmp("-h", argIn[1]) == 0 || strcmp("--help", argIn[1]) == 0))
    {
        usage(1);
    };

    time(&g_statsAll.timeStart);

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////// Parameters
    Parameters P; // all parameters
    P.inputParameters(argInN, argIn);

    *(P.inOut->logStdOut) << "\t" << P.commandLine << '\n';
    *(P.inOut->logStdOut) << "\tSTAR version: " << STAR_VERSION << "   compiled: " << COMPILATION_TIME_PLACE << '\n';
    *(P.inOut->logStdOut) << timeMonthDayTime(g_statsAll.timeStart) << " ..... started STAR run\n"
                          << flush;

    // runMode
    if (P.runMode == "alignReads" || P.runMode == "soloCellFiltering")
    {
        // continue
    }
    else if (P.runMode == "genomeGenerate")
    {
        { // normal genome generation
            Genome genomeMain(P, P.pGe);
            genomeMain.genomeGenerate();
        };

        if (P.pGe.transform.type > 0)
        { // generate original genome, in addition to the transfomed generated above
            P.pGe.transform.type = 0;
            P.pGe.transform.typeString = "None";
            P.pGe.transform.vcfFile = "-";
            P.pGe.gDir += "/OriginalGenome/";
            Genome genomeOrig(P, P.pGe);
            genomeOrig.genomeGenerate();
        };

        sysRemoveDir(P.outFileTmp);
        P.inOut->logMain << "DONE: Genome generation, EXITING\n"
                         << flush;
        exit(0);
    }
    else if (P.runMode == "liftOver")
    {
        for (uint ii = 0; ii < P.pGe.gChainFiles.size(); ii++)
        {
            Chain chain(P, P.pGe.gChainFiles.at(ii));
            chain.liftOverGTF(P.pGe.sjdbGTFfile, P.outFileNamePrefix + "GTFliftOver_" + to_string(ii + 1) + ".gtf");
            P.inOut->logMain << "DONE: lift-over of GTF file, EXITING\n"
                             << flush;
            exit(0);
        };
    }
    else
    {
        P.inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" << P.runMode << endl
                         << flush;
        exit(1);
    };

    // transcriptome placeholder (loaded only if P.quant.yes)
    Transcriptome *transcriptomeMain = nullptr;
    std::vector<double> vbGenePosterior;
    bool vbGenePosteriorReady = false;

    // --runMode soloCellFiltering executes and exits before genome/mapping
    if (P.runMode == "soloCellFiltering") {
        Transcriptome transcriptomeCellFilter(P);
        Solo soloCellFilter(P, transcriptomeCellFilter);
    }

    ////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// Genome
    Genome genomeMain(P, P.pGe);
    genomeMain.genomeLoad();

    if (P.pGe.transform.outYes) {
        genomeMain.Var = new Variation(P, genomeMain.chrStart, genomeMain.chrNameIndex, false);//no variation for mapGen, only for genOut
        genomeMain.genomeOut.g->Var = new Variation(P, genomeMain.genomeOut.g->chrStart, genomeMain.genomeOut.g->chrNameIndex, P.var.yes);
    } else {
        genomeMain.Var = new Variation(P, genomeMain.chrStart, genomeMain.chrNameIndex, P.var.yes);
    };

    SjdbClass sjdbLoci;

    if (P.sjdbInsert.pass1) {
        Genome genomeMain1 = genomeMain; // not sure if I need to create the copy - genomeMain1 below should not be changed
        sjdbInsertJunctions(P, genomeMain, genomeMain1, sjdbLoci);
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////START
    if (P.runThreadN > 1)
    {
        g_threadChunks.threadArray = new pthread_t[P.runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
        pthread_mutex_init(&g_threadChunks.mutexBAMsortBins, NULL);
        pthread_mutex_init(&g_threadChunks.mutexError, NULL);
        if (P.emitYNoYFastqyes) {
            for (uint32 imate = 0; imate < P.readNmates; imate++) {
                pthread_mutex_init(&g_threadChunks.mutexOutYFastq[imate], NULL);
                pthread_mutex_init(&g_threadChunks.mutexOutNoYFastq[imate], NULL);
            }
        }
    };

    g_statsAll.progressReportHeader(P.inOut->logProgress);

    /////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// 2-pass 1st pass
    twoPassRunPass1(P, genomeMain, transcriptomeMain, sjdbLoci);

    // Shared libem Transcriptome sequence cache for error model (read-only, shared across threads)
    std::unique_ptr<libem::Transcriptome> libem_transcriptome;
    
    if (P.quant.yes)
    { // load transcriptome
        transcriptomeMain = new Transcriptome(P);

        // SNP mask build pre-pass (if requested)
        bool hasMaskIn = !P.quant.slamSnpMask.maskIn.empty() && P.quant.slamSnpMask.maskIn != "-" && P.quant.slamSnpMask.maskIn != "None";
        bool hasVcfIn = !P.quant.slamSnpMask.vcfIn.empty() && P.quant.slamSnpMask.vcfIn != "-" && P.quant.slamSnpMask.vcfIn != "None";
        bool hasBuildFastqs = !P.quant.slamSnpMask.buildFastqsFofn.empty() && P.quant.slamSnpMask.buildFastqsFofn != "-" && P.quant.slamSnpMask.buildFastqsFofn != "None";
        bool hasBuildBam = !P.quant.slamSnpMask.buildBam.empty() && P.quant.slamSnpMask.buildBam != "-" && P.quant.slamSnpMask.buildBam != "None";
        
        if (hasBuildBam && !hasMaskIn && !hasBuildFastqs) {
            // Build mask from existing BAM
            // Note: Full integration pending - for now, users should call the external pileup_snp tool directly
            // See tools/pileup_snp/README.md for usage
            ostringstream errOut;
            errOut << "EXITING: --slamSnpMaskBuildBam is not yet fully integrated into STAR.\n"
                   << "SOLUTION: Use the external pileup_snp tool directly:\n"
                   << "  tools/pileup_snp/pileup_snp --bam <bam> --bed <candidates.bed> --ref <ref.fa> --output <mask.bed.gz>\n"
                   << "Then load the mask with --slamSnpMaskIn <mask.bed.gz>\n"
                   << "Future versions will support direct BAM-based mask building from within STAR.\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        if (hasVcfIn && !hasMaskIn) {
            P.inOut->logMain << timeMonthDayTime() << " ..... loading SNP mask from VCF\n" << flush;
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... loading SNP mask from VCF\n" << flush;

            SlamSnpMask maskFromVcf;
            SlamSnpMaskVcfOptions opts;
            opts.sample = P.quant.slamSnpMask.vcfSample;
            opts.mode = P.quant.slamSnpMask.vcfMode;
            opts.filter = P.quant.slamSnpMask.vcfFilter;
            opts.bedOut = P.quant.slamSnpMask.bedOut;
            opts.summaryOut = P.quant.slamSnpMask.summaryOut;
            SlamSnpMaskVcfStats stats;
            std::string err;
            if (!maskFromVcf.loadVcf(P.quant.slamSnpMask.vcfIn, *genomeMain.genomeOut.g, opts, &stats, &err)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal error loading VCF SNP mask: "
                       << P.quant.slamSnpMask.vcfIn << "\n";
                if (!err.empty()) {
                    errOut << "Details: " << err << "\n";
                }
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            P.inOut->logMain << timeMonthDayTime() << " ..... VCF mask loaded\n"
                            << "    records: " << stats.recordsTotal << "\n"
                            << "    filtered: " << stats.recordsFiltered << "\n"
                            << "    sites added: " << stats.sitesAdded << "\n";
            if (!P.quant.slamSnpMask.bedOut.empty()) {
                P.inOut->logMain << "    BED out: " << P.quant.slamSnpMask.bedOut << "\n";
            }
            if (!P.quant.slamSnpMask.summaryOut.empty()) {
                P.inOut->logMain << "    summary out: " << P.quant.slamSnpMask.summaryOut << "\n";
            }

            if (P.quant.slamSnpMask.buildOnly) {
                P.inOut->logMain << "Exiting after mask build (--slamSnpMaskOnly)\n" << flush;
                sysRemoveDir(P.outFileTmp);
                exit(0);
            }

            P.quant.slam.snpBed = P.quant.slamSnpMask.bedOut;
        } else if (hasBuildFastqs && !hasMaskIn) {
            // Build mask from FASTQs
            P.inOut->logMain << timeMonthDayTime() << " ..... starting SNP mask build pre-pass\n" << flush;
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... starting SNP mask build\n" << flush;
            
            SnpMaskBuild maskBuilder(P, *genomeMain.genomeOut.g);
            
            // Parse FOFN
            std::vector<std::pair<std::string, std::string>> fastqPairs;
            std::string err;
            if (!maskBuilder.parseFofn(P.quant.slamSnpMask.buildFastqsFofn, fastqPairs, &err)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: failed to parse FOFN "
                       << P.quant.slamSnpMask.buildFastqsFofn << "\n";
                if (!err.empty()) {
                    errOut << "Details: " << err << "\n";
                }
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            
            P.inOut->logMain << "Parsed " << fastqPairs.size() << " FASTQ pair(s) from FOFN\n";
            
            // Set up FASTQ files for mask build alignment
            // Save original readFilesNames
            std::vector<std::vector<std::string> > originalReadFilesNames = P.readFilesNames;
            P.readFilesNames.clear();
            P.readFilesNames.resize(1);  // One mate set
            
            // Convert fastqPairs to readFilesNames format
            for (const auto& pair : fastqPairs) {
                if (!pair.second.empty()) {
                    // Paired-end
                    P.readFilesNames[0].push_back(pair.first);
                    if (P.readFilesNames.size() < 2) {
                        P.readFilesNames.resize(2);
                    }
                    P.readFilesNames[1].push_back(pair.second);
                } else {
                    // Single-end
                    P.readFilesNames[0].push_back(pair.first);
                }
            }
            
            // Update readNends based on whether we have PE or SE
            P.readNends = (fastqPairs[0].second.empty()) ? 1 : 2;
            P.readNmates = P.readNends;
            
            // Re-initialize read files
            P.closeReadsFiles();
            P.readFilesInit();
            P.openReadsFiles();
            
            // Create ReadAlignChunk for mask build
            // Force-enable SLAM + SNP observation collection for this pre-pass, but do NOT swap SlamQuant
            // (swapping risks dangling pointers and also bypasses debug configuration done in ReadAlignChunk).
            bool originalSlamYes = P.quant.slam.yes;
            bool originalSnpDetect = P.quant.slam.snpDetect;
            P.quant.slam.yes = true;
            P.quant.slam.snpDetect = true;

            ReadAlignChunk* RAchunkMask = new ReadAlignChunk(P, genomeMain, transcriptomeMain, 0,
                                                            libem_transcriptome.get());
            
            // Run alignment to collect observations
            P.inOut->logMain << timeMonthDayTime() << " ..... aligning reads for mask build\n" << flush;
            RAchunkMask->processChunks();
            
            // Extract data from SlamQuant and run EM/filtering
            maskBuilder.extractFromSlamQuant(RAchunkMask->slamQuant);

            // If SNP-site debug was enabled, write debug outputs now (mask-only runs exit before SLAM quant writes).
            if (P.quant.slam.debugEnabled && RAchunkMask->slamQuant != nullptr) {
                RAchunkMask->slamQuant->writeDebug(*transcriptomeMain, P.quant.slam.errorRate, P.quant.slam.convRate);
            }
            
            // Restore original settings
            P.quant.slam.yes = originalSlamYes;
            P.quant.slam.snpDetect = originalSnpDetect;
            P.readFilesNames = originalReadFilesNames;
            // readNends is 1 for SE, 2 for PE (number of mate arrays, not file count)
            P.readNends = (originalReadFilesNames.size() > 1) ? 2 : 1;
            P.readNmates = P.readNends;
            
            // Close and re-open with original files
            P.closeReadsFiles();
            P.readFilesInit();
            if (P.runMode == "alignReads") {
                P.openReadsFiles();
            }
            
            // Clean up
            delete RAchunkMask;
            
            // Now run mask build processing (EM fit, filtering, etc.)
            // Log mask build parameters
            P.inOut->logMain << timeMonthDayTime() << " ..... building SNP mask using model: " 
                            << P.quant.slamSnpMask.model << "\n";
            P.inOut->logMain << "    mismatch counting mode (kMode): " << P.quant.slamSnpMask.kMode << "\n";
            if (P.quant.slamSnpMask.model == "binom") {
                double p_err = (P.quant.slamSnpMask.err >= 0.0) ? P.quant.slamSnpMask.err : P.quant.slam.snpErrUsed;
                if (p_err < P.quant.slam.snpErrMinThreshold) {
                    p_err = P.quant.slam.snpErrMinThreshold;
                }
                P.inOut->logMain << "    binomial model parameters:\n"
                                << "      p_err: " << std::fixed << std::setprecision(6) << p_err << "\n"
                                << "      pval threshold: " << std::fixed << std::setprecision(6) << P.quant.slamSnpMask.pval << "\n"
                                << "      minTcRatio: " << std::fixed << std::setprecision(3) << P.quant.slamSnpMask.minTcRatio << "\n"
                                << "      minCov: " << P.quant.slamSnpMask.minCov << "\n"
                                << "      minAlt: " << P.quant.slamSnpMask.minAlt << "\n";
            } else {
                P.inOut->logMain << "    EM model parameters:\n"
                                << "      posterior threshold: " << std::fixed << std::setprecision(3) << P.quant.slamSnpMask.posterior << "\n"
                                << "      minCov: " << P.quant.slamSnpMask.minCov << "\n"
                                << "      minAlt: " << P.quant.slamSnpMask.minAlt << "\n"
                                << "      maxIter: " << P.quant.slamSnpMask.maxIter << "\n";
            }
            
            SnpMaskBuildStats buildStats;
            if (!maskBuilder.buildMask(fastqPairs, transcriptomeMain, &buildStats, &err)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal error during mask build\n";
                if (!err.empty()) {
                    errOut << "Details: " << err << "\n";
                }
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            
            // Log mask build results
            P.inOut->logMain << timeMonthDayTime() << " ..... SNP mask build completed\n"
                            << "    total sites: " << buildStats.totalSites << "\n"
                            << "    candidate sites: " << buildStats.candidateSites << "\n"
                            << "    masked sites: " << buildStats.maskedSites << "\n"
                            << "    global baseline before: " << std::fixed << std::setprecision(6) << buildStats.globalBaselineBefore << "\n"
                            << "    global baseline after: " << std::fixed << std::setprecision(6) << buildStats.globalBaselineAfter << "\n";
            if (P.quant.slamSnpMask.model == "em") {
                const auto& emResult = maskBuilder.getEMResult();
                P.inOut->logMain << "    EM iterations: " << emResult.iterations << "\n"
                                << "    EM converged: " << (emResult.converged ? "yes" : "no") << "\n"
                                << "    p_ERR: " << std::fixed << std::setprecision(6) << emResult.p_ERR << "\n";
            }
            
            // Write outputs
            if (!maskBuilder.writeBed(P.quant.slamSnpMask.bedOut, *genomeMain.genomeOut.g, &err)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal error writing BED output\n";
                if (!err.empty()) {
                    errOut << "Details: " << err << "\n";
                }
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            
            if (!P.quant.slamSnpMask.summaryOut.empty()) {
                if (!maskBuilder.writeSummary(P.quant.slamSnpMask.summaryOut, buildStats, &err)) {
                    P.inOut->logMain << "WARNING: failed to write summary: " << err << "\n";
                }
            }
            
            P.inOut->logMain << timeMonthDayTime() << " ..... finished SNP mask build\n" << flush;
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished SNP mask build\n" << flush;
            
            // If --slamSnpMaskOnly, exit here
            if (P.quant.slamSnpMask.buildOnly) {
                P.inOut->logMain << "Exiting after mask build (--slamSnpMaskOnly)\n" << flush;
                sysRemoveDir(P.outFileTmp);
                exit(0);
            }
            
            // Otherwise, load the newly built mask for use in main SLAM run
            P.quant.slam.snpBed = P.quant.slamSnpMask.bedOut;
        }
        
        // Load SNP mask (either from --slamSnpMaskIn or from newly built mask)
        if (hasMaskIn) {
            P.quant.slam.snpBed = P.quant.slamSnpMask.maskIn;
        }
        
        if (P.quant.slam.yes && !P.quant.slam.snpBed.empty()) {
            P.quant.slam.snpMask = new SlamSnpMask();
            std::string err;
            if (!P.quant.slam.snpMask->loadBed(P.quant.slam.snpBed, *genomeMain.genomeOut.g, &err)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: failed to load slamSnpBed "
                       << P.quant.slam.snpBed << "\n";
                if (!err.empty()) {
                    errOut << "Details: " << err << "\n";
                }
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            P.inOut->logMain << "Loaded SLAM SNP BED (" << P.quant.slam.snpMask->size()
                             << " positions): " << P.quant.slam.snpBed << "\n";
        }
        
        // Validate --trimSource if provided
        if (P.quant.slam.yes && !P.quant.slam.trimSource.empty() && P.quant.slam.trimSource != "-") {
            // Check if trimSource file exists and is readable
            struct stat trimSourceStat;
            if (stat(P.quant.slam.trimSource.c_str(), &trimSourceStat) != 0) {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: --trimSource file does not exist or is not accessible: "
                       << P.quant.slam.trimSource << "\n"
                       << "SOLUTION: check the file path and permissions.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            
            // Check if trimScope is per-file (trimSource is ignored in that case)
            if (P.quant.slam.trimScope == "per-file") {
                P.inOut->logMain << "WARNING: --trimSource is ignored when --trimScope per-file is set. "
                                 << "Per-file mode computes trims independently for each file.\n";
                P.quant.slam.trimSource = "";  // Clear to avoid confusion
            } else {
                P.inOut->logMain << "SLAM auto-trim: using --trimSource file for trim computation: " 
                                 << P.quant.slam.trimSource << "\n";
            }
        }
        
        // Load transcript sequences for error model if enabled
        if (P.quant.transcriptVB.yes && P.quant.transcriptVB.errorModelMode != "off") {
            // Determine FASTA path: P.pGe.transcriptomeFasta else transcriptome.fa
            std::string fasta_path;
            if (!P.pGe.transcriptomeFasta.empty() && P.pGe.transcriptomeFasta != "-") {
                fasta_path = P.pGe.transcriptomeFasta;
            } else {
                // Try transcriptome.fa in genome directory
                fasta_path = P.pGe.gDir + "/transcriptome.fa";
            }
            
            libem_transcriptome.reset(new libem::Transcriptome());
            if (!libem_transcriptome->loadFromFasta(fasta_path)) {
                // Failed to load - disable error model
                P.inOut->logMain << "WARNING: Failed to load transcript sequences from " << fasta_path 
                                 << ". Error model will be disabled.\n";
                libem_transcriptome.reset();
            } else {
                // Reorder by STAR transcript names to match BAM header order
                std::vector<std::string> star_names(transcriptomeMain->nTr);
                for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
                    star_names[i] = transcriptomeMain->trID[i];
                }
                if (!libem_transcriptome->reorderByNames(star_names)) {
                    P.inOut->logMain << "WARNING: Failed to reorder transcript sequences to match STAR order. "
                                     << "Error model may use incorrect sequences.\n";
                }
                P.inOut->logMain << "Loaded " << libem_transcriptome->size() 
                                 << " transcript sequences for error model\n";
            }
        }
    }

    // Pre-initialize inline CB correction whitelist (needed before mapping)
    if (P.pSolo.inlineCBCorrection && P.pSolo.cbWLyes && !P.pSolo.cbWLstr.empty()) {
        InlineCBCorrection::initializeWhitelist(P.pSolo);
        P.inOut->logMain << "[INLINE-CB-INIT] size=" << P.pSolo.cbWLstr.size()
                         << " exact_map=" << InlineCBCorrection::exactMapSize()
                         << " variant_map=" << InlineCBCorrection::variantMapSize()
                         << " variant_collisions=" << InlineCBCorrection::variantCollisionSize()
                         << " collision_max_fanout=" << InlineCBCorrection::variantCollisionMaxFanout()
                         << "\n";
    }

    // Initialize global gene→probe cache before mapping (needed by inline resolver)
    if (transcriptomeMain != nullptr && P.pSolo.type != 0) {
        if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
            ProbeListIndex probeIdxInit;
            uint32_t deprecatedCount = 0;
            if (probeIdxInit.load(P.pSolo.probeListPath, P.pSolo.removeDeprecated, &deprecatedCount)) {
                SoloFeature::initGeneProbeIdx(*transcriptomeMain, &probeIdxInit);
                P.inOut->logMain << "[GENE-PROBE-INIT] Pre-mapping init done from " << P.pSolo.probeListPath << "\n";
                if (P.pSolo.removeDeprecated && deprecatedCount > 0) {
                    P.inOut->logMain << "[GENE-PROBE-INIT] Removed " << deprecatedCount << " deprecated entries from probe list\n";
                }
            } else {
                P.inOut->logMain << "[GENE-PROBE-INIT] WARNING: failed to load probe list at " << P.pSolo.probeListPath << "\n";
            }
        } else {
            P.inOut->logMain << "[GENE-PROBE-INIT] Skipped pre-mapping init (no probe list path)\n";
        }
    }

    // initialize Stats
    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started mapping\n"
                        << flush;
    g_statsAll.timeLastReport = g_statsAll.timeStartMap;

    // SAM headers
    samHeaders(P, *genomeMain.genomeOut.g, transcriptomeMain);

    // initialize chimeric parameters here - note that chimeric parameters require samHeader
    P.pCh.initialize(&P);

    // Apply limitBAMsortRAM fallback before creating SamtoolsSorter
    if (P.outBAMcoord && P.limitBAMsortRAM == 0) {
        // make it equal to the genome size
        P.limitBAMsortRAM = genomeMain.nGenome + genomeMain.SA.lengthByte + genomeMain.SAi.lengthByte;
    }

    // Initialize samtools sorter if needed
    if (P.outBAMsortMethod == "samtools" && P.outBAMcoord) {
        g_samtoolsSorter = new SamtoolsSorter(P.limitBAMsortRAM, 
                                               P.outBAMsortingThreadNactual,
                                               P.outBAMsortTmpDir,
                                               P);
    }

    // this does not seem to work at the moment
    // P.inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

    // === SLAM STATS COLLECTION: trimScope=first (single detection pass, then full mapping) ===
    // Always run detection pass when SLAM is enabled to collect variance stats and compute error rate
    // Trims are only computed/applied if --autoTrim variance is set
    volatile bool runSlamDetectionPass = true;  // keep main mapping flow independent of SLAM detection
    bool doSlamDetectionPass = (P.quant.slam.yes && P.quant.slam.trimScope == "first" && !P.quant.slam.autoTrimComputed);
    
    if (runSlamDetectionPass) {
        if (doSlamDetectionPass) {
        
        // Determine trim source: use --trimSource if provided, otherwise first input file
        bool usingTrimSource = !P.quant.slam.trimSource.empty() && P.quant.slam.trimSource != "-";
        std::string trimSourcePath = usingTrimSource ? P.quant.slam.trimSource : 
            (P.readFilesNames[0].size() > 0 ? P.readFilesNames[0][0] : "");
        
        if (usingTrimSource) {
            P.inOut->logMain << timeMonthDayTime() << " ..... starting SLAM stats collection (single-threaded, --trimSource)\n" << flush;
            P.inOut->logMain << "    trim_source_file=" << trimSourcePath << "\n";
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... starting SLAM stats collection (from --trimSource)\n" << flush;
        } else {
            P.inOut->logMain << timeMonthDayTime() << " ..... starting SLAM stats collection (single-threaded, trimScope=first)\n" << flush;
            P.inOut->logMain << "    trim_source_file=" << trimSourcePath << " (first input file)\n";
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... starting SLAM stats collection\n" << flush;
        }
        
        // Save original settings
        uint64_t originalReadMapNumber = P.readMapNumber;
        int originalRunThreadN = P.runThreadN;
        std::vector<std::vector<std::string>> originalReadFilesNames;
        
        // If using --trimSource, temporarily replace input files
        if (usingTrimSource) {
            originalReadFilesNames = P.readFilesNames;
            
            // Close current input files
            P.closeReadsFiles();
            
            // Set up single-file input for detection
            P.readFilesNames[0].clear();
            P.readFilesNames[0].push_back(P.quant.slam.trimSource);
            if (P.readNends > 1) {
                // For paired-end, we only use R1 for detection (trimSource should be R1)
                P.readFilesNames[1].clear();
                P.readFilesNames[1].push_back(P.quant.slam.trimSource);  // Same file for both mates
            }
            
            // Reset counters and reopen with trimSource
            P.iReadAll = 0;
            g_threadChunks.chunkInN = 0;
            g_threadChunks.chunkOutN = 0;
            P.readFilesIndex = 0;
            g_bamRecordIndex.store(0);
            P.openReadsFiles();
        }
        
        // Set detection mode: single thread, limited reads
        P.quant.slam.autoTrimDetectionPass = true;
        P.runThreadN = 1;  // Single-threaded for detection
        P.readMapNumber = static_cast<uint64_t>(P.quant.slam.autoTrimDetectionReads);
        if (P.readMapNumber > originalReadMapNumber && originalReadMapNumber > 0) {
            P.readMapNumber = originalReadMapNumber;  // Don't exceed original limit
        }
        
        // Create detection-only RAchunk (single thread)
        ReadAlignChunk *RAdetect = new ReadAlignChunk(P, genomeMain, transcriptomeMain, 0, 
                                                       libem_transcriptome.get());
        
        // Run detection pass - only collecting variance stats (no outputs)
        RAdetect->processChunks();
        
        // Check how many reads were actually processed
        uint64_t detectionReadsProcessed = P.iReadAll;
        
        // Compute global error rate and trims from variance stats
        if (RAdetect->slamQuant != nullptr && RAdetect->slamQuant->varianceAnalysisEnabled()) {
            // Estimate read length from processed reads
            uint32_t readLength = 100;
            if (RAdetect->RA != nullptr && RAdetect->RA->readLength[0] > 0) {
                readLength = static_cast<uint32_t>(RAdetect->RA->readLength[0] + RAdetect->RA->readLength[1]);
            }
            
            const SlamVarianceAnalyzer* analyzer = RAdetect->slamQuant->varianceAnalyzer();
            
            // Compute global T→C error rate (always, even if trimming disabled)
            uint64_t t_total = 0, tc_total = 0;
            double p_est = 0.0;
            // Initialize trim windows from manual trims (will be updated if auto-trim succeeds)
            int trim5p_for_err = P.quant.slam.compatTrim5p;
            int trim3p_for_err = P.quant.slam.compatTrim3p;
            
            if (analyzer != nullptr) {
                // Compute trims only if auto-trim mode is explicitly enabled
                SlamVarianceTrimResult trimResult;
                bool trimComputed = false;
                
                if (P.quant.slam.autoTrimMode == "variance") {
                    trimResult = RAdetect->slamQuant->computeVarianceTrim(readLength);

                    if (!trimResult.success || trimResult.readsAnalyzed < static_cast<uint64_t>(P.quant.slam.autoTrimMinReads)) {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL ERROR: --autoTrim variance requested but insufficient stdev data to compute trims.\n"
                               << "DETAILS: reads_analyzed=" << trimResult.readsAnalyzed
                               << " detection_reads_processed=" << detectionReadsProcessed
                               << " min_reads=" << P.quant.slam.autoTrimMinReads << "\n"
                               << "SOLUTION: increase --autoTrimDetectionReads/--autoTrimMinReads, "
                               << "or ensure the detection pass collects T->C stdev data.\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
                    }

                    P.quant.slam.autoTrim5p = trimResult.trim5p;
                    P.quant.slam.autoTrim3p = trimResult.trim3p;
                    P.quant.slam.autoTrimComputed = true;
                    P.quant.slam.autoTrimFileIndex = 0;

                    // Update manual trim settings so SlamCompat uses them
                    P.quant.slam.compatTrim5p = trimResult.trim5p;
                    P.quant.slam.compatTrim3p = trimResult.trim3p;

                    trim5p_for_err = trimResult.trim5p;
                    trim3p_for_err = trimResult.trim3p;
                    trimComputed = true;

                    P.inOut->logMain << "SLAM auto-trim (segmented regression) computed:\n"
                                     << "    trim5p=" << trimResult.trim5p
                                     << " trim3p=" << trimResult.trim3p
                                     << " mode=" << trimResult.mode << "\n"
                                     << "    breakpoints: b1=" << trimResult.kneeBin5p
                                     << " b2=" << trimResult.kneeBin3p
                                     << " total_sse=" << trimResult.totalSSE << "\n"
                                     << "    reads_analyzed=" << trimResult.readsAnalyzed
                                     << " detection_reads_processed=" << detectionReadsProcessed
                                     << " min_reads=" << P.quant.slam.autoTrimMinReads << "\n"
                                     << "    smooth_window=" << P.quant.slam.autoTrimSmoothWindow
                                     << " min_seg_len=" << P.quant.slam.autoTrimSegMinLen
                                     << " max_trim=" << P.quant.slam.autoTrimMaxTrim
                                     << " scope=" << P.quant.slam.trimScope << "\n"
                                     << "    trim_source=" << trimSourcePath
                                     << (usingTrimSource ? " (--trimSource)" : " (first input)") << "\n";
                }
                
                // Compute error rate: use trimmed window if trims were computed, otherwise full window
                std::tie(t_total, tc_total, p_est) = analyzer->computeGlobalTcErrorRate(
                    trim5p_for_err, trim3p_for_err, readLength);
                P.quant.slam.snpErrEst = p_est;
                
                // Apply fallback threshold
                if (p_est >= P.quant.slam.snpErrMinThreshold) {
                    P.quant.slam.snpErrUsed = p_est;
                    P.quant.slam.snpErrFallbackReason = "";
                } else {
                    P.quant.slam.snpErrUsed = P.quant.slam.snpErrMinThreshold;
                    P.quant.slam.snpErrFallbackReason = "p_est < threshold";
                }
                
                P.inOut->logMain << "SLAM global T→C error rate (from detection pass";
                if (trimComputed) {
                    P.inOut->logMain << ", trimmed window";
                }
                P.inOut->logMain << "):\n"
                                 << "    t_total=" << t_total << " tc_total=" << tc_total
                                 << " p_est=" << std::fixed << std::setprecision(6) << p_est
                                 << " p_used=" << std::fixed << std::setprecision(6) << P.quant.slam.snpErrUsed;
                if (!P.quant.slam.snpErrFallbackReason.empty()) {
                    P.inOut->logMain << " (fallback: " << P.quant.slam.snpErrFallbackReason << ")";
                }
                P.inOut->logMain << "\n";
                
                // Cache variance curve for comprehensive QC (detection pass)
                if (analyzer != nullptr) {
                    const auto& vstats = analyzer->getStats();
                    uint32_t maxPos = 0;
                    for (const auto& kv : vstats) {
                        if (kv.first > maxPos) {
                            maxPos = kv.first;
                        }
                    }
                    P.quant.slam.varianceStddevTcRate.assign(maxPos + 1, std::numeric_limits<double>::quiet_NaN());
                    for (const auto& kv : vstats) {
                        P.quant.slam.varianceStddevTcRate[kv.first] = kv.second.stddevTcRate();
                    }
                } else {
                    P.quant.slam.varianceStddevTcRate.clear();
                }

                // Write QC outputs (always write if analyzer exists, even when trim detection fails)
                if (analyzer != nullptr) {
                    std::string qcJsonPath = P.quant.slam.slamQcJson;
                    if (qcJsonPath.empty() || qcJsonPath == "-") {
                        qcJsonPath = P.outFileNamePrefix + "slam_qc.json";
                    }
                    // For stats-only mode or failed trim detection, pass null trimResult
                    const SlamVarianceTrimResult* trimResultPtr = trimComputed ? &trimResult : nullptr;
                    int trim5p_val = trimComputed ? trimResult.trim5p : trim5p_for_err;
                    int trim3p_val = trimComputed ? trimResult.trim3p : trim3p_for_err;
                    uint64_t readsAnalyzed_val = trimComputed ? trimResult.readsAnalyzed : analyzer->readsAnalyzed();
                    bool writeResult = writeSlamQcJson(*analyzer, qcJsonPath, 
                                        P.quant.slam.autoTrimFileIndex, P.quant.slam.trimScope,
                                        trim5p_val, trim3p_val, readsAnalyzed_val, trimResultPtr,
                                        trimSourcePath,
                                        P.quant.slam.snpErrEst, P.quant.slam.snpErrUsed, P.quant.slam.snpErrFallbackReason);
                    if (writeResult) {
                        P.inOut->logMain << "SLAM QC JSON written to: " << qcJsonPath << "\n";
                        
                        // Only write HTML if trim was computed (not for stats-only or failed trim)
                        if (trimComputed) {
                            std::string qcHtmlPath = P.quant.slam.slamQcHtml;
                            if (qcHtmlPath.empty() || qcHtmlPath == "-") {
                                qcHtmlPath = P.outFileNamePrefix + "slam_qc.html";
                            }
                            if (writeSlamQcHtml(qcJsonPath, qcHtmlPath, P.quant.slam.autoTrimFileIndex)) {
                                P.inOut->logMain << "SLAM QC HTML written to: " << qcHtmlPath << "\n";
                            }
                        }
                    }
                }
            }
        
        // Clean up detection RAchunk
        delete RAdetect;
        
        // End detection mode
        P.quant.slam.autoTrimDetectionPass = false;
        
        // Restore settings
        P.runThreadN = originalRunThreadN;
        P.readMapNumber = originalReadMapNumber;
        
        // REWIND: Close and reopen input files to start from beginning
        P.closeReadsFiles();
        
        // If we used --trimSource, restore original input files
        if (usingTrimSource) {
            P.readFilesNames = originalReadFilesNames;
            P.inOut->logMain << "SLAM auto-trim: restoring original input files for main mapping pass\n";
        } else {
            P.inOut->logMain << "SLAM auto-trim: rewinding input files for main mapping pass\n";
        }
        
        P.iReadAll = 0;  // Reset read counter
        g_threadChunks.chunkInN = 0;  // Reset chunk counter
        g_threadChunks.chunkOutN = 0;  // Reset output chunk counter
        P.readFilesIndex = 0;  // Reset file index
        g_bamRecordIndex.store(0);  // Reset BAM record index
        P.openReadsFiles();
        
        // Reset global stats after detection so final stats are clean
        g_statsAll.resetN();
        time(&g_statsAll.timeStartMap);
        g_statsAll.timeLastReport = g_statsAll.timeStartMap;
        
        if (P.quant.slam.autoTrimComputed) {
            P.inOut->logMain << timeMonthDayTime() << " ..... finished SLAM stats collection (trims computed)\n" << flush;
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished SLAM stats collection, "
                               << "trim5p=" << P.quant.slam.autoTrim5p << " trim3p=" << P.quant.slam.autoTrim3p
                               << " snp_err_used=" << std::fixed << std::setprecision(6) << P.quant.slam.snpErrUsed << "\n" << flush;
        } else {
            P.inOut->logMain << timeMonthDayTime() << " ..... finished SLAM stats collection\n" << flush;
            *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished SLAM stats collection, "
                               << "snp_err_used=" << std::fixed << std::setprecision(6) << P.quant.slam.snpErrUsed << "\n" << flush;
        }
        }
    }
    
    // Flag to track if per-file processing already ran mapping
    bool perFileMappingDone = false;
    
    // Declare RAchunk array early so it can be used by both per-file and normal paths
    ReadAlignChunk **RAchunk = new ReadAlignChunk*[P.runThreadN];
    for (int ii = 0; ii < P.runThreadN; ii++) {
        RAchunk[ii] = nullptr;
    }
    
    // === SLAM STATS COLLECTION: trimScope=per-file (detection + mapping per file) ===
    // Always run detection pass when SLAM is enabled to collect variance stats and compute error rate
    // Trims are only computed/applied if --autoTrim variance is set
    if (P.quant.slam.yes && P.quant.slam.trimScope == "per-file") {
        P.inOut->logMain << timeMonthDayTime() << " ..... starting SLAM stats collection with trimScope=per-file\n" << flush;
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... starting SLAM stats collection (per-file mode)\n" << flush;
        
        // Get total file count
        P.quant.slam.totalFileCount = static_cast<int>(P.readFilesNames[0].size());
        P.quant.slam.perFileProcessing = true;
        
        // Save original settings
        uint64_t originalReadMapNumber = P.readMapNumber;
        int originalRunThreadN = P.runThreadN;
        
        // Create RAchunks for mapping (will be reused for each file)
        for (int ii = 0; ii < P.runThreadN; ii++) {
            RAchunk[ii] = new ReadAlignChunk(P, genomeMain, transcriptomeMain, ii, 
                                              libem_transcriptome.get());
        }
        
        // Stats accumulator for tracking across files (detection reads excluded)
        Stats accumulatedStats;
        
        // Preserve the overall mapping start time for accurate speed/time reporting
        time_t overallTimeStartMap;
        time(&overallTimeStartMap);
        
        // Track cumulative read index offset for global read IDs across files
        uint64 cumulativeReadOffset = 0;
        
        // Process each file: detection + mapping
        for (int fileIdx = 0; fileIdx < P.quant.slam.totalFileCount; fileIdx++) {
            P.quant.slam.currentFileIndex = fileIdx;
            
            P.inOut->logMain << "\n=== SLAM per-file processing: file " << fileIdx 
                             << " of " << P.quant.slam.totalFileCount << " ===\n";
            
            // --- PHASE 1: Detection for this file ---
            P.inOut->logMain << timeMonthDayTime() << " ..... detection phase for file " << fileIdx << "\n";
            
            // Rewind and skip to the target file for detection
            P.closeReadsFiles();
            P.iReadAll = 0;
            g_threadChunks.chunkInN = 0;
            g_threadChunks.chunkOutN = 0;
            P.readFilesIndex = 0;
            P.openReadsFiles();
            
            P.quant.slam.skipToFileIndex = fileIdx;  // Skip to target file
            P.quant.slam.autoTrimDetectionPass = true;
            P.quant.slam.autoTrimComputed = false;  // Reset for this file
            P.runThreadN = 1;  // Single-threaded for detection
            P.readMapNumber = static_cast<uint64_t>(P.quant.slam.autoTrimDetectionReads);
            
            // Create detection-only RAchunk
            ReadAlignChunk *RAdetect = new ReadAlignChunk(P, genomeMain, transcriptomeMain, 0, 
                                                           libem_transcriptome.get());
            
            // Run detection pass for this file only
            RAdetect->processChunks();
            
            uint64_t detectionReadsProcessed = P.iReadAll;
            
            // Compute trims and global error rate from variance stats
            if (RAdetect->slamQuant != nullptr && RAdetect->slamQuant->varianceAnalysisEnabled()) {
                uint32_t readLength = 100;
                if (RAdetect->RA != nullptr && RAdetect->RA->readLength[0] > 0) {
                    readLength = static_cast<uint32_t>(RAdetect->RA->readLength[0] + RAdetect->RA->readLength[1]);
                }
                
                const SlamVarianceAnalyzer* analyzer = RAdetect->slamQuant->varianceAnalyzer();
                uint64_t t_total = 0, tc_total = 0;
                double p_est = 0.0;
                // Initialize trim windows from manual trims (will be updated if auto-trim succeeds)
                int trim5p_for_err = P.quant.slam.compatTrim5p;
                int trim3p_for_err = P.quant.slam.compatTrim3p;
            
            // Compute trims first (if auto-trim enabled, not stats-only)
            SlamVarianceTrimResult trimResult;
            bool trimComputed = false;
            
            if (P.quant.slam.autoTrimMode == "variance") {
                trimResult = RAdetect->slamQuant->computeVarianceTrim(readLength);

                if (!trimResult.success || trimResult.readsAnalyzed < static_cast<uint64_t>(P.quant.slam.autoTrimMinReads)) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR: --autoTrim variance requested but insufficient stdev data to compute trims.\n"
                           << "DETAILS: file_index=" << fileIdx
                           << " reads_analyzed=" << trimResult.readsAnalyzed
                           << " min_reads=" << P.quant.slam.autoTrimMinReads << "\n"
                           << "SOLUTION: increase --autoTrimDetectionReads/--autoTrimMinReads, "
                           << "or ensure the detection pass collects T->C stdev data.\n";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
                }

                P.quant.slam.autoTrim5p = trimResult.trim5p;
                P.quant.slam.autoTrim3p = trimResult.trim3p;
                P.quant.slam.autoTrimComputed = true;
                P.quant.slam.autoTrimFileIndex = fileIdx;
                P.quant.slam.compatTrim5p = trimResult.trim5p;
                P.quant.slam.compatTrim3p = trimResult.trim3p;

                trim5p_for_err = trimResult.trim5p;
                trim3p_for_err = trimResult.trim3p;
                trimComputed = true;

                P.inOut->logMain << "SLAM auto-trim (file " << fileIdx << ", segmented regression):\n"
                                 << "    trim5p=" << trimResult.trim5p
                                 << " trim3p=" << trimResult.trim3p
                                 << " mode=" << trimResult.mode << "\n"
                                 << "    breakpoints: b1=" << trimResult.kneeBin5p
                                 << " b2=" << trimResult.kneeBin3p
                                 << " total_sse=" << trimResult.totalSSE << "\n"
                                 << "    reads_analyzed=" << trimResult.readsAnalyzed << "\n";
            }
                
                // Compute error rate: use trimmed window if trims were computed, otherwise full window
                if (analyzer != nullptr) {
                    std::tie(t_total, tc_total, p_est) = analyzer->computeGlobalTcErrorRate(
                        trim5p_for_err, trim3p_for_err, readLength);
                    P.quant.slam.snpErrEst = p_est;
                    
                    // Apply fallback threshold
                    if (p_est >= P.quant.slam.snpErrMinThreshold) {
                        P.quant.slam.snpErrUsed = p_est;
                        P.quant.slam.snpErrFallbackReason = "";
                    } else {
                        P.quant.slam.snpErrUsed = P.quant.slam.snpErrMinThreshold;
                        P.quant.slam.snpErrFallbackReason = "p_est < threshold";
                    }
                    
                    P.inOut->logMain << "SLAM global T→C error rate (file " << fileIdx;
                    if (trimComputed) {
                        P.inOut->logMain << ", trimmed window";
                    }
                    P.inOut->logMain << "):\n"
                                     << "    t_total=" << t_total << " tc_total=" << tc_total
                                     << " p_est=" << std::fixed << std::setprecision(6) << p_est
                                     << " p_used=" << std::fixed << std::setprecision(6) << P.quant.slam.snpErrUsed;
                    if (!P.quant.slam.snpErrFallbackReason.empty()) {
                        P.inOut->logMain << " (fallback: " << P.quant.slam.snpErrFallbackReason << ")";
                    }
                    P.inOut->logMain << "\n";
                    
                    // Cache variance curve for comprehensive QC (per-file detection pass)
                    if (analyzer != nullptr) {
                        const auto& vstats = analyzer->getStats();
                        uint32_t maxPos = 0;
                        for (const auto& kv : vstats) {
                            if (kv.first > maxPos) {
                                maxPos = kv.first;
                            }
                        }
                        P.quant.slam.varianceStddevTcRate.assign(maxPos + 1, std::numeric_limits<double>::quiet_NaN());
                        for (const auto& kv : vstats) {
                            P.quant.slam.varianceStddevTcRate[kv.first] = kv.second.stddevTcRate();
                        }
                    } else {
                        P.quant.slam.varianceStddevTcRate.clear();
                    }

                    // Write per-file QC outputs (always write if analyzer exists)
                    if (analyzer != nullptr) {
                        std::string qcJsonPath = P.quant.slam.slamQcJson;
                        if (qcJsonPath.empty() || qcJsonPath == "-") {
                            qcJsonPath = P.outFileNamePrefix + "slam_qc_file" + std::to_string(fileIdx) + ".json";
                        } else if (P.quant.slam.totalFileCount > 1) {
                            // Append file index to user-specified path
                            size_t dotPos = qcJsonPath.rfind('.');
                            if (dotPos != std::string::npos) {
                                qcJsonPath = qcJsonPath.substr(0, dotPos) + "_file" + std::to_string(fileIdx) + qcJsonPath.substr(dotPos);
                            }
                        }
                        // For per-file mode, each file is its own trim source
                        std::string perFileTrimSource = P.readFilesNames[0].size() > static_cast<size_t>(fileIdx) ? 
                            P.readFilesNames[0][fileIdx] : "";
                        // For stats-only mode or failed trim detection, pass null trimResult
                        const SlamVarianceTrimResult* trimResultPtr = trimComputed ? &trimResult : nullptr;
                        int trim5p_val = trimComputed ? trimResult.trim5p : trim5p_for_err;
                        int trim3p_val = trimComputed ? trimResult.trim3p : trim3p_for_err;
                        uint64_t readsAnalyzed_val = trimComputed ? trimResult.readsAnalyzed : analyzer->readsAnalyzed();
                        if (writeSlamQcJson(*analyzer, qcJsonPath, fileIdx, P.quant.slam.trimScope,
                                            trim5p_val, trim3p_val, readsAnalyzed_val, trimResultPtr,
                                            perFileTrimSource,
                                            P.quant.slam.snpErrEst, P.quant.slam.snpErrUsed, P.quant.slam.snpErrFallbackReason)) {
                            P.inOut->logMain << "SLAM QC JSON written to: " << qcJsonPath << "\n";
                        }
                    }
                }
            }
            
            delete RAdetect;
            P.quant.slam.autoTrimDetectionPass = false;
            
            // Reset stats after each file's detection (detection reads should not show in final totals)
            // But preserve the overall timeStartMap for accurate speed/time reporting in Log.final.out
            g_statsAll.resetN();
            g_statsAll.timeStartMap = overallTimeStartMap;  // Use overall start time, not per-file
            time(&g_statsAll.timeLastReport);
            
            // Restore accumulated stats from previous files (saved at end of previous file's mapping)
            if (fileIdx > 0) {
                g_statsAll.addStats(accumulatedStats);
            }
            
            // --- PHASE 2: Rewind to start of this file and map ---
            P.inOut->logMain << timeMonthDayTime() << " ..... mapping phase for file " << fileIdx 
                             << " with trim5p=" << P.quant.slam.compatTrim5p 
                             << " trim3p=" << P.quant.slam.compatTrim3p << "\n";
            
            // Rewind to start of files and skip to target file
            P.closeReadsFiles();
            // Restore read indices from end of previous file (maintains global monotonic IDs)
            P.iReadAll = cumulativeReadOffset;
            // Note: g_bamRecordIndex is not reset - it continues monotonically
            g_threadChunks.chunkInN = 0;
            g_threadChunks.chunkOutN = 0;
            P.readFilesIndex = 0;
            P.openReadsFiles();
            
            P.quant.slam.skipToFileIndex = fileIdx;  // Skip to target file for mapping
            P.runThreadN = originalRunThreadN;
            P.readMapNumber = originalReadMapNumber;
            
            // Reinitialize SlamCompat for all threads with new trims
            for (int ii = 0; ii < P.runThreadN; ii++) {
                if (RAchunk[ii] != nullptr) {
                    RAchunk[ii]->reinitSlamCompat(P.quant.slam.compatTrim5p, P.quant.slam.compatTrim3p);
                }
            }
            
            // Run mapping for this file
            mapThreadsSpawn(P, RAchunk);
            
            P.inOut->logMain << timeMonthDayTime() << " ..... finished mapping file " << fileIdx << "\n";
            
            // Save accumulated stats at end of this file's mapping (for next file to restore)
            accumulatedStats = g_statsAll;
            
            // Save read indices for next file (to maintain global monotonic IDs)
            cumulativeReadOffset = P.iReadAll;
        }
        
        // DON'T delete RAchunk here - keep it for post-processing (SLAM quantification, etc.)
        // The RAchunk array will be deleted at the end along with the normal path
        
        P.quant.slam.perFileProcessing = false;
        P.runThreadN = originalRunThreadN;
        P.readMapNumber = originalReadMapNumber;
        
        P.inOut->logMain << timeMonthDayTime() << " ..... finished SLAM per-file processing\n" << flush;
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished SLAM per-file processing\n" << flush;
        
        perFileMappingDone = true;
    }

    // prepare chunks and spawn mapping threads (skip if per-file already handled mapping)
    if (!perFileMappingDone) {
        for (int ii = 0; ii < P.runThreadN; ii++)
        {
            RAchunk[ii] = new ReadAlignChunk(P, genomeMain, transcriptomeMain, ii, 
                                              libem_transcriptome.get());  // Pass shared transcriptome
        };
    }

    // === LIBRARY FORMAT DETECTION (single-threaded) ===
    if (P.quant.transcriptVB.yes && P.quant.transcriptVB.libType == "A") {
        P.inOut->logMain << "Starting library format auto-detection "
                         << "(first " << P.quant.transcriptVB.autoDetectWindow 
                         << " reads)...\n" << flush;
        
        // Create shared detector (will be accessed by TranscriptQuantEC during voting)
        P.quant.transcriptVB.libFormatDetector = new LibFormatDetector(
            P.quant.transcriptVB.autoDetectWindow);
        
        // Set detection mode BEFORE processing
        P.quant.transcriptVB.inDetectionMode = true;
        
        // Temporarily set read limit to detection window
        uint64_t originalLimit = P.readMapNumber;
        P.readMapNumber = P.quant.transcriptVB.autoDetectWindow;
        
        // Process first N reads using RAchunk[0] single-threaded
        // Voting happens INSIDE addReadAlignments() when inDetectionMode=true
        // NOTE: If total reads <= autoDetectWindow, RAchunk[0] will process all reads here
        // and then mapThreadsSpawn() will run again on RAchunk[0] with no remaining reads.
        // This is fine for datasets larger than autoDetectWindow, but avoid very small
        // read counts with auto-detect enabled to prevent potential duplicate output.
        RAchunk[0]->processChunks();
        
        // Restore original limit
        P.readMapNumber = originalLimit;
        
        // Finalize detection - HARD FAILURE if ambiguous (exits before returning)
        LibraryFormat detected = P.quant.transcriptVB.libFormatDetector->finalizeOrFail(
            P.inOut->logMain);
        
        // Store detected format as uint8_t
        P.quant.transcriptVB.detectedLibFormatId = detected.typeId();
        P.quant.transcriptVB.detectionComplete = true;
        P.quant.transcriptVB.inDetectionMode = false;  // Detection done
        
        // formatName() is defined in LibFormatDetection.h (see Section 9)
        P.inOut->logMain << "Detected library format: " << formatName(detected)
                         << " (formatID=" << static_cast<int>(detected.typeId()) << ")"
                         << " from " << P.iReadAll << " reads\n" << flush;
        
        // Clean up detector (no longer needed)
        delete P.quant.transcriptVB.libFormatDetector;
        P.quant.transcriptVB.libFormatDetector = nullptr;
        
        // Reset global counters after detection completes
        // This ensures pre-burn-in gating starts fresh for the main mapping pass (matches Salmon's behavior)
        extern std::atomic<uint64_t> global_processed_fragments;
        global_processed_fragments.store(0, std::memory_order_relaxed);
        Parameters::global_fld_obs_count.store(0, std::memory_order_relaxed);
    }

    if (P.runRestart.type != 1 && !perFileMappingDone)
        mapThreadsSpawn(P, RAchunk);

    if (P.outFilterBySJoutStage == 1 && !perFileMappingDone)
    { // completed stage 1, go to stage 2
        P.inOut->logMain << "Completed stage 1 mapping of outFilterBySJout mapping\n"
                         << flush;
        outputSJ(RAchunk, P); // collapse novel junctions
        P.readFilesIndex = -1;

        P.outFilterBySJoutStage = 2;
        if (P.outBAMcoord)
        {
            for (int it = 0; it < P.runThreadN; it++)
            { // prepare the unmapped bin
                RAchunk[it]->chunkOutBAMcoord->coordUnmappedPrepareBySJout();
            };
        };

        mapThreadsSpawn(P, RAchunk);
    };

    // close some BAM files
    if (P.inOut->outBAMfileUnsorted != NULL)
    {
        bgzf_flush(P.inOut->outBAMfileUnsorted);
        bgzf_close(P.inOut->outBAMfileUnsorted);
        if (P.emitNoYBAMyes) {
            if (P.inOut->outBAMfileY != NULL) bgzf_close(P.inOut->outBAMfileY);
            if (P.inOut->outBAMfileNoY != NULL) bgzf_close(P.inOut->outBAMfileNoY);
        }
    };
    if (P.inOut->outBAMfileUnsortedSoloTmp.is_open())
    {
        P.inOut->outBAMfileUnsortedSoloTmp.flush();
        P.inOut->outBAMfileUnsortedSoloTmp.close();
    };
    if (P.inOut->outQuantBAMfile != NULL)
    {
        bgzf_flush(P.inOut->outQuantBAMfile);
        bgzf_close(P.inOut->outQuantBAMfile);
    };


    time(&g_statsAll.timeFinishMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                        << flush;
    P.inOut->logMain << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                     << "RAM after mapping:\n"
                     << linuxProcMemory() << flush;

    // no need for genome anymore, free the memory
    genomeMain.freeMemory();

    // aggregate output junctions
    // collapse splice junctions from different threads/chunks, and output them
    if (P.runRestart.type != 1 && P.outSJ.yes)
        outputSJ(RAchunk, P);

    // solo counts
    Solo soloMain(RAchunk, P, *RAchunk[0]->chunkTr);
    
    // Skip Solo processing if inline replayer already ran (it produces MEX directly)
    if (!P.pSolo.useInlineReplayer) {
        soloMain.processAndOutput();
    } else {
        P.inOut->logMain << timeMonthDayTime() << " ..... skipping Solo processing (inline replayer already produced MEX)" << endl;
    }

    // Process Cell Ranger multi config if enabled
    if (!P.crMulti.crMultiConfig.empty()) {
        processCrMultiConfig(P);
    }

    // Note: Two-pass unsorted CB/UB tag injection removed - not used in inline flex path

    if (P.quant.geCount.yes)
    { // output gene quantifications
        for (int ichunk = 1; ichunk < P.runThreadN; ichunk++)
        { // sum counts from all chunks into 0th chunk
            RAchunk[0]->chunkTr->quants->addQuants(*(RAchunk[ichunk]->chunkTr->quants));
        };
        RAchunk[0]->chunkTr->quantsOutput();
    };

    if (P.runThreadN > 1 && P.outSAMorder == "PairedKeepInputOrder")
    { // concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };

    if (P.emitYReadNamesyes) {
        for (int ichunk = 0; ichunk < P.runThreadN; ichunk++) {
            if (RAchunk[ichunk] != nullptr &&
                RAchunk[ichunk]->RA != nullptr &&
                RAchunk[ichunk]->RA->chunkOutYReadNames.is_open()) {
                RAchunk[ichunk]->RA->chunkOutYReadNames.flush();
                RAchunk[ichunk]->RA->chunkOutYReadNames.close();
            }
        }

        ofstream yNamesOut(P.outYReadNamesFile.c_str());
        if (yNamesOut.fail()) {
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT ERROR: could not create Y read names file: "
                   << P.outYReadNamesFile << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        uint yNameChunk = 0;
        RAchunk[0]->chunkFilesCat(&yNamesOut, P.outFileTmp + "/YReadNames.out.thread", yNameChunk);
    }
    
    // Concatenate Y/noY FASTQ gzip files (if compression is enabled)
    if (P.emitYNoYFastqyes && P.emitYNoYFastqCompression == "gz") {
        for (uint32 imate = 0; imate < P.readNmates; imate++) {
            // Concatenate Y FASTQ gzip files
            ofstream yFastqOut(P.outYFastqFile[imate].c_str(), ios::binary);
            if (yFastqOut.fail()) {
                ostringstream errOut;
                errOut << "EXITING because of fatal OUTPUT ERROR: could not create Y FASTQ file: "
                       << P.outYFastqFile[imate] << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            for (int ichunk = 0; ichunk < P.runThreadN; ichunk++) {
                ostringstream chunkName;
                chunkName << P.outFileTmp << "/YFastq.mate" << imate << ".thread" << ichunk << ".gz";
                ifstream chunkIn(chunkName.str().c_str(), ios::binary);
                if (chunkIn.good()) {
                    yFastqOut << chunkIn.rdbuf();
                    chunkIn.close();
                    remove(chunkName.str().c_str());
                }
            }
            yFastqOut.close();
            
            // Concatenate noY FASTQ gzip files
            ofstream noYFastqOut(P.outNoYFastqFile[imate].c_str(), ios::binary);
            if (noYFastqOut.fail()) {
                ostringstream errOut;
                errOut << "EXITING because of fatal OUTPUT ERROR: could not create noY FASTQ file: "
                       << P.outNoYFastqFile[imate] << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
            for (int ichunk = 0; ichunk < P.runThreadN; ichunk++) {
                ostringstream chunkName;
                chunkName << P.outFileTmp << "/noYFastq.mate" << imate << ".thread" << ichunk << ".gz";
                ifstream chunkIn(chunkName.str().c_str(), ios::binary);
                if (chunkIn.good()) {
                    noYFastqOut << chunkIn.rdbuf();
                    chunkIn.close();
                    remove(chunkName.str().c_str());
                }
            }
            noYFastqOut.close();
        }
    }

    bamSortByCoordinate(P, RAchunk, *genomeMain.genomeOut.g, soloMain);
    
    // Transcript quantification (TranscriptVB mode)
    // TODO: Debug - check if transcriptomeMain is valid
    if (P.quant.transcriptVB.yes && transcriptomeMain != nullptr && transcriptomeMain->nTr > 0) {
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started transcript quantification\n"
                          << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started transcript quantification\n";
        
        // 1. Merge EC tables from all threads
        TranscriptQuantEC mergedEC(transcriptomeMain->nTr, 0, "", 0, P);
        for (int ichunk = 0; ichunk < P.runThreadN; ichunk++) {
            if (RAchunk[ichunk] != nullptr && 
                RAchunk[ichunk]->RA != nullptr && 
                RAchunk[ichunk]->RA->quantEC != nullptr) {
                mergedEC.merge(*RAchunk[ichunk]->RA->quantEC);
            }
        }
        
        mergedEC.finalize();
        
        *P.inOut->logStdOut << "Merged " << mergedEC.getECTable().n_ecs << " equivalence classes from " << P.runThreadN << " threads\n"
                          << flush;
        P.inOut->logMain << "Merged " << mergedEC.getECTable().n_ecs << " equivalence classes from " << P.runThreadN << " threads\n";
        
        // Log drop statistics (per trace plan Step 3)
        P.inOut->logMain << "EC building drop statistics:\n"
                         << "  dropped_incompat: " << mergedEC.getDroppedIncompat() << "\n"
                         << "  dropped_missing_mate_fields: " << mergedEC.getDroppedMissingMateFields() << "\n"
                         << "  dropped_unknown_obs_fmt: " << mergedEC.getDroppedUnknownObsFmt() << "\n";
        
        // Note: Transcript lengths are already set per-thread in ReadAlign.cpp constructor
        // and carried through merge(), so no need to set them again here.
        
        // 2. Initialize transcript state
        TranscriptState state;
        state.resize(transcriptomeMain->nTr);
        
        // Compute effective lengths using FLD (if available) or fallback to mean fragment length
        const FLDAccumulator& observedFLD = mergedEC.getObservedFLD();
        double meanFragLen = 200.0;  // Default fallback
        bool use_fld = observedFLD.isValid();
        
        if (use_fld) {
            meanFragLen = observedFLD.getMean();
            double fragLenStdDev = observedFLD.getStdDev();
            *P.inOut->logStdOut << "Fragment length distribution: mean=" << meanFragLen 
                              << ", stddev=" << fragLenStdDev 
                              << ", fragments=" << observedFLD.getTotalFragments() << "\n"
                              << flush;
            P.inOut->logMain << "Fragment length distribution: mean=" << meanFragLen 
                           << ", stddev=" << fragLenStdDev 
                           << ", fragments=" << observedFLD.getTotalFragments() << "\n";
        }
        
        // Log GC observations if GC bias is enabled
        if (P.quant.transcriptVB.gcBias) {
            const GCFragModel& observedGC = mergedEC.getObservedGC();
            const auto& gcCounts = observedGC.getCounts();
            double totalObs = 0.0;
            for (int i = 0; i < 101; i++) {
                totalObs += gcCounts[i];
            }
            if (totalObs > 100) {
                *P.inOut->logStdOut << "GC bias: collected " << static_cast<uint64_t>(totalObs) << " fragment observations\n"
                                  << flush;
                P.inOut->logMain << "GC bias: collected " << static_cast<uint64_t>(totalObs) << " fragment observations\n";
            }
        }
        
        // Compute effective lengths using FLD PMF (if available) or simple mean-based calculation
        std::vector<double> fld_pmf;
        if (use_fld) {
            fld_pmf = observedFLD.getPMF();
        }
        
        // Build raw_lengths_int vector
        std::vector<int32_t> raw_lengths_int(transcriptomeMain->nTr);
        for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
            raw_lengths_int[i] = static_cast<int32_t>(transcriptomeMain->trLen[i]);
        }
        
        // Compute effective lengths
        std::vector<double> eff_lengths;
        if (use_fld && !fld_pmf.empty()) {
            // Use EffectiveLengthCalculator via wrapper (avoids Transcriptome conflict)
            eff_lengths = computeEffectiveLengthsFromPMFWrapper(fld_pmf, raw_lengths_int);
        } else {
            // Fallback: simple mean-based calculation
            eff_lengths.resize(transcriptomeMain->nTr);
            for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
                double rawLen = static_cast<double>(transcriptomeMain->trLen[i]);
                double effLen = rawLen - meanFragLen + 1.0;
                if (effLen < 1.0) effLen = 1.0;
                if (effLen > rawLen) effLen = rawLen;
                eff_lengths[i] = effLen;
            }
        }
        
        // Populate transcript state
        for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
            state.names[i] = transcriptomeMain->trID[i];
            double rawLen = static_cast<double>(transcriptomeMain->trLen[i]);
            state.lengths[i] = rawLen;
            state.eff_lengths[i] = eff_lengths[i];
        }
        
        // Note: GC-corrected effective lengths would require transcript sequences
        // For now, FLD-based effective lengths are computed above
        if (P.quant.transcriptVB.gcBias) {
            *P.inOut->logStdOut << "GC bias: GC correction requires transcript sequences (not yet implemented)\n"
                              << flush;
            P.inOut->logMain << "GC bias: GC correction requires transcript sequences (not yet implemented)\n";
        }
        
        // 4. Run VB/EM quantification
        EMParams params;
        params.use_vb = P.quant.transcriptVB.vb;
        params.vb_prior = P.quant.transcriptVB.vbPrior;
        // Use defaults from EMParams (max_iters=10000, min_iters=100, tolerance=0.01)
        // Do NOT override for VB - let VB use same defaults as EM for Salmon parity
        // Thread count: use OMP default unless explicitly set (we'll pass --runThreadN externally for parity tests)
        params.threads = 0;  // 0 = use OMP default (multi-thread capable)
        
        EMResult result;
        if (params.use_vb) {
            result = run_vb(mergedEC.getECTable(), state, params);
        } else {
            result = run_em(mergedEC.getECTable(), state, params);
        }

        if (P.quant.slam.dumpWeightsMode == 1 && transcriptomeMain->nGe > 0) {
            vbGenePosterior.assign(transcriptomeMain->nGe, 0.0);
            for (uint32_t tr = 0; tr < transcriptomeMain->nTr; ++tr) {
                uint32_t gene = transcriptomeMain->trGene[tr];
                if (gene < transcriptomeMain->nGe && tr < result.counts.size()) {
                    vbGenePosterior[gene] += result.counts[tr];
                }
            }
            vbGenePosteriorReady = true;
        }
        
        *P.inOut->logStdOut << "Quantification converged: " << (result.converged ? "yes" : "no")
                          << ", iterations: " << result.iterations << "\n"
                          << flush;
        P.inOut->logMain << "Quantification converged: " << (result.converged ? "yes" : "no")
                         << ", iterations: " << result.iterations << "\n";
        
        // 5. Write quant.sf
        writeQuantSF(result, state, P.quant.transcriptVB.outFile);
        
        // 6. Write quant.genes.sf (gene-level aggregation, Legacy mode)
        if (P.quant.transcriptVB.geneOutput) {
            int geneResult = writeQuantGeneSF(result, state, *transcriptomeMain, 
                                               P.quant.transcriptVB.outFileGene);
            
            if (geneResult == 1) {
                // File open error
                P.inOut->logMain << "ERROR: Failed to open gene output file: " 
                                 << P.quant.transcriptVB.outFileGene << "\n";
            } else {
                P.inOut->logMain << "Gene quantification written to: " 
                                 << P.quant.transcriptVB.outFileGene << "\n";
                
                if (geneResult == 2) {
                    // MissingGeneID warning
                    P.inOut->logMain << "WARNING: transcripts with MissingGeneID were aggregated "
                                     << "to a single gene entry in quant.genes.sf\n";
                }
            }
            
            // 6b. Write tximport-style gene output if enabled
            if (P.quant.transcriptVB.genesTximport) {
                int tximportResult = writeQuantGeneSFTximport(result, state, *transcriptomeMain,
                                                              P.quant.transcriptVB.outFileGeneTximport);
                
                if (tximportResult == 1) {
                    P.inOut->logMain << "ERROR: Failed to open tximport gene output file: "
                                     << P.quant.transcriptVB.outFileGeneTximport << "\n";
                } else {
                    P.inOut->logMain << "Gene quantification (tximport mode) written to: "
                                     << P.quant.transcriptVB.outFileGeneTximport << "\n";
                }
            }
        }
        
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished transcript quantification\n"
                          << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... finished transcript quantification\n";
    }

    // SLAM quantification (gene-level)
    if (P.quant.slam.yes && transcriptomeMain != nullptr && transcriptomeMain->nGe > 0) {
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started SLAM quantification\n"
                            << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started SLAM quantification\n";

        // Merge counts from all threads
        // Note: With rewind-based auto-trim, trims are already applied during main mapping pass
        SlamQuant mergedSlam(transcriptomeMain->nGe, P.quant.slam.snpDetect);
        if (P.quant.slam.debugEnabled) {
            mergedSlam.initDebug(*transcriptomeMain, P.quant.slam.debugGenes, P.quant.slam.debugReads,
                                 static_cast<size_t>(P.quant.slam.debugMaxReads),
                                 P.quant.slam.debugOutPrefix);
        }
        
        for (int ichunk = 0; ichunk < P.runThreadN; ++ichunk) {
            if (RAchunk[ichunk] != nullptr && RAchunk[ichunk]->slamQuant != nullptr) {
                mergedSlam.merge(*RAchunk[ichunk]->slamQuant);
            }
        }
        
        SlamSnpBufferStats snpStats;
        mergedSlam.finalizeSnpMask(&snpStats);
        if (mergedSlam.snpDetectEnabled()) {
            P.inOut->logMain << "SLAM SNP detect: buffered_reads=" << snpStats.bufferedReads
                             << " avg_mismatches=" << snpStats.avgMismatches
                             << " avg_mismatches_kept=" << snpStats.avgMismatchesKept
                             << " buffer_bytes=" << snpStats.bufferBytes
                             << " snp_sites=" << snpStats.maskEntries
                             << " snp_blacklist=" << snpStats.blacklistEntries
                             << " mismatch_frac=" << snpStats.mismatchFracUsed
                             << " mode=" << snpStats.mismatchFracMode;
            if (snpStats.mismatchFracAuto > 0.0) {
                P.inOut->logMain << " auto_frac=" << snpStats.mismatchFracAuto
                                 << " knee_bin=" << snpStats.kneeBin
                                 << " eligible_sites=" << snpStats.eligibleSites;
            }
            P.inOut->logMain << "\n";
        }
        mergedSlam.write(*transcriptomeMain, P.quant.slam.outFile,
                         P.quant.slam.errorRate, P.quant.slam.convRate);
        if (P.quant.slam.grandSlamOut != 0) {
            std::string gsOut = P.outFileNamePrefix + "SlamQuant.grandslam.tsv";
            mergedSlam.writeGrandSlam(*transcriptomeMain, gsOut, P.outFileNamePrefix,
                                      P.quant.slam.errorRate, P.quant.slam.convRate);
        }

        // Optional: write dump for external re-quant
        if (!P.quant.slam.dumpBinary.empty() && P.quant.slam.dumpBinary != "-" &&
            P.quant.slam.dumpBinary != "None") {
            SlamDumpMetadata meta;
            meta.errorRate = P.quant.slam.errorRate;
            meta.convRate = P.quant.slam.convRate;
            meta.weightMode = P.quant.slam.weightMode;
            meta.geneIds = transcriptomeMain->geID;
            meta.geneNames = transcriptomeMain->geName;
            if (genomeMain.genomeOut.g != nullptr) {
                meta.chrNames = genomeMain.genomeOut.g->chrName;
                meta.chrStart.clear();
                meta.chrStart.reserve(genomeMain.genomeOut.g->chrStart.size());
                for (auto v : genomeMain.genomeOut.g->chrStart) {
                    meta.chrStart.push_back(static_cast<uint64_t>(v));
                }
            }
            std::vector<const SlamReadBuffer*> buffers;
            buffers.reserve(P.runThreadN);
            for (int ichunk = 0; ichunk < P.runThreadN; ++ichunk) {
                if (RAchunk[ichunk] != nullptr && RAchunk[ichunk]->slamQuant != nullptr) {
                    buffers.push_back(RAchunk[ichunk]->slamQuant->dumpBuffer());
                }
            }
            std::string dumpErr;
            if (writeSlamDump(P.quant.slam.dumpBinary, meta, buffers, P.quant.slam.dumpMaxReads, &dumpErr)) {
                P.inOut->logMain << "SLAM dump written to: " << P.quant.slam.dumpBinary << "\n";
            } else {
                P.inOut->logMain << "WARNING: failed to write SLAM dump: " << dumpErr << "\n";
            }
        }
        // Optional: write weight sidecar for external re-quant (uses same buffers/order).
        if (!P.quant.slam.dumpWeights.empty() && P.quant.slam.dumpWeights != "-" &&
            P.quant.slam.dumpWeights != "None") {
            SlamDumpMetadata meta;
            meta.errorRate = P.quant.slam.errorRate;
            meta.convRate = P.quant.slam.convRate;
            meta.weightMode = P.quant.slam.weightMode;
            std::vector<const SlamReadBuffer*> buffers;
            buffers.reserve(P.runThreadN);
            for (int ichunk = 0; ichunk < P.runThreadN; ++ichunk) {
                if (RAchunk[ichunk] != nullptr && RAchunk[ichunk]->slamQuant != nullptr) {
                    buffers.push_back(RAchunk[ichunk]->slamQuant->dumpBuffer());
                }
            }
            std::string wErr;
            std::vector<double> overrideWeights;
            const std::vector<double>* overridePtr = nullptr;
            if (P.quant.slam.dumpWeightsMode == 1) {
                if (!vbGenePosteriorReady) {
                    P.inOut->logMain << "WARNING: vbGene weights requested but TranscriptVB not available; "
                                     << "falling back to dump weights\n";
                } else {
                    if (!buildVbGeneWeights(buffers, vbGenePosterior, P.quant.slam.dumpMaxReads,
                                            &overrideWeights, &wErr)) {
                        P.inOut->logMain << "WARNING: failed to build vbGene weights: " << wErr
                                         << "; falling back to dump weights\n";
                    } else {
                        overridePtr = &overrideWeights;
                    }
                }
            }
            if (writeSlamWeights(P.quant.slam.dumpWeights, meta, buffers, P.quant.slam.dumpMaxReads,
                                 overridePtr, &wErr)) {
                P.inOut->logMain << "SLAM weights written to: " << P.quant.slam.dumpWeights << "\n";
            } else {
                P.inOut->logMain << "WARNING: failed to write SLAM weights: " << wErr << "\n";
            }
        }

        // Write diagnostics
        std::string diagFile = P.quant.slam.outFile + ".diagnostics";
        mergedSlam.writeDiagnostics(diagFile);

        // Write transition summaries (GEDI-style and compact)
        std::string transitionsFile = P.quant.slam.outFile + ".transitions.tsv";
        mergedSlam.writeTransitions(transitionsFile);
        std::string mismatchFile = P.quant.slam.outFile + ".mismatches.tsv";
        mergedSlam.writeMismatches(mismatchFile, P.outFileNamePrefix);
        std::string mismatchDetailsFile = P.quant.slam.outFile + ".mismatchdetails.tsv";
        mergedSlam.writeMismatchDetails(mismatchDetailsFile);
        if (P.quant.slam.debugEnabled) {
            mergedSlam.writeDebug(*transcriptomeMain, P.quant.slam.errorRate, P.quant.slam.convRate);
        }
        
        // Write top mismatches if reference file exists
        std::string refFile = P.pGe.gDir + "/../expected/fixture_ref_human.tsv.gz";
        std::ifstream refCheck(refFile.c_str());
        if (!refCheck.good()) {
            refFile = "test/fixtures/slam/expected/fixture_ref_human.tsv.gz";
            refCheck.open(refFile.c_str());
        }
        if (refCheck.good()) {
            refCheck.close();
            std::string mismatchFile = P.quant.slam.outFile + ".top_mismatches";
            mergedSlam.writeTopMismatches(*transcriptomeMain, refFile, mismatchFile, 20);
        }

        P.inOut->logMain << "SLAM quantification written to: "
                         << P.quant.slam.outFile << "\n";
        P.inOut->logMain << "SLAM diagnostics written to: "
                         << diagFile << "\n";
        
        // Log auto-trim summary
        if (P.quant.slam.autoTrimComputed) {
            P.inOut->logMain << "SLAM auto-trim applied: trim5p=" << P.quant.slam.autoTrim5p
                             << " trim3p=" << P.quant.slam.autoTrim3p
                             << " scope=" << P.quant.slam.trimScope
                             << " file_index=" << P.quant.slam.autoTrimFileIndex
                             << " mode=rewind\n";
        }
        
        // Log compat mode summary if enabled
        if (P.quant.slam.compatIntronic || P.quant.slam.compatLenientOverlap ||
            P.quant.slam.compatOverlapWeight || P.quant.slam.compatIgnoreOverlap ||
            P.quant.slam.compatTrim5p != 0 || P.quant.slam.compatTrim3p != 0 ||
            P.quant.slam.autoTrimComputed) {
            P.inOut->logMain << "SLAM compat(" << P.quant.slam.compatModeStr << "): "
                             << "alignsIntronic=" << mergedSlam.diagnostics().compatAlignsReclassifiedIntronic
                             << " alignsLenient=" << mergedSlam.diagnostics().compatAlignsLenientAccepted
                             << " alignsWeightAdj=" << mergedSlam.diagnostics().compatAlignsOverlapWeightApplied
                             << " posSkipOvlp=" << mergedSlam.diagnostics().compatPositionsSkippedOverlap
                             << " posSkipTrim=" << mergedSlam.diagnostics().compatPositionsSkippedTrim << "\n";
        }
        P.inOut->logMain << "SLAM transition summary written to: "
                         << transitionsFile << "\n";
        P.inOut->logMain << "SLAM mismatch summary written to: "
                         << mismatchFile << "\n";
        P.inOut->logMain << "SLAM mismatch details written to: "
                         << mismatchDetailsFile << "\n";
        if (P.quant.slam.debugEnabled) {
            P.inOut->logMain << "SLAM debug outputs written to: "
                             << P.quant.slam.debugOutPrefix << ".gene.tsv and "
                             << P.quant.slam.debugOutPrefix << ".reads.tsv\n";
        }
        
        // Generate comprehensive QC report if requested
        if (!P.quant.slam.slamQcReport.empty()) {
            std::string qcJsonPath = P.quant.slam.slamQcReport + ".slam_qc.json";
            std::string qcHtmlPath = P.quant.slam.slamQcReport + ".slam_qc.html";
            
            // Get trim result if available (recompute from merged SlamQuant if variance analysis enabled)
            SlamVarianceTrimResult* trimResultPtr = nullptr;
            SlamVarianceTrimResult trimResult;
            if (mergedSlam.varianceAnalysisEnabled() && P.quant.slam.autoTrimComputed) {
                // Use default read length (100) - trim result is optional for QC report
                uint32_t readLength = 100;
                trimResult = mergedSlam.computeVarianceTrim(readLength);
                if (trimResult.success) {
                    trimResultPtr = &trimResult;
                }
            }
            
            int trim5p = P.quant.slam.autoTrimComputed ? P.quant.slam.autoTrim5p : 0;
            int trim3p = P.quant.slam.autoTrimComputed ? P.quant.slam.autoTrim3p : 0;
            
            const std::vector<double>* varianceCurve = nullptr;
            if (!P.quant.slam.varianceStddevTcRate.empty()) {
                varianceCurve = &P.quant.slam.varianceStddevTcRate;
            }
            if (writeSlamQcComprehensiveJson(mergedSlam, qcJsonPath, trim5p, trim3p, trimResultPtr, varianceCurve)) {
                P.inOut->logMain << "SLAM comprehensive QC JSON written to: " << qcJsonPath << "\n";
                
                if (writeSlamQcComprehensiveHtml(qcJsonPath, qcHtmlPath)) {
                    P.inOut->logMain << "SLAM comprehensive QC HTML written to: " << qcHtmlPath << "\n";
                } else {
                    P.inOut->logMain << "WARNING: Failed to write SLAM comprehensive QC HTML\n";
                }
            } else {
                P.inOut->logMain << "WARNING: Failed to write SLAM comprehensive QC JSON\n";
            }
        }
        
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished SLAM quantification\n"
                            << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... finished SLAM quantification\n";
    }

    // Trim QC report (FastQC-like), optional and independent of trimming
    if (P.trimQcEnabled) {
        TrimQcCollector mergedTrimQc;
        bool hasTrimQc = false;
        for (int ichunk = 0; ichunk < P.runThreadN; ++ichunk) {
            if (RAchunk[ichunk] != nullptr && RAchunk[ichunk]->RA != nullptr &&
                RAchunk[ichunk]->RA->trimQc.enabled()) {
                if (!hasTrimQc) {
                    mergedTrimQc = RAchunk[ichunk]->RA->trimQc;
                    hasTrimQc = true;
                } else {
                    mergedTrimQc.merge(RAchunk[ichunk]->RA->trimQc);
                }
            }
        }
        if (hasTrimQc) {
            bool wantHtml = (!P.trimQcReport.empty() && P.trimQcReport != "-") ||
                            (!P.trimQcHtml.empty() && P.trimQcHtml != "-");
            bool wantJson = (!P.trimQcReport.empty() && P.trimQcReport != "-") ||
                            (!P.trimQcJson.empty() && P.trimQcJson != "-") || wantHtml;

            std::string qcJsonPath = P.trimQcJson;
            std::string qcHtmlPath = P.trimQcHtml;
            if (!P.trimQcReport.empty() && P.trimQcReport != "-") {
                qcJsonPath = P.trimQcReport + ".trim_qc.json";
                qcHtmlPath = P.trimQcReport + ".trim_qc.html";
            }
            if (qcJsonPath.empty() || qcJsonPath == "-") {
                qcJsonPath = P.outFileNamePrefix + "trim_qc.json";
            }
            if (qcHtmlPath.empty() || qcHtmlPath == "-") {
                qcHtmlPath = P.outFileNamePrefix + "trim_qc.html";
            }

            std::string stageLabel = (P.trimCutadapt == "Yes") ? "trimmed" : "raw";
            if (wantJson) {
                if (writeTrimQcJson(mergedTrimQc, g_statsAll, qcJsonPath, stageLabel)) {
                    P.inOut->logMain << "Trim QC JSON written to: " << qcJsonPath << "\n";
                } else {
                    P.inOut->logMain << "WARNING: Failed to write trim QC JSON\n";
                }
            }
            if (wantHtml) {
                if (writeTrimQcHtml(qcJsonPath, qcHtmlPath)) {
                    P.inOut->logMain << "Trim QC HTML written to: " << qcHtmlPath << "\n";
                } else {
                    P.inOut->logMain << "WARNING: Failed to write trim QC HTML\n";
                }
            }
        }
    }

    if (RAchunk != nullptr) {
        for (int ichunk = 0; ichunk < P.runThreadN; ++ichunk) {
            delete RAchunk[ichunk];
            RAchunk[ichunk] = nullptr;
        }
        delete[] RAchunk;
        RAchunk = nullptr;
    }

    if (transcriptomeMain != nullptr) {
        delete transcriptomeMain;
        transcriptomeMain = nullptr;
    }

    // Variation objects are heap-allocated and stored as raw pointers in Genome.
    // Genome is copied in several places (e.g., 2-pass), so we avoid a Genome
    // destructor and instead explicitly delete the Variation objects at shutdown.
    if (genomeMain.Var != nullptr) {
        delete genomeMain.Var;
        genomeMain.Var = nullptr;
    }
    if (genomeMain.genomeOut.g != nullptr && genomeMain.genomeOut.g != &genomeMain) {
        if (genomeMain.genomeOut.g->Var != nullptr) {
            delete genomeMain.genomeOut.g->Var;
            genomeMain.genomeOut.g->Var = nullptr;
        }
    }

    // Free auxiliary genome index buffers after all outputs that require them (e.g., SJ).
    delete[] genomeMain.chrBin;
    genomeMain.chrBin = nullptr;
    delete[] genomeMain.genomeSAindexStart;
    genomeMain.genomeSAindexStart = nullptr;

    // wiggle output
    if (P.outWigFlags.yes)
    {
        *(P.inOut->logStdOut) << timeMonthDayTime() << " ..... started wiggle output\n"
                              << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started wiggle output\n"
                         << flush;
        string wigOutFileNamePrefix = P.outFileNamePrefix + "Signal";
        signalFromBAM(P.outBAMfileCoordName, wigOutFileNamePrefix, P);
    }

    g_statsAll.writeLines(P.inOut->outChimJunction, P.pCh.outJunctionFormat, "#", STAR_VERSION + string("   ") + P.commandLine);

    g_statsAll.progressReport(P.inOut->logProgress);
    P.inOut->logProgress << "ALL DONE!\n"
                         << flush;
    P.inOut->logFinal.open((P.outFileNamePrefix + "Log.final.out").c_str());
    g_statsAll.reportFinal(P.inOut->logFinal);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... finished successfully\n"
                        << flush;

    P.inOut->logMain << "ALL DONE!\n"
                     << flush;
    if (P.outTmpKeep == "None")
    {
        sysRemoveDir(P.outFileTmp);
    };

    P.closeReadsFiles(); // this will kill the readFilesCommand processes if necessary
    // genomeMain.~Genome(); //need explicit call because of the 'delete P.inOut' below, which will destroy P.inOut->logStdOut
    if (genomeMain.sharedMemory != NULL)
    { // need explicit call because this destructor will write to files which are deleted by 'delete P.inOut' below
        delete genomeMain.sharedMemory;
        genomeMain.sharedMemory = NULL;
    };

    // Cleanup samtools sorter if it exists
    if (g_samtoolsSorter != nullptr) {
        delete g_samtoolsSorter;
        g_samtoolsSorter = nullptr;
    }

    delete P.inOut; // to close files

    return 0;
}
}
