#include "bamSortByCoordinate.h"
#include "BAMfunctions.h"
#include "BAMbinSortByCoordinate.h"
#include "BAMbinSortUnmapped.h"
#include "ErrorWarning.h"
#include "bam_cat.h"
#include "SamtoolsSorter.h"
#include "GlobalVariables.h"
#include "OcmMultiConfig.h"
#include "OcmMultiMaterialize.h"
#include "SoloBamParsing.h"
#include "streamFuns.h"
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

static string lowerStringLocal(string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool unsetTokenLocal(const string& value) {
    return value.empty() || value == "-";
}

static string resolveOcmConfigPath(const Parameters& P) {
    if (!unsetTokenLocal(P.pfMulti.ocmMultiConfig)) {
        return P.pfMulti.ocmMultiConfig;
    }
    if (!unsetTokenLocal(P.pfMulti.pfMultiConfig)) {
        return P.pfMulti.pfMultiConfig;
    }
    return string();
}

static string resolveOcmProjectRoot(const string& outPrefix) {
    string prefix = outPrefix;
    if (!prefix.empty() && prefix.back() != '/') {
        prefix.push_back('/');
    }
    const string runSuffix = "/run/";
    size_t runPos = prefix.rfind(runSuffix);
    if (runPos != string::npos) {
        const string samplesToken = "/samples/";
        size_t samplesPos = prefix.rfind(samplesToken, runPos);
        if (samplesPos != string::npos) {
            return prefix.substr(0, samplesPos + 1);
        }
        if (runPos == 0) {
            return "./";
        }
        return prefix.substr(0, runPos + 1);
    }
    size_t lastSlash = prefix.find_last_of('/');
    if (lastSlash == string::npos || lastSlash == 0) {
        return "./";
    }
    return prefix.substr(0, lastSlash + 1);
}

struct OcmBamSampleHandle {
    string sampleId;
    BGZF* bam = nullptr;
    uint64_t records = 0;
};

struct OcmBamRouter {
    bool enabled = false;
    vector<OcmBamSampleHandle> samples;
    map<string, vector<size_t>> tagToSampleIndices;
    BGZF* unassigned = nullptr;
    uint64_t unassignedRecords = 0;
    uint64_t assignedRecords = 0;
};

static void closeOcmBamRouter(OcmBamRouter& router) {
    for (auto& sample : router.samples) {
        if (sample.bam != nullptr) {
            bgzf_close(sample.bam);
            sample.bam = nullptr;
        }
    }
    if (router.unassigned != nullptr) {
        bgzf_close(router.unassigned);
        router.unassigned = nullptr;
    }
}

static BGZF* openOcmBamFile(const string& path, Parameters& P) {
    BGZF* bam = bgzf_open(path.c_str(), ("w" + to_string((long long)P.outBAMcompression)).c_str());
    if (bam == nullptr) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open OCM BAM output: " << path << "\n";
        errOut << "SOLUTION: check that the disk is not full and the output directory is writable";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    return bam;
}

static void initializeOcmBamRouter(Parameters& P, Genome& genome, OcmBamRouter& router) {
    const string mode = lowerStringLocal(P.pfMulti.ocmMultiBamSplit);
    if (mode == "no" || mode.empty()) {
        return;
    }
    if (!P.outBAMunsorted) {
        if (mode == "yes") {
            P.inOut->logMain << "WARNING: --ocmMultiBamSplit yes requires unsorted BAM output; skipping OCM BAM split.\n";
        }
        return;
    }
    const string configPath = resolveOcmConfigPath(P);
    if (configPath.empty()) {
        if (mode == "yes") {
            P.inOut->logMain << "WARNING: --ocmMultiBamSplit yes requires --ocmMultiConfig or --pfMultiConfig; skipping OCM BAM split.\n";
        }
        return;
    }

    PfMultiConfig::Config config;
    try {
        config = OcmMultiConfig::parseAndValidate(configPath, P.inOut->logMain);
    } catch (const std::exception& ex) {
        if (mode == "auto") {
            P.inOut->logMain << "WARNING: --ocmMultiBamSplit auto could not parse OCM config; skipping: "
                             << ex.what() << "\n";
            return;
        }
        exitWithError("EXITING because of fatal OCM BAM split config error: " + string(ex.what()) + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    if (config.samples.empty()) {
        if (mode == "yes") {
            P.inOut->logMain << "WARNING: OCM BAM split requested but config has no [samples]; skipping.\n";
        }
        return;
    }

    const string projectRoot = resolveOcmProjectRoot(P.outFileNamePrefix);
    const string multiCountDir = projectRoot + "outs/multi/count/";
    createDirectory(multiCountDir, P.runDirPerm, "OCM multi BAM output", P);
    router.unassigned = openOcmBamFile(multiCountDir + "unassigned_alignments.bam", P);
    outBAMwriteHeader(router.unassigned, P.samHeader, genome.chrNameAll, genome.chrLengthAll);

    for (const auto& sample : config.samples) {
        const string sampleCountDir = projectRoot + "outs/per_sample_outs/" + sample.sample_id + "/count/";
        createDirectory(sampleCountDir, P.runDirPerm, "OCM per-sample BAM output", P);
        OcmBamSampleHandle handle;
        handle.sampleId = sample.sample_id;
        handle.bam = openOcmBamFile(sampleCountDir + "sample_alignments.bam", P);
        outBAMwriteHeader(handle.bam, P.samHeader, genome.chrNameAll, genome.chrLengthAll);
        const size_t index = router.samples.size();
        router.samples.push_back(handle);
        for (const auto& ocmId : sample.resolvedOcmIds()) {
            router.tagToSampleIndices[ocmId].push_back(index);
        }
    }

    router.enabled = true;
    P.inOut->logMain << "OCM BAM split enabled during tagged unsorted BAM replay: "
                     << router.samples.size() << " samples, projectRoot=" << projectRoot << "\n";
}

static void routeOcmBamRecord(OcmBamRouter& router, const char* bamData, uint32_t bamSize) {
    if (!router.enabled) {
        return;
    }
    const string cb = SoloBamParsing::parseCB(bamData, bamSize);
    const string ocmId = cb.empty() ? string() : OcmMultiMaterialize::classifyBarcodeTag(cb);
    auto it = router.tagToSampleIndices.find(ocmId);
    if (it == router.tagToSampleIndices.end() || it->second.empty()) {
        if (router.unassigned != nullptr) {
            bgzf_write(router.unassigned, bamData, bamSize);
            ++router.unassignedRecords;
        }
        return;
    }
    for (size_t sampleIndex : it->second) {
        if (sampleIndex >= router.samples.size() || router.samples[sampleIndex].bam == nullptr) {
            continue;
        }
        bgzf_write(router.samples[sampleIndex].bam, bamData, bamSize);
        ++router.samples[sampleIndex].records;
    }
    ++router.assignedRecords;
}

static void logOcmBamRouterSummary(const OcmBamRouter& router, Parameters& P) {
    if (!router.enabled) {
        return;
    }
    P.inOut->logMain << "OCM BAM split completed: assigned_records=" << router.assignedRecords
                     << " unassigned_records=" << router.unassignedRecords << "\n";
    for (const auto& sample : router.samples) {
        P.inOut->logMain << "  OCM BAM sample " << sample.sampleId
                         << " records=" << sample.records << "\n";
    }
}

} // namespace

// Helper function to check if a chromosome is Y chromosome using genome.yTids
static bool isYChromosome(const char* bamData, const Genome& genome) {
    const uint32_t* bam32 = reinterpret_cast<const uint32_t*>(bamData);
    int32_t refID = static_cast<int32_t>(bam32[1]);
    if (refID < 0) {
        return false;
    }
    // Use genome.yTids which contains the proper Y chromosome tid set
    return genome.yTids.count(refID) > 0;
}

// Finalize samtools sorting and write output
static void bamSortSamtoolsFinalize(Parameters& P, Genome& genome, Solo& solo) {
    if (g_samtoolsSorter == nullptr) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: SamtoolsSorter not initialized but --outBAMsortMethod samtools was specified\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... started samtools BAM sorting\n" << flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... started samtools BAM sorting\n" << flush;
    
    // Finalize sorting
    g_samtoolsSorter->finalize();
    
    // Open output handles based on settings
    BGZF* bgzfPrimary = nullptr;
    BGZF* bgzfY = nullptr;
    BGZF* bgzfNoY = nullptr;
    
    if (!P.emitNoYBAMyes || P.keepBAMyes) {
        bgzfPrimary = bgzf_open(P.outBAMfileCoordName.c_str(), 
                                ("w" + to_string((long long)P.outBAMcompression)).c_str());
        if (bgzfPrimary == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not open output bam file: " << P.outBAMfileCoordName << "\n";
            errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        outBAMwriteHeader(bgzfPrimary, P.samHeaderSortedCoord, 
                          genome.chrNameAll, genome.chrLengthAll);
    }
    
    if (P.emitNoYBAMyes) {
        bgzfY = bgzf_open(P.outBAMfileYName.c_str(), 
                          ("w" + to_string((long long)P.outBAMcompression)).c_str());
        if (bgzfY == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not open output Y bam file: " << P.outBAMfileYName << "\n";
            errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        outBAMwriteHeader(bgzfY, P.samHeaderSortedCoord, 
                          genome.chrNameAll, genome.chrLengthAll);
        
        bgzfNoY = bgzf_open(P.outBAMfileNoYName.c_str(), 
                            ("w" + to_string((long long)P.outBAMcompression)).c_str());
        if (bgzfNoY == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not open output noY bam file: " << P.outBAMfileNoYName << "\n";
            errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        outBAMwriteHeader(bgzfNoY, P.samHeaderSortedCoord, 
                          genome.chrNameAll, genome.chrLengthAll);
    }
    
    // Add temp buffer for tag injection
    char bam1[BAM_ATTR_MaxSize];
    
    // Stream sorted records via k-way merge
    const char* bamData;
    uint32_t bamSize;
    uint32_t readId;
    bool hasY;
    uint64_t recordCount = 0;
    
    while (g_samtoolsSorter->nextRecord(&bamData, &bamSize, &readId, &hasY)) {
        char* bam0 = const_cast<char*>(bamData);
        uint32 size0 = bamSize;
        
        // Inject CB/UB tags if requested
        if (solo.pSolo.samAttrYes) {
            // Pass readId directly (no shifting needed)
            solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]
                ->addBAMtags(bam0, size0, bam1, readId);
        }
        
        // Determine if this is a Y chromosome alignment
        bool isYChrom = hasY || isYChromosome(bam0, genome);
        
        // Use bam0/size0 for bgzf_write (may point to bam1 after tag injection)
        if (bgzfPrimary) {
            bgzf_write(bgzfPrimary, bam0, size0);
        }
        if (P.emitNoYBAMyes) {
            bgzf_write(isYChrom ? bgzfY : bgzfNoY, bam0, size0);
        }
        recordCount++;
    }
    
    // Close handles
    if (bgzfPrimary) bgzf_close(bgzfPrimary);
    if (bgzfY) bgzf_close(bgzfY);
    if (bgzfNoY) bgzf_close(bgzfNoY);
    
    P.inOut->logMain << "samtools sorting completed: " << recordCount << " records sorted\n";
    
    // Cleanup
    delete g_samtoolsSorter;
    g_samtoolsSorter = nullptr;
}

void bamSortByCoordinate (Parameters &P, ReadAlignChunk **RAchunk, Genome &genome, Solo &solo) {
    if (P.outBAMcoord) {//sort BAM if needed
        // Branch to samtools backend if enabled
        if (P.outBAMsortMethod == "samtools") {
            bamSortSamtoolsFinalize(P, genome, solo);
            return;
        }
        
        // Continue with legacy STAR bin sorter
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started sorting BAM\n" <<flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started sorting BAM\n" <<flush;
        uint32 nBins=P.outBAMcoordNbins;

        //check max size needed for sorting
        uint maxMem=0;
        for (uint32 ibin=0; ibin<nBins-1; ibin++) {//check all bins
            uint binS=0;
            for (int it=0; it<P.runThreadN; it++) {//collect sizes from threads
                binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin]+24*RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
            };
            if (binS>maxMem) maxMem=binS;
        };

        uint64 unmappedReadsN = 0;
        for (int it=0; it<P.runThreadN; it++)
            unmappedReadsN += RAchunk[it]->chunkOutBAMcoord->binTotalN[nBins-1];

        P.inOut->logMain << "Max memory needed for sorting = "<<maxMem<<endl;
        if (maxMem>P.limitBAMsortRAM) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: not enough memory for BAM sorting: \n";
            errOut <<"SOLUTION: re-run STAR with at least --limitBAMsortRAM " <<maxMem+1000000000;
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        } else if(maxMem==0 && unmappedReadsN==0) {//both mapped and unmapped reads are absent
            P.inOut->logMain << "WARNING: nothing to sort - no output alignments" <<endl;
            if (P.emitNoYBAMyes) {
                // Create empty Y/noY BAM files with headers
                BGZF *bgzfOut_Y = bgzf_open(P.outBAMfileYName.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
                if (bgzfOut_Y==NULL) {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal ERROR: could not open output Y bam file: " << P.outBAMfileYName << "\n";
                    errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };
                outBAMwriteHeader(bgzfOut_Y,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
                bgzf_close(bgzfOut_Y);
                
                BGZF *bgzfOut_noY = bgzf_open(P.outBAMfileNoYName.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
                if (bgzfOut_noY==NULL) {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal ERROR: could not open output noY bam file: " << P.outBAMfileNoYName << "\n";
                    errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };
                outBAMwriteHeader(bgzfOut_noY,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
                bgzf_close(bgzfOut_noY);
                
                // Create primary BAM if keepBAM is enabled
                if (P.keepBAMyes) {
                    BGZF *bgzfOut;
                    bgzfOut=bgzf_open(P.outBAMfileCoordName.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
                    if (bgzfOut==NULL) {
                        ostringstream errOut;
                        errOut <<"EXITING because of fatal ERROR: could not open output bam file: " << P.outBAMfileCoordName << "\n";
                        errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                    };
                    outBAMwriteHeader(bgzfOut,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
                    bgzf_close(bgzfOut);
                }
            } else {
                // Normal mode: create primary BAM
                BGZF *bgzfOut;
                bgzfOut=bgzf_open(P.outBAMfileCoordName.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
                if (bgzfOut==NULL) {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal ERROR: could not open output bam file: " << P.outBAMfileCoordName << "\n";
                    errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };
                outBAMwriteHeader(bgzfOut,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
                bgzf_close(bgzfOut);
            }
        } else {//sort
            uint totalMem=0;
            #pragma omp parallel num_threads(P.outBAMsortingThreadNactual)
            #pragma omp for schedule (dynamic,1)
            for (uint32 ibin1=0; ibin1<nBins; ibin1++) {
                uint32 ibin=nBins-1-ibin1;//reverse order to start with the last bin - unmapped reads

                uint binN=0, binS=0;
                for (int it=0; it<P.runThreadN; it++) {//collect sizes from threads
                    binN += RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
                    binS += RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin];
                };

                if (binS==0) continue; //empty bin

                if (ibin == nBins-1) {//last bin for unmapped reads
                    BAMbinSortUnmapped(ibin,P.runThreadN,P.outBAMsortTmpDir, P, genome, solo);
                } else {
                    uint newMem=binS+binN*24;
                    bool boolWait=true;
                    while (boolWait) {
                        #pragma omp critical
                        if (totalMem+newMem < P.limitBAMsortRAM) {
                            boolWait=false;
                            totalMem+=newMem;
                        };
                        sleep(0.1);
                    };
                    BAMbinSortByCoordinate(ibin,binN,binS,P.runThreadN,P.outBAMsortTmpDir, P, genome, solo);
                    #pragma omp critical
                    totalMem-=newMem;//"release" RAM
                };
            };

            //concatenate all BAM files, using bam_cat
            if (P.emitNoYBAMyes) {
                // Concatenate Y and noY bin files separately
                vector <string> bamBinNamesV_Y, bamBinNamesV_noY;
                for (uint32 ibin=0; ibin<nBins; ibin++) {
                    string yFile = P.outBAMsortTmpDir+"/b"+std::to_string((uint) ibin)+"_Y";
                    string noYFile = P.outBAMsortTmpDir+"/b"+std::to_string((uint) ibin)+"_noY";
                    struct stat buffer;
                    if (stat(yFile.c_str(), &buffer) == 0) {
                        bamBinNamesV_Y.push_back(yFile);
                    }
                    if (stat(noYFile.c_str(), &buffer) == 0) {
                        bamBinNamesV_noY.push_back(noYFile);
                    }
                };
                
                if (!bamBinNamesV_Y.empty()) {
                    char **bamBinNames_Y = new char* [bamBinNamesV_Y.size()];
                    for (uint32 ibin=0; ibin<bamBinNamesV_Y.size(); ibin++) {
                        bamBinNames_Y[ibin] = (char*) bamBinNamesV_Y.at(ibin).c_str();
                    };
                    bam_cat(bamBinNamesV_Y.size(), bamBinNames_Y, 0, P.outBAMfileYName.c_str());
                    delete [] bamBinNames_Y;
                }
                
                if (!bamBinNamesV_noY.empty()) {
                    char **bamBinNames_noY = new char* [bamBinNamesV_noY.size()];
                    for (uint32 ibin=0; ibin<bamBinNamesV_noY.size(); ibin++) {
                        bamBinNames_noY[ibin] = (char*) bamBinNamesV_noY.at(ibin).c_str();
                    };
                    bam_cat(bamBinNamesV_noY.size(), bamBinNames_noY, 0, P.outBAMfileNoYName.c_str());
                    delete [] bamBinNames_noY;
                }
                
                // Concatenate primary BAM if keepBAM is enabled
                if (P.keepBAMyes) {
                    char **bamBinNames = new char* [nBins];
                    vector <string> bamBinNamesV;
                    for (uint32 ibin=0; ibin<nBins; ibin++) {
                        bamBinNamesV.push_back(P.outBAMsortTmpDir+"/b"+std::to_string((uint) ibin));
                        struct stat buffer;
                        if (stat (bamBinNamesV.back().c_str(), &buffer) != 0) {//check if file exists
                            bamBinNamesV.pop_back();
                        };
                    };
                    for (uint32 ibin=0; ibin<bamBinNamesV.size(); ibin++) {
                        bamBinNames[ibin] = (char*) bamBinNamesV.at(ibin).c_str();
                    };
                    bam_cat(bamBinNamesV.size(), bamBinNames, 0, P.outBAMfileCoordName.c_str());
                    delete [] bamBinNames;
                }
            } else {
                // Normal mode: concatenate primary BAM
                char **bamBinNames = new char* [nBins];
                vector <string> bamBinNamesV;
                for (uint32 ibin=0; ibin<nBins; ibin++) {
                    bamBinNamesV.push_back(P.outBAMsortTmpDir+"/b"+std::to_string((uint) ibin));
                    struct stat buffer;
                    if (stat (bamBinNamesV.back().c_str(), &buffer) != 0) {//check if file exists
                        bamBinNamesV.pop_back();
                    };
                };
                for (uint32 ibin=0; ibin<bamBinNamesV.size(); ibin++) {
                    bamBinNames[ibin] = (char*) bamBinNamesV.at(ibin).c_str();
                };
                bam_cat(bamBinNamesV.size(), bamBinNames, 0, P.outBAMfileCoordName.c_str());
                delete [] bamBinNames;
            }
        };
    };    
};

// Finalize unsorted BAM with CB/UB tag injection
// Uses g_unsortedTagBuffer which buffers records in noSort mode during mapping
void bamUnsortedWithTags(Parameters &P, Genome &genome, Solo &solo) {
    if (g_unsortedTagBuffer == nullptr) {
        return;  // Not in buffered mode
    }
    
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... writing unsorted BAM with CB/UB tags\n" << flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... writing unsorted BAM with CB/UB tags\n" << flush;
    
    // Finalize the buffer (no sorting due to noSort mode)
    g_unsortedTagBuffer->finalize();
    
    // Open output BAM file
    BGZF* bgzfOut = bgzf_open(P.outBAMfileUnsortedName.c_str(),
                              ("w" + to_string((long long)P.outBAMcompression)).c_str());
    if (bgzfOut == nullptr) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open output unsorted bam file: " << P.outBAMfileUnsortedName << "\n";
        errOut << "SOLUTION: check that the disk is not full";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    
    // Write BAM header
    outBAMwriteHeader(bgzfOut, P.samHeader, genome.chrNameAll, genome.chrLengthAll);
    
    // Optional Y/noY split files
    BGZF* bgzfY = nullptr;
    BGZF* bgzfNoY = nullptr;
    if (P.emitNoYBAMyes) {
        bgzfY = bgzf_open(P.outBAMfileYName.c_str(),
                          ("w" + to_string((long long)P.outBAMcompression)).c_str());
        if (bgzfY == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not open Y bam file\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        outBAMwriteHeader(bgzfY, P.samHeader, genome.chrNameAll, genome.chrLengthAll);
        
        bgzfNoY = bgzf_open(P.outBAMfileNoYName.c_str(),
                            ("w" + to_string((long long)P.outBAMcompression)).c_str());
        if (bgzfNoY == nullptr) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not open noY bam file\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        outBAMwriteHeader(bgzfNoY, P.samHeader, genome.chrNameAll, genome.chrLengthAll);
    }

    OcmBamRouter ocmRouter;
    initializeOcmBamRouter(P, genome, ocmRouter);
    
    // Temp buffer for tag injection
    char bam1[BAM_ATTR_MaxSize];
    
    // Stream records with tag injection
    const char* bamData;
    uint32_t bamSize;
    uint32_t readId;
    bool hasY;
    uint64_t recordCount = 0;
    
    while (g_unsortedTagBuffer->nextRecord(&bamData, &bamSize, &readId, &hasY)) {
        char* bam0 = const_cast<char*>(bamData);
        uint32 size0 = bamSize;
        
        // Inject CB/UB tags using the same path as sorted BAM
        if (solo.pSolo.samAttrYes) {
            solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]
                ->addBAMtags(bam0, size0, bam1, readId);
        }
        
        // Write to output (bam0/size0 may point to bam1 after tag injection)
        bgzf_write(bgzfOut, bam0, size0);
        routeOcmBamRecord(ocmRouter, bam0, size0);
        
        // Y/noY split
        if (P.emitNoYBAMyes) {
            // Check if this alignment is Y-chromosome
            bool isYChrom = hasY || isYChromosome(bam0, genome);
            bgzf_write(isYChrom ? bgzfY : bgzfNoY, bam0, size0);
        }
        
        recordCount++;
    }
    
    // Close output files
    bgzf_close(bgzfOut);
    if (bgzfY) bgzf_close(bgzfY);
    if (bgzfNoY) bgzf_close(bgzfNoY);
    closeOcmBamRouter(ocmRouter);
    
    P.inOut->logMain << "Unsorted BAM with CB/UB tags completed: " << recordCount << " records\n";
    logOcmBamRouterSummary(ocmRouter, P);
    
    // Cleanup
    delete g_unsortedTagBuffer;
    g_unsortedTagBuffer = nullptr;
};
