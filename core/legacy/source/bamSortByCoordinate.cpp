#include "bamSortByCoordinate.h"
#include "BAMfunctions.h"
#include "BAMbinSortByCoordinate.h"
#include "BAMbinSortUnmapped.h"
#include "ErrorWarning.h"
#include "bam_cat.h"
#include "SamtoolsSorter.h"
#include "GlobalVariables.h"
#include <sys/stat.h>

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
