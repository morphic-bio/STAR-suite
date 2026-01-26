#include "ReadAlignChunk.h"
#include "alignment_model.h"  // For Transcriptome
#include "SlamQuant.h"
#include "SlamCompat.h"
#include <pthread.h>
#include "ErrorWarning.h"
#include "ProbeListIndex.h"
#include "Solo.h"
#include <fstream>
#include <cstdlib>
#include <unordered_set>
#include <zlib.h>

namespace {
std::vector<uint8_t> buildSlamAllowedGenes(const Transcriptome& tr) {
    std::vector<uint8_t> allowed;
    if (tr.geBiotype.empty()) {
        return allowed;
    }
    bool anyBiotype = false;
    for (const auto& bt : tr.geBiotype) {
        if (!bt.empty()) {
            anyBiotype = true;
            break;
        }
    }
    if (!anyBiotype) {
        return allowed;
    }
    static const std::unordered_set<std::string> allowedBiotypes = {
        "protein_coding",
        "lincRNA",
        "antisense",
        "IG_LV_gene",
        "IG_V_gene",
        "IG_V_pseudogene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_J_pseudogene",
        "IG_C_gene",
        "IG_C_pseudogene",
        "TR_V_gene",
        "TR_V_pseudogene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_J_pseudogene",
        "TR_C_gene",
        "synthetic"
    };
    allowed.assign(tr.nGe, 0);
    size_t n = tr.geBiotype.size();
    if (n > tr.nGe) {
        n = tr.nGe;
    }
    for (size_t i = 0; i < n; ++i) {
        if (allowedBiotypes.count(tr.geBiotype[i]) > 0) {
            allowed[i] = 1;
        }
    }
    return allowed;
}
}

ReadAlignChunk::ReadAlignChunk(Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk,
                                const libem::Transcriptome* libemTr) : P(Pin), mapGen(genomeIn) {//initialize chunk

    iThread=iChunk;
    chunkTr = nullptr;
    slamQuant = nullptr;
    slamCompat = nullptr;
    RA = nullptr;
    chunkIn = nullptr;
    readInStream = nullptr;
    chunkOutBAM = nullptr;
    chunkOutBAM1 = nullptr;
    chunkOutSJ = nullptr;
    chunkOutSJ1 = nullptr;
    chunkOutBAMcoord = nullptr;
    chunkOutBAMunsorted = nullptr;
    chunkOutBAMquant = nullptr;
    chunkQuants = nullptr;
    chunkOutBAMstream = nullptr;

    if ( P.quant.yes ) {//allocate transcriptome structures
        chunkTr=new Transcriptome(*TrIn);
        chunkTr->quantsAllocate();
    } else {
        chunkTr=NULL;
    };

    if (P.quant.slam.yes && chunkTr != nullptr) {
        // For SNP mask build pre-pass, allow "alt" to mean any mismatch (GEDI-like) vs conversions only.
        bool snpObsAnyMismatch = false;
        bool hasBuildFastqs = !P.quant.slamSnpMask.buildFastqsFofn.empty() && P.quant.slamSnpMask.buildFastqsFofn != "-" &&
                              P.quant.slamSnpMask.buildFastqsFofn != "None" && P.quant.slamSnpMask.buildFastqsFofn != "none";
        bool hasBuildBam = !P.quant.slamSnpMask.buildBam.empty() && P.quant.slamSnpMask.buildBam != "-" &&
                           P.quant.slamSnpMask.buildBam != "None" && P.quant.slamSnpMask.buildBam != "none";
        if (hasBuildFastqs || hasBuildBam) {
            snpObsAnyMismatch = (P.quant.slamSnpMask.kMode == "any");
        }

        slamQuant = new SlamQuant(chunkTr->nGe, buildSlamAllowedGenes(*chunkTr), P.quant.slam.snpDetect, P.quant.slam.snpDetectFrac, snpObsAnyMismatch);
        // Enable dump buffer for external re-quant (skip auto-trim detection pass).
        bool wantDump = !P.quant.slam.dumpBinary.empty() && P.quant.slam.dumpBinary != "-" &&
                        P.quant.slam.dumpBinary != "None";
        bool wantWeights = !P.quant.slam.dumpWeights.empty() && P.quant.slam.dumpWeights != "-" &&
                           P.quant.slam.dumpWeights != "None";
        if ((wantDump || wantWeights) && !P.quant.slam.autoTrimDetectionPass) {
            slamQuant->enableDumpBuffer(P.quant.slam.dumpMaxReads);
        }
        if (P.quant.slam.debugEnabled) {
            slamQuant->initDebug(*chunkTr, P.quant.slam.debugGenes, P.quant.slam.debugReads,
                                 static_cast<size_t>(P.quant.slam.debugMaxReads),
                                 P.quant.slam.debugOutPrefix);
        }
        if (!P.quant.slam.debugSnpLoc.empty() && P.quant.slam.debugSnpLoc != "-" &&
            P.quant.slam.debugSnpLoc != "None" && P.quant.slam.debugSnpLoc != "none") {
            // Parse <chrom>:<pos1> (1-based) and convert to STAR absolute genome coordinate.
            std::string loc = P.quant.slam.debugSnpLoc;
            if (loc.rfind("chr", 0) == 0) {
                // allow chr prefix; genome names typically already include it
            }
            size_t colon = loc.find(':');
            if (colon != std::string::npos && colon + 1 < loc.size()) {
                std::string chr = loc.substr(0, colon);
                std::string posStr = loc.substr(colon + 1);
                uint64_t pos1 = 0;
                try {
                    pos1 = static_cast<uint64_t>(std::stoull(posStr));
                } catch (...) {
                    pos1 = 0;
                }
                if (pos1 > 0) {
                    // Find chromosome index
                    int chrIdx = -1;
                    for (size_t i = 0; i < mapGen.chrName.size(); ++i) {
                        if (mapGen.chrName[i] == chr) {
                            chrIdx = static_cast<int>(i);
                            break;
                        }
                    }
                    // Also try stripping/adding "chr" if needed
                    if (chrIdx < 0 && chr.rfind("chr", 0) == 0) {
                        std::string chr2 = chr.substr(3);
                        for (size_t i = 0; i < mapGen.chrName.size(); ++i) {
                            if (mapGen.chrName[i] == chr2) {
                                chrIdx = static_cast<int>(i);
                                break;
                            }
                        }
                    } else if (chrIdx < 0) {
                        std::string chr2 = "chr" + chr;
                        for (size_t i = 0; i < mapGen.chrName.size(); ++i) {
                            if (mapGen.chrName[i] == chr2) {
                                chrIdx = static_cast<int>(i);
                                break;
                            }
                        }
                    }
                    if (chrIdx >= 0 && static_cast<size_t>(chrIdx) < mapGen.chrStart.size()) {
                        uint64_t pos0 = pos1 - 1; // 1-based -> 0-based
                        uint64_t absPos = mapGen.chrStart[chrIdx] + pos0;
                        slamQuant->enableSnpSiteDebug(absPos, P.quant.slam.debugSnpWindow, loc);
                    }
                }
            }
        }
        
        // Enable variance analysis during detection pass (single-threaded)
        // With rewind approach, detection pass collects variance stats, then files are rewound
        // and main mapping pass uses computed trims from the start
        // Always enabled during detection pass when SLAM is active
        if (P.quant.slam.autoTrimDetectionPass) {
            slamQuant->enableVarianceAnalysis(
                P.quant.slam.autoTrimMaxReads, 
                P.quant.slam.autoTrimMinReads,
                P.quant.slam.autoTrimSmoothWindow,
                P.quant.slam.autoTrimSegMinLen,
                P.quant.slam.autoTrimMaxTrim);
        }
        
        // Create SlamCompat if any compat mode is enabled or auto-trim is active
        bool needsCompat = P.quant.slam.compatIntronic || P.quant.slam.compatLenientOverlap ||
            P.quant.slam.compatOverlapWeight || P.quant.slam.compatIgnoreOverlap ||
            P.quant.slam.compatTrim5p != 0 || P.quant.slam.compatTrim3p != 0 ||
            (P.quant.slam.autoTrimMode == "variance" && P.quant.slam.autoTrimComputed);
        if (needsCompat) {
            SlamCompatConfig cfg;
            cfg.intronic = P.quant.slam.compatIntronic;
            cfg.lenientOverlap = P.quant.slam.compatLenientOverlap;
            cfg.overlapWeight = P.quant.slam.compatOverlapWeight;
            cfg.ignoreOverlap = P.quant.slam.compatIgnoreOverlap;
            // Use auto-trim values if computed, otherwise use manual trims
            if (P.quant.slam.autoTrimComputed) {
                cfg.trim5p = P.quant.slam.autoTrim5p;
                cfg.trim3p = P.quant.slam.autoTrim3p;
            } else {
                cfg.trim5p = P.quant.slam.compatTrim5p;
                cfg.trim3p = P.quant.slam.compatTrim3p;
            }
            slamCompat = new SlamCompat(*chunkTr, cfg);
        }
    }

    RA = new ReadAlign(P, mapGen, chunkTr, iChunk, libemTr);//new local copy of RA for each chunk

    RA->iRead=0;
    RA->slamQuant = slamQuant;
    RA->slamSnpMask = P.quant.slam.snpMask;
    RA->slamCompat = slamCompat;

    chunkIn=new char* [P.readNends];
    readInStream=new istringstream* [P.readNends];
    
    for (uint ii=0;ii<P.readNends;ii++) {
       chunkIn[ii]=new char[P.chunkInSizeBytesArray];//reserve more space to finish loading one read
       memset(chunkIn[ii],'\n',P.chunkInSizeBytesArray);
       readInStream[ii] = new istringstream;
       readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii],P.chunkInSizeBytesArray);
       RA->readInStream[ii]=readInStream[ii];
    };


    if (P.outSAMbool) {
        chunkOutBAM=new char [P.chunkOutBAMsizeBytes];
        RA->outBAMarray=chunkOutBAM;
        chunkOutBAMstream=new ostringstream;
        chunkOutBAMstream->rdbuf()->pubsetbuf(chunkOutBAM,P.chunkOutBAMsizeBytes);
        RA->outSAMstream=chunkOutBAMstream;
        RA->outSAMstream->seekp(0,ios::beg);
        chunkOutBAMtotal=0;
    };

    if (P.outBAMunsorted) {
        // Direct mode: use regular BAM output
        chunkOutBAMunsorted = new BAMoutput (P.inOut->outBAMfileUnsorted, P);
        RA->outBAMunsorted = chunkOutBAMunsorted;
    } else {
        chunkOutBAMunsorted=NULL;
        RA->outBAMunsorted=NULL;
    };

    if (P.outBAMcoord) {
        chunkOutBAMcoord = new BAMoutput (iChunk, P.outBAMsortTmpDir, P);
        RA->outBAMcoord = chunkOutBAMcoord;
    } else {
        chunkOutBAMcoord=NULL;
        RA->outBAMcoord=NULL;
    };

    if ( P.quant.trSAM.bamYes ) {
        chunkOutBAMquant = new BAMoutput (P.inOut->outQuantBAMfile,P);
        // Transcriptome BAM goes to separate file, skip g_unsortedTagBuffer routing
        chunkOutBAMquant->setSkipGlobalBuffer(true);
        RA->outBAMquant = chunkOutBAMquant;
    } else {
        chunkOutBAMquant=NULL;
        RA->outBAMquant=NULL;
    };

    if (P.outSJ.yes) {
        chunkOutSJ  = new OutSJ (P.limitOutSJcollapsed, P, mapGen);
        RA->chunkOutSJ  = chunkOutSJ;
    } else {
        RA->chunkOutSJ  = NULL;
    };

    if (P.outFilterBySJoutStage == 1) {
        chunkOutSJ1 = new OutSJ (P.limitOutSJcollapsed, P, mapGen);
        RA->chunkOutSJ1 = chunkOutSJ1;
    } else {
        RA->chunkOutSJ1  = NULL;
    };

    
    

    if (P.pCh.segmentMin>0) {
       if (P.pCh.out.samOld) {
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.sam.thread", iChunk, RA->chunkOutChimSAM);
       };
       if (P.pCh.out.junctions) {
            chunkFstreamOpen(P.outFileTmp + "/Chimeric.out.junction.thread", iChunk, *RA->chunkOutChimJunction);
       };
    };

    if (P.outReadsUnmapped=="Fastx" ) {
        for (uint32 imate=0; imate < P.readNends; imate++) 
            chunkFstreamOpen(P.outFileTmp + "/Unmapped.out.mate"+ to_string(imate) +".thread",iChunk, RA->chunkOutUnmappedReadsStream[imate]);
    };

    if (P.outFilterType=="BySJout") {
        chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate1.thread",iChunk, RA->chunkOutFilterBySJoutFiles[0]);
        if (P.readNends==2) chunkFstreamOpen(P.outFileTmp + "/FilterBySJoutFiles.mate2.thread",iChunk, RA->chunkOutFilterBySJoutFiles[1]); //here we do not output barcode read
    };

    if (P.emitYReadNamesyes) {
        chunkFstreamOpen(P.outFileTmp + "/YReadNames.out.thread", iChunk, RA->chunkOutYReadNames);
    }
    
    if (P.emitYNoYFastqyes) {
        for (uint32 imate = 0; imate < P.readNmates; imate++) {
            if (P.emitYNoYFastqCompression == "gz") {
                // Open gzip-compressed streams
                ostringstream yName, noYName;
                yName << P.outFileTmp << "/YFastq.mate" << imate << ".thread" << iChunk << ".gz";
                noYName << P.outFileTmp << "/noYFastq.mate" << imate << ".thread" << iChunk << ".gz";
                RA->chunkOutYFastqGz[imate] = gzopen(yName.str().c_str(), "wb");
                RA->chunkOutNoYFastqGz[imate] = gzopen(noYName.str().c_str(), "wb");
                if (RA->chunkOutYFastqGz[imate] == nullptr || RA->chunkOutNoYFastqGz[imate] == nullptr) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR: could not create Y/noY FASTQ output files\n";
                    errOut << "Solution: check that you have permission to write and disk space\n";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                }
                RA->chunkOutYFastqStream[imate].setstate(ios::badbit); // Mark fstream as unused
                RA->chunkOutNoYFastqStream[imate].setstate(ios::badbit);
            } else {
                // Open uncompressed streams
                chunkFstreamOpen(P.outFileTmp + "/YFastq.mate" + to_string(imate) + ".thread", iChunk, RA->chunkOutYFastqStream[imate]);
                chunkFstreamOpen(P.outFileTmp + "/noYFastq.mate" + to_string(imate) + ".thread", iChunk, RA->chunkOutNoYFastqStream[imate]);
                RA->chunkOutYFastqGz[imate] = nullptr;
                RA->chunkOutNoYFastqGz[imate] = nullptr;
            }
        }
    }

    if (P.wasp.yes) {
        RA->waspRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
    };
    if (P.peOverlap.yes) {
        RA->peMergeRA= new ReadAlign(Pin,genomeIn,TrIn,iChunk);
        delete RA->peMergeRA->chunkOutChimJunction;
        RA->peMergeRA->chunkOutChimJunction=RA->chunkOutChimJunction;//point to the same out-stream
        RA->peMergeRA->ownsChunkOutChimJunction_ = false;
        RA->peMergeRA->chimDet->ostreamChimJunction=RA->peMergeRA->chunkOutChimJunction;
        RA->peMergeRA->outBAMunsorted=RA->outBAMunsorted;
        RA->peMergeRA->outBAMcoord=RA->outBAMcoord;
    };
};

ReadAlignChunk::~ReadAlignChunk() {
    // Clean up owned resources
    delete RA;
    RA = nullptr;

    if (readInStream != nullptr) {
        for (uint ii = 0; ii < P.readNends; ++ii) {
            delete readInStream[ii];
        }
        delete[] readInStream;
        readInStream = nullptr;
    }

    if (chunkIn != nullptr) {
        for (uint ii = 0; ii < P.readNends; ++ii) {
            delete[] chunkIn[ii];
        }
        delete[] chunkIn;
        chunkIn = nullptr;
    }

    delete chunkOutBAMstream;
    chunkOutBAMstream = nullptr;

    delete[] chunkOutBAM;
    chunkOutBAM = nullptr;

    delete[] chunkOutBAM1;
    chunkOutBAM1 = nullptr;

    delete chunkOutBAMcoord;
    chunkOutBAMcoord = nullptr;

    delete chunkOutBAMunsorted;
    chunkOutBAMunsorted = nullptr;

    delete chunkOutBAMquant;
    chunkOutBAMquant = nullptr;

    delete chunkOutSJ;
    chunkOutSJ = nullptr;

    delete chunkOutSJ1;
    chunkOutSJ1 = nullptr;

    delete chunkTr;
    chunkTr = nullptr;

    delete slamCompat;
    slamCompat = nullptr;
    // Note: slamQuant is merged into global stats before destruction in STAR.cpp,
    // but we should still clean up the per-chunk instance
    delete slamQuant;
    slamQuant = nullptr;
};

void ReadAlignChunk::reinitSlamCompat(int trim5p, int trim3p) {
    if (slamCompat != nullptr) {
        // Update existing SlamCompat with new trim values
        slamCompat->updateTrims(trim5p, trim3p);
    } else if (P.quant.slam.yes && chunkTr != nullptr) {
        // Create SlamCompat if it doesn't exist yet (trims were computed but no other compat features)
        SlamCompatConfig cfg;
        cfg.intronic = P.quant.slam.compatIntronic;
        cfg.lenientOverlap = P.quant.slam.compatLenientOverlap;
        cfg.overlapWeight = P.quant.slam.compatOverlapWeight;
        cfg.ignoreOverlap = P.quant.slam.compatIgnoreOverlap;
        cfg.trim5p = trim5p;
        cfg.trim3p = trim3p;
        slamCompat = new SlamCompat(*chunkTr, cfg);
        RA->slamCompat = slamCompat;
    }
}

///////////////
void ReadAlignChunk::chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut) {//open fstreams for chunks
    ostringstream fNameStream1;
    fNameStream1 << filePrefix << iChunk;
    string fName1=fNameStream1.str();
    P.inOut->logMain << "Opening the file: " << fName1 << " ... " <<flush;

    remove(fName1.c_str()); //remove the file
    fstreamOut.open(fName1.c_str(),ios::out); //create empty file
    fstreamOut.close();
    fstreamOut.open(fName1.c_str(), ios::in | ios::out); //re-open the file in in/out mode

    if (fstreamOut.fail()) {
        P.inOut->logMain << "failed!\n";
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create output file "<< fName1 << "\n";
        errOut << "Solution: check that you have permission to write this file\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    P.inOut->logMain << "ok" <<endl;
};

void ReadAlignChunk::chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal){
    chunkOut.flush();
    chunkOut.seekg(0,ios::beg);
    if (mutexFlag) pthread_mutex_lock(&mutexVal);
    allOut << chunkOut.rdbuf();
    allOut.clear();
    allOut.flush();
    allOut.clear();
    if (mutexFlag) pthread_mutex_unlock(&mutexVal);
    chunkOut.clear();
    chunkOut.seekp(0,ios::beg); //set put pointer at the beginning
};


void ReadAlignChunk::chunkFilesCat(ostream *allOut, string filePrefix, uint &iC) {//concatenates a file into main output
            while (true) {
                ostringstream name1("");
                name1 << filePrefix <<iC;
                ifstream fileChunkIn(name1.str().c_str());
                if (fileChunkIn.good()) {
                    *allOut << fileChunkIn.rdbuf();
                    allOut->flush();
                    allOut->clear();
                    fileChunkIn.close();
                    fileChunkIn.clear();
                    remove(name1.str().c_str());
                    iC++;
                } else {
                    fileChunkIn.close();
                    break;
                };
            };
};
