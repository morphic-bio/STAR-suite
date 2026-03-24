#include "SoloReadFeature.h"
#include "SoloFeature.h"
#include "Transcriptome.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ReadAnnotations.h"
#include "SoloReadBarcode.h"
#include "hash_shims_cpp_compat.h"
#include "ReadAlign.h"
#include "SoloBinarySpool.h"
#include "SoloReadFeature_record_shared.h"
#include "ErrorWarning.h"
#include <unordered_set>
#include <sstream>
#include <string>
#include <cstdlib>

#ifdef DEBUG_CB_UB_PARITY
// Optional tracing of specific readIds via STAR_DEBUG_TRACE_READS=1,2,3
static std::unordered_set<uint32_t> buildTraceReadSetWriter() {
    std::unordered_set<uint32_t> out;
    const char* env = std::getenv("STAR_DEBUG_TRACE_READS");
    if (!env || env[0]=='\0') return out;
    std::string s(env);
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try {
            out.insert(static_cast<uint32_t>(std::stoul(tok)));
        } catch (...) {
            // ignore malformed entries
        }
    }
    return out;
}
static const std::unordered_set<uint32_t> g_traceReadsWriter = buildTraceReadSetWriter();
#endif

void record_base(SoloReadFeature *soloReadFeat, SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (soloReadFeat->pSolo.type==0)
        return;

    if (soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]) {//readFlag

        if (nTr==1) {
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.genomeU);
        } else if (nTr>1) {
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.genomeM);
        };

        for (uint64 itr=0; itr<nTr; itr++) {
            if (soloReadFeat->P.pGe.chrSet.mito.count(alignOut[itr]->Chr) == 1) {
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.mito);
            };
        };

        switch (readAnnot.annotFeatures[soloReadFeat->featureType].ovType) {
            case ReadAnnotFeature::overlapTypes::exonic : 
            case ReadAnnotFeature::overlapTypes::exonic50p :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.exonic);
                break;
            case ReadAnnotFeature::overlapTypes::intronic : 
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.intronic);
                break;
            case ReadAnnotFeature::overlapTypes::exonicAS : 
            case ReadAnnotFeature::overlapTypes::exonic50pAS :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.exonicAS);            
                break;
            case ReadAnnotFeature::overlapTypes::intronicAS :
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.intronicAS);            
        };

        if (soloBar.cbMatch<0 && soloReadFeat->pSolo.cbWLyes) {//no CB match in the WL
            if (readAnnot.annotFeatures[soloReadFeat->featureType].fSet.size()==1) {
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
            } else if (readAnnot.annotFeatures[soloReadFeat->featureType].fSet.size()>1){
                soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureM);
            };
            soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.cbMatch);//this will counts reads with no CB match
            soloReadFeat->readFlag.countsAddNoCB();
        };
    };

    if (soloBar.cbMatch<0) {
        return;
    }

       
    ReadSoloFeatures reFe;
    reFe.alignOut=alignOut;
    reFe.indAnnotTr = 0;    

    uint32 nFeat=0; //number of features in this read (could be >1 for SJs)
    if (nTr==0) {//unmapped
        soloReadFeat->stats.V[soloReadFeat->stats.noUnmapped]++;
        
    } else {
        switch (soloReadFeat->featureType) {
            case SoloFeatureTypes::Gene :
            case SoloFeatureTypes::GeneFull :
            case SoloFeatureTypes::GeneFull_Ex50pAS :
            case SoloFeatureTypes::GeneFull_ExonOverIntron :
                {
                    auto *readGe = &readAnnot.annotFeatures[soloReadFeat->featureType].fSet;

                    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq) {
                        for (int32 itr=nTr-1; itr>=0; itr--) {
                            if (readAnnot.annotFeatures[soloReadFeat->featureType].fAlign[itr].size() > 0) {
                                reFe.indAnnotTr = itr;
                                break;//indAnnotTr is the first genic align, and is used for Smart-seq deduplication, below in outputReadCB()
                                      //TODO this is ambiguous for reads that map to the same gene multiple times. *Last* align is chosen to match the previous version.
                            };
                        };
                    };
                        
                    if (readGe->size()==0) {//check genes
                        soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;//no gene
                    } else if (readGe->size()>1) {
                        soloReadFeat->stats.V[soloReadFeat->stats.MultiFeature]++;//multi-gene reads
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureM);
                        if (nTr>1)
                            soloReadFeat->stats.V[soloReadFeat->stats.subMultiFeatureMultiGenomic]++;//multigene caused by multimapper
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsWriter.empty() && iRead!=(uint64)-1 && g_traceReadsWriter.count((uint32_t)iRead)) {
                            std::string genes;
                            genes.reserve(readGe->size()*8);
                            bool first=true;
                            for (auto gVal : *readGe) {
                                if (!first) genes.push_back(',');
                                first=false;
                                genes += std::to_string(gVal);
                            }
                            fprintf(stderr, "[TRACE annotate] read=%llu genes=%s nTr=%u cbMatch=%d cbMatchString=%s\n",
                                    (unsigned long long)iRead, genes.c_str(), (unsigned)nTr,
                                    soloBar.cbMatch, soloBar.cbMatchString.c_str());
                        }
#endif
                            
                        if (soloReadFeat->pSolo.multiMap.yes.multi) {//output multimappers
                            reFe.geneMult.resize(readGe->size());
                            uint32 ii=0;
                            for (auto &g: *readGe) {//set high bit to mark multimappers
                                reFe.geneMult[ii] = g | geneMultMark;
                                ++ii;
                            };
                                
                            nFeat = outputReadCB_base(soloReadFeat->streamReads, iRead, soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                        };
                    } else {//unique-gene reads
                        reFe.gene = *readGe->begin();
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
                        nFeat = outputReadCB_base(soloReadFeat->streamReads, (soloReadFeat->readIndexYes ? iRead : (uint64)-1), soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                    };
                };
                break;
        
            case SoloFeatureTypes::SJ : 
                if (nTr>1) {//reject all multimapping reads
                    soloReadFeat->stats.V[soloReadFeat->stats.subMultiFeatureMultiGenomic]++;
                    soloReadFeat->stats.V[soloReadFeat->stats.MultiFeature]++;
                } else {//one gene or no gene
                    alignOut[0]->extractSpliceJunctions(reFe.sj, reFe.sjAnnot);
                    if ( reFe.sj.empty() ) {
                        soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                    } else {//good junction
                        soloReadFeat->readFlag.setBit(soloReadFeat->readFlag.featureU);
                        nFeat = outputReadCB_base(soloReadFeat->streamReads, (soloReadFeat->readIndexYes ? iRead : (uint64)-1), soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                    };
                };                  
                break;
        
            case SoloFeatureTypes::Transcript3p :
                if (readAnnot.transcriptConcordant.size()==0 || soloBar.cbMatch>1) {//do not record ambiguous CB  
                    soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                } else {
                    nFeat = outputReadCB_base(soloReadFeat->streamReads, iRead, soloReadFeat->featureType, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
                };
                if (readAnnot.transcriptConcordant.size()==1 && readAnnot.transcriptConcordant[0][1] < soloReadFeat->transcriptDistCount.size()) {
                    //read maps to one transcript - use for distTTS distribution function
                    soloReadFeat->transcriptDistCount[readAnnot.transcriptConcordant[0][1]]++;
                };
                
                break;

            case SoloFeatureTypes::Velocyto :
                //different record: iRead, nTr, tr1, type1, tr2, type2 ...
                if (readAnnot.trVelocytoType.size()>0) {//otherwise, no gene
                    
                    sort(readAnnot.trVelocytoType.begin(), readAnnot.trVelocytoType.end(),
                         [](const trTypeStruct &t1, const trTypeStruct &t2) {return t1.tr < t2.tr;});

                    if (soloReadFeat->streamReads) {
                        *soloReadFeat->streamReads << iRead <<' '<< readAnnot.trVelocytoType.size();
                        for (auto &tt: readAnnot.trVelocytoType)
                             *soloReadFeat->streamReads <<' '<< tt.tr <<' '<< (uint32) tt.type;
                        *soloReadFeat->streamReads <<'\n';
                    }
                    nFeat=1;
                } else {
                    soloReadFeat->stats.V[soloReadFeat->stats.noNoFeature]++;
                };
                break; //no need to go with downstream processing                
                
        };//switch (featureType)
    };//if (nTr==0)
    
    if ( nFeat==0 && (soloReadFeat->readInfoYes | soloReadFeat->pSolo.readStatsYes[soloReadFeat->featureType]) ) {//no feature, but readInfo requested
        outputReadCB_base(soloReadFeat->streamReads, iRead, (uint32)-1, soloBar, reFe, readAnnot, soloReadFeat->readFlag, soloReadFeat);
    };
    
    if (nFeat==0)
        return; //no need to record the number of reads per CB
    
    if (soloReadFeat->pSolo.cbWLyes) {//WL
        for (auto &cbi : soloBar.cbMatchInd)
            soloReadFeat->cbReadCount[cbi] += nFeat;
    } else {//no WL
        soloReadFeat->cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    };
    
    return;
};

uint32 outputReadCB_base(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
                         const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                         SoloReadFeature *soloReadFeat)
{   
    /*format of the temp output file
     * UMI [iRead] type feature* cbMatchString
     *             0=exact match, 1=one non-exact match, 2=multipe non-exact matches
     *                   gene or sj[0] sj[1]
     *                         CB or nCB {CB Qual, ...}
     */
    
    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq && featureType!=-1) {//need to calculate "UMI" from align start/end
        soloBar.umiB=reFe.alignOut[reFe.indAnnotTr]->chrStartLengthExtended();
    };

    const bool binarySpool = soloReadFeat != nullptr && soloReadFeat->binarySpool;
    std::vector<std::pair<uint32_t, uint8_t>> ambiguousCbCandidates;
    if (binarySpool && soloBar.cbMatch > 1) {
        ambiguousCbCandidates.reserve(static_cast<size_t>(soloBar.cbMatch));
        std::istringstream cbMatchStream(soloBar.cbMatchString);
        uint32_t cbIdx = 0;
        char qualChar = 0;
        while (cbMatchStream >> cbIdx >> qualChar) {
            ambiguousCbCandidates.emplace_back(cbIdx, static_cast<uint8_t>(qualChar));
        }
        if (static_cast<int32>(ambiguousCbCandidates.size()) != soloBar.cbMatch) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not encode ambiguous CB payload into experimental binary Solo spool\n"
                   << "cbMatch=" << soloBar.cbMatch << " parsedCandidates=" << ambiguousCbCandidates.size() << "\n";
            exitWithError(errOut.str(), std::cerr, soloReadFeat->P.inOut->logMain, EXIT_CODE_PARAMETER, soloReadFeat->P);
        }
    }

    auto writeCbPayload = [&]() {
        const bool binarySpoolInMemoryNow = binarySpool && soloReadFeat->binarySpoolInMemory;
        fstream *activeStream = binarySpool ? soloReadFeat->streamReads : streamOut;
        if (activeStream == nullptr && !binarySpoolInMemoryNow) {
            return;
        }
        if (binarySpool) {
            if (soloBar.cbMatch <= 1) {
                uint64_t cbValue = 0;
                if (!soloBar.cbMatchInd.empty()) {
                    cbValue = soloBar.cbMatchInd[0];
                } else if (!soloBar.cbMatchString.empty()) {
                    cbValue = static_cast<uint64_t>(std::stoull(soloBar.cbMatchString));
                } else {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal ERROR: missing resolved CB payload for experimental binary Solo spool\n";
                    exitWithError(errOut.str(), std::cerr, soloReadFeat->P.inOut->logMain, EXIT_CODE_PARAMETER, soloReadFeat->P);
                }
                if (binarySpoolInMemoryNow) {
                    SoloBinarySpool::writeResolvedCb(soloReadFeat->binarySpoolBuffer, cbValue);
                } else {
                    SoloBinarySpool::writeResolvedCb(*activeStream, cbValue);
                }
            } else {
                for (const auto &candidate : ambiguousCbCandidates) {
                    if (binarySpoolInMemoryNow) {
                        SoloBinarySpool::writeAmbiguousCandidate(soloReadFeat->binarySpoolBuffer, candidate.first, candidate.second);
                    } else {
                        SoloBinarySpool::writeAmbiguousCandidate(*activeStream, candidate.first, candidate.second);
                    }
                }
            }
        } else {
            *activeStream << soloBar.cbMatchString;
        }
    };
    auto binaryHeaderBytes = [&](bool writeReadIndex) -> size_t {
        return sizeof(uint64_t) + sizeof(uint32_t) + sizeof(int32_t) +
               (writeReadIndex ? (sizeof(uint64_t) + sizeof(uint32_t)) : 0);
    };
    auto binaryPayloadBytes = [&]() -> size_t {
        if (soloBar.cbMatch <= 1) {
            return sizeof(uint64_t);
        }
        return static_cast<size_t>(soloBar.cbMatch) * (sizeof(uint32_t) + sizeof(uint8_t));
    };
    auto maybeSpillBeforeBinaryWrite = [&](bool writeReadIndex) {
        if (binarySpool && soloReadFeat->binarySpoolInMemory) {
            soloReadFeat->maybeSpillBinarySpool(binaryHeaderBytes(writeReadIndex) + binaryPayloadBytes());
        }
    };
    
    uint64 nout=1;
    
    switch (featureType) {
        case -1 : {
            // no feature, output for readInfo
            if (streamOut || (binarySpool && soloReadFeat->binarySpoolInMemory)) {
                if (binarySpool) {
                    maybeSpillBeforeBinaryWrite(true);
                    fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                    if (soloReadFeat->binarySpoolInMemory) {
                        SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, true, soloBar.umiB, iRead, readFlag.flag, static_cast<uint32_t>(-1), soloBar.cbMatch);
                    } else {
                        SoloBinarySpool::writeRecordHeader(*activeStream, true, soloBar.umiB, iRead, readFlag.flag, static_cast<uint32_t>(-1), soloBar.cbMatch);
                    }
                    writeCbPayload();
                } else {
                    *streamOut << soloBar.umiB <<' '<< iRead <<' '<< readFlag.flag <<' '<< -1 <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                }
            }
            break;
        }
            
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron : {
            if (streamOut || (binarySpool && soloReadFeat->binarySpoolInMemory)) {
                if (!reFe.geneMult.empty()) {
                    for (uint32_t geneIdx : reFe.geneMult) {
                        if (binarySpool) {
                            maybeSpillBeforeBinaryWrite(iRead != (uint64)-1);
                            fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                            if (soloReadFeat->binarySpoolInMemory) {
                                SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, geneIdx, soloBar.cbMatch);
                            } else {
                                SoloBinarySpool::writeRecordHeader(*activeStream, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, geneIdx, soloBar.cbMatch);
                            }
                            writeCbPayload();
                        } else {
                            *streamOut << soloBar.umiB <<' ';
                            if (iRead != (uint64)-1) {
                                *streamOut << iRead <<' '<< readFlag.flag <<' ';
                            }
                            *streamOut << geneIdx <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                        }
                    }
                    nout = reFe.geneMult.size();
                } else {
                    if (binarySpool) {
                        maybeSpillBeforeBinaryWrite(iRead != (uint64)-1);
                        fstream *activeStream = soloReadFeat->binarySpoolInMemory ? nullptr : soloReadFeat->streamReads;
                        if (soloReadFeat->binarySpoolInMemory) {
                            SoloBinarySpool::writeRecordHeader(soloReadFeat->binarySpoolBuffer, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, reFe.gene, soloBar.cbMatch);
                        } else {
                            SoloBinarySpool::writeRecordHeader(*activeStream, iRead != (uint64)-1, soloBar.umiB, iRead, readFlag.flag, reFe.gene, soloBar.cbMatch);
                        }
                        writeCbPayload();
                    } else {
                        *streamOut << soloBar.umiB <<' ';
                        if (iRead != (uint64)-1) {
                            *streamOut << iRead <<' '<< readFlag.flag <<' ';
                        }
                        *streamOut << reFe.gene <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
                    }
                    nout = 1;
                }
            }
            break;
        }

        case SoloFeatureTypes::SJ : {
            //sj - two numbers, multiple sjs per read
            if (streamOut) {
                for (auto &sj : reFe.sj) {
                    *streamOut << soloBar.umiB <<' ';//UMI
                    if ( iRead != (uint64)-1 )
                        *streamOut << iRead <<' '<< readFlag.flag <<' ';//iRead            
                    *streamOut << sj[0] <<' '<< sj[1] <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n' << flush;
                }
            }
            nout=reFe.sj.size();
            break;
        }

        case SoloFeatureTypes::Transcript3p : {
            //transcript,distToTTS structure
            if (streamOut) {
                *streamOut << soloBar.cbMatchString <<' ';            
                *streamOut << soloBar.umiB <<' ';
                *streamOut << readAnnot.transcriptConcordant.size();
                for (auto &tt: readAnnot.transcriptConcordant) {
                    *streamOut <<' '<< tt[0] <<' '<< tt[1];
                };
                if ( iRead != (uint64)-1 )
                    *streamOut  <<' '<< iRead;//iRead
                *streamOut  <<'\n';
            }
            nout=1;

            break;
        }
        default:
            break;
    }; //switch (featureType)
    
    return nout;
};
