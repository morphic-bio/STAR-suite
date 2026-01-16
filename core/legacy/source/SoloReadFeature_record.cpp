#include "SoloReadFeature.h"
#include "SoloFeature.h"
#include "GeneResolver.h"
#include "ProbeListIndex.h"
#include "SampleDetector.h"
#include "Transcriptome.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ReadAnnotations.h"
#include "SoloReadBarcode.h"
#include "hash_shims_cpp_compat.h"
#include "ReadAlign.h"
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <string>
#include <cstdio>
#include <mutex>
#include <cstdlib>
#include <algorithm>
#include "FlexDebugCounters.h"

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

// Debug counters for probe/genomic resolution - now use thread-local counters via FlexDebugCounters.h
// Enable with -DDEBUG_FLEX_COUNTERS at compile time
#ifdef DEBUG_FLEX_COUNTERS
extern "C" uint64_t solo_probe_align_count() { return flexCountersSumAll().probeAlignCount; }
extern "C" uint64_t solo_genomic_align_count() { return flexCountersSumAll().genomicAlignCount; }
extern "C" uint64_t solo_probe_resolved_count() { return flexCountersSumAll().probeResolvedCount; }
extern "C" uint64_t solo_genomic_resolved_count() { return flexCountersSumAll().genomicResolvedCount; }
extern "C" uint64_t solo_resolver_dropped_count() { return flexCountersSumAll().resolverDropped; }
extern "C" uint64_t solo_genomic_dropped_mapq() { return flexCountersSumAll().genomicDroppedMapq; }
extern "C" uint64_t solo_genomic_dropped_nm() { return flexCountersSumAll().genomicDroppedNM; }
extern "C" uint64_t solo_genomic_dropped_mmrate() { return flexCountersSumAll().genomicDroppedMMrate; }
extern "C" uint64_t solo_probe_dropped_nm() { return flexCountersSumAll().probeDroppedNM; }
extern "C" uint64_t solo_probe_dropped_mmrate() { return flexCountersSumAll().probeDroppedMMrate; }
extern "C" uint64_t solo_probe_missing_idx_count() { return flexCountersSumAll().probeMissingIdx; }
extern "C" uint64_t solo_resolver_drop_probe_disagree_count() { return flexCountersSumAll().resolverDropProbeDisagree; }
extern "C" uint64_t solo_resolver_drop_genomic_disagree_count() { return flexCountersSumAll().resolverDropGenomicDisagree; }
extern "C" uint64_t solo_resolver_drop_mixed_count() { return flexCountersSumAll().resolverDropMixed; }
extern "C" uint64_t solo_resolver_drop_no_candidates_count() { return flexCountersSumAll().resolverDropNoCandidates; }
extern "C" uint64_t solo_resolver_keep_probe_count() { return flexCountersSumAll().resolverKeepProbe; }
extern "C" uint64_t solo_resolver_keep_genomic_count() { return flexCountersSumAll().resolverKeepGenomic; }
extern "C" uint64_t solo_genomic_align_with_probe_genes_count() { return flexCountersSumAll().genomicAlignWithProbeGenes; }
extern "C" uint64_t solo_genomic_only_reads_with_probe_genes_count() { return flexCountersSumAll().genomicOnlyReadsWithProbeGenes; }
extern "C" uint64_t solo_genomic_only_probe_gene_count() { return flexCountersSumAll().genomicOnlyProbeGeneCount; }
#else
extern "C" uint64_t solo_probe_align_count() { return 0; }
extern "C" uint64_t solo_genomic_align_count() { return 0; }
extern "C" uint64_t solo_probe_resolved_count() { return 0; }
extern "C" uint64_t solo_genomic_resolved_count() { return 0; }
extern "C" uint64_t solo_resolver_dropped_count() { return 0; }
extern "C" uint64_t solo_genomic_dropped_mapq() { return 0; }
extern "C" uint64_t solo_genomic_dropped_nm() { return 0; }
extern "C" uint64_t solo_genomic_dropped_mmrate() { return 0; }
extern "C" uint64_t solo_probe_dropped_nm() { return 0; }
extern "C" uint64_t solo_probe_dropped_mmrate() { return 0; }
extern "C" uint64_t solo_probe_missing_idx_count() { return 0; }
extern "C" uint64_t solo_resolver_drop_probe_disagree_count() { return 0; }
extern "C" uint64_t solo_resolver_drop_genomic_disagree_count() { return 0; }
extern "C" uint64_t solo_resolver_drop_mixed_count() { return 0; }
extern "C" uint64_t solo_resolver_drop_no_candidates_count() { return 0; }
extern "C" uint64_t solo_resolver_keep_probe_count() { return 0; }
extern "C" uint64_t solo_resolver_keep_genomic_count() { return 0; }
extern "C" uint64_t solo_genomic_align_with_probe_genes_count() { return 0; }
extern "C" uint64_t solo_genomic_only_reads_with_probe_genes_count() { return 0; }
extern "C" uint64_t solo_genomic_only_probe_gene_count() { return 0; }
#endif

// Global ProbeListIndex loader for inline path (independent of SoloFeature owner)
static ProbeListIndex* gProbeIndex = nullptr;
static bool gProbeIndexLoaded = false;

// Reject logging infrastructure for trace-drops debugging
static FILE* g_rejectLogFile = nullptr;
static std::mutex g_rejectLogMutex;
static bool g_rejectLogInitialized = false;
static bool g_rejectLogEnabled = false;
static bool g_rejectLogTraceQname = false;
static std::unordered_map<uint64_t, std::string> g_iReadToQname;
static std::mutex g_qnameMapMutex;

// Initialize reject logging from environment variable
static void initRejectLogging() {
    if (g_rejectLogInitialized) return;
    g_rejectLogInitialized = true;
    
    const char* logPath = std::getenv("STAR_INLINE_REJECT_LOG");
    if (!logPath || logPath[0] == '\0') {
        g_rejectLogEnabled = false;
        return;
    }
    
    g_rejectLogFile = fopen(logPath, "a");
    if (!g_rejectLogFile) {
        fprintf(stderr, "[WARNING] Failed to open reject log file: %s\n", logPath);
        g_rejectLogEnabled = false;
        return;
    }
    
    g_rejectLogEnabled = true;
    
    // Check for optional QNAME tracing
    const char* traceQname = std::getenv("STAR_INLINE_TRACE_QNAME");
    if (traceQname && traceQname[0] != '\0' && (traceQname[0] == '1' || traceQname[0] == 'y' || traceQname[0] == 'Y')) {
        g_rejectLogTraceQname = true;
    }
    
    // Write TSV header
    fprintf(g_rejectLogFile, "sample_tag\tiRead\tqname\tcb16\tumi24_hex\tfeatureType\tisProbe\tgeneIdx15\treason\textra\n");
    fflush(g_rejectLogFile);
    
    // Register cleanup
    std::atexit([]() {
        std::lock_guard<std::mutex> lock(g_rejectLogMutex);
        if (g_rejectLogFile) {
            fclose(g_rejectLogFile);
            g_rejectLogFile = nullptr;
        }
    });
}

// Store iRead -> qname mapping (called from ReadAlign::outputAlignments)
// Note: This function is called from ReadAlign_outputAlignments.cpp
void storeQnameMapping(uint64_t iRead, const char* qname) {
    if (!g_rejectLogTraceQname || iRead == (uint64_t)-1 || !qname) {
        return;
    }
    std::lock_guard<std::mutex> lock(g_qnameMapMutex);
    g_iReadToQname[iRead] = std::string(qname);
}

// Get qname for iRead (thread-safe)
static std::string getQnameForIRead(uint64_t iRead) {
    if (!g_rejectLogTraceQname || iRead == (uint64_t)-1) {
        return std::string();
    }
    std::lock_guard<std::mutex> lock(g_qnameMapMutex);
    auto it = g_iReadToQname.find(iRead);
    if (it != g_iReadToQname.end()) {
        return it->second;
    }
    return std::string();
}

// Helper function to format sample tag as numeric_index:label
static std::string formatSampleTag(uint8_t detectedToken, const ParametersSolo& pSolo) {
    if (detectedToken == 0xFF) {
        return std::string(); // No tag detected
    }
    
    uint16_t sampleIdx = SampleDetector::sampleIndexForToken(detectedToken);
    if (sampleIdx == 0) {
        return std::string(); // Invalid token
    }
    
    // Get user-provided label first (supports arbitrary labels)
    std::string label = SampleDetector::labelForIndexStatic(sampleIdx);
    
    // If no user label, fall back to canonical sequence (8-mer)
    if (label.empty()) {
        label = SampleDetector::canonicalForIndexStatic(sampleIdx);
    }
    
    // Final fallback: BC format with sample index
    if (label.empty()) {
        char buf[32];
        snprintf(buf, sizeof(buf), "BC%u", sampleIdx);
        label = buf;
    }
    
    char result[64];
    snprintf(result, sizeof(result), "%u:%s", sampleIdx, label.c_str());
    return std::string(result);
}

// Log reject reason (thread-safe)
static void logRejectReason(const SoloReadBarcode& soloBar, uint64_t iRead, int32_t featureType,
                            uint8_t isProbe, uint16_t geneIdx15, const char* reason, const char* extra,
                            const ParametersSolo& pSolo) {
    // Only log if inline hash mode is enabled
    if (!pSolo.inlineHashMode) {
        return;
    }
    
    // Initialize logging on first call
    if (!g_rejectLogInitialized) {
        initRejectLogging();
    }
    
    if (!g_rejectLogEnabled || !g_rejectLogFile) {
        return;
    }
    
    // Format sample tag once (needed for filtering and logging)
    std::string sampleTag = formatSampleTag(soloBar.detectedSampleToken, pSolo);
    
    // Optional filters at write time
    const char* filterSample = std::getenv("STAR_INLINE_REJECT_FILTER_SAMPLE");
    if (filterSample && filterSample[0] != '\0') {
        if (sampleTag.find(filterSample) == std::string::npos && !sampleTag.empty()) {
            return; // Filter out this sample
        }
    }
    
    const char* filterCB = std::getenv("STAR_INLINE_REJECT_FILTER_CB");
    if (filterCB && filterCB[0] != '\0') {
        if (soloBar.cbMatchString != filterCB) {
            return; // Filter out this CB
        }
    }
    
    std::lock_guard<std::mutex> lock(g_rejectLogMutex);
    
    // Format iRead (empty if not available)
    std::string iReadStr;
    if (iRead != (uint64_t)-1) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%llu", (unsigned long long)iRead);
        iReadStr = buf;
    }
    
    // Get qname if tracing enabled
    std::string qnameStr = getQnameForIRead(iRead);
    
    // Format CB16 (empty if no CB match)
    std::string cb16Str = soloBar.cbMatchString.empty() ? std::string() : soloBar.cbMatchString;
    
    // Format UMI24 hex
    uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
    char umi24Hex[16];
    snprintf(umi24Hex, sizeof(umi24Hex), "0x%06X", umi24);
    
    // Get feature type name
    std::string featTypeName;
    if (featureType >= 0 && featureType < (int32_t)SoloFeatureTypes::Names.size()) {
        featTypeName = SoloFeatureTypes::Names[featureType];
    } else if (featureType == -1) {
        featTypeName = "NoFeature";
    } else {
        featTypeName = "Unknown";
    }
    
    // Write TSV line
    fprintf(g_rejectLogFile, "%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\n",
            sampleTag.c_str(),
            iReadStr.c_str(),
            qnameStr.c_str(),
            cb16Str.c_str(),
            umi24Hex,
            featTypeName.c_str(),
            (unsigned)isProbe,
            (unsigned)geneIdx15,
            reason,
            extra ? extra : "");
    fflush(g_rejectLogFile);
}

const ProbeListIndex* getGlobalProbeIndex(const SoloReadFeature* rf) {
    if (gProbeIndexLoaded) return gProbeIndex;
    gProbeIndexLoaded = true;
    if (rf == nullptr) return nullptr;
    const std::string& path = rf->pSolo.probeListPath;
    if (path.empty() || path == "-") return nullptr;
    ProbeListIndex* idx = new ProbeListIndex();
    uint32_t deprecatedCount = 0;
    if (!idx->load(path, rf->pSolo.removeDeprecated, &deprecatedCount)) {
        delete idx;
        return nullptr;
    }
    if (rf->pSolo.removeDeprecated && deprecatedCount > 0) {
        // Note: Cannot log here as we don't have access to logMain in this context
        // Logging will happen in STAR.cpp initialization
    }
    gProbeIndex = idx;
    return gProbeIndex;
}

class ReadSoloFeatures {
public:
    uint32 gene;
    vector<uint32> geneMult;
    vector<array<uint64,2>> sj;
    bool sjAnnot;
    uint32 indAnnotTr; //index of the annotated transcript
    Transcript **alignOut;
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
                    const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot, const SoloReadFlagClass &readFlag,
                    SoloReadFeature *soloReadFeat = nullptr); // Optional pointer for inline hash capture

void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    // Check for NO_GENE debug mode
    const char* dbgNoGene = std::getenv("STAR_INLINE_DEBUG_NOGENE");
    bool debugNoGene = dbgNoGene && dbgNoGene[0] != '\0';
    
    if (pSolo.type==0)
        return;

    if (pSolo.readStatsYes[featureType]) {//readFlag

        if (nTr==1) {
            readFlag.setBit(readFlag.genomeU);
        } else if (nTr>1) {
            readFlag.setBit(readFlag.genomeM);
        };

        for (uint64 itr=0; itr<nTr; itr++) {
            if (P.pGe.chrSet.mito.count(alignOut[itr]->Chr) == 1) {
                readFlag.setBit(readFlag.mito);
            };
        };

        switch (readAnnot.annotFeatures[featureType].ovType) {
            case ReadAnnotFeature::overlapTypes::exonic : 
            case ReadAnnotFeature::overlapTypes::exonic50p :
                readFlag.setBit(readFlag.exonic);
                break;
            case ReadAnnotFeature::overlapTypes::intronic : 
                readFlag.setBit(readFlag.intronic);
                break;
            case ReadAnnotFeature::overlapTypes::exonicAS : 
            case ReadAnnotFeature::overlapTypes::exonic50pAS :
                readFlag.setBit(readFlag.exonicAS);            
                break;
            case ReadAnnotFeature::overlapTypes::intronicAS :
                readFlag.setBit(readFlag.intronicAS);            
        };

        if (soloBar.cbMatch<0 && pSolo.cbWLyes) {//no CB match in the WL
            if (readAnnot.annotFeatures[featureType].fSet.size()==1) {
                readFlag.setBit(readFlag.featureU);
            } else if (readAnnot.annotFeatures[featureType].fSet.size()>1){
                readFlag.setBit(readFlag.featureM);
            };
            readFlag.setBit(readFlag.cbMatch);//this will counts reads with no CB match
            readFlag.countsAddNoCB();
        };
    };

    if (soloBar.cbMatch<0) {
        // Log NO_CB_MATCH before returning
        if (pSolo.inlineHashMode) {
            logRejectReason(soloBar, iRead, featureType, 0, 0, "NO_CB_MATCH", "", pSolo);
        }
        return;
    }

       
    ReadSoloFeatures reFe;
    reFe.alignOut=alignOut;
    reFe.indAnnotTr = 0;    

    uint32 nFeat=0; //number of features in this read (could be >1 for SJs)
    if (nTr==0) {//unmapped
        stats.V[stats.noUnmapped]++;
        
    } else {
        switch (featureType) {
            case SoloFeatureTypes::Gene :
            case SoloFeatureTypes::GeneFull :
            case SoloFeatureTypes::GeneFull_Ex50pAS :
            case SoloFeatureTypes::GeneFull_ExonOverIntron :
                {
                    auto *readGe = &readAnnot.annotFeatures[featureType].fSet;

                    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq) {
                        for (int32 itr=nTr-1; itr>=0; itr--) {
                            if (readAnnot.annotFeatures[featureType].fAlign[itr].size() > 0) {
                                reFe.indAnnotTr = itr;
                                break;//indAnnotTr is the first genic align, and is used for Smart-seq deduplication, below in outputReadCB()
                                      //TODO this is ambiguous for reads that map to the same gene multiple times. *Last* align is chosen to match the previous version.
                            };
                        };
                    };
                        
                    if (readGe->size()==0) {//check genes
                        stats.V[stats.noNoFeature]++;//no gene
                        
                        // Debug instrumentation for NO_GENE reads
                        if (debugNoGene) {
                            std::string qn = getQnameForIRead(iRead);
                            std::string cbStr = soloBar.cbMatchString.empty() ? "(no_cb)" : soloBar.cbMatchString;
                            
                            std::fprintf(stderr,
                                "[NO_GENE_DBG] qname=%s cb=%s nTr=%u feature=%s\n",
                                qn.empty() ? "(unknown)" : qn.c_str(),
                                cbStr.c_str(),
                                (unsigned)nTr,
                                SoloFeatureTypes::Names[featureType].c_str());
                            
                            const auto &featAnnot = readAnnot.annotFeatures[featureType];
                            for (uint32_t ia = 0; ia < featAnnot.fAlign.size(); ++ia) {
                                Transcript *tr = (alignOut && alignOut[ia]) ? alignOut[ia] : nullptr;
                                const char* chrName = tr && !tr->chrName.empty() ? tr->chrName.c_str() : 
                                                      (tr ? "(chr_idx)" : "(null)");
                                
                                std::fprintf(stderr, "  [NO_GENE_DBG] ia=%u chr=%s genes={", ia, chrName);
                                bool first = true;
                                for (uint32_t g : featAnnot.fAlign[ia]) {
                                    if (!first) std::fprintf(stderr, ",");
                                    first = false;
                                    std::fprintf(stderr, "%u", g);
                                }
                                std::fprintf(stderr, "}\n");
                            }
                        }
                        
                        // Log NO_GENE
                        if (pSolo.inlineHashMode) {
                            logRejectReason(soloBar, iRead, featureType, 0, 0, "NO_GENE", "", pSolo);
                        }
                    } else if (readGe->size()>1) {
                        stats.V[stats.MultiFeature]++;//multi-gene reads
                        readFlag.setBit(readFlag.featureM);
                        if (nTr>1)
                            stats.V[stats.subMultiFeatureMultiGenomic]++;//multigene caused by multimapper
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
                            
                        if (pSolo.multiMap.yes.multi) {//output multimappers
                            reFe.geneMult.resize(readGe->size());
                            uint32 ii=0;
                            for (auto &g: *readGe) {//set high bit to mark multimappers
                                reFe.geneMult[ii] = g | geneMultMark;
                                ++ii;
                            };
                                
                            nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot, readFlag, this);
                        } else {
                            // Multi-gene but multimappers disabled - drop the read
                            if (pSolo.inlineHashMode) {
                                logRejectReason(soloBar, iRead, featureType, 0, 0, "MULTIGENE_DROPPED", "", pSolo);
                            }
                        };
                    } else {//unique-gene reads
                        reFe.gene = *readGe->begin();
                        readFlag.setBit(readFlag.featureU);
                        nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (uint64)-1), featureType, soloBar, reFe, readAnnot, readFlag, this);
                    };

                    //debug
                    //{
                    //    if (SoloFeatureTypes::Gene==featureType && ((readFlag.checkBit(readFlag.featureM)==1|readFlag.checkBit(readFlag.featureU)==1) != readFlag.checkBit(readFlag.exonic) ))
                    //        cout << iRead<<" "<<readFlag.checkBit(readFlag.featureM)<<" "<<readFlag.checkBit(readFlag.featureU)<<" "<<readFlag.checkBit(readFlag.exonic)<<endl;
                    //};                    
                };
                break;
        
            case SoloFeatureTypes::SJ : 
                if (nTr>1) {//reject all multimapping reads
                    stats.V[stats.subMultiFeatureMultiGenomic]++;
                    stats.V[stats.MultiFeature]++;
                //} else if (readAnnot.geneConcordant.size()>1){//for SJs, still check genes, no feature if multi-gene
                //    stats.V[stats.MultiFeature]++;
                } else {//one gene or no gene
                    alignOut[0]->extractSpliceJunctions(reFe.sj, reFe.sjAnnot);
                    //if ( reFe.sj.empty() || (reFe.sjAnnot && readAnnot.geneConcordant.size()==0) ) {//no junctions, or annotated junction but no gene (i.e. read does not fully match transcript)
                    if ( reFe.sj.empty() ) {
                        stats.V[stats.noNoFeature]++;
                    } else {//good junction
                        readFlag.setBit(readFlag.featureU);
                        nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (uint64)-1), featureType, soloBar, reFe, readAnnot, readFlag, this);
                    };
                };                  
                break;
        
            case SoloFeatureTypes::Transcript3p :
                if (readAnnot.transcriptConcordant.size()==0 || soloBar.cbMatch>1) {//do not record ambiguous CB  
                    stats.V[stats.noNoFeature]++;
                } else {
                    nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot, readFlag, this);
                };
                if (readAnnot.transcriptConcordant.size()==1 && readAnnot.transcriptConcordant[0][1] < transcriptDistCount.size()) {
                    //read maps to one transcript - use for distTTS distribution function
                    transcriptDistCount[readAnnot.transcriptConcordant[0][1]]++;
                };
                
                break;

            case SoloFeatureTypes::Velocyto :
                //different record: iRead, nTr, tr1, type1, tr2, type2 ...
                if (readAnnot.trVelocytoType.size()>0) {//otherwise, no gene
                    
                    sort(readAnnot.trVelocytoType.begin(), readAnnot.trVelocytoType.end(),
                         [](const trTypeStruct &t1, const trTypeStruct &t2) {return t1.tr < t2.tr;});

                    if (!pSolo.inlineHashMode && streamReads) {
                        *streamReads << iRead <<' '<< readAnnot.trVelocytoType.size();
                        for (auto &tt: readAnnot.trVelocytoType)
                             *streamReads <<' '<< tt.tr <<' '<< (uint32) tt.type;
                        *streamReads <<'\n';
                    }
                    nFeat=1;
                } else {
                    stats.V[stats.noNoFeature]++;
                };
                break; //no need to go with downstream processing                
                
        };//switch (featureType)
    };//if (nTr==0)
    
    if ( nFeat==0 && (readInfoYes | pSolo.readStatsYes[featureType]) ) {//no feature, but readInfo requested
        outputReadCB(streamReads, iRead, (uint32)-1, soloBar, reFe, readAnnot, readFlag, this);
    };
    
    if (nFeat==0)
        return; //no need to record the number of reads per CB
    
    if (pSolo.cbWLyes) {//WL
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nFeat;
    } else {//no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    };
    
    /*//if we wanted to record all barcodes, even those without features
        if (!pSolo.cbWLyes) {//no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    } else if (nFeat>0){//WL && nFeat>0
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nFeat;
    };
    */

    return;
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, 
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
    
    // Helper lambda to extract tag whitelist index from detected sample token
    // Tag extraction: Uses SampleDetector::sampleIndexForToken() which returns 1-based sample index
    // Default path: tagIdx=0 if sample detection disabled or no tag detected (detectedSampleToken==0xFF)
    // Bounds: Clamped to 5 bits (0-31) to fit in packed key format [CB20][UMI24][GENE15][TAG5]
    auto extractTagIdx = [&]() -> uint8_t {
        uint8_t tagIdx = 0; // Default: no tag
        if (soloBar.detectedSampleToken != 0xFF) {
            uint16_t sampleIdx = SampleDetector::sampleIndexForToken(soloBar.detectedSampleToken);
            if (sampleIdx > 0) {
                tagIdx = static_cast<uint8_t>(sampleIdx & 0x1F); // Clamp to 5 bits (max 31)
            }
        }
        return tagIdx;
    };

    // Compute tag once per read. If sample tagging is enabled (whitelist provided or require-match set)
    // and no tag was detected, drop the read to avoid propagating tagIdx=0 into the hash/MEX.
    const uint8_t tagIdx = extractTagIdx();
    const bool dropUnmatchedTag = ( (!soloBar.pSolo.sampleWhitelistPath.empty()) || soloBar.pSolo.sampleRequireMatch ) && (tagIdx == 0);
    if (dropUnmatchedTag) {
        return 0;
    }
    
    // Helper lambda to handle ambiguous CB accumulation with gene/tag info
    // Ambiguous CB: multiple whitelist candidates (cbMatchInd.size() > 1) or cbMatch > 1 (multiple matches)
    auto accumulateAmbiguousCB = [&](uint16_t geneIdx, uint8_t tagIdx) {
        bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
        if (!soloReadFeat || !soloReadFeat->inlineHash_ || !isAmbiguous || soloBar.cbMatchString.empty() || soloBar.cbMatchInd.empty()) {
            return;
        }
        // Ambiguous CB: accumulate for resolution with gene/tag/umi info
        uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
        ReadAlign::AmbigKey ambigKey = ReadAlign::hashCbSeq(soloBar.cbMatchString);
        auto &entry = soloReadFeat->pendingAmbiguous_[ambigKey];
        
        if (entry.candidateIdx.empty()) {
            // First time seeing this ambiguous CB: initialize candidates
            entry.candidateIdx.reserve(soloBar.cbMatchInd.size());
            for (auto idx : soloBar.cbMatchInd) {
                entry.candidateIdx.push_back(static_cast<uint32_t>(idx + 1)); // 1-based for resolver
            }
            entry.cbSeq = soloBar.cbMatchString;
            entry.cbQual = soloBar.cbQual;
            entry.umiCounts.reserve(32);
        }
        
        // Accumulate UMI count (24-bit packed UMI -> count)
        entry.umiCounts[umi24]++;
        
        // Store gene/tag/umi observation for later hash creation after resolution
        SoloReadFeature::ExtendedAmbiguousEntry::AmbiguousObservation obs;
        obs.geneIdx = geneIdx;
        obs.tagIdx = tagIdx;
        obs.umi24 = umi24;
        obs.count = 1;
        entry.observations.push_back(obs);
    };
    
    uint64 nout=1;
    
    switch (featureType) {
        case -1 : {
            // no feature, output for readInfo
            const uint16_t geneIdx = 0; // No feature
            const bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
            if (isAmbiguous && !soloBar.cbMatchInd.empty()) {
                accumulateAmbiguousCB(geneIdx, tagIdx);
            } else if (soloReadFeat && soloReadFeat->inlineHash_ != nullptr && soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                uint32_t cbIdx = soloBar.cbMatchInd[0];
                uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
                uint64_t key = packCgAggKey(cbIdx, umi24, geneIdx, tagIdx);
                int absent;
                khiter_t iter = kh_put(cg_agg, soloReadFeat->inlineHash_, key, &absent);
                if (absent) {
                    kh_val(soloReadFeat->inlineHash_, iter) = 1;
                } else {
                    kh_val(soloReadFeat->inlineHash_, iter)++;
                }
            }
            if (!soloBar.pSolo.inlineHashMode && streamOut) {
                *streamOut << soloBar.umiB <<' '<< iRead <<' '<< readFlag.flag <<' '<< -1 <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
            }
            break;
        }
            
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron : {
            const bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);

            // Build CandidateView per alignment (probe/genomic) for probe-aware resolution
            std::vector<CandidateView> candidates;
            const auto &featAnnot = readAnnot.annotFeatures[featureType];
            const auto &fAlign = featAnnot.fAlign;
            uint32_t nAlign = fAlign.size();

            auto lookupProbeIdx = [&](uint32_t geneIdx) -> uint16_t {
                return SoloFeature::getProbeIdxForGene(geneIdx);
            };

            // One-time sanity check removed after validation (previously dumped gene→probe table)
            // If future debugging is needed, reintroduce lightweight logging here.

            bool hasProbeCandidates = false;
            bool hasGenomicProbeGenes = false;

            for (uint32_t ia = 0; ia < nAlign; ++ia) {
                if (!reFe.alignOut || reFe.alignOut[ia] == nullptr) continue;
                Transcript *tr = reFe.alignOut[ia];
                const std::string &chrName = tr->chrName;
                bool isProbeChr = !chrName.empty() && chrName.rfind("ENSG", 0) == 0;

                auto isCanonicalProbeCigar = [](const std::string &cig) -> bool {
                    // Mirror bam_to_counts probe CIGAR whitelist: 40S50M, 50M40S, 50M
                    // This is a conservative filter; genomic alignments are not filtered here.
                    return cig == "40S50M" || cig == "50M40S" || cig == "50M";
                };

                CandidateView cv;
                cv.mapq = tr->mapq;
                cv.asScore = tr->asScore;
                cv.nm = tr->nm;
                cv.probeCigarOk = true;
                cv.zgGeneIdx15.clear();

                if (isProbeChr) {
                    cv.isGenomic = false;
                    cv.probeCigarOk = isCanonicalProbeCigar(tr->cigarString);
                    FLEX_COUNT_INC(probeAlignCount);
                    
                    // NM filter for probe alignments (match bam_to_counts behavior)
                    // soloNMmax: max allowed mismatches; -1 disables
                    if (soloBar.pSolo.nmMax >= 0 && cv.nm >= 0) {
                        if (cv.nm > soloBar.pSolo.nmMax) {
                            FLEX_COUNT_INC(probeDroppedNM);
                            continue; // skip high-mismatch probe alignment
                        }
                    }
                    
                    // MM rate filter for probe alignments (match bam_to_counts behavior)
                    // soloMMrateMax: max allowed mismatch rate (NM / aligned_len)
                    // Only apply if nmMax is disabled (-1) and we have valid NM
                    if (soloBar.pSolo.nmMax < 0 && cv.nm >= 0 && soloBar.pSolo.mmRateMax > 0) {
                        // For probes, use mapped length from transcript
                        int alignedLen = static_cast<int>(tr->exons[tr->nExons-1][EX_G] - tr->exons[0][EX_G] + tr->exons[tr->nExons-1][EX_L]);
                        if (alignedLen > 0) {
                            double mmRate = static_cast<double>(cv.nm) / alignedLen;
                            if (mmRate > soloBar.pSolo.mmRateMax) {
                                FLEX_COUNT_INC(probeDroppedMMrate);
                                continue; // skip high-mismatch-rate probe alignment
                            }
                        }
                    }
                    
                    uint16_t probeIdx = 0;
                    const ProbeListIndex* probeIdxTable = getGlobalProbeIndex(soloReadFeat);
                    if (probeIdxTable) {
                        std::string ensgId = chrName.substr(0, std::min<size_t>(15, chrName.size()));
                        probeIdx = probeIdxTable->geneIndex15(ensgId);
                    }
                    if (probeIdx != 0) {
                        cv.geneIdx15 = probeIdx;
                        cv.zgGeneIdx15.push_back(probeIdx);
                        hasProbeCandidates = true;
                    } else {
                        FLEX_COUNT_INC(probeMissingIdx);
                        // Log NO_PROBE_IDX
                        if (soloBar.pSolo.inlineHashMode) {
                            char extraBuf[256];
                            snprintf(extraBuf, sizeof(extraBuf), "chr=%s", chrName.c_str());
                            logRejectReason(soloBar, iRead, featureType, 1, 0, "NO_PROBE_IDX", extraBuf, soloBar.pSolo);
                        }
                        continue; // no usable gene for this alignment
                    }
                } else {
                    cv.isGenomic = true;
                    FLEX_COUNT_INC(genomicAlignCount);
                    
                    // MAPQ filter for genomic alignments (match bam_to_counts behavior)
                    // Only applies when soloMapqMode is "genomic" or "all"
                    if (soloBar.pSolo.mapqMode == ParametersSolo::MapqGenomicOnly ||
                        soloBar.pSolo.mapqMode == ParametersSolo::MapqAll) {
                        if (cv.mapq >= 0 && cv.mapq < soloBar.pSolo.mapqThreshold) {
                            FLEX_COUNT_INC(genomicDroppedMapq);
                            continue; // skip low-MAPQ genomic alignment
                        }
                    }
                    
                    // NM filter for genomic alignments (match bam_to_counts behavior)
                    // soloNMmax: max allowed mismatches; -1 disables
                    if (soloBar.pSolo.nmMax >= 0 && cv.nm >= 0) {
                        if (cv.nm > soloBar.pSolo.nmMax) {
                            FLEX_COUNT_INC(genomicDroppedNM);
                            continue; // skip high-mismatch genomic alignment
                        }
                    }
                    
                    // MM rate filter for genomic alignments (match bam_to_counts behavior)
                    // soloMMrateMax: max allowed mismatch rate (NM / read_len)
                    // Only apply if nmMax is disabled (-1) and we have valid NM
                    if (soloBar.pSolo.nmMax < 0 && cv.nm >= 0 && soloBar.pSolo.mmRateMax > 0) {
                        // Use nM (mismatches) and aligned length from transcript
                        int alignedLen = static_cast<int>(tr->exons[tr->nExons-1][EX_G] - tr->exons[0][EX_G] + tr->exons[tr->nExons-1][EX_L]);
                        if (alignedLen > 0) {
                            double mmRate = static_cast<double>(cv.nm) / alignedLen;
                            if (mmRate > soloBar.pSolo.mmRateMax) {
                                FLEX_COUNT_INC(genomicDroppedMMrate);
                                continue; // skip high-mismatch-rate genomic alignment
                            }
                        }
                    }
                    
                    // Map per-alignment genes to probe indices (if available)
                    if (ia < fAlign.size()) {
                        for (uint32_t g : fAlign[ia]) {
                            uint16_t probeIdx = lookupProbeIdx(g);
                            if (probeIdx != 0) {
                                cv.zgGeneIdx15.push_back(probeIdx);
                            }
                        }
                    }
                    if (cv.zgGeneIdx15.empty()) {
                        continue; // genomic align with no probe-mapped gene
                    }
                    // Genomic alignment that yielded >=1 probe gene
                    FLEX_COUNT_INC(genomicAlignWithProbeGenes);
                    hasGenomicProbeGenes = true;
                    if (cv.zgGeneIdx15.size() == 1) {
                        cv.geneIdx15 = cv.zgGeneIdx15[0];
                    } else {
                        cv.geneIdx15 = 0;
                    }
                }

                candidates.push_back(std::move(cv));
            }

            // Track reads where we have NO probe-chromosome candidates but DO
            // have genomic alignments that map to probe genes via ZG/lookup.
            if (!hasProbeCandidates && hasGenomicProbeGenes) {
                FLEX_COUNT_INC(genomicOnlyReadsWithProbeGenes);
                std::unordered_set<uint16_t> geneSet;
                for (const auto &cv : candidates) {
                    if (!cv.isGenomic) continue;
                    for (uint16_t g : cv.zgGeneIdx15) {
                        if (g != 0) geneSet.insert(g);
                    }
                }
                FLEX_COUNT_ADD(genomicOnlyProbeGeneCount, geneSet.size());
            }

            // Helper lambda to build candidate summary
            auto buildCandidateSummary = [](const std::vector<CandidateView>& candidates) {
                    std::vector<uint16_t> probeGenes;
                    std::vector<uint16_t> genomicGenes;
                size_t probeCount = 0;
                size_t genomicCount = 0;
                
                    for (const auto &cv : candidates) {
                        if (!cv.isGenomic) {
                        probeCount++;
                            if (cv.geneIdx15 != 0) {
                                probeGenes.push_back(cv.geneIdx15);
                            }
                            for (uint16_t g : cv.zgGeneIdx15) {
                                if (g != 0) probeGenes.push_back(g);
                            }
                        } else {
                        genomicCount++;
                            if (cv.geneIdx15 != 0) {
                                genomicGenes.push_back(cv.geneIdx15);
                            }
                            for (uint16_t g : cv.zgGeneIdx15) {
                                if (g != 0) genomicGenes.push_back(g);
                            }
                        }
                    }
                
                    // Remove duplicates
                    std::sort(probeGenes.begin(), probeGenes.end());
                    probeGenes.erase(std::unique(probeGenes.begin(), probeGenes.end()), probeGenes.end());
                    std::sort(genomicGenes.begin(), genomicGenes.end());
                    genomicGenes.erase(std::unique(genomicGenes.begin(), genomicGenes.end()), genomicGenes.end());
                    
                return std::make_tuple(probeGenes, genomicGenes, probeCount, genomicCount);
            };

            if (candidates.empty()) {
                FLEX_COUNT_INC(resolverDropped);
                FLEX_COUNT_INC(resolverDropNoCandidates);
                // Log RESOLVER_DROP with no candidates
                if (soloBar.pSolo.inlineHashMode) {
                    logRejectReason(soloBar, iRead, featureType, 0, 0, "RESOLVER_DROP", "reason=NO_CANDIDATES", soloBar.pSolo);
                }
                break;
            }

            uint16_t resolvedGeneIdx = resolveGeneFromCandidates(candidates);
            
            // If resolver returns 0, drop the read (don't insert into hash)
            // Note: This only happens for:
            //   - Probe gene conflicts (disagreeing probe genes)
            //   - Empty/no genes
            if (resolvedGeneIdx == 0) {
                FLEX_COUNT_INC(resolverDropped);
                // Log RESOLVER_DROP with candidate summary
                if (soloBar.pSolo.inlineHashMode) {
                    // Build candidate summary using helper
                    auto [probeGenes, genomicGenes, probeCount, genomicCount] = buildCandidateSummary(candidates);
                    
                    // Determine reason
                    std::string reason;
                    if (probeGenes.size() > 1 && genomicGenes.empty()) {
                        reason = "PROBE_DISAGREE";
                        FLEX_COUNT_INC(resolverDropProbeDisagree);
                    } else if (genomicGenes.size() > 1 && probeGenes.empty()) {
                        reason = "GENOMIC_DISAGREE";
                        FLEX_COUNT_INC(resolverDropGenomicDisagree);
                    } else if (probeGenes.size() > 0 && genomicGenes.size() > 0) {
                        reason = "MIXED";
                        FLEX_COUNT_INC(resolverDropMixed);
                    } else {
                        // Fallback: one gene in only one side (rare)
                        reason = probeGenes.size() > 0 ? "PROBE_DISAGREE" : "GENOMIC_DISAGREE";
                        if (probeGenes.size() > 0) {
                            FLEX_COUNT_INC(resolverDropProbeDisagree);
                        } else {
                            FLEX_COUNT_INC(resolverDropGenomicDisagree);
                        }
                    }
                    
                    // Format extra field: reason=...;probe={p1,p2};genomic={g1}
                    std::string extraStr = "reason=" + reason + ";probe={";
                    for (size_t i = 0; i < probeGenes.size(); ++i) {
                        if (i > 0) extraStr += ",";
                        char buf[16];
                        snprintf(buf, sizeof(buf), "%u", (unsigned)probeGenes[i]);
                        extraStr += buf;
                    }
                    extraStr += "};genomic={";
                    for (size_t i = 0; i < genomicGenes.size(); ++i) {
                        if (i > 0) extraStr += ",";
                        char buf[16];
                        snprintf(buf, sizeof(buf), "%u", (unsigned)genomicGenes[i]);
                        extraStr += buf;
                    }
                    extraStr += "}";
                    
                    // Determine isProbe: if all candidates are probe, set to 1; if all genomic, set to 0; if mixed, use majority
                    uint8_t isProbeDrop = 0;
                    if (probeCount > 0 && genomicCount == 0) {
                        isProbeDrop = 1; // All probe
                    } else if (probeCount == 0 && genomicCount > 0) {
                        isProbeDrop = 0; // All genomic
                    } else {
                        // Mixed: use majority or last candidate type
                        isProbeDrop = (probeCount >= genomicCount) ? 1 : 0;
                    }
                    
                    logRejectReason(soloBar, iRead, featureType, isProbeDrop, 0, "RESOLVER_DROP", extraStr.c_str(), soloBar.pSolo);
                }
                break; // drop read
            }
            // Count resolution type based on winning candidate
            bool resolvedGenomic = false;
            const CandidateView* winningCandidate = nullptr;
            for (const auto &cv : candidates) {
                if (cv.geneIdx15 == resolvedGeneIdx || 
                    std::find(cv.zgGeneIdx15.begin(), cv.zgGeneIdx15.end(), resolvedGeneIdx) != cv.zgGeneIdx15.end()) {
                    resolvedGenomic = cv.isGenomic;
                    winningCandidate = &cv;
                    break;
                }
            }
            if (resolvedGenomic) {
                FLEX_COUNT_INC(genomicResolvedCount);
                FLEX_COUNT_INC(resolverKeepGenomic);
            } else {
                FLEX_COUNT_INC(probeResolvedCount);
                FLEX_COUNT_INC(resolverKeepProbe);
            }
            
            // Log KEEP_HASH with winning candidate info
            if (soloBar.pSolo.inlineHashMode) {
                uint8_t isProbeKeep = winningCandidate ? (!winningCandidate->isGenomic ? 1 : 0) : 0;
                
                // Build candidate summary using helper
                auto [probeGenes, genomicGenes, probeCount, genomicCount] = buildCandidateSummary(candidates);
                
                // Determine win
                std::string win = winningCandidate && !winningCandidate->isGenomic ? "PROBE" : "GENOMIC";
                
                // Format extra field: win=...;gid=...;probe={...};genomic={...}
                std::string extraStr = "win=" + win + ";gid=";
                char gidBuf[16];
                snprintf(gidBuf, sizeof(gidBuf), "%u", (unsigned)resolvedGeneIdx);
                extraStr += gidBuf;
                extraStr += ";probe={";
                for (size_t i = 0; i < probeGenes.size(); ++i) {
                    if (i > 0) extraStr += ",";
                    char buf[16];
                    snprintf(buf, sizeof(buf), "%u", (unsigned)probeGenes[i]);
                    extraStr += buf;
                }
                extraStr += "};genomic={";
                for (size_t i = 0; i < genomicGenes.size(); ++i) {
                    if (i > 0) extraStr += ",";
                    char buf[16];
                    snprintf(buf, sizeof(buf), "%u", (unsigned)genomicGenes[i]);
                    extraStr += buf;
                }
                extraStr += "}";
                
                logRejectReason(soloBar, iRead, featureType, isProbeKeep, resolvedGeneIdx, "KEEP_HASH", extraStr.c_str(), soloBar.pSolo);
            }
            
            // Insert exactly one quartet key with resolved gene (no per-gene fanout)
            if (soloReadFeat && soloReadFeat->inlineHash_ != nullptr && soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                if (isAmbiguous) {
                    // Log AMBIG_CB when accumulating (will be resolved later)
                    if (soloBar.pSolo.inlineHashMode) {
                        char extraBuf[128];
                        snprintf(extraBuf, sizeof(extraBuf), "candidates=%zu", soloBar.cbMatchInd.size());
                        logRejectReason(soloBar, iRead, featureType, 
                                       (!winningCandidate || !winningCandidate->isGenomic ? 1 : 0), 
                                       resolvedGeneIdx, "AMBIG_CB", extraBuf, soloBar.pSolo);
                    }
                    accumulateAmbiguousCB(resolvedGeneIdx, tagIdx);
                } else {
                    uint32_t cbIdx = soloBar.cbMatchInd[0];
                    uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
                    uint64_t key = packCgAggKey(cbIdx, umi24, resolvedGeneIdx, tagIdx);
                    int absent;
                    khiter_t iter = kh_put(cg_agg, soloReadFeat->inlineHash_, key, &absent);
                    if (absent) {
                        kh_val(soloReadFeat->inlineHash_, iter) = 1;
                    } else {
                        kh_val(soloReadFeat->inlineHash_, iter)++;
                    }
                }
            }
            
            // Write to stream (legacy mode)
            if (!soloBar.pSolo.inlineHashMode && streamOut) {
                *streamOut << soloBar.umiB <<' ';//UMI
                if ( iRead != (uint64)-1 )
                    *streamOut << iRead <<' '<< readFlag.flag <<' ';//iRead
                *streamOut << resolvedGeneIdx <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
            }
            
            nout = 1; // Single resolved gene
            break;
        }

        case SoloFeatureTypes::SJ : {
            //sj - two numbers, multiple sjs per read
            const bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
            for (auto &sj : reFe.sj) {
                // Inline hash capture for SJ (use sj[0] as "gene" identifier)
                uint16_t geneIdx = sj[0]; // Use first SJ coordinate as gene identifier
                // Check for ambiguous CB (multiple candidates) vs non-ambiguous (single candidate)
                if (isAmbiguous && !soloBar.cbMatchInd.empty()) {
                    accumulateAmbiguousCB(geneIdx, tagIdx);
                } else if (soloReadFeat && soloReadFeat->inlineHash_ != nullptr && soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                    uint32_t cbIdx = soloBar.cbMatchInd[0];
                    uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
                    uint64_t key = packCgAggKey(cbIdx, umi24, geneIdx, tagIdx);
                    int absent;
                    khiter_t iter = kh_put(cg_agg, soloReadFeat->inlineHash_, key, &absent);
                    if (absent) {
                        kh_val(soloReadFeat->inlineHash_, iter) = 1;
                    } else {
                        kh_val(soloReadFeat->inlineHash_, iter)++;
                    }
                }
                if (!soloBar.pSolo.inlineHashMode && streamOut) {
                    *streamOut << soloBar.umiB <<' ';//UMI
                    if ( iRead != (uint64)-1 )
                        *streamOut << iRead <<' '<< readFlag.flag <<' ';//iRead            
                    *streamOut << sj[0] <<' '<< sj[1] <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n' << flush;
                }
            };
            nout=reFe.sj.size();
            break;
        }

        case SoloFeatureTypes::Transcript3p : {
            //transcript,distToTTS structure
            // Inline hash capture for Transcript3p (use transcript ID as gene identifier)
            if (!readAnnot.transcriptConcordant.empty()) {
                uint16_t geneIdx = readAnnot.transcriptConcordant[0][0]; // Use first transcript ID
                // Check for ambiguous CB (multiple candidates) vs non-ambiguous (single candidate)
                const bool isAmbiguous = (soloBar.cbMatchInd.size() > 1) || (soloBar.cbMatch > 1);
                if (isAmbiguous && !soloBar.cbMatchInd.empty()) {
                    accumulateAmbiguousCB(geneIdx, tagIdx);
                } else if (soloReadFeat && soloReadFeat->inlineHash_ != nullptr && soloBar.cbMatch >= 0 && soloBar.cbMatch <= 1 && !soloBar.cbMatchInd.empty()) {
                    uint32_t cbIdx = soloBar.cbMatchInd[0];
                    uint32_t umi24 = soloBar.umiB & 0xFFFFFF;
                    uint64_t key = packCgAggKey(cbIdx, umi24, geneIdx, tagIdx);
                    int absent;
                    khiter_t iter = kh_put(cg_agg, soloReadFeat->inlineHash_, key, &absent);
                    if (absent) {
                        kh_val(soloReadFeat->inlineHash_, iter) = 1;
                    } else {
                        kh_val(soloReadFeat->inlineHash_, iter)++;
                    }
                }
            }
            if (!soloBar.pSolo.inlineHashMode && streamOut) {
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
