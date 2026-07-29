#include "ReadAlign.h"
#include "readLoad.h"
#include "readBarcodeLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#include "SampleDetector.h"
#include "SoloReadFeature_record_shared.h"
#include "FlexPipeline.h"
#include "FlexHashScreen.h"
#include "SpatialGex.h"
#include "libtrim/trim.h"
#include <cstdlib>
#include <atomic>
#include <chrono>
#include <mutex>

// --- Hash screen dump (env-gated: STAR_DUMP_HASH_SCREEN=<path>) ---
// Binary format: 8-byte magic "HSCRN001", 8-byte nReads placeholder,
// then per-read: {uint32 readLen, uint16 sampleIdx, uint8 action,
//   uint8 cacheClass, uint16 geneIdx15, int8 offset, uint8 pad,
//   char readSeq[readLen]}.
namespace {
struct HashScreenDumpState {
    FILE* fp = nullptr;
    std::mutex mu;
    uint64_t nReads = 0;
};
static HashScreenDumpState g_hsDump;
static std::once_flag g_hsDumpOnce;

void hashScreenDumpFinalize() {
    if (!g_hsDump.fp) return;
    std::fseek(g_hsDump.fp, 8, SEEK_SET);
    std::fwrite(&g_hsDump.nReads, 8, 1, g_hsDump.fp);
    std::fclose(g_hsDump.fp);
    g_hsDump.fp = nullptr;
}

void hashScreenDumpInit() {
    const char* path = std::getenv("STAR_DUMP_HASH_SCREEN");
    if (!path || path[0] == '\0') return;
    g_hsDump.fp = std::fopen(path, "wb");
    if (!g_hsDump.fp) return;
    const char magic[8] = {'H','S','C','R','N','0','0','1'};
    std::fwrite(magic, 1, 8, g_hsDump.fp);
    uint64_t placeholder = 0;
    std::fwrite(&placeholder, 8, 1, g_hsDump.fp);
    std::atexit(hashScreenDumpFinalize);
}

void hashScreenDumpWrite(const char* readSeq, uint32_t readLen,
                         uint16_t sampleIdx, const FlexHashScreenDecision& d) {
    std::call_once(g_hsDumpOnce, hashScreenDumpInit);
    if (!g_hsDump.fp) return;
    std::lock_guard<std::mutex> lock(g_hsDump.mu);
    std::fwrite(&readLen, 4, 1, g_hsDump.fp);
    std::fwrite(&sampleIdx, 2, 1, g_hsDump.fp);
    uint8_t action = static_cast<uint8_t>(d.action);
    std::fwrite(&action, 1, 1, g_hsDump.fp);
    std::fwrite(&d.cacheClass, 1, 1, g_hsDump.fp);
    std::fwrite(&d.geneIdx15, 2, 1, g_hsDump.fp);
    std::fwrite(&d.offset, 1, 1, g_hsDump.fp);
    std::fwrite(&d.negativeCode, 1, 1, g_hsDump.fp);
    std::fwrite(readSeq, 1, readLen, g_hsDump.fp);
    ++g_hsDump.nReads;
}
} // namespace

// Static counter for debug logging (guarded by STAR_TRIM_DEBUG_N env var)
// Use atomic for thread-safety when debug logging is enabled
static atomic<uint64_t> g_trimDebugCount(0);
static atomic<int64_t> g_trimDebugMax(-1);  // -1 means not initialized

int ReadAlign::oneRead() {//process one read: load, map, write

    //load read name, sequence, quality from the streams into internal arrays
    int readStatus[P.readNends];

    for (uint32 im=0; im<P.readNends; im++) {
        readStatus[im] = readLoad(*(readInStream[im]), P, readLength[im], readLengthOriginal[im], readNameMates[im], Read0[im], Read1[im], Qual0[im], clipMates[im], iReadAll, readFilesIndex, readFilter, readNameExtra[im]);
        if (readStatus[im] != readStatus[0]) {//check if the end of file was reached or not for all files
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: read files are not consistent, reached the end of the one before the other one\n";
            errOut << "SOLUTION: Check you your input files: they may be corrupted\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
    };

    if (readStatus[0]==-1) {//finished with the stream
        return -1;
    };

    // Integrated spatial mode is restricted to the ordinary paired FASTQ
    // reader. Keep its per-read payload in the gated Pipeline so default and
    // Flex packet paths retain the original ReadAlign object layout and work.
    if (P.spatialGexPipeline != nullptr) {
        if (P.readNends != 2 || Read0 == nullptr || Qual0 == nullptr
            || readLengthOriginal[1] == 0 || iReadAll == 0) {
            exitWithError("EXITING because integrated spatial mode did not receive paired raw R1 state\n",
                          std::cerr, P.inOut->logMain,
                          EXIT_CODE_INCONSISTENT_DATA, P);
        }
        std::string spatialError;
        if (!P.spatialGexPipeline->decodeCurrentThread(
                Read0[1], readLengthOriginal[1], Qual0[1],
                readLengthOriginal[1], iReadAll - 1, spatialError)) {
            exitWithError("EXITING because integrated spatial R1 decoding failed: "
                              + spatialError + "\n",
                          std::cerr, P.inOut->logMain,
                          EXIT_CODE_INCONSISTENT_DATA, P);
        }
    }
    
    return oneReadLoaded(readStatus[0]);

};

int ReadAlign::oneReadLoaded(const int readStatus0) {
    hasYAlignment_ = false;

    // Increment read counters BEFORE trimming (so dropped reads are counted)
    statsRA.readN++;
    statsRA.readBases += readLength[0] + (P.readNmates == 2 ? readLength[1] : 0);

    const auto completeSpatialWithoutFeature = [&]() {
        if (P.spatialGexPipeline == nullptr) return;
        if (iReadAll == 0) {
            exitWithError(
                "EXITING because integrated spatial mode lost a filtered read ordinal\n",
                std::cerr, P.inOut->logMain,
                EXIT_CODE_INCONSISTENT_DATA, P);
        }
        std::string spatialError;
        const spatial_gex::FeatureEvidenceClass source =
            P.soloSpatialFlexIntegratedEnabled
                ? spatial_gex::FeatureEvidenceClass::FlexAlignment
                : spatial_gex::FeatureEvidenceClass::Gex;
        if (!P.spatialGexPipeline->completeCurrentThread(
                source, false, 0, iReadAll - 1, spatialError)) {
            exitWithError(
                "EXITING because integrated spatial filtered-read completion failed: "
                    + spatialError + "\n",
                std::cerr, P.inOut->logMain,
                EXIT_CODE_INCONSISTENT_DATA, P);
        }
    };
    
    // Quality encoding debug check (BEFORE trimming) - guarded by STAR_TRIM_DEBUG_N env var
    int64_t debugMax = g_trimDebugMax.load();
    if (debugMax == -1) {
        const char* debugEnv = getenv("STAR_TRIM_DEBUG_N");
        int64_t newMax = debugEnv ? atol(debugEnv) : 0;
        int64_t expected = -1;
        if (g_trimDebugMax.compare_exchange_strong(expected, newMax)) {
            debugMax = newMax;
        } else {
            debugMax = g_trimDebugMax.load();
        }
    }
    uint64_t debugCount = g_trimDebugCount.fetch_add(1);
    if (debugMax > 0 && debugCount < (uint64_t)debugMax && P.trimCutadapt == "Yes") {
        // Log quality encoding for first N pairs before trimming
        uint32_t checkLen1 = min((uint32_t)readLength[0], 10u);
        uint32_t checkLen2 = (P.readNmates == 2) ? min((uint32_t)readLength[1], 10u) : 0;
        
        // Find min/max quality bytes
        uint8_t minQual1 = 255, maxQual1 = 0;
        uint8_t minQual2 = 255, maxQual2 = 0;
        for (uint32_t i = 0; i < readLength[0]; i++) {
            uint8_t q = (uint8_t)Qual0[0][i];
            if (q < minQual1) minQual1 = q;
            if (q > maxQual1) maxQual1 = q;
        }
        if (P.readNmates == 2) {
            for (uint32_t i = 0; i < readLength[1]; i++) {
                uint8_t q = (uint8_t)Qual0[1][i];
                if (q < minQual2) minQual2 = q;
                if (q > maxQual2) maxQual2 = q;
            }
        }
        
        P.inOut->logMain << "QUAL_CHECK: " << readNameMates[0]
                         << " len=" << readLength[0] << "," << readLength[1]
                         << " qual1_first10=\"";
        for (uint32_t i = 0; i < checkLen1; i++) {
            P.inOut->logMain << Qual0[0][i];
        }
        P.inOut->logMain << "\" qual1_bytes=";
        for (uint32_t i = 0; i < checkLen1; i++) {
            P.inOut->logMain << (uint32_t)(uint8_t)Qual0[0][i];
            if (i < checkLen1 - 1) P.inOut->logMain << ",";
        }
        P.inOut->logMain << " qual1_phred=";
        for (uint32_t i = 0; i < checkLen1; i++) {
            int phred = (int)(uint8_t)Qual0[0][i] - 33;
            P.inOut->logMain << phred;
            if (i < checkLen1 - 1) P.inOut->logMain << ",";
        }
        P.inOut->logMain << " qual1_min=" << (uint32_t)minQual1 << " qual1_max=" << (uint32_t)maxQual1;
        if (P.readNmates == 2) {
            P.inOut->logMain << " qual2_first10=\"";
            for (uint32_t i = 0; i < checkLen2; i++) {
                P.inOut->logMain << Qual0[1][i];
            }
            P.inOut->logMain << "\" qual2_bytes=";
            for (uint32_t i = 0; i < checkLen2; i++) {
                P.inOut->logMain << (uint32_t)(uint8_t)Qual0[1][i];
                if (i < checkLen2 - 1) P.inOut->logMain << ",";
            }
            P.inOut->logMain << " qual2_phred=";
            for (uint32_t i = 0; i < checkLen2; i++) {
                int phred = (int)(uint8_t)Qual0[1][i] - 33;
                P.inOut->logMain << phred;
                if (i < checkLen2 - 1) P.inOut->logMain << ",";
            }
            P.inOut->logMain << " qual2_min=" << (uint32_t)minQual2 << " qual2_max=" << (uint32_t)maxQual2;
        }
        P.inOut->logMain << endl;
    }
    
    // Cutadapt-style trimming (if enabled)
    if (P.trimCutadapt == "Yes") {
        if (P.readNmates == 2) {
        struct TrimParams params;
        trim_params_init(&params);
        params.quality_cutoff = P.trimCutadaptQuality;
        params.min_length = P.trimCutadaptMinLength;
        // Set compatibility mode
        if (P.trimCutadaptCompat == "Cutadapt3") {
            params.compat_mode = TRIM_COMPAT_CUTADAPT3;
        } else {
            params.compat_mode = TRIM_COMPAT_OFF;  // Default: "-" or "Off"
        }
        // Set adapters if custom (validate exactly 2 adapters)
        if (P.trimCutadaptAdapter.size() == 2 && P.trimCutadaptAdapter[0] != "-" && P.trimCutadaptAdapter[1] != "-") {
            params.adapter_r1 = P.trimCutadaptAdapter[0].c_str();
            params.adapter_r2 = P.trimCutadaptAdapter[1].c_str();
        } else if (P.trimCutadaptAdapter.size() > 0 && (P.trimCutadaptAdapter[0] != "-" || (P.trimCutadaptAdapter.size() > 1 && P.trimCutadaptAdapter[1] != "-"))) {
            // Invalid adapter specification
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: --trimCutadaptAdapter requires exactly 2 adapter sequences (R1 and R2), separated by space\n";
            errOut << "Provided: " << P.trimCutadaptAdapter.size() << " adapter(s)\n";
            errOut << "SOLUTION: Provide both R1 and R2 adapters, or use '-' to use default TruSeq adapters\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        struct TrimResult result1, result2;
        uint32_t len1 = (uint32_t)readLength[0];
        uint32_t len2 = (uint32_t)readLength[1];
        uint32_t origLen1 = len1;  // Save original lengths for debug logging
        uint32_t origLen2 = len2;
        trim_pair(Read0[0], Qual0[0], &len1,
                  Read0[1], Qual0[1], &len2,
                  &params, &result1, &result2);
        readLength[0] = len1;
        readLength[1] = len2;
        
        // Accumulate stats using centralized helper (counting reads, not pairs)
        struct TrimStats trimStats = {0, 0, 0, 0, 0};
        trim_stats_add(&trimStats, &result1);  // R1
        trim_stats_add(&trimStats, &result2);  // R2
        
        // Copy accumulated stats to STAR's Stats class
        statsRA.trimReadsProcessed += trimStats.reads_processed;
        statsRA.trimReadsTrimmed += trimStats.reads_trimmed;
        statsRA.trimReadsTooShort += trimStats.reads_too_short;
        statsRA.trimBasesQualityTrimmed += trimStats.bases_quality_trimmed;
        statsRA.trimBasesAdapterTrimmed += trimStats.bases_adapter_trimmed;
        
        // Pair-level counters (for direct comparison with Trim Galore)
        statsRA.trimPairsProcessed++;
        
            // Debug logging (guarded by STAR_TRIM_DEBUG_N env var)
        // Use same counter value from start of function (already incremented)
        if (debugMax > 0 && debugCount < (uint64_t)debugMax) {
            P.inOut->logMain << "TRIM_DEBUG: " << readNameMates[0]
                             << " params{quality_cutoff=" << static_cast<int>(params.quality_cutoff)
                             << " min_length=" << params.min_length << "}"
                             << " origLen=" << origLen1 << "," << origLen2
                             << " postLen=" << len1 << "," << len2
                             << " R1{dropped=" << result1.dropped
                             << " new_length=" << result1.new_length
                             << " qual3p=" << result1.qual_trimmed_3p 
                             << " qual5p=" << result1.qual_trimmed_5p 
                             << " adapter=" << result1.adapter_trimmed << "}"
                             << " R2{dropped=" << result2.dropped
                             << " new_length=" << result2.new_length
                             << " qual3p=" << result2.qual_trimmed_3p 
                             << " qual5p=" << result2.qual_trimmed_5p 
                             << " adapter=" << result2.adapter_trimmed << "}"
                             << endl;
            // Counter already incremented at start of function
        }
        
        // Handle dropped reads
        if (result1.dropped || result2.dropped) {
            readFilter = 'Y';  // Fail QC - same as other filters
            statsRA.trimPairsDropped++;
            
            // Debug log dropped pairs (use same counter value from start of function)
            if (debugMax > 0 && debugCount < (uint64_t)debugMax) {
                P.inOut->logMain << "TRIM_DROP: " << readNameMates[0]
                                 << " reason=minLength postLen=" << len1 << "," << len2
                                 << " minLen=" << params.min_length << endl;
            }
            
            // Skip mapping for this pair (readN already incremented above)
            completeSpatialWithoutFeature();
            return 0;
        }
        
        // Pair kept after trimming
        statsRA.trimPairsKept++;
        
        // Update original lengths after trimming
        readLengthOriginal[0] = readLength[0];
        readLengthOriginal[1] = readLength[1];
        
        // Re-convert to numeric after trimming
        convertNucleotidesToNumbers(Read0[0], Read1[0], readLength[0]);
        convertNucleotidesToNumbers(Read0[1], Read1[1], readLength[1]);
        } else if (P.readNmates == 1) {
            // Single-end trimming
            struct TrimParams params;
            trim_params_init(&params);
            params.quality_cutoff = P.trimCutadaptQuality;
            params.min_length = P.trimCutadaptMinLength;
            // Use R1 adapter for single-end
            const char* adapter = TRUSEQ_ADAPTER_R1;
            if (P.trimCutadaptAdapter.size() >= 1 && P.trimCutadaptAdapter[0] != "-") {
                adapter = P.trimCutadaptAdapter[0].c_str();
            }
            
            struct TrimResult result;
            uint32_t len1 = (uint32_t)readLength[0];
            result = trim_read(Read0[0], Qual0[0], len1, adapter, &params);
            readLength[0] = result.new_length;
            
            // Accumulate stats using centralized helper
            struct TrimStats trimStats = {0, 0, 0, 0, 0};
            trim_stats_add(&trimStats, &result);
            
            // Copy accumulated stats to STAR's Stats class
            statsRA.trimReadsProcessed += trimStats.reads_processed;
            statsRA.trimReadsTrimmed += trimStats.reads_trimmed;
            statsRA.trimReadsTooShort += trimStats.reads_too_short;
            statsRA.trimBasesQualityTrimmed += trimStats.bases_quality_trimmed;
            statsRA.trimBasesAdapterTrimmed += trimStats.bases_adapter_trimmed;
            
            // Handle dropped reads
            if (result.dropped) {
                readFilter = 'Y';
                completeSpatialWithoutFeature();
                return 0;
            }
            
            // Update length and re-convert
            readLengthOriginal[0] = readLength[0];
            convertNucleotidesToNumbers(Read0[0], Read1[0], readLength[0]);
        }
    }
    
    if (P.outFilterBySJoutStage != 2) {
        for (uint32 im=0; im<P.readNmates; im++) {//not readNends: the barcode quality will be calculated separately
            for (uint64 ix=clipMates[im][0].clippedN; ix<readLengthOriginal[im]-clipMates[im][1].clippedN; ix++) {
                qualHist[im][(uint8)Qual0[im][ix]]++;
            };
        };
    };
    
    if (P.readNmates==2) {//combine two mates together
        Lread=readLength[0]+readLength[1]+1;
        readLengthPairOriginal=readLengthOriginal[0]+readLengthOriginal[1]+1;
        if (Lread>DEF_readSeqLengthMax) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: Lread of the pair = " << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<endl;
            errOut << "Read Name="<<readNameMates[0]<<endl;
            errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR"<<endl<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        //marker for spacer base
        Read1[0][readLength[0]]=MARK_FRAG_SPACER_BASE;
        
        //copy 2nd mate into Read1[0] & reverse-complement
        complementSeqNumbers(Read1[1],Read1[0]+readLength[0]+1,readLength[1]);//complement. Here Read1[1] is still the 2nd mate's numeric-sequence. Later Read1[1] will be reverse complement of the combined read.
        for (uint ii=0;ii<readLength[1]/2;ii++) {
            swap(Read1[0][Lread-ii-1],Read1[0][ii+readLength[0]+1]); //reverse
        };

    } else {//1 mate

        if (readStatus0==-1) {//finished with the stream
            return -1;
        };

        Lread=readLength[0];
        readLengthPairOriginal=readLengthOriginal[0];
        readLength[1]=0;

    };
      
    readFileType=readStatus0;

    hashScreenDecision_ = FlexHashScreenDecision();
    if (P.pSolo.hashScreenEnabled && soloRead != nullptr && soloRead->readBar != nullptr &&
        soloRead->readFeat != nullptr && P.pSolo.featureYes[SoloFeatureTypes::Gene] &&
        P.readNmates > 0) {
        const bool spatialFlex = P.soloSpatialFlexIntegratedEnabled;
        uint16_t hashScreenSampleIdx = 0;
        // Read0[0] is still ASCII-encoded at this point (numeric conversion
        // happens later at complementSeqNumbers). The hash screen encodes
        // A/C/G/T characters; moving this call after convertNucleotidesToNumbers
        // would silently break classification.
        if (spatialFlex) {
            hashScreenDecision_ =
                FlexHashScreenCache::instance().classifyReadH0H1Offset0(
                    Read0[0], readLengthOriginal[0]);
        } else {
            soloRead->readBar->getCBandUMI(
                Read0, Qual0, readLengthOriginal, readNameExtra[0],
                readFilesIndex, readName);
            const auto sampleDetectStart = std::chrono::steady_clock::now();
            detectSampleFromRawR2();
            const auto sampleDetectEnd = std::chrono::steady_clock::now();
            statsRA.sampleDetectPreAlignCalls++;
            statsRA.sampleDetectPreAlignNs += static_cast<uint64>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    sampleDetectEnd - sampleDetectStart).count());
            soloRead->readBar->detectedSampleToken = detectedSampleByte_;
            hashScreenSampleIdx =
                SampleDetector::sampleIndexForToken(detectedSampleByte_);
            if (sampleDetReady_) {
                hashScreenDecision_ = FlexHashScreenCache::instance().classifyRead(
                    Read0[0], readLengthOriginal[0], hashScreenSampleIdx);
            } else {
                // Single-sample/no-tag Flex libraries have no runtime sample
                // index. Reuse the sample-free H0/H1 classifier.
                hashScreenDecision_ =
                    FlexHashScreenCache::instance().classifyReadH0H1Offset0(
                        Read0[0], readLengthOriginal[0]);
            }
        }
        hashScreenDumpWrite(Read0[0], readLengthOriginal[0], hashScreenSampleIdx, hashScreenDecision_);
        if (hashScreenDecision_.action == FlexHashScreenDecision::Keep) {
            if (spatialFlex) {
                if (iReadAll == 0 || hashScreenDecision_.geneIdx15 == 0
                    || (hashScreenDecision_.cacheClass != 0
                        && hashScreenDecision_.cacheClass != 1)) {
                    exitWithError(
                        "EXITING because native spatial Flex received an invalid "
                        "H0/H1 cache keep decision\n",
                        std::cerr, P.inOut->logMain,
                        EXIT_CODE_INCONSISTENT_DATA, P);
                }
                std::string spatialError;
                const spatial_gex::FeatureEvidenceClass source =
                    hashScreenDecision_.cacheClass == 0
                        ? spatial_gex::FeatureEvidenceClass::FlexH0
                        : spatial_gex::FeatureEvidenceClass::FlexH1;
                if (!P.spatialGexPipeline->completeCurrentThread(
                        source, true, hashScreenDecision_.geneIdx15 - 1,
                        iReadAll - 1, spatialError)) {
                    exitWithError(
                        "EXITING because native spatial Flex cache evidence "
                        "completion failed: " + spatialError + "\n",
                        std::cerr, P.inOut->logMain,
                        EXIT_CODE_INCONSISTENT_DATA, P);
                }
                ++statsRA.hashScreenKeep;
                return 0;
            }
            const bool noBarcode = (soloRead->readBar->cbMatch < 0);
            soloRead->readFlagReset();
            SoloReadFeature *geneFeat = soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
            bool handled = record_flex_hash_screen_keep(geneFeat, *soloRead->readBar, iReadAll,
                                                        hashScreenDecision_.geneIdx15,
                                                        hashScreenDecision_.cacheClass,
                                                        hashScreenDecision_.probeRegion);
            if (handled) {
                statsRA.hashScreenKeep++;
                if (noBarcode) {
                    statsRA.hashScreenKeepNoBarcode++;
                }
                return 0;
            }
            statsRA.hashScreenPass++;
        } else if (hashScreenDecision_.action == FlexHashScreenDecision::Deny) {
            statsRA.hashScreenDeny++;
            if (spatialFlex) {
                if (iReadAll == 0) {
                    exitWithError(
                        "EXITING because native spatial Flex lost the hash-deny "
                        "read ordinal\n",
                        std::cerr, P.inOut->logMain,
                        EXIT_CODE_INCONSISTENT_DATA, P);
                }
                std::string spatialError;
                if (!P.spatialGexPipeline->completeCurrentThread(
                        spatial_gex::FeatureEvidenceClass::FlexHashDeny,
                        false, 0, iReadAll - 1, spatialError)) {
                    exitWithError(
                        "EXITING because native spatial Flex hash-deny completion "
                        "failed: " + spatialError + "\n",
                        std::cerr, P.inOut->logMain,
                        EXIT_CODE_INCONSISTENT_DATA, P);
                }
                return 0;
            }
            soloRead->readFlagReset();
            SoloReadFeature *geneFeat = soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
            record_flex_hash_screen_deny(geneFeat, *soloRead->readBar, iReadAll, "NEG_PROBE_AMBIG");
            return 0;
        } else {
            statsRA.hashScreenPass++;
        }
    } else {
        // Read skipped the hash screen entirely (disabled, missing soloRead,
        // wrong feature config, etc.). Dump with action=Disabled so the replay
        // tool can count and verify these reads are correctly routed.
        // sampleIdx=0 because sample detection did not run.
        hashScreenDumpWrite(Read0[0], readLengthOriginal[0], 0, hashScreenDecision_);
    }

    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
    };

    //max number of mismatches allowed for this read
    outFilterMismatchNmaxTotal=min(P.outFilterMismatchNmax, (uint) (P.outFilterMismatchNoverReadLmax*(readLength[0]+readLength[1])));

    //map the read
    const auto alignCoreStart = std::chrono::steady_clock::now();
    if (P.pGe.gType==101) {//SpliceGraph
        mapOneReadSpliceGraph();
    } else {//all other cases - standard alignment algorithm
        mapOneRead();
    };
    const auto alignCoreEnd = std::chrono::steady_clock::now();
    statsRA.alignCoreCalls++;
    statsRA.alignCoreNs += static_cast<uint64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(alignCoreEnd - alignCoreStart).count());

    peOverlapMergeMap();
    
    multMapSelect();
    
    mappedFilter();  
    
    transformGenome();//for now genome transformation happens after multimapper selection, and mapping filter

    if (!peOv.yes) {//if the alignment was not mates merged - otherwise the chimeric detection was already done
        chimericDetection();
    };

    if (P.pCh.out.bam && chimRecord) {//chimeric alignment was recorded in main BAM files, and it contains the representative portion, so non-chimeric aligmnent is not output
        completeSpatialWithoutFeature();
        return 0;
    };

    waspMap();

    #ifdef OFF_BEFORE_OUTPUT
        #warning OFF_BEFORE_OUTPUT
        completeSpatialWithoutFeature();
        return 0;
    #endif

    //write out alignments
    outputAlignments();

    {
    #ifdef DEBUG_OutputLastRead
        lastReadStream.seekp(ios::beg);
        lastReadStream << iReadAll <<" "<< readName <<endl;
    #endif
    };

    return 0;

};

int ReadAlign::oneReadFromCbqView(const star::input::CbqReadView& view) {
    const int readStatus = loadCbqReadView(view);
    return oneReadLoaded(readStatus);
}

int ReadAlign::oneReadFromPacket(EnrichedPacket &pkt) {
    // Copy read data from packet into ReadAlign buffers
    std::strncpy(readName, pkt.name, DEF_readNameLengthMax - 1);
    readName[DEF_readNameLengthMax - 1] = '\0';
    readNameMates[0][0] = '\0';
    readNameMates[1][0] = '\0';

    const uint64_t pipelineSourceOrdinal = pkt.iReadAll;
    // The fused pipeline uses a zero-based logical ordinal. The ordinary
    // ReadAlign/output path uses one-based iReadAll and subtracts one when it
    // completes integrated spatial evidence.
    iReadAll = P.soloSpatialFlexIntegratedEnabled
        ? pipelineSourceOrdinal + 1
        : pipelineSourceOrdinal;
    readFilesIndex = pkt.readFilesIndex;
    readFilter = pkt.readFilter;

    // Copy R2 (mate 0) and R1 (mate 1) sequences and qualities
    std::memcpy(Read0[0], pkt.seq[0], pkt.readLen[0]);
    Read0[0][pkt.readLen[0]] = '\0';
    std::memcpy(Qual0[0], pkt.qual[0], pkt.readLen[0]);
    Qual0[0][pkt.readLen[0]] = '\0';
    readLength[0] = pkt.readLen[0];
    readLengthOriginal[0] = pkt.readLen[0];

    std::memcpy(Read0[1], pkt.seq[1], pkt.readLen[1]);
    Read0[1][pkt.readLen[1]] = '\0';
    std::memcpy(Qual0[1], pkt.qual[1], pkt.readLen[1]);
    Qual0[1][pkt.readLen[1]] = '\0';
    readLength[1] = pkt.readLen[1];
    readLengthOriginal[1] = pkt.readLen[1];

    if (P.soloSpatialFlexIntegratedEnabled) {
        if (P.spatialGexPipeline == nullptr) {
            exitWithError(
                "EXITING because native spatial Flex pipeline lost its accumulator\n",
                std::cerr, P.inOut->logMain,
                EXIT_CODE_INCONSISTENT_DATA, P);
        }
        std::string spatialError;
        if (!P.spatialGexPipeline->decodeCurrentThread(
                Read0[1], readLengthOriginal[1], Qual0[1],
                readLengthOriginal[1], pipelineSourceOrdinal, spatialError)) {
            exitWithError(
                "EXITING because native spatial Flex pipeline raw-R1 decoding failed: "
                    + spatialError + "\n",
                std::cerr, P.inOut->logMain,
                EXIT_CODE_INCONSISTENT_DATA, P);
        }
    }

    statsRA.readN++;
    statsRA.readBases += readLength[0] + readLength[1];

    readFileType = 2; // FASTQ

    // Restore CB/UMI state from triage's pre-extracted fields
    if (soloRead != nullptr && soloRead->readBar != nullptr) {
        soloRead->readBar->cbMatch = pkt.cbMatch;
        soloRead->readBar->cbMatchInd.resize(pkt.cbMatchIndN);
        for (uint32_t i = 0; i < pkt.cbMatchIndN; ++i)
            soloRead->readBar->cbMatchInd[i] = pkt.cbMatchInd[i];
        soloRead->readBar->umiB = pkt.umiB;
        soloRead->readBar->detectedSampleToken = pkt.detectedSampleToken;
    }

    // Convert ASCII sequences to numeric encoding (must happen before hash screen
    // since classifyRead operates on ASCII Read0, but alignment needs numeric Read1)
    convertNucleotidesToNumbers(Read0[0], Read1[0], readLength[0]);
    convertNucleotidesToNumbers(Read0[1], Read1[1], readLength[1]);

    // Apply 5' and 3' clipping on numeric sequence
    clipMates[0][0].clip(readLength[0], Read1[0]);
    clipMates[0][1].clip(readLength[0], Read1[0]);
    clipMates[1][0].clip(readLength[1], Read1[1]);
    clipMates[1][1].clip(readLength[1], Read1[1]);

    // H0+H1 hash screen is now fully handled pre-queue by the fused lane thread.
    // Alignment workers receive only true MISS reads — proceed directly to alignment.
    statsRA.hashScreenPass++;

    // Combine PE mates with spacer base (mirrors oneRead path)
    if (P.readNmates == 2) {
        Lread = readLength[0] + readLength[1] + 1;
        readLengthPairOriginal = readLengthOriginal[0] + readLengthOriginal[1] + 1;
        Read1[0][readLength[0]] = MARK_FRAG_SPACER_BASE;
        complementSeqNumbers(Read1[1], Read1[0] + readLength[0] + 1, readLength[1]);
        for (uint ii = 0; ii < readLength[1] / 2; ii++) {
            swap(Read1[0][Lread - ii - 1], Read1[0][ii + readLength[0] + 1]);
        }
    } else {
        Lread = readLength[0];
        readLengthPairOriginal = readLengthOriginal[0];
        readLength[1] = 0;
    }

    // Full-read complement and reverse complement for alignment
    complementSeqNumbers(Read1[0], Read1[1], Lread);
    for (uint ii = 0; ii < Lread; ii++) {
        Read1[2][Lread - ii - 1] = Read1[1][ii];
    }

    outFilterMismatchNmaxTotal = min(P.outFilterMismatchNmax,
        (uint)(P.outFilterMismatchNoverReadLmax * (readLength[0] + readLength[1])));

    const auto alignCoreStart = std::chrono::steady_clock::now();
    if (P.pGe.gType == 101) {
        mapOneReadSpliceGraph();
    } else {
        mapOneRead();
    }
    const auto alignCoreEnd = std::chrono::steady_clock::now();
    statsRA.alignCoreCalls++;
    statsRA.alignCoreNs += static_cast<uint64>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(alignCoreEnd - alignCoreStart).count());

    peOverlapMergeMap();
    multMapSelect();
    mappedFilter();
    transformGenome();

    if (!peOv.yes) {
        chimericDetection();
    }

    if (P.pCh.out.bam && chimRecord) {
        return 0;
    }

    waspMap();

    outputAlignments();

    return 0;
}
