#include "ReadAlign.h"
#include "readLoad.h"
#include "readBarcodeLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#include "libtrim/trim.h"
#include <cstdlib>
#include <atomic>

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
    
    // Increment read counters BEFORE trimming (so dropped reads are counted)
    statsRA.readN++;
    statsRA.readBases += readLength[0] + (P.readNmates == 2 ? readLength[1] : 0);
    
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

        if (readStatus[0]==-1) {//finished with the stream
            return -1;
        };

        Lread=readLength[0];
        readLengthPairOriginal=readLengthOriginal[0];
        readLength[1]=0;

    };
      
    readFileType=readStatus[0];

    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
    };

    //max number of mismatches allowed for this read
    outFilterMismatchNmaxTotal=min(P.outFilterMismatchNmax, (uint) (P.outFilterMismatchNoverReadLmax*(readLength[0]+readLength[1])));

    //map the read
    if (P.pGe.gType==101) {//SpliceGraph
        mapOneReadSpliceGraph();
    } else {//all other cases - standard alignment algorithm
        mapOneRead();
    };

    peOverlapMergeMap();
    
    multMapSelect();
    
    mappedFilter();  
    
    transformGenome();//for now genome transformation happens after multimapper selection, and mapping filter

    if (!peOv.yes) {//if the alignment was not mates merged - otherwise the chimeric detection was already done
        chimericDetection();
    };

    if (P.pCh.out.bam && chimRecord) {//chimeric alignment was recorded in main BAM files, and it contains the representative portion, so non-chimeric aligmnent is not output
        return 0;
    };

    waspMap();

    #ifdef OFF_BEFORE_OUTPUT
        #warning OFF_BEFORE_OUTPUT
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


