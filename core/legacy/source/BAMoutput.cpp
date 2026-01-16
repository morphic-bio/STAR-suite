#include "BAMoutput.h"
#include "SampleDetector.h"
#include "UmiCodec.h"
#include <sys/stat.h>
#include "GlobalVariables.h"
#include "SamtoolsSorter.h"
#include <pthread.h>
#include "serviceFuns.cpp"
#include "ThreadControl.h"
#include "streamFuns.h"
#include "SoloBamParsing.h"
#include "ProbeListIndex.h"
#include "PackedReadInfo.h"
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <algorithm>
#include <atomic>
#include "FlexDebugCounters.h"

// Debug flag for Task 2 staging verification
static const bool g_debugTask2 = (std::getenv("STAR_DEBUG_TASK2") != nullptr);
// Debug flag for sample detection (STAR_DEBUG_SAMPLE)
static const bool g_debugSample = (std::getenv("STAR_DEBUG_SAMPLE") != nullptr);
// Keep debugSampleDetectionCount as atomic since it's used for conditional debug limiting (not a hot path)
static std::atomic<int> g_debugSampleDetectionCount{0};
// recordsWithoutSample and ambiguousCbCount now use thread-local FlexDebugCounters

namespace {

std::string normalizeCbTag(const std::string &cbRaw) {
    size_t dash = cbRaw.find('-');
    if (dash != std::string::npos) {
        return cbRaw.substr(0, dash);
    }
    return cbRaw;
}

uint32_t lookupCbIndex(const std::string &cbRaw, const Parameters &P) {
    if (cbRaw.empty() || P.pSolo.cbWLstr.empty()) {
        return 0u;
    }
    std::string cb = normalizeCbTag(cbRaw);
    if (cb.empty()) {
        return 0u;
    }
    const auto &wl = P.pSolo.cbWLstr;
    auto it = std::lower_bound(wl.begin(), wl.end(), cb);
    if (it != wl.end() && *it == cb) {
        return static_cast<uint32_t>(std::distance(wl.begin(), it) + 1);
    }
    return 0u;
}

} // namespace

BAMoutput::BAMoutput (int iChunk, string tmpDir, Parameters &Pin) : P(Pin){//allocate bam array

    nBins=P.outBAMcoordNbins;
    binSize=P.chunkOutBAMsizeBytes/nBins;
    bamArraySize=binSize*nBins;
    bamArray = new char [bamArraySize];

    bamDir=tmpDir+to_string((uint) iChunk);//local directory for this thread (iChunk)

    mkdir(bamDir.c_str(),P.runDirPerm);
    binStart=new char* [nBins];
    binBytes=new uint64 [nBins];
    binStream=new ofstream* [nBins];
    binTotalN=new uint [nBins];
    binTotalBytes=new uint [nBins];
    for (uint ii=0;ii<nBins;ii++) {
        binStart[ii]=bamArray+bamArraySize/nBins*ii;
        binBytes[ii]=0;
        binStream[ii]=&ofstrOpen((bamDir +"/"+to_string(ii)).c_str(), ERROR_OUT, P);    //open temporary files
        binTotalN[ii]=0;
        binTotalBytes[ii]=0;
    };

    binSize1=binStart[nBins-1]-binStart[0];
    nBins=1;//start with one bin to estimate genomic bin sizes

    sampleDet_ = nullptr;
    sampleDetReady_ = false;
    ambiguousCbSpilloverEnabled_ = false;
    
    if (!P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
        !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
        sampleDet_ = new SampleDetector(P.pSolo);
        if (sampleDet_->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
            sampleDet_->loadProbes(P.pSolo.sampleProbesPath)) {
            sampleDetReady_ = true;
        } else {
            delete sampleDet_;
            sampleDet_ = nullptr;
        }
    }
    
    // Note: ambiguousCbSpilloverFile_ will be opened in setKeyAggregator() when aggregator is set
};

BAMoutput::BAMoutput (BGZF *bgzfBAMin, Parameters &Pin) : P(Pin){//allocate BAM array with one bin, streamed directly into bgzf file

    bamArraySize=P.chunkOutBAMsizeBytes;
    bamArray = new char [bamArraySize];
    binBytes1=0;
    bgzfBAM=bgzfBAMin;
    //not used
    binSize=0;
    binStream=NULL;
    binStart=NULL;
    binBytes=NULL;
    binTotalBytes=NULL;
    binTotalN=NULL;
    nBins=0;
    
    // Initialize staging queues for flush-synchronized ledger
    pendingReadIds_.clear();
    pendingAux_.clear();
    
    // Task 2: Initialize solo meta staging
    pendingSoloMeta_.clear();
    lastQname_.clear();
    sampleDet_ = nullptr;
    sampleDetReady_ = false;
    
    // Initialize Y-chromosome split handles
    if (P.emitNoYBAMyes) {
        bgzfBAM_Y = P.inOut->outBAMfileY;
        bgzfBAM_noY = P.inOut->outBAMfileNoY;
        suppressPrimary_ = !P.keepBAMyes;  // suppress primary unless --keepBAM
    } else {
        bgzfBAM_Y = nullptr;
        bgzfBAM_noY = nullptr;
        suppressPrimary_ = false;
    }
    
    if (!P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
        !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
        sampleDet_ = new SampleDetector(P.pSolo);
        if (sampleDet_->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
            sampleDet_->loadProbes(P.pSolo.sampleProbesPath)) {
            sampleDetReady_ = true;
            if (g_debugSample) {
                P.inOut->logMain << "BAMoutput: SampleDetector initialized successfully (whitelist=" 
                                 << P.pSolo.sampleWhitelistPath << ", probes=" 
                                 << P.pSolo.sampleProbesPath << ")" << std::endl;
            }
        } else {
            P.inOut->logMain << "WARNING: BAMoutput SampleDetector initialization failed (whitelist=" 
                             << P.pSolo.sampleWhitelistPath << ", probes=" 
                             << P.pSolo.sampleProbesPath << ")" << std::endl;
            delete sampleDet_;
            sampleDet_ = nullptr;
        }
    } else {
        if (g_debugSample) {
            P.inOut->logMain << "BAMoutput: SampleDetector not initialized (whitelist=" 
                             << (P.pSolo.sampleWhitelistPath.empty() ? "empty" : P.pSolo.sampleWhitelistPath)
                             << ", probes=" 
                             << (P.pSolo.sampleProbesPath.empty() ? "empty" : P.pSolo.sampleProbesPath) << ")" << std::endl;
        }
    }
    
};

BAMoutput::~BAMoutput() {
    if (sampleDet_ != nullptr) {
        delete sampleDet_;
        sampleDet_ = nullptr;
    }
    
    // Close ambiguous CB spillover file if open
    if (ambiguousCbSpilloverFile_.is_open()) {
#ifdef DEBUG_FLEX_COUNTERS
        uint64_t ambigTotal = flexCountersSumAll().ambiguousCbCount;
        if (ambigTotal > 0) {
            P.inOut->logMain << "[CB-RESOLVER] Final ambiguous CB count: " << ambigTotal 
                             << " (written to spillover file)" << std::endl;
        }
#endif
        ambiguousCbSpilloverFile_.close();
    }
}

void BAMoutput::unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2, uint64_t iReadAll, uint8_t sampleByte,
                                  uint32_t cbIdxPlus1, uint32_t umi24, const std::string &cbSeq, bool hasY) {
    if (bamSize==0) return; //no output, could happen if one of the mates is not mapped
    
    // Skip output during auto-trim detection pass
    if (P.quant.slam.autoTrimDetectionPass) return;

    if (binBytes1+bamSize2 > bamArraySize) {//write out this buffer
        flushPendingToLedgerAndDisk();
    }

    // Stage metadata before copying BAM data
    pendingReadIds_.push_back(static_cast<uint32_t>(iReadAll));
    
    // Task 2/3: Stage solo metadata for direct-mode keys emission
    PendingSoloMeta meta;
    meta.iReadIdx = static_cast<uint32_t>(iReadAll);
    
    // sampleByte convention: mate 0 (R2/sample read) gets detected token or 0u, mate 1 (R1/genomic) gets 0xFFu
    // allowSampleDetection = true means this is the sample read (mate 0)
    bool allowSampleDetection = (sampleByte != 0xFFu);
    uint8_t sampleToken = allowSampleDetection ? (sampleByte & 0x1F) : 0u;
    uint32_t detectedIdx = 0u;
    
    // If sampleByte already contains a detected token (non-zero), use it directly
    // Otherwise, try to detect from BAM as fallback (for cases where upstream detection wasn't done)
    if (allowSampleDetection) {
        if (sampleToken != 0) {
            // Use token passed from upstream detection (ReadAlign detected from raw sequence)
            // Resolve sequential index from token LUT
            detectedIdx = SampleDetector::sampleIndexForToken(sampleToken);
        } else if (sampleDetReady_) {
            // Fallback: try to detect from BAM record (for two-pass mode or edge cases)
            detectedIdx = sampleDet_->detectSampleIndexFromBam(bamIn, bamSize);
            if (detectedIdx > 0) {
                // Update token from detected index (token = detectedIdx & 0x1F)
                // Token registration happens automatically inside detectSampleIndex()
                sampleToken = static_cast<uint8_t>(detectedIdx & 0x1F);
            }
        }
        
        // Debug instrumentation: log first few detections
        if (g_debugSample) {
            int count = g_debugSampleDetectionCount.fetch_add(1);
            if (count < 10) {
                P.inOut->logMain << "BAMoutput::unsortedOneAlign: sampleByte=" << (unsigned int)sampleByte
                                 << ", allowSampleDetection=" << allowSampleDetection
                                 << ", sampleDetReady_=" << sampleDetReady_
                                 << ", detectedIdx=" << detectedIdx 
                                 << ", sampleToken=" << (unsigned int)sampleToken << std::endl;
            }
        }
    } else if (g_debugSample && allowSampleDetection) {
        // Log when detection should run but doesn't
        static int skipCount = 0;
        if (skipCount++ < 5) {
            P.inOut->logMain << "BAMoutput::unsortedOneAlign: SKIP detection - allowSampleDetection=" 
                             << allowSampleDetection << ", sampleDetReady_=" << sampleDetReady_ << std::endl;
        }
    }
    
    // Use detected index directly; fallback to LUT only if detection didn't run or returned 0
    uint32_t sampleIdxFull = detectedIdx;
    if (sampleIdxFull == 0 && sampleToken != 0) {
        // Fallback: try to resolve from token LUT (useful for two-pass replay)
        sampleIdxFull = SampleDetector::sampleIndexForToken(sampleToken);
    }
    if (sampleIdxFull > std::numeric_limits<uint16_t>::max()) {
        sampleIdxFull = std::numeric_limits<uint16_t>::max();
    }
    meta.sampleIdx = static_cast<uint16_t>(sampleIdxFull);
    bool hasSample = (meta.sampleIdx != 0u);
    
    // Track records without sample indices for warning
    if (!hasSample && allowSampleDetection) {
        FLEX_COUNT_INC(recordsWithoutSample);
    }
    
    // Debug instrumentation: log staging result
    if (g_debugSample && g_debugSampleDetectionCount <= 5 && allowSampleDetection) {
        P.inOut->logMain << "BAMoutput::unsortedOneAlign: staged meta.sampleIdx=" << meta.sampleIdx 
                         << ", hasSample=" << hasSample << std::endl;
    }
    
    // Use extracted UMI from upstream detection if available, otherwise fall back to parsing BAM tags
    if (umi24 != 0) {
        // Use upstream-extracted UMI (from raw FASTQ sequence)
        meta.umi24 = umi24;
        meta.urValid = true; // Assume valid if extracted upstream
    } else {
        // Fallback: parse UR tag from BAM record
        std::string urStr = SoloBamParsing::parseUR(bamIn, bamSize);
        uint32_t urPacked = UINT32_MAX;
        bool urValid = false;
        if (!urStr.empty() && urStr.length() == 12) {
            // Check for N characters (invalid UR)
            bool hasN = (urStr.find('N') != std::string::npos || urStr.find('n') != std::string::npos);
            if (!hasN) {
                urPacked = encodeUMI12(urStr);
                urValid = (urPacked != UINT32_MAX);
            } else {
                urValid = false; // UR contains N, invalid
            }
        }
        meta.umi24 = urValid ? (urPacked & 0xFFFFFFu) : 0u;
        meta.urValid = urValid;
    }
    
    // Build auxWord: [drop(1) | tag_valid(1) | ur_valid(1) | tag_code(5) | packed_ur(24)]
    // drop flag will be set later during cells_allow/assignment filtering
    uint32_t auxWord = 0;
    uint8_t tagCode = hasSample ? static_cast<uint8_t>(sampleToken & 0x1Fu) : 0u;
    bool urValid = (meta.umi24 != 0u);
    if (hasSample) {
        auxWord |= (1 << 30);  // tag_valid
    }
    if (urValid) {
        auxWord |= (1 << 29);  // ur_valid
        auxWord |= (meta.umi24 & 0xFFFFFF);  // packed_ur (24 bits)
    }
    auxWord |= (static_cast<uint32_t>(tagCode) << 24);  // tag_code (5 bits)
    // drop flag (bit 31) remains 0 for now, will be set during filtering
    
    pendingAux_.push_back(auxWord);
    
    // Extract QNAME and core fields for drop logic
    std::string qname;
    int mapq = 0, nm = -1, nh = 1, readLen = 0;
    uint16_t flag = 0;
    std::string readBases; 
    SoloBamParsing::parseCoreAndNM(bamIn, bamSize, qname, mapq, nm, nh, readLen, flag, readBases);
    
    // Use extracted CB whitelist index from upstream detection if available, otherwise fall back to parsing BAM tags
    if (cbIdxPlus1 != 0) {
        // Use upstream-extracted CB index (from raw FASTQ sequence)
        meta.cbIdxPlus1 = cbIdxPlus1;
        meta.cbSeq.clear(); // Not needed for resolved CBs
    } else {
        // Fallback: use provided CB sequence (from ReadAlign) or parse from BAM record
        if (!cbSeq.empty()) {
            // Use CB sequence passed from ReadAlign (for ambiguous CBs)
            meta.cbSeq = cbSeq;
            meta.cbIdxPlus1 = lookupCbIndex(cbSeq, P);
        } else {
            // Parse CB tag from BAM record
            std::string cbTag = SoloBamParsing::parseCB(bamIn, bamSize);
            meta.cbSeq = cbTag; // Store CB sequence for Phase 2 resolution lookup
            meta.cbIdxPlus1 = lookupCbIndex(cbTag, P);
        }
        // If still 0, may be ambiguous - will be resolved by Bayesian resolver
    }
    
    // Set newGroup bit if QNAME changed
    meta.newGroup = (qname != lastQname_) ? 1 : 0;
    if (meta.newGroup) {
        lastQname_ = qname;
    }
    
    // Gene index resolution - set to 0 (no gene info in this path)
    meta.geneIdx15 = 0;
    meta.zgGeneIdx15.clear();
    
    // Task 3: Complete drop flag logic
    uint8_t flags = 0;
    // Bit 0: MAPQ threshold (simplified - use common threshold of 30)
    int mapqThreshold = 30; // Common MAPQ threshold for single-cell
    if (mapq < mapqThreshold) flags |= 0x01;
    
    // Bit 1: NH (multi-mapper check)
    if (nh > 1) flags |= 0x02;
    
    // Bit 2: NM (mismatch rate check - placeholder, full calc needs CIGAR)
    if (nm >= 0 && readLen > 0) {
        double mmRate = (double)nm / readLen;
        if (mmRate > 0.05) flags |= 0x04; // 5% threshold placeholder
    }
    
    // Bit 3: Gene validity
    if (meta.geneIdx15 == 0) flags |= 0x08;
    
    // Bit 4: Sample validity (check hasSample flag, not bit pattern)
    if (!hasSample) flags |= 0x10;
    
    // Bit 5: Reserved for future use
    
    meta.dropFlags = flags;
    meta.hasY = hasY;  // Store Y flag for routing
    
    pendingSoloMeta_.push_back(meta);
    
    memcpy(bamArray+binBytes1, bamIn, bamSize);
    binBytes1 += bamSize;
};

void BAMoutput::flushPendingToLedgerAndDisk() {
    if (g_threadChunks.threadBool) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    
    // Task 2: Debug logging to verify staging counts match
    if (g_debugTask2 && !pendingReadIds_.empty()) {
        size_t ledgerCount = pendingReadIds_.size();
        size_t soloMetaCount = pendingSoloMeta_.size();
        if (ledgerCount != soloMetaCount) {
            std::cerr << "[TASK2_ERROR] Staging count mismatch: ledger=" << ledgerCount 
                      << " soloMeta=" << soloMetaCount << std::endl;
        } else {
            std::cerr << "[TASK2_OK] Flushing " << ledgerCount << " records" << std::endl;
        }
        // Verify newGroup distribution
        size_t groupCount = 0;
        for (const auto& meta : pendingSoloMeta_) {
            if (meta.newGroup) groupCount++;
        }
        std::cerr << "[TASK2_OK] Groups in batch: " << groupCount << std::endl;
    }
    
    // Note: SoloTagLedger append removed - ledger not used in inline flex path
    
    // Write BAM buffer to BGZF stream(s)
    if (P.emitNoYBAMyes && bgzfBAM_Y != nullptr && bgzfBAM_noY != nullptr) {
        // Y-split mode: route records individually based on hasY flag
        uint64_t offset = 0;
        size_t metaIdx = 0;
        while (offset < binBytes1 && metaIdx < pendingSoloMeta_.size()) {
            // Extract BAM record size (first 4 bytes)
            if (offset + 4 > binBytes1) break;
            uint32_t recSize = *((uint32_t*)(bamArray + offset));
            if (offset + recSize + 4 > binBytes1) break;  // Safety check
            
            // Route to appropriate handle based on hasY flag
            BGZF *targetHandle = pendingSoloMeta_[metaIdx].hasY ? bgzfBAM_Y : bgzfBAM_noY;
            bgzf_write(targetHandle, bamArray + offset, recSize + 4);
            
            // Also write to primary if not suppressed
            if (!suppressPrimary_ && bgzfBAM != nullptr) {
                bgzf_write(bgzfBAM, bamArray + offset, recSize + 4);
            }
            
            offset += recSize + 4;
            metaIdx++;
        }
    } else {
        // Normal mode: write entire buffer to primary BAM
        bgzf_write(bgzfBAM, bamArray, binBytes1);
    }
    
    if (g_threadChunks.threadBool) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
    
    // Log ambiguous CB summary periodically (every 100K records or at end)
#ifdef DEBUG_FLEX_COUNTERS
    static uint64_t logCounter = 0;
    logCounter += pendingSoloMeta_.size();
    uint64_t ambigTotal = flexCountersSumAll().ambiguousCbCount;
    if (ambigTotal > 0 && (logCounter % 100000 == 0 || pendingSoloMeta_.empty())) {
        P.inOut->logMain << "[CB-RESOLVER] Ambiguous CBs encountered: " << ambigTotal 
                         << " (spillover file: " << (ambiguousCbSpilloverEnabled_ ? "enabled" : "disabled") << ")" << std::endl;
    }
#endif
    
    // Clear staging and reset buffer position
    pendingReadIds_.clear();
    pendingAux_.clear();
    pendingSoloMeta_.clear(); // Task 2: Clear solo meta staging
    binBytes1 = 0;
}

void BAMoutput::unsortedFlush () {//flush all alignments
    flushPendingToLedgerAndDisk();
    
    // Warn if many records were processed without sample indices
#ifdef DEBUG_FLEX_COUNTERS
    uint64_t noSampleCount = flexCountersSumAll().recordsWithoutSample;
    if (noSampleCount > 1000000 && sampleDetReady_) {
        P.inOut->logMain << "WARNING: " << noSampleCount 
                         << " sample-read records were processed without detecting sample indices. "
                         << "Check that --soloSampleWhitelist and --soloSampleProbes are correctly specified." << std::endl;
    }
#endif
};

void BAMoutput::coordOneAlign (char *bamIn, uint bamSize, uint iRead, bool hasY) {
    // Skip output during auto-trim detection pass
    if (P.quant.slam.autoTrimDetectionPass) return;
    
    // Branch to samtools sorter if enabled
    if (P.outBAMsortMethod == "samtools" && g_samtoolsSorter != nullptr) {
        // Extract readId from iRead (bits[63:32])
        uint32_t readId = static_cast<uint32_t>(iRead >> 32);
        g_samtoolsSorter->addRecord(bamIn, bamSize, readId, hasY);
        return;
    }

    uint32 *bamIn32;
    uint alignG;
    uint32 iBin=0;

    if (bamSize==0) {
        return; //no output, could happen if one of the mates is not mapped
    } else {
        //determine which bin this alignment belongs to
        bamIn32=(uint32*) bamIn;
        alignG=( ((uint) bamIn32[1]) << 32 ) | ( (uint)bamIn32[2] );
        if (bamIn32[1] == ((uint32) -1) ) {//unmapped
            iBin=P.outBAMcoordNbins-1;
        } else if (nBins>1) {//bin starts have already been determined
            iBin=binarySearch1a <uint64> (alignG, P.outBAMsortingBinStart, (int32) (nBins-1));
        };
    };

    //write buffer is filled
    if (binBytes[iBin]+bamSize+sizeof(uint64) > ( (iBin>0 || nBins>1) ? binSize : binSize1) ) {//write out this buffer
        if ( nBins>1 || iBin==(P.outBAMcoordNbins-1) ) {//normal writing, bins have already been determined
            binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
            binBytes[iBin]=0;//rewind the buffer
        } else {//the first chunk of reads was written in one bin, need to determine bin sizes, and re-distribute reads into bins
            coordBins();
            coordOneAlign (bamIn, bamSize, iRead, hasY);//record the current align into the new bins
            return;
        };
    };

    // Encode Y flag in bit 63 of iRead (convert to uint64)
    uint64 iReadWithY = static_cast<uint64>(iRead);
    if (hasY) {
        iReadWithY |= (1ULL << 63);  // set Y bit
    }

    //record this alignment in its bin
    memcpy(binStart[iBin]+binBytes[iBin], bamIn, bamSize);
    binBytes[iBin] += bamSize;
    memcpy(binStart[iBin]+binBytes[iBin], &iReadWithY, sizeof(uint64));
    binBytes[iBin] += sizeof(uint64);
    binTotalBytes[iBin] += bamSize+sizeof(uint64);
    binTotalN[iBin] += 1;
    return;
};

void BAMoutput::coordBins() {//define genomic starts for bins
    nBins=P.outBAMcoordNbins;//this is the true number of bins

    //mutex here
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexBAMsortBins);
    if (P.outBAMsortingBinStart[0]!=0) {//it's set to 0 only after the bin sizes are determined
        //extract coordinates and sort
        uint64 *startPos = new uint64 [binTotalN[0]+1];//array of aligns start positions
        for (uint ib=0,ia=0;ia<binTotalN[0];ia++) {
            uint32 *bamIn32=(uint32*) (binStart[0]+ib);
            startPos[ia]  =( ((uint64) bamIn32[1]) << 32) | ( (uint64)bamIn32[2] );
            ib+=bamIn32[0]+sizeof(uint32)+sizeof(uint64);//note: iRead is now uint64 with Y-bit
        };
        qsort((void*) startPos, binTotalN[0], sizeof(uint64), funCompareUint1);

        //determine genomic starts of the bins
        P.inOut->logMain << "BAM sorting: "<<binTotalN[0]<< " mapped reads\n";
        P.inOut->logMain << "BAM sorting bins genomic start loci:\n";

        P.outBAMsortingBinStart[0]=0;
        for (uint32 ib=1; ib<(nBins-1); ib++) {
            P.outBAMsortingBinStart[ib]=static_cast<uint64>(startPos[binTotalN[0]/(nBins-1)*ib]);
            P.inOut->logMain << ib <<"\t"<< (P.outBAMsortingBinStart[ib]>>32) << "\t" << ((P.outBAMsortingBinStart[ib]<<32)>>32) <<endl;
            //how to deal with equal boundaries???
        };
        delete [] startPos;
    };
    //mutex here
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexBAMsortBins);

    //re-allocate binStart
    uint binTotalNold=binTotalN[0];
    char *binStartOld=new char [binSize1];
    memcpy(binStartOld,binStart[0],binBytes[0]);

    binBytes[0]=0;
    binTotalN[0]=0;
    binTotalBytes[0]=0;

    //re-bin all aligns
    for (uint ib=0,ia=0;ia<binTotalNold;ia++) {
        uint32 *bamIn32=(uint32*) (binStartOld+ib);
        uint ib1=ib+bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        uint64 iReadWithY = *((uint64*) (binStartOld+ib1));  // Read as uint64 (includes Y-bit)
        bool hasY = (iReadWithY >> 63) & 1;  // Extract Y-bit
        uint iRead = static_cast<uint>(iReadWithY & 0x7FFFFFFF);  // Mask out Y-bit for original iRead
        coordOneAlign (binStartOld+ib, (uint) (bamIn32[0]+sizeof(uint32)), iRead, hasY);
        ib=ib1+sizeof(uint64);//iReadWithY at the end of the BAM record (now 8 bytes)
    };
    delete [] binStartOld;
    return;
};

void BAMoutput::coordFlush () {//flush all alignments
    if (nBins==1) {
        coordBins();
    };
    for (uint32 iBin=0; iBin<nBins; iBin++) {
        binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
        binStream[iBin]->flush();
        binBytes[iBin]=0;//rewind the buffer
    };
};

void BAMoutput::coordUnmappedPrepareBySJout () {//flush all alignments
    uint iBin=P.outBAMcoordNbins-1;
    binStream[iBin]->write(binStart[iBin],binBytes[iBin]);
    binStream[iBin]->flush();
    binBytes[iBin]=0;//rewind the buffer
    binStream[iBin]->close();
    binStream[iBin]->open((bamDir +"/"+to_string(iBin)+".BySJout").c_str());
};
