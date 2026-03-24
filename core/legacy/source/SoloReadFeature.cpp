#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "SoloBinarySpool.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"

SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk)
             : featureType(feTy), P(Pin), pSolo(P.pSolo), binarySpool(false), binarySpoolInMemory(false), binarySpoolMemoryLimitBytes(0), streamReads(nullptr), binarySpoolReadPos(0), inlineHash_(nullptr), readIdTracker_(nullptr)
{
    if (pSolo.type==0)
        return;
//     if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
//         return;
    
    readInfoYes = pSolo.readInfoYes[featureType];
    readIndexYes = pSolo.readIndexYes[featureType];
    
    if (pSolo.cbWLyes) {
        cbReadCount.resize(pSolo.cbWLsize,0);
    };

    const bool wantBinarySpool = SoloBinarySpool::envEnabled();
    const bool wantBinarySpoolInMemory = SoloBinarySpool::envInMemoryEnabled();
    const uint64_t binarySpoolMemoryLimitBytesEnv = SoloBinarySpool::envInMemoryLimitBytes();
    if (wantBinarySpool && P.runRestart.type == 1) {
        ostringstream errOut;
        errOut << "EXITING because STAR_SOLO_BINARY_SPOOL is not compatible with restart mode in this experiment.\n"
               << "SOLUTION: rerun without restart files or disable STAR_SOLO_BINARY_SPOOL.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    if (pSolo.inlineHashMode) {
        // Initialize inline hash instead of opening temp stream file
        inlineHash_ = kh_init(cg_agg);
        streamReads = nullptr; // Do NOT open stream file in inline hash mode
        
        // Initialize parallel readId tracker for sorted BAM CB/UB tag injection
        if (pSolo.trackReadIdsForTags) {
            readIdTracker_ = kh_init(readid_cbumi);
        }
    } else if (iChunk>=0) {
        binarySpool = wantBinarySpool && SoloBinarySpool::supportsFeature(featureType);
        binarySpoolInMemory = binarySpool && wantBinarySpoolInMemory;
        binarySpoolFileName = P.outFileTmp+"/solo"+SoloFeatureTypes::Names[featureType]+'_'+std::to_string(iChunk);
        binarySpoolMemoryLimitBytes = binarySpoolInMemory ? binarySpoolMemoryLimitBytesEnv : 0;
        if (binarySpoolInMemory) {
            streamReads = nullptr;
            SoloBinarySpool::writeFileHeader(binarySpoolBuffer, featureType, readIndexYes);
            if (iChunk == 0) {
                P.inOut->logMain << "Using experimental in-memory binary Solo temp spool for feature "
                                 << SoloFeatureTypes::Names[featureType] << endl;
                if (binarySpoolMemoryLimitBytes > 0) {
                    P.inOut->logMain << "Experimental in-memory Solo spool limit: "
                                     << (binarySpoolMemoryLimitBytes / (1024ull * 1024ull))
                                     << " MiB per thread; spilling overflow to disk" << endl;
                }
            }
        } else if (binarySpool) {
            streamReads = &fstrOpenBinary(binarySpoolFileName, ERROR_OUT, P, true);
            SoloBinarySpool::writeFileHeader(*streamReads, featureType, readIndexYes);
            if (iChunk == 0) {
                P.inOut->logMain << "Using experimental binary Solo temp spool for feature "
                                 << SoloFeatureTypes::Names[featureType] << endl;
            }
        } else {
            //open with flagDelete=false, i.e. try to keep file if it exists
            streamReads = &fstrOpen(binarySpoolFileName, ERROR_OUT, P, false);
            if (wantBinarySpool && iChunk == 0 && !SoloBinarySpool::supportsFeature(featureType)) {
                P.inOut->logMain << "WARNING: STAR_SOLO_BINARY_SPOOL is enabled, but feature "
                                 << SoloFeatureTypes::Names[featureType]
                                 << " is not supported by the experimental binary spool; using legacy text temp stream"
                                 << endl;
            }
        }
    };
    
    if (featureType==SoloFeatureTypes::Transcript3p)
        transcriptDistCount.resize(10000,0);
};

SoloReadFeature::~SoloReadFeature() {
    if (inlineHash_) {
        kh_destroy(cg_agg, inlineHash_);
        inlineHash_ = nullptr;
    }
    if (readIdTracker_) {
        kh_destroy(readid_cbumi, readIdTracker_);
        readIdTracker_ = nullptr;
    }
}

void SoloReadFeature::addCounts(const SoloReadFeature &rfIn)
{
    if (pSolo.cbWLyes) {//WL
        for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
            cbReadCount[ii] += rfIn.cbReadCount[ii];
        };
    } else {
        for (auto ii=rfIn.cbReadCountMap.cbegin(); ii!=rfIn.cbReadCountMap.cend(); ++ii) {
            cbReadCountMap[ii->first] += ii->second;
        };
    };
    
    if (transcriptDistCount.size()>0) {
        for (uint32 ii=0; ii<transcriptDistCount.size(); ii++)
            transcriptDistCount[ii] += rfIn.transcriptDistCount[ii];
    };
};

void SoloReadFeature::addStats(const SoloReadFeature &rfIn)
{
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];

    for (uint32 ii=0; ii<readFlag.nBits; ii++)
        readFlag.flagCountsNoCB[ii] += rfIn.readFlag.flagCountsNoCB[ii];
};

void SoloReadFeature::statsOut(ofstream &streamOut)
{
    //streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(50) << stats.names[ii] << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

void SoloReadFeature::mergeInlineHash(SoloReadFeature &other)
{
    if (!inlineHash_ || !other.inlineHash_) {
        return;
    }
    
    // Merge hash tables: iterate over source hash, add counts
    for (khiter_t iter = kh_begin(other.inlineHash_); iter != kh_end(other.inlineHash_); ++iter) {
        if (!kh_exist(other.inlineHash_, iter)) continue;
        
        uint64_t key = kh_key(other.inlineHash_, iter);
        uint32_t count = kh_val(other.inlineHash_, iter);
        
        int absent;
        khiter_t dest_iter = kh_put(cg_agg, inlineHash_, key, &absent);
        if (absent) {
            kh_val(inlineHash_, dest_iter) = count;
        } else {
            kh_val(inlineHash_, dest_iter) += count;
        }
    }
    
    // Merge readIdTracker_ if both have it
    // Note: For readId tracking, we keep ALL entries (no collision resolution needed)
    // Each readId should only appear in one thread's tracker
    if (readIdTracker_ && other.readIdTracker_) {
        for (khiter_t iter = kh_begin(other.readIdTracker_); iter != kh_end(other.readIdTracker_); ++iter) {
            if (!kh_exist(other.readIdTracker_, iter)) continue;
            
            uint32_t readId = kh_key(other.readIdTracker_, iter);
            uint64_t val = kh_val(other.readIdTracker_, iter);
            
            int absent;
            khiter_t dest_iter = kh_put(readid_cbumi, readIdTracker_, readId, &absent);
            // Should always be absent (each readId processed by one thread)
            kh_val(readIdTracker_, dest_iter) = val;
        }
    }
    
    // Merge ambiguous CB structs: combine UMI counts and observations on key collision
    for (const auto &kv : other.pendingAmbiguous_) {
        ReadAlign::AmbigKey key = kv.first;
        const ExtendedAmbiguousEntry &otherEntry = kv.second;
        
        auto &entry = pendingAmbiguous_[key];
        if (entry.candidateIdx.empty()) {
            // First time seeing this ambiguous CB: copy entire entry
            entry.candidateIdx = otherEntry.candidateIdx;
            entry.cbSeq = otherEntry.cbSeq;
            entry.cbQual = otherEntry.cbQual;
            entry.umiCounts = otherEntry.umiCounts;
            entry.observations = otherEntry.observations;
        } else {
            // Merge UMI counts
            for (const auto &umiCount : otherEntry.umiCounts) {
                entry.umiCounts[umiCount.first] += umiCount.second;
            }
            // Merge observations (gene/tag/umi combinations)
            entry.observations.insert(entry.observations.end(), 
                                     otherEntry.observations.begin(), 
                                     otherEntry.observations.end());
        }
    }
}

void SoloReadFeature::maybeSpillBinarySpool(size_t extraBytes)
{
    if (!binarySpoolInMemory || binarySpoolMemoryLimitBytes == 0) {
        return;
    }
    if (binarySpoolBuffer.totalBytes + extraBytes <= binarySpoolMemoryLimitBytes) {
        return;
    }

    streamReads = &fstrOpenBinary(binarySpoolFileName, ERROR_OUT, P, true);
    SoloBinarySpool::flushToStream(binarySpoolBuffer, *streamReads);
    streamReads->flush();
    binarySpoolBuffer.clear();
    binarySpoolReadPos = 0;
    binarySpoolInMemory = false;

    P.inOut->logMain << "Spilling experimental in-memory binary Solo temp spool to disk for feature "
                     << SoloFeatureTypes::Names[featureType]
                     << " after reaching "
                     << (binarySpoolMemoryLimitBytes / (1024ull * 1024ull))
                     << " MiB per-thread limit" << endl;
}
