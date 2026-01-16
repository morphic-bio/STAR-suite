#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "ProbeListIndex.h"
#include "ReadAlignChunk.h"
#include "Genome.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstdio>

// Static cache for gene→probe index (15-bit), 0 = not a probe gene
static std::vector<uint16_t> gGeneToProbeIdx;
// Static cache for chromosome→genomic flag
static std::vector<uint8_t> gChrIsGenomic;
static size_t gChrProbeCount = 0, gChrGenomicCount = 0, gChrENSGCount = 0;
static std::vector<std::string> gChrProbeExamples;
static std::vector<std::string> gChrGenomicExamples;
static bool gChrFlagsInitialized = false;

SoloFeature::SoloFeature(Parameters &Pin, ReadAlignChunk **RAchunk, Transcriptome &inTrans, int32 feTy, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll)
            : P(Pin), RAchunk(RAchunk), Trans(inTrans), featureType(feTy), soloFeatAll(soloFeatAll), pSolo(P.pSolo), readBarSum(readBarSumIn),
              umisBeforeTotal(0), umisAfterTotal(0), readsBeforeTotal(0), readsAfterTotal(0),
              mergesTotal(0), componentsTotal(0), componentsCappedTotal(0), componentsBelowThresholdTotal(0),
              maxComponentSeen(0), readsURGrouped(0), readsURMissing(0),
              cellsInput(0), cellsKept(0), ambiguousCellsDropped(0)
{
    // Initialize gene→probe cache once if probe index is available
    initGeneProbeIdx(Trans, getProbeListIndex());
    // Debug: summarize gene→probe mapping (helps diagnose genomic rescue issues)
    if (!gGeneToProbeIdx.empty()) {
        size_t nonZero = 0;
        for (size_t i = 0; i < gGeneToProbeIdx.size(); ++i) {
            if (gGeneToProbeIdx[i] != 0) {
                nonZero++;
            }
        }
        P.inOut->logMain << "[GENE-PROBE] cache_size=" << gGeneToProbeIdx.size()
                         << " mapped_genes=" << nonZero << "\n";
        if (nonZero == 0) {
            P.inOut->logMain << "[GENE-PROBE] WARNING: no gene IDs mapped to probe list; gGeneToProbeIdx is all zeros\n";
        } else {
            const auto& geID = Trans.geIDCanonical.empty() ? Trans.geID : Trans.geIDCanonical;
            size_t printed = 0;
            for (size_t i = 0; i < gGeneToProbeIdx.size() && printed < 5; ++i) {
                if (gGeneToProbeIdx[i] != 0 && i < geID.size()) {
                    P.inOut->logMain << "  [GENE-PROBE] geneIdx=" << i
                                     << " geID=" << geID[i]
                                     << " -> probeIdx=" << gGeneToProbeIdx[i] << "\n";
                    printed++;
                }
            }
        }
    }
    // Initialize chr→genomic flags (probe vs genomic) once
    initChrGenomicFlags(RAchunk);
    // Silent: chr flag info no longer logged
    // Initialize component size histogram
    for (int i = 0; i < 5; i++) {
        componentSizeHist[i] = 0;
    }
    
    if (featureType>=0) {//otherwise we do not need these arrays - e.g. with --runMode soloCellFiltering 
        readFeatSum = new SoloReadFeature(featureType,P,-1);
        readFeatAll = new SoloReadFeature*[P.runThreadN];
        readFeatSum->setOwner(this);
    };
    
    //number of features
    switch (featureType) {
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron :
        case SoloFeatureTypes::Velocyto :
            featuresNumber=Trans.nGe;
            break;
        case SoloFeatureTypes::SJ :
            featuresNumber=P.sjAll[0].size();
            break;
        default:
            featuresNumber = -1; //undefined
    };
    
    // Load assignment files if UMI correction is enabled
    if (pSolo.umiCorrectionMode > 0 && featureType == SoloFeatureTypes::Gene) {
        if (!pSolo.cellsAllowPath.empty() && pSolo.cellsAllowPath != "-") {
            if (!loadCellsAllowList(pSolo.cellsAllowPath)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: Failed to load cells allow-list from " << pSolo.cellsAllowPath << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
        }
        if (!pSolo.sampleAssignmentsPath.empty() && pSolo.sampleAssignmentsPath != "-") {
            if (!loadSampleAssignments(pSolo.sampleAssignmentsPath)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: Failed to load sample assignments from " << pSolo.sampleAssignmentsPath << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
        }

    }
};

const ProbeListIndex* SoloFeature::getProbeListIndex() const {
    // Lazily load from --soloProbeList once per process
    static ProbeListIndex *fallbackProbeIndex = nullptr;
    static bool attemptedLoad = false;
    if (!attemptedLoad) {
        attemptedLoad = true;
        if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
            ProbeListIndex *idx = new ProbeListIndex();
            uint32_t deprecatedCount = 0;
            if (idx->load(P.pSolo.probeListPath, P.pSolo.removeDeprecated, &deprecatedCount)) {
                fallbackProbeIndex = idx;
                if (P.pSolo.removeDeprecated && deprecatedCount > 0) {
                    P.inOut->logMain << "[PROBE-LIST] Removed " << deprecatedCount << " deprecated entries from probe list\n";
                }
            } else {
                delete idx;
            }
        }
    }
    return fallbackProbeIndex;
}

void SoloFeature::initGeneProbeIdx(const Transcriptome& trans, const ProbeListIndex* probeIdx) {
    if (!gGeneToProbeIdx.empty()) {
        return; // already initialized
    }
    gGeneToProbeIdx.assign(trans.nGe, 0);
    if (!probeIdx || probeIdx->empty()) {
        // Cannot log via P here – this is a static context.
        return;
    }
    const auto& geID = trans.geIDCanonical.empty() ? trans.geID : trans.geIDCanonical;
    for (size_t i = 0; i < geID.size(); ++i) {
        if (!geID[i].empty()) {
            uint16_t pIdx = probeIdx->geneIndex15(geID[i]);
            if (pIdx != 0 && pIdx <= 0x7FFFu) {
                gGeneToProbeIdx[i] = pIdx;
            }
        }
    }
}

uint16_t SoloFeature::getProbeIdxForGene(uint32_t geneIdx) {
    // Inline path uses gene indices coming from Transcriptome::geneFull.g[gi1]
    // (stored in ReadAnnotFeature::fAlign). In most Solo code these gene IDs
    // are treated as 1-based when indexing features, but gGeneToProbeIdx is
    // filled with 0-based indices into geID/geIDCanonical. To be robust, try
    // both conventions:
    //
    //  1) direct 0-based index (geneIdx)
    //  2) 1-based to 0-based (geneIdx-1), if the first lookup is empty.
    //
    // This preserves existing behavior for correctly aligned 0-based indices,
    // while allowing genomic rescue when fAlign carries 1-based gene IDs.

    if (geneIdx < gGeneToProbeIdx.size()) {
        uint16_t p = gGeneToProbeIdx[geneIdx];
        if (p != 0) return p;
    }
    if (geneIdx > 0) {
        uint32_t idx1 = geneIdx - 1;
        if (idx1 < gGeneToProbeIdx.size()) {
            uint16_t p = gGeneToProbeIdx[idx1];
            if (p != 0) return p;
        }
    }
    return 0;
}

void SoloFeature::initChrGenomicFlags(ReadAlignChunk **RAchunkIn) {
    if (!gChrIsGenomic.empty()) {
        return; // already initialized
    }
    if (!RAchunkIn || !RAchunkIn[0]) {
        return;
    }
    // Use full chromosome list to cover probe pseudo-chromosomes
    const auto &chrNames = RAchunkIn[0]->mapGen.chrNameAll.empty() ? RAchunkIn[0]->mapGen.chrName : RAchunkIn[0]->mapGen.chrNameAll;
    gChrIsGenomic.resize(chrNames.size(), 1);
    size_t probeCount = 0, genomicCount = 0, ensgCount = 0;
    gChrProbeExamples.clear();
    gChrGenomicExamples.clear();
    for (size_t i = 0; i < chrNames.size(); ++i) {
        const std::string &nm = chrNames[i];
        bool startsWithENSG = nm.rfind("ENSG", 0) == 0;
        bool isGenomic = !startsWithENSG;
        gChrIsGenomic[i] = isGenomic ? 1 : 0;
        if (startsWithENSG) {
            ensgCount++;
            if (gChrProbeExamples.size() < 5) gChrProbeExamples.push_back(nm);
        } else {
            if (gChrGenomicExamples.size() < 5) gChrGenomicExamples.push_back(nm);
        }
        if (isGenomic) genomicCount++; else probeCount++;
    }
    gChrProbeCount = probeCount;
    gChrGenomicCount = genomicCount;
    gChrENSGCount = ensgCount;
    gChrFlagsInitialized = true;
}

bool SoloFeature::isChrGenomic(uint32_t chrIdx) {
    if (chrIdx < gChrIsGenomic.size()) {
        return gChrIsGenomic[chrIdx] != 0;
    }
    return true; // default to genomic if unknown
}

size_t SoloFeature::geneProbeCacheSize() {
    return gGeneToProbeIdx.size();
}

// Shared helper implementations for readInfo management
void SoloFeature::resetPackedStorage(uint32_t nReads)
{
    // Skip allocation when minimal memory flag is on
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        return;
    }
    // Skip allocation when inline CB correction is active (Solo structures not used)
    if (pSolo.inlineCBCorrection) {
        // Assert that packedReadInfo stays empty
        assert(packedReadInfo.data.empty());
        return;
    }
    packedReadInfo.init(nReads, pSolo.cbWLstr.size(), pSolo.umiL);
}

void SoloFeature::recordReadInfo(uint32_t readId, uint32_t cbIdx, uint32_t umiPacked, uint8_t status)
{
    // Skip entirely when minimal memory flag is on
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        return;
    }
    if (packedReadInfo.data.empty()) {
        return;
    }
#ifndef DEBUG_CB_UB_PARITY
    static const bool allowCbConflict = (std::getenv("STAR_ALLOW_CB_CONFLICT") != nullptr);
#else
    static const bool allowCbConflict = (std::getenv("STAR_ALLOW_CB_CONFLICT") != nullptr);
#endif
#ifdef DEBUG_CB_UB_PARITY
    if (parityEnabled) {
        dbgWriteTotal++;
        if (status==0) dbgWriteStatus0++;
        else if (status==1) dbgWriteStatus1++;
        else if (status==2) dbgWriteStatus2++;
    }
#endif
    // Packed path stores validity in status; map sentinels to in-range placeholders.
    // status: 0=missing CB, 1=good, 2=invalid UMI
    uint32_t cbStore = cbIdx;
    uint32_t umiStore = umiPacked;
#ifdef DEBUG_CB_UB_PARITY
    // For parity validation, do not drop CB/UB on status==0; keep them to observe divergence
    if (parityEnabled && status==0) {
        status = 1;
        if (cbStore == (uint32)-1) cbStore = 0u;
        if (umiStore == (uint32)-1) umiStore = 0u;
    }
#endif
    if (status==0) { // missing CB ⇒ clear both
        cbStore = 0u;
        umiStore = 0u;
    } else if (status==2) { // invalid UMI ⇒ clear UMI only
        if (umiStore == (uint32)-1) umiStore = 0u;
        // keep cbStore within bounds; if sentinel, clear it
        if (cbStore == (uint32)-1) cbStore = 0u;
    } else {
        // good path: sanitize accidental sentinels
        if (cbStore == (uint32)-1) cbStore = 0u;
        if (umiStore == (uint32)-1) umiStore = 0u;
    }

    // Guard: same readId must not receive conflicting CB/UMI/status assignments
    uint8_t existingStatus = packedReadInfo.getStatus(readId);
    if (existingStatus != 0) { // already written once
        uint32_t existingCB = packedReadInfo.getCB(readId);
        uint32_t existingUMI = packedReadInfo.getUMI(readId);
        if ((existingCB != cbStore || existingUMI != umiStore || existingStatus != status) && !allowCbConflict) {
            fprintf(stderr, "[ERROR] Conflicting CB/UMI/status for readId=%u existing(cb=%u,umi=%u,status=%u) new(cb=%u,umi=%u,status=%u)\n",
                    readId, existingCB, existingUMI, existingStatus, cbStore, umiStore, status);
            std::exit(1);
        }
    }
    packedReadInfo.set(readId, cbStore, umiStore, status);
}

uint32_t SoloFeature::getPackedCB(uint32_t readId) const
{
    return packedReadInfo.getCB(readId);
}

uint32_t SoloFeature::getPackedUMI(uint32_t readId) const
{
    return packedReadInfo.getUMI(readId);
}

uint8_t SoloFeature::getPackedStatus(uint32_t readId) const
{
    return packedReadInfo.getStatus(readId);
}

#ifdef DEBUG_CB_UB_PARITY
void SoloFeature::resetDebugStatusCounters() {
    debugStatusCounters.clear();
    dbgBufferedRecords = 0;
    dbgBufferedCBs = 0;
    dbgWriteTotal = 0;
    dbgWriteSiteUnique = 0;
    dbgWriteSiteCliqueGood = 0;
    dbgWriteSiteCliqueDrop = 0;
    dbgWriteSiteMulti = 0;
    dbgRAWPass = 0;
    dbgRAWFail = 0;
    dbgWriteStatus0 = 0;
    dbgWriteStatus1 = 0;
    dbgWriteStatus2 = 0;
    dbgReadStatus0 = 0;
}

void SoloFeature::noteDebugStatus(uint8_t status, const char* reason) {
    std::string key = std::string(reason ? reason : "unknown");
    key += "|s=" + std::to_string(status);
    debugStatusCounters[key]++;
}

void SoloFeature::logDebugStatusCounters() {
    if (!parityEnabled) return;
    if (debugStatusCounters.empty()) return;
    P.inOut->logMain << "[PARITY-STATUS] counts\n";
    for (const auto &kv : debugStatusCounters) {
        P.inOut->logMain << "  " << kv.first << " -> " << kv.second << "\n";
    }
}

void SoloFeature::logDebugStageCounters() {
    if (!parityEnabled) return;
    P.inOut->logMain << "[PARITY-STAGE] loader_status_good=" << debugStatusCounters["good|s=1"]
                     << " loader_status_other=" << (debugStatusCounters["noTooManyWLmatches|s=0"] + debugStatusCounters["noMMtoWLwithoutExact|s=0"])
                     << " buffered_records=" << dbgBufferedRecords
                     << " buffered_cbs=" << dbgBufferedCBs
                     << " packed_writes_total=" << dbgWriteTotal
                     << " site_unique=" << dbgWriteSiteUnique
                     << " site_clique_good=" << dbgWriteSiteCliqueGood
                     << " site_clique_drop=" << dbgWriteSiteCliqueDrop
                     << " site_multi=" << dbgWriteSiteMulti
                     << " packed_writes_s0=" << dbgWriteStatus0
                     << " packed_writes_s1=" << dbgWriteStatus1
                     << " packed_writes_s2=" << dbgWriteStatus2
                     << " raw_pass=" << dbgRAWPass
                     << " raw_fail=" << dbgRAWFail
                     << " readback_status0=" << dbgReadStatus0
                     << " parity_inline_reads=" << parityReadsInline
                     << " parity_mismatches=" << parityMismatches
                     << "\n";
}
#endif

void SoloFeature::clearLarge()
{
    cbFeatureUMImap.clear();
    cbFeatureUMImap.shrink_to_fit();
    countCellGeneUMI.clear();
    countCellGeneUMI.shrink_to_fit();
    countCellGeneUMIindex.clear();
    countCellGeneUMIindex.shrink_to_fit();
    countMatMult.i.clear();
    countMatMult.i.shrink_to_fit();
    countMatMult.m.clear();
    countMatMult.m.shrink_to_fit();
    //indCB.clear(); //needed for Velocyto
    //indCB.shrink_to_fit();
    indCBwl.clear();
    indCBwl.shrink_to_fit();
    nGenePerCB.clear();
    nGenePerCB.shrink_to_fit();
    nGenePerCBmulti.clear();
    nGenePerCBmulti.shrink_to_fit();
    nReadPerCB.clear();
    nReadPerCB.shrink_to_fit();
    nReadPerCBtotal.clear();
    nReadPerCBtotal.shrink_to_fit();
    nReadPerCBunique.clear();
    nReadPerCBunique.shrink_to_fit();
    nUMIperCB.clear();
    nUMIperCB.shrink_to_fit();
    nUMIperCBmulti.clear();
    nUMIperCBmulti.shrink_to_fit();
    nUMIperCBsorted.clear();
    nUMIperCBsorted.shrink_to_fit();
    sjAll[0].clear();
    sjAll[0].shrink_to_fit();
    sjAll[1].clear();
    sjAll[1].shrink_to_fit();
    
    // Clear UMI correction data structures
    cellsAllowSet.clear();
    assignmentsMap.clear();
    umiGroupHistograms.clear();
    
#ifdef DEBUG_CB_UB_PARITY
    // Cleanup parity validation resources
    if (parityEnabled && mismatchLog.is_open()) {
        mismatchLog.close();
    }
    legacyCBUBByRead.clear();
#endif

    // Defensive cleanup for minimal memory mode
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        // Clear packed read info if somehow allocated (defensive - should never be allocated)
        packedReadInfo.data.clear();
        packedReadInfo.data.shrink_to_fit();
        
        // Ensure inline hash is destroyed (defensive - should already be destroyed)
        if (readFeatSum && readFeatSum->inlineHash_) {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
    }
};

// Load cells allow-list (post-EmptyDrops CB16s or full composites CB24)
// Supports both formats: 16bp entries act as wildcards (any tag), 24bp entries require exact match
bool SoloFeature::loadCellsAllowList(const std::string& path) {
    P.inOut->logMain << "Loading cells allow-list from " << path << "..." << endl;
    
    std::ifstream fp(path);
    if (!fp.is_open()) {
        P.inOut->logMain << "ERROR: Cannot open cells allow-list file " << path << endl;
        return false;
    }
    
    cellsAllowSet.clear();
    int count = 0;
    int composite_count = 0;
    int cb16_count = 0;
    
    std::string line;
    while (std::getline(fp, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty() || line[0] == '#') continue;
        
        size_t len = line.length();
        std::string key;
        
        if (len >= 24) {
            // Full composite (CB16+TAG8) - exact match required
            key = line.substr(0, 24);
            composite_count++;
        } else if (len >= 16) {
            // CB16 only - acts as wildcard allowing any tag
            key = line.substr(0, 16);
            cb16_count++;
        } else {
            // Too short, skip
            continue;
        }
        
        cellsAllowSet.insert(key);
        count++;
    }
    
    if (composite_count > 0) {
        P.inOut->logMain << "Loaded " << count << " entries from allow-list (" 
                         << composite_count << " composites, " << cb16_count << " CB16 wildcards)" << endl;
    } else {
        P.inOut->logMain << "Loaded " << count << " CB16s from allow-list" << endl;
    }
    
    return count > 0;
}

// Load sample assignments TSV (cb16\tsample_idx\tsample_tag\tstatus)
// Format: full composite barcode (CB16+TAG8) or CB16, sample_idx, optional sample_tag, optional status
bool SoloFeature::loadSampleAssignments(const std::string& path) {
    P.inOut->logMain << "Loading sample assignments from " << path << "..." << endl;
    
    std::ifstream fp(path);
    if (!fp.is_open()) {
        P.inOut->logMain << "ERROR: Cannot open sample assignments file " << path << endl;
        return false;
    }
    
    assignmentsMap.clear();
    int line_num = 0;
    int count = 0;
    int ambiguous_count = 0;
    
    std::string line;
    while (std::getline(fp, line)) {
        line_num++;
        
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty() || line[0] == '#') continue;
        
        // Parse TSV: cb16\tsample_idx[\tsample_tag\tstatus]
        std::istringstream iss(line);
        std::string cb_input, sample_idx_str, sample_tag, status;
        
        if (!(iss >> cb_input >> sample_idx_str)) {
            P.inOut->logMain << "Warning: Sample assignments line " << line_num 
                            << " has no tab separator, skipping" << endl;
            continue;
        }
        
        // Parse optional columns
        if (iss >> sample_tag) {
            iss >> status;  // May be empty
        }
        
        // Extract CB16 (first 16 chars)
        std::string cb16;
        if (cb_input.length() >= 16) {
            cb16 = cb_input.substr(0, 16);
        } else {
            continue;  // Too short
        }
        
        // Parse sample_idx
        int sample_idx = std::atoi(sample_idx_str.c_str());
        
        // Handle ambiguous/no_call status
        if (!status.empty() && (status == "ambiguous" || status == "no_call")) {
            AssignmentInfo info;
            info.sampleIdx = 0;
            info.sampleTag = sample_tag;
            info.status = status;
            assignmentsMap[cb16] = info;
            ambiguous_count++;
            count++;
            continue;
        }
        
        // Validate sample_idx (1-16)
        if (sample_idx > 0 && sample_idx <= 16) {
            AssignmentInfo info;
            info.sampleIdx = static_cast<uint8_t>(sample_idx);
            info.sampleTag = sample_tag;
            // Normalize status: "assigned" -> "CONFIDENT", empty -> "CONFIDENT"
            if (status.empty() || status == "assigned") {
                info.status = "CONFIDENT";
            } else {
                info.status = status;
            }
            assignmentsMap[cb16] = info;
            count++;
        } else {
            P.inOut->logMain << "Warning: Invalid sample_idx " << sample_idx 
                            << " on line " << line_num << ", skipping" << endl;
        }
    }
    
    P.inOut->logMain << "Loaded " << count << " assignments from file (" 
                     << ambiguous_count << " ambiguous)" << endl;
    
    return count > 0;
}
