#include "IncludeDefine.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "OutSJ.h"
#include "sysRemoveDir.h"
#include "stringSubstituteAll.h"
#include SAMTOOLS_BGZF_H
#include "GlobalVariables.h"
#include "signalFromBAM.h"
#include "bamRemoveDuplicates.h"
#include "streamFuns.h"
#include <fstream>
#include <cctype>

// Define global atomic counter for processed read groups (matches Salmon's processedReads)
// Used for pre-burn-in gating: aux params are enabled when this count >= numPreBurninFrags (5000)
// Each read group (nAlignT > 0) counts as one, regardless of compat/FLD validity
std::atomic<uint64_t> global_processed_fragments{0};

// Define global atomic counter for FLD observations (fragments that contribute to FLD)
// Used for FLD statistics only, NOT for pre-burn-in gating
std::atomic<uint64_t> Parameters::global_fld_obs_count{0};

namespace {
string pathBasename(const string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos == string::npos) {
        return path;
    }
    return path.substr(pos + 1);
}

string outputDirFromPrefix(const string& prefix) {
    size_t pos = prefix.find_last_of("/\\");
    if (pos == string::npos) {
        return "";
    }
    return prefix.substr(0, pos + 1);
}

bool endsWith(const string& value, const string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool hasReadToken(const string& name) {
    return name.rfind("_R1") != string::npos ||
           name.rfind("_R2") != string::npos ||
           name.rfind("_r1") != string::npos ||
           name.rfind("_r2") != string::npos;
}

string insertTagBeforeReadToken(const string& name, const string& tag) {
    size_t pos = string::npos;
    auto updatePos = [&](size_t candidate) {
        if (candidate != string::npos) {
            pos = (pos == string::npos) ? candidate : max(pos, candidate);
        }
    };
    updatePos(name.rfind("_R1"));
    updatePos(name.rfind("_R2"));
    updatePos(name.rfind("_r1"));
    updatePos(name.rfind("_r2"));
    if (pos == string::npos) {
        return string();
    }
    return name.substr(0, pos) + tag + name.substr(pos);
}

string adjustCompressionExt(const string& name, const string& compression) {
    if (compression == "gz") {
        if (!endsWith(name, ".gz")) {
            return name + ".gz";
        }
        return name;
    }
    if (endsWith(name, ".gz")) {
        return name.substr(0, name.size() - 3);
    }
    return name;
}
}

//for mkfifo
#include <sys/stat.h>
#include <cstdlib>
#include <algorithm>

#define PAR_NAME_PRINT_WIDTH 30

Parameters::Parameters() {//initalize parameters info

    inOut = new InOutStreams;

    //versions
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "versionGenome", &versionGenome));

    //parameters
    parArray.push_back(new ParameterInfoVector <string> (-1, 2, "parametersFiles", &parametersFiles));

    //system
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sysShell", &sysShell));

    //run
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "runMode", &runModeIn));
    parArray.push_back(new ParameterInfoScalar <int> (-1, -1, "runThreadN", &runThreadN));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "runDirPerm", &runDirPermIn));
    parArray.push_back(new ParameterInfoScalar <int> (-1, -1, "runRNGseed", &runRNGseed));

    //genome
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeType", &pGe.gTypeString));    
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeDir", &pGe.gDir));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeLoad", &pGe.gLoad));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeFastaFiles", &pGe.gFastaFiles));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeChainFiles", &pGe.gChainFiles));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAindexNbases", &pGe.gSAindexNbases));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeChrBinNbits", &pGe.gChrBinNbits));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAsparseD", &pGe.gSAsparseD));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSuffixLengthMax", &pGe.gSuffixLengthMax));
    parArray.push_back(new ParameterInfoVector <uint> (-1, -1, "genomeFileSizes", &pGe.gFileSizes));
    //parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeConsensusFile", &pGe.gConsensusFile)); DEPRECATED
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeTransformType", &pGe.transform.typeString));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeTransformVCF", &pGe.transform.vcfFile));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeTransformOutput", &pGe.transform.output));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeChrSetMitochondrial", &pGe.chrSet.mitoStrings));

    // Flex gene probe parameters (50bp gene probes for Flex workflow)
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "flexGeneProbeSet", &pGe.flexGeneProbe.csvFile));
    parArray.push_back(new ParameterInfoScalar <uint32> (-1, -1, "flexGeneProbeLength", &pGe.flexGeneProbe.enforceLength));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "removeDeprecated", &pGe.flexGeneProbe.removeDeprecated));

    // CellRanger-style reference formatting parameters
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "cellrangerStyleIndex", &pGe.cellrangerStyle.indexEnabled));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "cellrangerStyleDownloadOnly", &pGe.cellrangerStyle.downloadOnly));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "cellrangerStyleCacheDir", &pGe.cellrangerStyle.cacheDir));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "allUntrustedUrl", &pGe.cellrangerStyle.allUntrustedUrl));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "faUrl", &pGe.cellrangerStyle.faUrl));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "gtfUrl", &pGe.cellrangerStyle.gtfUrl));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "autoCksumUpdate", &pGe.cellrangerStyle.autoCksumUpdate));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "replaceUnverifiableFiles", &pGe.cellrangerStyle.replaceUnverifiableFiles));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "cellrangerRefRelease", &pGe.cellrangerStyle.refRelease));
    
    // Auto-index workflow parameters
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "autoIndex", &pGe.autoIndexWorkflow.autoIndex));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "forceIndex", &pGe.autoIndexWorkflow.forceIndex));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "forceAllIndex", &pGe.autoIndexWorkflow.forceAllIndex));
    
    // Transcriptome FASTA generation parameters
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeGenerateTranscriptome", &pGe.transcriptomeGen.generateTranscriptome));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeGenerateTranscriptomeFasta", &pGe.transcriptomeGen.transcriptomeFastaPath));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeGenerateTranscriptomeOverwrite", &pGe.transcriptomeGen.overwrite));

    //read
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesType", &readFilesType));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesIn", &readFilesIn));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readFilesPrefix", &readFilesPrefix));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesCommand", &readFilesCommand));

    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readMatesLengthsIn", &readMatesLengthsIn));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "readMapNumber", &readMapNumber));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readNameSeparator", &readNameSeparator));
    parArray.push_back(new ParameterInfoScalar <uint32> (-1, -1, "readQualityScoreBase", &readQualityScoreBase));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesManifest", &readFilesManifest));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesSAMattrKeep", &readFiles.samAttrKeepIn));

    //parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readStrand", &pReads.strandString));


    //input from BAM
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "inputBAMfile", &inputBAMfile));

    //BAM processing
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "bamRemoveDuplicatesType", &removeDuplicates.mode));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "bamRemoveDuplicatesMate2basesN", &removeDuplicates.mate2basesN));

    //limits
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitGenomeGenerateRAM", &limitGenomeGenerateRAM));
    parArray.push_back(new ParameterInfoVector <uint64>   (-1, -1, "limitIObufferSize", &limitIObufferSize));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSAMoneReadBytes", &limitOutSAMoneReadBytes));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJcollapsed", &limitOutSJcollapsed));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJoneRead", &limitOutSJoneRead));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitBAMsortRAM", &limitBAMsortRAM));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitSjdbInsertNsj", &limitSjdbInsertNsj));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitNreadsSoft", &limitNreadsSoft));

    //output
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outFileNamePrefix", &outFileNamePrefix));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outTmpDir", &outTmpDir));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outTmpKeep", &outTmpKeep));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outStd", &outStd));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outReadsUnmapped", &outReadsUnmapped));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outQSconversionAdd", &outQSconversionAdd));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outMultimapperOrder", &outMultimapperOrder.mode));

    //outSAM
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMtype", &outSAMtype));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMmode", &outSAMmode));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMstrandField", &outSAMstrandField.in));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMattributes", &outSAMattributes));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMunmapped", &outSAMunmapped.mode));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMorder", &outSAMorder));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMprimaryFlag", &outSAMprimaryFlag));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMreadID", &outSAMreadID));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outSAMmapqUnique", &outSAMmapqUnique));
    parArray.push_back(new ParameterInfoScalar <uint16>        (-1, -1, "outSAMflagOR", &outSAMflagOR));
    parArray.push_back(new ParameterInfoScalar <uint16>        (-1, -1, "outSAMflagAND", &outSAMflagAND));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMattrRGline", &outSAMattrRGline));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMheaderHD", &outSAMheaderHD));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMheaderPG", &outSAMheaderPG));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMheaderCommentFile", &outSAMheaderCommentFile));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outBAMcompression", &outBAMcompression));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outBAMsortingThreadN", &outBAMsortingThreadN));
    parArray.push_back(new ParameterInfoScalar <uint32>        (-1, -1, "outBAMsortingBinsN", &outBAMsortingBinsN));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outBAMsortMethod", &outBAMsortMethod));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "emitNoYBAM", &emitNoYBAM));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "emitYReadNames", &emitYReadNames));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "keepBAM", &keepBAM));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "noYOutput", &noYOutput));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "YOutput", &YOutput));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "YReadNamesOutput", &YReadNamesOutput));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "emitYNoYFastq", &emitYNoYFastq));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "emitYNoYFastqCompression", &emitYNoYFastqCompression));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "YFastqOutputPrefix", &YFastqOutputPrefix));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "noYFastqOutputPrefix", &noYFastqOutputPrefix));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSAMfilter", &outSAMfilter.mode));
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outSAMmultNmax", &outSAMmultNmax));
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outSAMattrIHstart", &outSAMattrIHstart));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outSAMtlen", &outSAMtlen));

    //outSJ
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "outSJtype", &outSJ.type));
    
    //output SJ filtering
    parArray.push_back(new ParameterInfoScalar <string>  (-1, -1, "outSJfilterReads", &outSJfilterReads));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountUniqueMin", &outSJfilterCountUniqueMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountTotalMin", &outSJfilterCountTotalMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterOverhangMin", &outSJfilterOverhangMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterDistToOtherSJmin", &outSJfilterDistToOtherSJmin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterIntronMaxVsReadN", &outSJfilterIntronMaxVsReadN));

    //output wiggle
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigType", &outWigType));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigStrand", &outWigStrand));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outWigReferencesPrefix", &outWigReferencesPrefix));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "outWigNorm", &outWigNorm));

    //output filtering
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterType", &outFilterType) );

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMultimapNmax", &outFilterMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterMultimapScoreRange", &outFilterMultimapScoreRange));

    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterScoreMin", &outFilterScoreMin));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterScoreMinOverLread", &outFilterScoreMinOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMatchNmin", &outFilterMatchNmin));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMatchNminOverLread", &outFilterMatchNminOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMismatchNmax", &outFilterMismatchNmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMismatchNoverLmax", &outFilterMismatchNoverLmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMismatchNoverReadLmax", &outFilterMismatchNoverReadLmax));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterIntronMotifs", &outFilterIntronMotifs));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterIntronStrands", &outFilterIntronStrands));

    //clipping
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "clipAdapterType", &pClip.adapterType));
    parArray.push_back(new ParameterInfoVector <uint32>   (-1, -1, "clip5pNbases", &pClip.in[0].N));
    parArray.push_back(new ParameterInfoVector <uint32>   (-1, -1, "clip3pNbases", &pClip.in[1].N));
    parArray.push_back(new ParameterInfoVector <uint32>   (-1, -1, "clip5pAfterAdapterNbases", &pClip.in[0].NafterAd));
    parArray.push_back(new ParameterInfoVector <uint32>   (-1, -1, "clip3pAfterAdapterNbases", &pClip.in[1].NafterAd));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "clip5pAdapterSeq", &pClip.in[0].adSeq));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "clip3pAdapterSeq", &pClip.in[1].adSeq));
    parArray.push_back(new ParameterInfoVector <double> (-1, -1, "clip5pAdapterMMp", &pClip.in[0].adMMp));
    parArray.push_back(new ParameterInfoVector <double> (-1, -1, "clip3pAdapterMMp", &pClip.in[1].adMMp));

    //cutadapt-style trimming
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "trimCutadapt", &trimCutadapt));
    parArray.push_back(new ParameterInfoScalar <uint8> (-1, -1, "trimCutadaptQuality", &trimCutadaptQuality));
    parArray.push_back(new ParameterInfoScalar <uint32> (-1, -1, "trimCutadaptMinLength", &trimCutadaptMinLength));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "trimCutadaptAdapter", &trimCutadaptAdapter));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "trimCutadaptCompat", &trimCutadaptCompat));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "trimQcReport", &trimQcReport));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "trimQcJson", &trimQcJson));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "trimQcHtml", &trimQcHtml));
    parArray.push_back(new ParameterInfoScalar <uint64> (-1, -1, "trimQcMaxReads", &trimQcMaxReads));

    //binning, anchors, windows
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winBinNbits", &winBinNbits));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorDistNbins", &winAnchorDistNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winFlankNbins", &winFlankNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorMultimapNmax", &winAnchorMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "winReadCoverageRelativeMin", &winReadCoverageRelativeMin));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winReadCoverageBasesMin", &winReadCoverageBasesMin));

    //scoring
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGap", &scoreGap));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapNoncan", &scoreGapNoncan));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapGCAG", &scoreGapGCAG));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapATAC", &scoreGapATAC));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreStitchSJshift", &scoreStitchSJshift));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "scoreGenomicLengthLog2scale", &scoreGenomicLengthLog2scale));

    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelBase", &scoreDelBase));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelOpen", &scoreDelOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsOpen", &scoreInsOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsBase", &scoreInsBase));

    //alignment
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchLmax", &seedSearchLmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchStartLmax", &seedSearchStartLmax));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "seedSearchStartLmaxOverLread", &seedSearchStartLmaxOverLread));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerReadNmax", &seedPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerWindowNmax", &seedPerWindowNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedNoneLociPerWindow", &seedNoneLociPerWindow));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedMultimapNmax", &seedMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSplitMin", &seedSplitMin));
    parArray.push_back(new ParameterInfoScalar <uint64>       (-1, -1, "seedMapMin", &seedMapMin));
    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMin", &alignIntronMin));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMax", &alignIntronMax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignMatesGapMax", &alignMatesGapMax));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerReadNmax", &alignTranscriptsPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJoverhangMin", &alignSJoverhangMin));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJDBoverhangMin", &alignSJDBoverhangMin));
    parArray.push_back(new ParameterInfoVector <int32>      (-1, -1, "alignSJstitchMismatchNmax", &alignSJstitchMismatchNmax));

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSplicedMateMapLmin", &alignSplicedMateMapLmin));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "alignSplicedMateMapLminOverLmate", &alignSplicedMateMapLminOverLmate));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignWindowsPerReadNmax", &alignWindowsPerReadNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerWindowNmax", &alignTranscriptsPerWindowNmax));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignEndsType", &alignEndsType.in));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignSoftClipAtReferenceEnds", &alignSoftClipAtReferenceEnds.in));

    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "alignEndsProtrude", &alignEndsProtrude.in));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "alignInsertionFlush", &alignInsertionFlush.in));

    //peOverlap
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "peOverlapNbasesMin", &peOverlap.NbasesMin));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "peOverlapMMp", &peOverlap.MMp));

    //chimeric
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimSegmentMin", &pCh.segmentMin));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreMin", &pCh.scoreMin));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreDropMax", &pCh.scoreDropMax));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreSeparation", &pCh.scoreSeparation));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "chimScoreJunctionNonGTAG", &pCh.scoreJunctionNonGTAG));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMainSegmentMultNmax", &pCh.mainSegmentMultNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimJunctionOverhangMin", &pCh.junctionOverhangMin));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "chimOutType", &pCh.out.type));
    parArray.push_back(new ParameterInfoVector <string>     (-1, -1, "chimFilter", &pCh.filter.stringIn));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimSegmentReadGapMax", &pCh.segmentReadGapMax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMultimapNmax", &pCh.multimapNmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimMultimapScoreRange", &pCh.multimapScoreRange));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimNonchimScoreDropMin", &pCh.nonchimScoreDropMin));
    parArray.push_back(new ParameterInfoVector <int>        (-1, -1, "chimOutJunctionFormat", &pCh.outJunctionFormat));

    //sjdb
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "sjdbFileChrStartEnd", &pGe.sjdbFileChrStartEnd));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfile", &pGe.sjdbGTFfile));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFchrPrefix", &pGe.sjdbGTFchrPrefix));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "transcriptomeFasta", &pGe.transcriptomeFasta));
    
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfeatureExon", &pGe.sjdbGTFfeatureExon));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFtagExonParentTranscript", &pGe.sjdbGTFtagExonParentTranscript));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFtagExonParentGene", &pGe.sjdbGTFtagExonParentGene));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "sjdbGTFtagExonParentGeneName", &pGe.sjdbGTFtagExonParentGeneName));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "sjdbGTFtagExonParentGeneType", &pGe.sjdbGTFtagExonParentGeneType));

    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "sjdbOverhang", &pGe.sjdbOverhang));
    pGe.sjdbOverhang_par=parArray.size()-1;
    parArray.push_back(new ParameterInfoScalar <int>    (-1, -1, "sjdbScore", &pGe.sjdbScore));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbInsertSave", &pGe.sjdbInsertSave));

    //variation
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "varVCFfile", &var.vcfFile));

    //WASP
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "waspOutputMode", &wasp.outputMode));

    //quant
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "quantMode", &quant.mode));
    parArray.push_back(new ParameterInfoScalar <int>     (-1, -1, "quantTranscriptomeBAMcompression", &quant.trSAM.bamCompression));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "quantTranscriptomeSAMoutput", &quant.trSAM.output));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "quantVBgcBias", &quant.transcriptVB.gcBiasInt));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "quantVBprior", &quant.transcriptVB.vbPrior));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "quantVBem", &quant.transcriptVB.quantVBemInt)); // If true, use EM instead of VB
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "quantVBgenes", &quant.transcriptVB.geneOutputInt));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "quantVBgenesMode", &quant.transcriptVB.genesModeStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "quantVBLibType", &quant.transcriptVB.libType));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "quantVBAutoDetectWindow", &quant.transcriptVB.autoDetectWindow));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "quantVBTrace", &quant.transcriptVB.traceFile));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "quantVBTraceLimit", &quant.transcriptVB.traceLimit));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "quantVBErrorModel", &quant.transcriptVB.errorModelMode));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamQuantMode", &quant.slam.modeInt));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpBed", &quant.slam.snpBed));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamSnpDetect", &quant.slam.snpDetectInt));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpDetectFrac", &quant.slam.snpDetectFrac));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamStrandness", &quant.slam.strandnessStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamWeightMode", &quant.slam.weightModeStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDebugGeneList", &quant.slam.debugGeneList));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDebugReadList", &quant.slam.debugReadList));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDebugOutPrefix", &quant.slam.debugOutPrefix));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamDebugMaxReads", &quant.slam.debugMaxReads));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDebugSnpLoc", &quant.slam.debugSnpLoc));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamDebugSnpWindow", &quant.slam.debugSnpWindow));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamErrorRate", &quant.slam.errorRate));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamConvRate", &quant.slam.convRate));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamOutFile", &quant.slam.outFile));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamCompatMode", &quant.slam.compatModeStr));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatIntronic", &quant.slam.compatIntronicInt));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatLenientOverlap", &quant.slam.compatLenientOverlapInt));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatOverlapWeight", &quant.slam.compatOverlapWeightInt));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatIgnoreOverlap", &quant.slam.compatIgnoreOverlapInt));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatTrim5p", &quant.slam.compatTrim5p));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamCompatTrim3p", &quant.slam.compatTrim3p));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "autoTrim", &quant.slam.autoTrimMode));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "trimScope", &quant.slam.trimScope));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "trimSource", &quant.slam.trimSource));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimMaxReads", &quant.slam.autoTrimMaxReads));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimMinReads", &quant.slam.autoTrimMinReads));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimSmoothWindow", &quant.slam.autoTrimSmoothWindow));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimSegMinLen", &quant.slam.autoTrimSegMinLen));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimMaxTrim", &quant.slam.autoTrimMaxTrim));
    parArray.push_back(new ParameterInfoScalar <uint64_t> (-1, -1, "autoTrimBufferReads", &quant.slam.autoTrimBufferReads));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "autoTrimDetectionReads", &quant.slam.autoTrimDetectionReads));
    parArray.push_back(new ParameterInfoScalar <double>  (-1, -1, "snpErrMinThreshold", &quant.slam.snpErrMinThreshold));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamQcJson", &quant.slam.slamQcJson));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamQcHtml", &quant.slam.slamQcHtml));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamQcReport", &quant.slam.slamQcReport));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamGrandSlamOut", &quant.slam.grandSlamOut));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDumpBinary", &quant.slam.dumpBinary));
    parArray.push_back(new ParameterInfoScalar <uint64_t> (-1, -1, "slamDumpMaxReads", &quant.slam.dumpMaxReads));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDumpWeights", &quant.slam.dumpWeights));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamDumpWeightsMode", &quant.slam.dumpWeightsModeStr));
    
    // SLAM SNP mask build parameters
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskIn", &quant.slamSnpMask.maskIn));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskVcfIn", &quant.slamSnpMask.vcfIn));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskVcfSample", &quant.slamSnpMask.vcfSample));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskVcfMode", &quant.slamSnpMask.vcfMode));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskVcfFilter", &quant.slamSnpMask.vcfFilter));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskBuildFastqs", &quant.slamSnpMask.buildFastqsFofn));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskBuildBam", &quant.slamSnpMask.buildBam));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "slamSnpMaskOnly", &quant.slamSnpMask.buildOnlyInt));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskBedOut", &quant.slamSnpMask.bedOut));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskSummaryOut", &quant.slamSnpMask.summaryOut));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskBamOut", &quant.slamSnpMask.bamOut));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskCompat", &quant.slamSnpMask.compat));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskKMode", &quant.slamSnpMask.kMode));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "slamSnpMaskModel", &quant.slamSnpMask.model));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpMaskPval", &quant.slamSnpMask.pval));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpMaskMinTcRatio", &quant.slamSnpMask.minTcRatio));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpMaskErr", &quant.slamSnpMask.err));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskMinCov", &quant.slamSnpMask.minCov));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskMinAlt", &quant.slamSnpMask.minAlt));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpMaskPosterior", &quant.slamSnpMask.posterior));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "slamSnpMaskMaxIter", &quant.slamSnpMask.maxIter));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "slamSnpMaskConvergeRelLL", &quant.slamSnpMask.convergeRelLL));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskJunctionFlank", &quant.slamSnpMask.junctionFlank));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskIndelFlank", &quant.slamSnpMask.indelFlank));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskMinMapQ", &quant.slamSnpMask.minMapQ));
    parArray.push_back(new ParameterInfoScalar <int32>    (-1, -1, "slamSnpMaskMinBaseQ", &quant.slamSnpMask.minBaseQ));

    //2-pass
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "twopass1readsN", &twoPass.pass1readsN));
    twoPass.pass1readsN_par=parArray.size()-1;
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "twopassMode", &twoPass.mode));

    //solo
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloType", &pSolo.typeStr));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloCBstart", &pSolo.cbS));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloUMIstart", &pSolo.umiS));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloCBlen", &pSolo.cbL));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloUMIlen", &pSolo.umiL));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloBarcodeReadLength", &pSolo.bL));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloBarcodeMate", &pSolo.barcodeReadIn));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloCBwhitelist", &pSolo.soloCBwhitelist));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloStrand", &pSolo.strandStr));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloOutFileNames", &pSolo.outFileNames));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloFeatures", &pSolo.featureIn));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloUMIdedup", &pSolo.umiDedup.typesIn));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloAdapterSequence",&pSolo.adapterSeq));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloAdapterMismatchesNmax", &pSolo.adapterMismatchesNmax));
    parArray.push_back(new ParameterInfoScalar <string>    (-1, -1, "soloCBmatchWLtype", &pSolo.CBmatchWL.type));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloCBposition",&pSolo.cbPositionStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloUMIposition",&pSolo.umiPositionStr));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloCellFilter",&pSolo.cellFilter.type));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloUMIfiltering",&pSolo.umiFiltering.type));
    
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloMultiMappers", &pSolo.multiMap.typesIn));
    
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloClusterCBfile",&pSolo.clusterCBfile));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloOutFormatFeaturesGeneField3",&pSolo.outFormat.featuresGeneField3));
    
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloInputSAMattrBarcodeSeq",&pSolo.samAtrrBarcodeSeq));
    parArray.push_back(new ParameterInfoVector <string>   (-1, -1, "soloInputSAMattrBarcodeQual",&pSolo.samAtrrBarcodeQual));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloCellReadStats",&pSolo.readStats.type));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloCBtype",&pSolo.CBtype.typeString));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloAddTagsToUnsorted",&pSolo.addTagsToUnsortedStr));
    // legacy sidecar flag removed; tag-table export is always enabled via unified writer
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloWriteKeysBin",&pSolo.writeKeysBinStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSkipProcessing",&pSolo.skipProcessingStr));

    // CR-compatible keys mode (handoff)
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloKeysCompat", &pSolo.keysCompatStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloProbeList", &pSolo.probeListPath));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloRemoveDeprecated", &pSolo.removeDeprecatedStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSampleWhitelist", &pSolo.sampleWhitelistPath));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSampleProbes", &pSolo.sampleProbesPath));
    parArray.push_back(new ParameterInfoScalar <uint32>   (-1, -1, "soloSampleProbeOffset", &pSolo.sampleProbeOffset));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSampleSearchNearby", &pSolo.sampleSearchNearbyStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSampleStrictMatch", &pSolo.sampleStrictMatchStr));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "soloMapqThreshold", &pSolo.mapqThreshold));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloMapqMode", &pSolo.mapqModeStr));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "soloNMmax", &pSolo.nmMax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "soloMMrateMax", &pSolo.mmRateMax));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloUniqueHarmonize", &pSolo.uniqueHarmonizeStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloBarcodesObservedOnly", &pSolo.barcodesObservedOnlyStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloCellsAllow", &pSolo.cellsAllowPath));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloSampleAssignments", &pSolo.sampleAssignmentsPath));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloUMICorrection", &pSolo.umiCorrectionModeStr));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "soloUMIMinCount", &pSolo.umiMinCount));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "soloUMIRatioThresh", &pSolo.umiRatioThresh));
    parArray.push_back(new ParameterInfoScalar <int>      (-1, -1, "soloMaxComponentSize", &pSolo.maxComponentSize));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloUMICorrectionUseTags", &pSolo.umiCorrectionUseTagsStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloAssignAmbiguous", &pSolo.assignAmbiguousStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloUseInlineReplayer", &pSolo.useInlineReplayerStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloInlineCBCorrection", &pSolo.inlineCBCorrectionStr));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "soloInlineHashMode", &pSolo.inlineHashModeStr));

    // Flex omnibus flag
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "flex", &pSolo.flexModeStr));
    
    // FlexFilter inline integration
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloRunFlexFilter", &pSolo.runFlexFilterStr));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexFatalOnError", &pSolo.flexFilterFatalOnErrorStr));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexTotalExpected", &pSolo.flexFilterTotalExpected));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexExpectedCellsTotal", &pSolo.flexFilterExpectedCellsTotal));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexExpectedCellsPerTag", &pSolo.flexFilterExpectedCellsPerTag));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexAllowedTags", &pSolo.flexFilterAllowedTagsPath));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexOutputPrefix", &pSolo.flexFilterOutputPrefix));

    // OrdMag parameters
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexOrdmagNsamples", &pSolo.flexFilterOrdmagNsamples));
    parArray.push_back(new ParameterInfoScalar<double>(-1, -1, "soloFlexOrdmagUmiMin", &pSolo.flexFilterOrdmagUmiMin));
    parArray.push_back(new ParameterInfoScalar<double>(-1, -1, "soloFlexOrdmagTargetPct", &pSolo.flexFilterOrdmagTargetPct));

    // EmptyDrops parameters
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexEdLower", &pSolo.flexFilterEdLower));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexEdNiters", &pSolo.flexFilterEdNiters));
    parArray.push_back(new ParameterInfoScalar<double>(-1, -1, "soloFlexEdFdrThreshold", &pSolo.flexFilterEdFdrThreshold));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexEdMaxTotalBuckets", &pSolo.flexFilterEdMaxTotalBuckets));

    // Occupancy parameters
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexTotalPartitions", &pSolo.flexFilterTotalPartitions));
    parArray.push_back(new ParameterInfoScalar<double>(-1, -1, "soloFlexRecoveryFactor", &pSolo.flexFilterRecoveryFactor));
    parArray.push_back(new ParameterInfoScalar<double>(-1, -1, "soloFlexOccupancyPercentile", &pSolo.flexFilterOccupancyPercentile));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexLowUmiThreshold", &pSolo.flexFilterLowUmiThreshold));

    // Simple EmptyDrops (fallback filter) flags
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexUseSimpleED", &pSolo.flexFilterUseSimpleEDStr));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexSimpleEDMinRescues", &pSolo.flexFilterSimpleEDMinRescues));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexSimpleEDMinAmbient", &pSolo.flexFilterSimpleEDMinAmbient));
    parArray.push_back(new ParameterInfoScalar<uint32>(-1, -1, "soloFlexSimpleEDMinCandidates", &pSolo.flexFilterSimpleEDMinCandidates));

    // Debug flags
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexDebugOutputDir", &pSolo.flexFilterDebugOutputDir));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexDebugTagLog", &pSolo.flexFilterDebugTagLogStr));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexDisableOccupancy", &pSolo.flexFilterDisableOccupancyStr));
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexInvariantChecks", &pSolo.flexFilterInvariantChecksStr));

    // Output options
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexKeepCBTag", &pSolo.flexFilterKeepCBTagStr));

    // Minimal memory mode
    parArray.push_back(new ParameterInfoScalar<string>(-1, -1, "soloFlexMinimalMemory", &pSolo.soloFlexMinimalMemoryStr));

    parameterInputName.push_back("Default");
    parameterInputName.push_back("Command-Line-Initial");
    parameterInputName.push_back("Command-Line");
    parameterInputName.push_back("genomeParameters.txt");

};


void Parameters::inputParameters (int argInN, char* argIn[]) {//input parameters: default, from files, from command line
    
    //hard-coded parameters
    runRestart.type=0;

///////// Default parameters

    #include "parametersDefault.xxd"
    string parString( (const char*) parametersDefault,parametersDefault_len);
    stringstream parStream (parString);
    
    scanAllLines(parStream, 0, -1);
    // Safety: if new parameters are missing from the baked-in defaults, seed them here
    for (auto *p : parArray) {
        if (p->nameString == "soloUMICorrectionUseTags" && p->inputLevel < 0) {
            // Conservative default mirrors the legacy behavior (ignore tags)
            pSolo.umiCorrectionUseTagsStr = "no";
            p->inputLevel = 0;
        }
    }
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel<0) {
            ostringstream errOut;
            errOut <<"BUG: DEFAULT parameter value not defined: "<<parArray[ii]->nameString;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

///////// Initial parameters from Command Line

    commandLine="";
    string commandLineFile="";

    if (argInN>1) {//scan parameters from command line
        commandLine += string(argIn[0]);
        for (int iarg=1; iarg<argInN; iarg++) {
            string oneArg=string(argIn[iarg]);

            if (oneArg=="--version") {//print version and exit
                std::cout << STAR_VERSION <<std::endl;
                exit(0);
            };

            size_t found = oneArg.find("=");
            if (found!=string::npos && oneArg.substr(0,2)=="--") {// --parameter=value
                string key = oneArg.substr(2, found - 2);
                string val = oneArg.substr(found + 1);
                if (val.find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
                    val ='\"' + val + '\"';
                };
                commandLineFile += '\n' + key + ' ' + val;
            } else if (oneArg.substr(0,2)=="--") {//parameter name, cut --
                commandLineFile +='\n' + oneArg.substr(2);
            } else {//parameter value
                if (oneArg.find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
                    oneArg ='\"'  + oneArg +'\"';
                };
                commandLineFile +=' ' + oneArg;
            };
            commandLine += ' ' + oneArg;
        };
        istringstream parStreamCommandLine(commandLineFile);
        scanAllLines(parStreamCommandLine, 1, 2); //read only initial Command Line parameters
    };

	createDirectory(outFileNamePrefix, runDirPerm, "--outFileNamePrefix", *this);

    outLogFileName=outFileNamePrefix + "Log.out";
    inOut->logMain.open(outLogFileName.c_str());
    if (inOut->logMain.fail()) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL ERROR: could not create output file: "<<outFileNamePrefix + "Log.out"<<"\n";
        errOut <<"SOLUTION: check if the path " << outFileNamePrefix << " exists and you have permissions to write there\n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    inOut->logMain << "STAR version=" << STAR_VERSION << "\n";
    inOut->logMain << "STAR compilation time,server,dir=" << COMPILATION_TIME_PLACE << "\n";
    inOut->logMain << "STAR git: " << GIT_BRANCH_COMMIT_DIFF << "\n";
    #ifdef COMPILE_FOR_LONG_READS
           inOut->logMain << "Compiled for LONG reads" << "\n";
    #endif

    //define what goes to cout
    if (outStd=="Log") {
        inOut->logStdOut=& std::cout;
    } else if (outStd=="SAM" || outStd=="BAM_Unsorted" || outStd=="BAM_SortedByCoordinate" || outStd=="BAM_Quant") {
        inOut->logStdOutFile.open((outFileNamePrefix + "Log.std.out").c_str());
        inOut->logStdOut= & inOut->logStdOutFile;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL PARAMETER error: outStd="<<outStd <<" is not a valid value of the parameter\n";
        errOut <<"SOLUTION: provide a valid value fot outStd: Log / SAM / BAM_Unsorted / BAM_SortedByCoordinate";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    /*
    inOut->logMain << "##### DEFAULT parameters:\n" <<flush;
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };
    */
    
    inOut->logMain <<"##### Command Line:\n"<<commandLine <<endl ;

    inOut->logMain << "##### Initial USER parameters from Command Line:\n";
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==1) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };

///////// Parameters files

    if (parametersFiles.at(0) != "-") {//read parameters from a user-defined file
        for (uint ii=0; ii<parametersFiles.size(); ii++) {
            parameterInputName.push_back(parametersFiles.at(ii));
            ifstream parFile(parametersFiles.at(ii).c_str());
            if (parFile.good()) {
                inOut->logMain << "##### USER parameters from user-defined parameters file " <<parametersFiles.at(ii)<< ":\n" <<flush;
                scanAllLines(parFile, parameterInputName.size()-1, -1);
                parFile.close();
            } else {
                ostringstream errOut;
                errOut <<"EXITING because of fatal input ERROR: could not open user-defined parameters file " <<parametersFiles.at(ii)<< "\n" <<flush;
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };

///////// Command Line Final

    if (argInN>1) {//scan all parameters from command line and override previous values
        inOut->logMain << "###### All USER parameters from Command Line:\n" <<flush;
        istringstream parStreamCommandLine(commandLineFile);
        scanAllLines(parStreamCommandLine, 2, -1);
    };

    inOut->logMain << "##### Finished reading parameters from all sources\n\n" << flush;

    inOut->logMain << "##### Final user re-defined parameters-----------------:\n" << flush;

    ostringstream clFull;
    clFull << argIn[0];
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel>0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
            if (parArray[ii]->nameString != "parametersFiles" ) {
                clFull << "   --" << parArray[ii]->nameString << " " << *(parArray[ii]);
            };
        };
    };
    commandLineFull=clFull.str();
    inOut->logMain << "\n-------------------------------\n##### Final effective command line:\n" <<  clFull.str() << "\n";

    /*
    //     parOut.close();
    inOut->logMain << "\n##### Final parameters after user input--------------------------------:\n" << flush;
    //     parOut.open("Parameters.all.out");
    for (uint ii=0; ii<parArray.size(); ii++) {
        inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
    };
    //     parOut.close();
    */
    inOut->logMain << "----------------------------------------\n\n" << flush;


    ///////////////////////////////////////// Old variables
    //splitting
    maxNsplit=10;


////////////////////////////////////////////////////// Calculate and check parameters
    iReadAll=0;
    g_bamRecordIndex = 0;
    
    pGe.initialize(this);

    //directory permissions TODO: this needs to be done before outPrefixFileName is created
    if (runDirPermIn=="User_RWX") {
        runDirPerm=S_IRWXU;
    } else if (runDirPermIn=="All_RWX") {
        runDirPerm= S_IRWXU | S_IRWXG | S_IRWXO;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --runDirPerm=" << runDirPerm << "\n";
        errOut << "SOLUTION: use one of the allowed values of --runDirPerm : 'User_RWX' or 'All_RWX' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outTmpDir=="-") {
        outFileTmp=outFileNamePrefix +"_STARtmp/";
        if (runRestart.type!=1)
            sysRemoveDir (outFileTmp);
    } else {
        outFileTmp=outTmpDir + "/";
    };

    if (mkdir (outFileTmp.c_str(),runDirPerm)!=0 && runRestart.type!=1) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal ERROR: could not make temporary directory: "<< outFileTmp<<"\n";
        errOut <<"SOLUTION: (i) please check the path and writing permissions \n (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //threaded or not
    g_threadChunks.threadBool=(runThreadN>1);

    //wigOut parameters
    if (outWigType.at(0)=="None") {
        outWigFlags.yes=false;
    } else if (outWigType.at(0)=="bedGraph") {
        outWigFlags.yes=true;
        outWigFlags.format=0;
    } else if (outWigType.at(0)=="wiggle") {
        outWigFlags.yes=true;
        outWigFlags.format=1;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigType=" << outWigType.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigType : 'None' or 'bedGraph' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    if (outWigStrand.at(0)=="Stranded") {
        outWigFlags.strand=true;
    } else if (outWigStrand.at(0)=="Unstranded") {
        outWigFlags.strand=false;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --outWigStrand=" << outWigStrand.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigStrand : 'Stranded' or 'Unstranded' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outWigType.size()==1) {//simple bedGraph
        outWigFlags.type=0;
    } else {
        if (outWigType.at(1)=="read1_5p") {
            outWigFlags.type=1;
        } else if (outWigType.at(1)=="read2") {
            outWigFlags.type=2;
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: unrecognized second option in --outWigType=" << outWigType.at(1) << "\n";
            errOut << "SOLUTION: use one of the allowed values of --outWigType : 'read1_5p' \n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    //wigOut parameters
    if (outWigNorm.at(0)=="None") {
        outWigFlags.norm=0;
    } else if (outWigNorm.at(0)=="RPM") {
        outWigFlags.norm=1;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal parameter ERROR: unrecognized option in --outWigNorm=" << outWigNorm.at(0) << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigNorm : 'None' or 'RPM' \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };


    //remove duplicates parameters
    if (removeDuplicates.mode=="UniqueIdentical")
    {
        removeDuplicates.yes=true;
        removeDuplicates.markMulti=true;
    } else if (removeDuplicates.mode=="UniqueIdenticalNotMulti")
    {
        removeDuplicates.yes=true;
        removeDuplicates.markMulti=false;
    } else if (removeDuplicates.mode!="-")
    {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --bamRemoveDuplicatesType="<<removeDuplicates.mode<<"\n";
            errOut << "SOLUTION: use allowed option: - or UniqueIdentical or UniqueIdenticalNotMulti";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    runMode=runModeIn[0];
    if (runMode=="alignReads") {
        inOut->logProgress.open((outFileNamePrefix + "Log.progress.out").c_str());
    } else if (runMode=="inputAlignmentsFromBAM") {
        //at the moment, only wiggle output is implemented
        if (outWigFlags.yes) {
            *inOut->logStdOut << timeMonthDayTime() << " ..... reading from BAM, output wiggle\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... reading from BAM, output wiggle\n" <<flush;
            string wigOutFileNamePrefix=outFileNamePrefix + "Signal";
            signalFromBAM(inputBAMfile, wigOutFileNamePrefix, *this);
            *inOut->logStdOut << timeMonthDayTime() << " ..... done\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... done\n" <<flush;
        } else if (removeDuplicates.mode!="-") {
            *inOut->logStdOut << timeMonthDayTime() << " ..... reading from BAM, remove duplicates, output BAM\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... reading from BAM, remove duplicates, output BAM\n" <<flush;
            bamRemoveDuplicates(inputBAMfile, (outFileNamePrefix+"Processed.out.bam").c_str(), *this);
            *inOut->logStdOut << timeMonthDayTime() << " ..... done\n" <<flush;
            inOut->logMain << timeMonthDayTime()    << " ..... done\n" <<flush;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of fatal INPUT ERROR: at the moment --runMode inputFromBAM only works with --outWigType bedGraph OR --bamRemoveDuplicatesType Identical"<<"\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        sysRemoveDir (outFileTmp);
        exit(0);
    };

    outSAMbool=false;
    outBAMunsorted=false;
    outBAMunsortedUseSoloTmp=false;
    outBAMcoord=false;
    emitNoYBAMyes=false;
    emitYReadNamesyes=false;
    keepBAMyes=false;
    emitYNoYFastqyes=false;
    if (runMode=="alignReads" && outSAMmode != "None") {//open SAM file and write header
        if (outSAMtype.at(0)=="BAM") {
            if (outSAMtype.size()<2) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETER error: missing BAM option\n";
                errOut <<"SOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted OR SortedByCoordinate OR both\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            for (uint32 ii=1; ii<outSAMtype.size(); ii++) {
                if (outSAMtype.at(ii)=="Unsorted") {
                    outBAMunsorted=true;
                } else if (outSAMtype.at(ii)=="SortedByCoordinate") {
                    outBAMcoord=true;
                } else {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal input ERROR: unknown value for the word " <<ii+1<<" of outSAMtype: "<< outSAMtype.at(ii) <<"\n";
                    errOut <<"SOLUTION: re-run STAR with one of the allowed values of --outSAMtype BAM Unsorted or SortedByCoordinate or both\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                };
            };
            //TODO check for conflicts
            if (outBAMunsorted) {
                if (outStd=="BAM_Unsorted") {
                    outBAMfileUnsortedName="-";
                } else {
                    outBAMfileUnsortedName=outFileNamePrefix + "Aligned.out.bam";
                };
                
                // BAM opening logic moved after Solo initialization (see after pSolo.initialize)
            };
            if (outBAMcoord) {
                if (outStd=="BAM_SortedByCoordinate") {
                    outBAMfileCoordName="-";
                } else {
                    outBAMfileCoordName=outFileNamePrefix + "Aligned.sortedByCoord.out.bam";
                };
                inOut->outBAMfileCoord = bgzf_open(outBAMfileCoordName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
                if (outBAMsortingThreadN==0) {
                    outBAMsortingThreadNactual=min(6, runThreadN);
                } else {
                    outBAMsortingThreadNactual=outBAMsortingThreadN;
                };
                outBAMcoordNbins=max((uint32)outBAMsortingThreadNactual*3,outBAMsortingBinsN);
                outBAMsortingBinStart= new uint64 [outBAMcoordNbins];
                outBAMsortingBinStart[0]=1;//this initial value means that the bin sizes have not been determined yet

                outBAMsortTmpDir=outFileTmp+"/BAMsort/";
                mkdir(outBAMsortTmpDir.c_str(),runDirPerm);
                
                // Validate outBAMsortMethod
                {
                    string t = outBAMsortMethod; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
                    if (t.empty() || t == "star") {
                        outBAMsortMethod = "star";
                    } else if (t == "samtools") {
                        outBAMsortMethod = "samtools";
                    } else {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL PARAMETER ERROR: --outBAMsortMethod must be 'star' or 'samtools', got: " << outBAMsortMethod << "\n";
                        errOut << "SOLUTION: use --outBAMsortMethod star or --outBAMsortMethod samtools";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    }
                }
            };
            
            // Parse Y-chromosome BAM split parameters
            {
                string t = emitNoYBAM; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
                if (t == "yes") emitNoYBAMyes = true;
                else if (t == "no" || t.empty()) emitNoYBAMyes = false;
                else {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --emitNoYBAM=" << emitNoYBAM << "\n";
                    errOut << "SOLUTION: use allowed option: yes OR no\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }
            {
                string t = keepBAM; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
                if (t == "yes") keepBAMyes = true;
                else if (t == "no" || t.empty()) keepBAMyes = false;
                else {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --keepBAM=" << keepBAM << "\n";
                    errOut << "SOLUTION: use allowed option: yes OR no\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }

            // Derive Y/noY BAM output paths
            if (emitNoYBAMyes) {
                // Handle user-specified override paths
                if (!noYOutput.empty() && noYOutput != "-") {
                    outBAMfileNoYName = noYOutput;
                } else if (!YOutput.empty() && YOutput != "-") {
                    // If only YOutput specified, derive noY from it
                    size_t pos = YOutput.find("_Y.bam");
                    if (pos != string::npos) {
                        outBAMfileNoYName = YOutput.substr(0, pos) + "_noY.bam";
                    } else {
                        outBAMfileNoYName = outFileNamePrefix + "Aligned.out_noY.bam";
                    }
                } else {
                    // Default: derive from primary BAM name
                    if (outBAMunsorted) {
                        if (outBAMfileUnsortedName == "-" || outBAMfileUnsortedName == "/dev/stdout") {
                            outBAMfileNoYName = "star_output_noY.bam";
                            outBAMfileYName = "star_output_Y.bam";
                        } else {
                            size_t pos = outBAMfileUnsortedName.find(".bam");
                            if (pos != string::npos) {
                                outBAMfileNoYName = outBAMfileUnsortedName.substr(0, pos) + "_noY.bam";
                                outBAMfileYName = outBAMfileUnsortedName.substr(0, pos) + "_Y.bam";
                            } else {
                                outBAMfileNoYName = outBAMfileUnsortedName + "_noY.bam";
                                outBAMfileYName = outBAMfileUnsortedName + "_Y.bam";
                            }
                        }
                    } else if (outBAMcoord) {
                        if (outBAMfileCoordName == "-" || outBAMfileCoordName == "/dev/stdout") {
                            outBAMfileNoYName = "star_output_noY.bam";
                            outBAMfileYName = "star_output_Y.bam";
                        } else {
                            size_t pos = outBAMfileCoordName.find(".bam");
                            if (pos != string::npos) {
                                outBAMfileNoYName = outBAMfileCoordName.substr(0, pos) + "_noY.bam";
                                outBAMfileYName = outBAMfileCoordName.substr(0, pos) + "_Y.bam";
                            } else {
                                outBAMfileNoYName = outBAMfileCoordName + "_noY.bam";
                                outBAMfileYName = outBAMfileCoordName + "_Y.bam";
                            }
                        }
                    } else {
                        // Fallback if neither unsorted nor coord is set
                        outBAMfileNoYName = outFileNamePrefix + "Aligned.out_noY.bam";
                        outBAMfileYName = outFileNamePrefix + "Aligned.out_Y.bam";
                    }
                }
                
                // Handle YOutput override
                if (!YOutput.empty() && YOutput != "-") {
                    outBAMfileYName = YOutput;
                }
            }
        } else if (outSAMtype.at(0)=="SAM") {
            if (outSAMtype.size()>1)
            {
                ostringstream errOut;
                errOut <<"EXITING because of fatal PARAMETER error: --outSAMtype SAM can cannot be combined with "<<outSAMtype.at(1)<<" or any other options\n";
                errOut <<"SOLUTION: re-run STAR with with '--outSAMtype SAM' only, or with --outSAMtype BAM Unsorted|SortedByCoordinate\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            outSAMbool=true;
            if (outStd=="SAM") {
                inOut->outSAM = & std::cout;
            } else {
                inOut->outSAMfile.open((outFileNamePrefix + "Aligned.out.sam").c_str());
                inOut->outSAM = & inOut->outSAMfile;
            };
        } else if (outSAMtype.at(0)=="None") {
            //nothing to do, all flags are already false
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of fatal input ERROR: unknown value for the first word of outSAMtype: "<< outSAMtype.at(0) <<"\n";
            errOut <<"SOLUTION: re-run STAR with one of the allowed values of outSAMtype: BAM or SAM \n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    {
        string t = emitYReadNames; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        if (t == "yes") emitYReadNamesyes = true;
        else if (t == "no" || t.empty()) emitYReadNamesyes = false;
        else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --emitYReadNames=" << emitYReadNames << "\n";
            errOut << "SOLUTION: use allowed option: yes OR no\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    }

    // Derive Y read names output path
    if (emitYReadNamesyes) {
        if (!YReadNamesOutput.empty() && YReadNamesOutput != "-") {
            outYReadNamesFile = YReadNamesOutput;
        } else {
            outYReadNamesFile = outFileNamePrefix + "Aligned.out_Y.names.txt";
        }
    }
    
    // Parse Y-chromosome FASTQ emission parameters
    {
        string t = emitYNoYFastq; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        if (t == "yes") emitYNoYFastqyes = true;
        else if (t == "no" || t.empty()) emitYNoYFastqyes = false;
        else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --emitYNoYFastq=" << emitYNoYFastq << "\n";
            errOut << "SOLUTION: use allowed option: yes OR no\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    }
    {
        string t = emitYNoYFastqCompression; std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        if (t.empty() || t == "gz") {
            emitYNoYFastqCompression = "gz";
        } else if (t == "none") {
            emitYNoYFastqCompression = "none";
        } else {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --emitYNoYFastqCompression=" << emitYNoYFastqCompression << "\n";
            errOut << "SOLUTION: use allowed option: gz OR none\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    }
    
    if (!outBAMcoord && outWigFlags.yes && runMode=="alignReads") {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: generating signal with --outWigType requires sorted BAM\n";
        errOut <<"SOLUTION: re-run STAR with with --outSAMtype BAM SortedByCoordinate, or, id you also need unsroted BAM, with --outSAMtype BAM SortedByCoordinate Unsorted\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //versions
    for (uint ii=0;ii<1;ii++) {
        if (parArray[ii]->inputLevel>0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal input ERROR: the version parameter "<< parArray[ii]->nameString << " cannot be re-defined by the user\n";
            errOut <<"SOLUTION: please remove this parameter from the command line or input files and re-start STAR\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    //run
    if (runThreadN<=0) {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: runThreadN must be >0, user-defined runThreadN="<<runThreadN<<"\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //
    if (outFilterType=="BySJout" && outSAMorder=="PairedKeepInputOrder") {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: --outFilterType=BySJout is not presently compatible with --outSAMorder=PairedKeepInputOrder\n";
        errOut <<"SOLUTION: re-run STAR without setting one of those parameters.\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    if (!outSAMbool && outSAMorder=="PairedKeepInputOrder") {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: --outSAMorder=PairedKeepInputOrder is presently only compatible with SAM output, i.e. default --outSMAtype SAM\n";
        errOut <<"SOLUTION: re-run STAR without --outSAMorder=PairedKeepInputOrder, or with --outSAMorder=PairedKeepInputOrder --outSMAtype SAM .\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    //SJ filtering
    for (int ii=0;ii<4;ii++) {
        if (outSJfilterOverhangMin.at(ii)<0) outSJfilterOverhangMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterCountUniqueMin.at(ii)<0) outSJfilterCountUniqueMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterCountTotalMin.at(ii)<0) outSJfilterCountTotalMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterDistToOtherSJmin.at(ii)<0) outSJfilterDistToOtherSJmin.at(ii)=numeric_limits<int32>::max();

        if (alignSJstitchMismatchNmax.at(ii)<0) alignSJstitchMismatchNmax.at(ii)=numeric_limits<int32>::max();
    };

    if (limitGenomeGenerateRAM==0) {//must be >0
        inOut->logMain <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=0\n";
        inOut->logMain <<"SOLUTION: please specify a >0 value for limitGenomeGenerateRAM\n"<<flush;
        exit(1);
    } else if (limitGenomeGenerateRAM>1000000000000) {//
        inOut->logMain <<"WARNING: specified limitGenomeGenerateRAM="<<limitGenomeGenerateRAM<<" bytes appears to be too large, if you do not have enough memory the code will crash!\n"<<flush;
    };

    outSAMfilter.KeepOnlyAddedReferences=false;
    outSAMfilter.KeepAllAddedReferences=false;
    outSAMfilter.yes=true;
    if (outSAMfilter.mode.at(0)=="KeepOnlyAddedReferences") {
        outSAMfilter.KeepOnlyAddedReferences=true;
    } else if (outSAMfilter.mode.at(0)=="KeepAllAddedReferences") {
        outSAMfilter.KeepAllAddedReferences=true;
    } else if (outSAMfilter.mode.at(0)=="None") {
      outSAMfilter.yes=false;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outSAMfilter: "<<outSAMfilter.mode.at(0) <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: KeepOnlyAddedReferences or None\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if ( (outSAMfilter.KeepOnlyAddedReferences || outSAMfilter.KeepAllAddedReferences) && pGe.gFastaFiles.at(0)=="-" ) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMfilter KeepOnlyAddedReferences OR KeepAllAddedReferences options can only be used if references are added on-the-fly with --genomeFastaFiles" <<"\n";
        errOut <<"SOLUTION: use default --outSAMfilter None, OR add references with --genomeFataFiles\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };


    if (outMultimapperOrder.mode=="Old_2.4") {
        outMultimapperOrder.random=false;
    } else if (outMultimapperOrder.mode=="Random") {
        outMultimapperOrder.random=true;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --outMultimapperOrder: "<<outMultimapperOrder.mode <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Old_2.4 or SortedByCoordinate or Random\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //read parameters
    readFilesInit();

    //two-pass
    if (parArray.at(twoPass.pass1readsN_par)->inputLevel>0  && twoPass.mode=="None") {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: --twopass1readsN is defined, but --twoPassMode is not defined\n";
        errOut << "SOLUTION: to activate the 2-pass mode, use --twopassMode Basic";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    twoPass.yes=false;
    twoPass.pass2=false;
    if (twoPass.mode!="None") {//2-pass parameters
        if (runMode!="alignReads") {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: 2-pass mapping option  can only be used with --runMode alignReads\n";
            errOut << "SOLUTION: remove --twopassMode option";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (twoPass.mode!="Basic") {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --twopassMode="<<twoPass.mode<<"\n";
            errOut << "SOLUTION: for the 2-pass mode, use allowed values --twopassMode: Basic";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (twoPass.pass1readsN==0) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: --twopass1readsN = 0 in the 2-pass mode\n";
            errOut << "SOLUTION: for the 2-pass mode, specify --twopass1readsN > 0. Use a very large number or -1 to map all reads in the 1st pass.\n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if (pGe.gLoad!="NoSharedMemory") {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: 2-pass method is not compatible with --genomeLoad "<<pGe.gLoad<<"\n";
            errOut << "SOLUTION: re-run STAR with --genomeLoad NoSharedMemory ; this is the only option compatible with --twopassMode Basic .\n";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        twoPass.yes=true;
        twoPass.dir=outFileNamePrefix+"_STARpass1/";
        sysRemoveDir (twoPass.dir);
        if (mkdir (twoPass.dir.c_str(),runDirPerm)!=0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not make pass1 directory: "<< twoPass.dir<<"\n";
            errOut <<"SOLUTION: please check the path and writing permissions \n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    // openReadFiles depends on twoPass for reading SAM header
    if (runMode=="alignReads" && pGe.gLoad!="Remove" && pGe.gLoad!="LoadAndExit") {//open reads files to check if they are present
        openReadsFiles();

        if (readNends > 2 && pSolo.typeStr=="None") {//could have >2 mates only for Solo
            ostringstream errOut;
            errOut <<"EXITING: because of fatal input ERROR: number of read mates files > 2: " <<readNends << "\n";
            errOut <<"SOLUTION:specify only one or two files in the --readFilesIn option. If file names contain spaces, use quotes: \"file name\"\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };

        if ( runMode=="alignReads" && outReadsUnmapped=="Fastx" ) {//open unmapped reads file
            for (uint imate=0;imate<readNends;imate++) {
                ostringstream ff;
                ff << outFileNamePrefix << "Unmapped.out.mate" << imate+1;
                inOut->outUnmappedReadsStream[imate].open(ff.str().c_str());
            };
        };
    };

    if (outSAMmapqUnique<0 || outSAMmapqUnique>255) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL input ERROR: out of range value for outSAMmapqUnique=" << outSAMmapqUnique <<"\n";
            errOut <<"SOLUTION: specify outSAMmapqUnique within the range of 0 to 255\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
        
    //variation
    var.yes=false;
    if (var.vcfFile!="-") {
        var.yes=true;
    };

    //WASP
    wasp.yes=false;
    wasp.SAMtag=false;
    if (wasp.outputMode=="SAMtag") {
        wasp.yes=true;
        wasp.SAMtag=true;
        var.heteroOnly=true;
    } else if (wasp.outputMode=="None") {
        //nothing to do
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented --waspOutputMode option: "<<wasp.outputMode <<"\n";
        errOut <<"SOLUTION: re-run STAR with allowed --waspOutputMode options: None or SAMtag\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (wasp.yes && !var.yes) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --waspOutputMode option requires VCF file: "<<wasp.outputMode <<"\n";
        errOut <<"SOLUTION: re-run STAR with --waspOutputMode ... and --varVCFfile /path/to/file.vcf\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

     if (wasp.yes && outSAMtype.at(0)!="BAM") {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --waspOutputMode requires output to BAM file\n";
        errOut <<"SOLUTION: re-run STAR with --waspOutputMode ... and --outSAMtype BAM ... \n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    //quantification parameters
    quant.yes=false;
    quant.geCount.yes=false;
    quant.trSAM.yes=false;
    quant.trSAM.bamYes=false;
    quant.trSAM.indel=false;
    quant.trSAM.softClip=false;
    quant.trSAM.singleEnd=false; 
    if (quant.mode.at(0) != "-") {
        quant.yes=true;
        for (uint32 ii=0; ii<quant.mode.size(); ii++) {
            if (quant.mode.at(ii)=="TranscriptomeSAM") {
                quant.trSAM.yes=true;

                if (quant.trSAM.bamCompression>-2)
                    quant.trSAM.bamYes=true;

                if (quant.trSAM.bamYes) {
                    if (outStd=="BAM_Quant") {
                        outFileNamePrefix="-";
                    } else {
                        outQuantBAMfileName=outFileNamePrefix + "Aligned.toTranscriptome.out.bam";
                    };
                    inOut->outQuantBAMfile=bgzf_open(outQuantBAMfileName.c_str(),("w"+to_string((long long) quant.trSAM.bamCompression)).c_str());
                };
                if (quant.trSAM.output=="BanSingleEnd_BanIndels_ExtendSoftclip") {
                    quant.trSAM.indel=false;
                    quant.trSAM.softClip=false;
                    quant.trSAM.singleEnd=false;
                } else if (quant.trSAM.output=="BanSingleEnd") {
                    quant.trSAM.indel=true;
                    quant.trSAM.softClip=true;
                    quant.trSAM.singleEnd=false;
                } else if (quant.trSAM.output=="BanSingleEnd_ExtendSoftclip") {
                    quant.trSAM.indel=true;
                    quant.trSAM.softClip=false;
                    quant.trSAM.singleEnd=false;
                };
            } else if  (quant.mode.at(ii)=="GeneCounts") {
                quant.geCount.yes=true;
                quant.geCount.outFile=outFileNamePrefix + "ReadsPerGene.out.tab";
            } else if (quant.mode.at(ii)=="TranscriptVB") {
                quant.transcriptVB.yes=true;
                quant.transcriptVB.outFile=outFileNamePrefix + "quant.sf";
            } else {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: unrecognized option in --quantMode=" << quant.mode.at(ii) << "\n";
                errOut << "SOLUTION: use one of the allowed values of --quantMode : TranscriptomeSAM, GeneCounts, TranscriptVB, or - .\n";
                exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };
    //these may be set in STARsolo or in SAM attributes
    quant.geneFull.yes=false;
    quant.gene.yes=false;

    // SLAM quantification (independent of quantMode)
    quant.slam.yes = (quant.slam.modeInt != 0);
    if (quant.slam.yes) {
        quant.yes = true;
        if (quant.slam.outFile.empty() || quant.slam.outFile == "-") {
            quant.slam.outFile = outFileNamePrefix + "SlamQuant.out";
        }
        if (quant.slam.errorRate <= 0.0 || quant.slam.errorRate >= 1.0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamErrorRate must be in (0,1)\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        if (quant.slam.convRate <= 0.0 || quant.slam.convRate >= 1.0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamConvRate must be in (0,1)\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        if (quant.slam.snpBed == "-" || quant.slam.snpBed == "None") {
            quant.slam.snpBed.clear();
        }
        quant.slam.snpDetect = (quant.slam.snpDetectInt != 0);
        if (quant.slam.snpDetect && !quant.slam.snpBed.empty()) {
            inOut->logMain << "WARNING: --slamSnpDetect ignored because --slamSnpBed is set.\n";
            quant.slam.snpDetect = false;
        }
        string strandLower = quant.slam.strandnessStr;
        if (strandLower.empty() || strandLower == "-") {
            strandLower = "unspecific";
        }
        for (auto& c : strandLower) c = std::tolower(c);
        if (strandLower == "sense" || strandLower == "forward") {
            quant.slam.strandness = 1;
            quant.slam.strandnessStr = "Sense";
        } else if (strandLower == "antisense" || strandLower == "reverse") {
            quant.slam.strandness = 2;
            quant.slam.strandnessStr = "Antisense";
        } else if (strandLower == "unspecific" || strandLower == "unstranded" ||
                   strandLower == "auto" || strandLower == "autodetect") {
            quant.slam.strandness = 0;
            quant.slam.strandnessStr = "Unspecific";
            if (strandLower == "auto" || strandLower == "autodetect") {
                inOut->logMain << "WARNING: --slamStrandness AutoDetect treated as Unspecific; "
                               << "STAR-SLAM does not auto-detect strandness.\n";
            }
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamStrandness must be Unspecific, Sense, or Antisense\n"
                   << "Got: " << quant.slam.strandnessStr << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        string weightModeLower = quant.slam.weightModeStr;
        if (weightModeLower.empty() || weightModeLower == "-") {
            weightModeLower = "alignments";
        }
        for (auto& c : weightModeLower) c = std::tolower(c);
        if (weightModeLower == "alignments" || weightModeLower == "nh" || weightModeLower == "ntr" ||
            weightModeLower == "weight") {
            quant.slam.weightMode = 0;
            quant.slam.weightModeStr = "Alignments";
        } else if (weightModeLower == "uniform" || weightModeLower == "none" || weightModeLower == "all") {
            quant.slam.weightMode = 1;
            quant.slam.weightModeStr = "Uniform";
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamWeightMode must be Alignments or Uniform\n"
                   << "Got: " << quant.slam.weightModeStr << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        string dumpWeightsModeLower = quant.slam.dumpWeightsModeStr;
        if (dumpWeightsModeLower.empty() || dumpWeightsModeLower == "-") {
            dumpWeightsModeLower = "dump";
        }
        for (auto& c : dumpWeightsModeLower) c = std::tolower(c);
        if (dumpWeightsModeLower == "dump" || dumpWeightsModeLower == "default") {
            quant.slam.dumpWeightsMode = 0;
            quant.slam.dumpWeightsModeStr = "dump";
        } else if (dumpWeightsModeLower == "vbgene" || dumpWeightsModeLower == "vb" ||
                   dumpWeightsModeLower == "gene") {
            quant.slam.dumpWeightsMode = 1;
            quant.slam.dumpWeightsModeStr = "vbGene";
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamDumpWeightsMode must be dump or vbGene\n"
                   << "Got: " << quant.slam.dumpWeightsModeStr << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        auto trimLine = [](const std::string& input) -> std::string {
            size_t start = input.find_first_not_of(" \t\r\n");
            if (start == std::string::npos) {
                return "";
            }
            size_t end = input.find_last_not_of(" \t\r\n");
            return input.substr(start, end - start + 1);
        };
        auto loadListFile = [&](const std::string& path, std::unordered_set<std::string>& out, bool stripAt) {
            std::ifstream in(path.c_str());
            if (!in.good()) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "cannot open " << path << "\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            std::string line;
            while (std::getline(in, line)) {
                line = trimLine(line);
                if (line.empty() || line[0] == '#') {
                    continue;
                }
                if (stripAt && line[0] == '@') {
                    line.erase(0, 1);
                }
                if (!line.empty()) {
                    out.insert(line);
                }
            }
        };
        if (quant.slam.debugMaxReads < 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamDebugMaxReads must be >= 0\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        if (!quant.slam.debugGeneList.empty() && quant.slam.debugGeneList != "-" &&
            quant.slam.debugGeneList != "None" && quant.slam.debugGeneList != "none") {
            loadListFile(quant.slam.debugGeneList, quant.slam.debugGenes, false);
        }
        if (!quant.slam.debugReadList.empty() && quant.slam.debugReadList != "-" &&
            quant.slam.debugReadList != "None" && quant.slam.debugReadList != "none") {
            loadListFile(quant.slam.debugReadList, quant.slam.debugReads, true);
        }
        if (quant.slam.debugSnpWindow < 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamDebugSnpWindow must be >= 0\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }

        bool debugSnpEnabled = !quant.slam.debugSnpLoc.empty() && quant.slam.debugSnpLoc != "-" &&
                               quant.slam.debugSnpLoc != "None" && quant.slam.debugSnpLoc != "none";

        quant.slam.debugEnabled = !quant.slam.debugGenes.empty() || !quant.slam.debugReads.empty() || debugSnpEnabled;
        if (quant.slam.debugEnabled) {
            if (quant.slam.debugOutPrefix.empty() || quant.slam.debugOutPrefix == "-") {
                quant.slam.debugOutPrefix = outFileNamePrefix + "SlamQuant.debug";
            }
        } else {
            quant.slam.debugOutPrefix.clear();
        }
        
        // Parse compatibility mode
        string compatModeLower = quant.slam.compatModeStr;
        if (compatModeLower.empty() || compatModeLower == "-") {
            compatModeLower = "none";
        }
        for (auto& c : compatModeLower) c = std::tolower(c);
        if (compatModeLower == "gedi") {
            quant.slam.compatIntronic = true;
            quant.slam.compatLenientOverlap = true;
            quant.slam.compatOverlapWeight = true;
            quant.slam.compatIgnoreOverlap = false;  // NOT enabled by default in gedi mode
        } else if (compatModeLower != "none") {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--slamCompatMode must be none or gedi\n"
                   << "Got: " << quant.slam.compatModeStr << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        
        // Apply SNP mask defaults and compatibility mode (--slamSnpMaskCompat gedi)
        // Sentinel value (-1 for int, -2.0 for double) indicates the parameter was not set by user.
        // In GEDI compat mode, unset parameters get GEDI-like values.
        // Otherwise, unset parameters get STAR stringent defaults.
        {
            bool isGediCompat = false;
            if (!quant.slamSnpMask.compat.empty()) {
                string maskCompatLower = quant.slamSnpMask.compat;
                for (auto& c : maskCompatLower) c = std::tolower(c);
                if (maskCompatLower == "gedi") {
                    isGediCompat = true;
                } else if (maskCompatLower != "none" && maskCompatLower != "-") {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL PARAMETER ERROR: "
                           << "--slamSnpMaskCompat must be 'gedi' or empty\n"
                           << "Got: " << quant.slamSnpMask.compat << "\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }
            
            // Track settings for logging
            vector<string> appliedSettings;
            vector<string> userOverrides;
            
            // STAR stringent defaults vs GEDI permissive defaults
            // pval: STAR=0.001, GEDI=0.001 (same) - no sentinel needed
            if (quant.slamSnpMask.pval < 0.0) {  // sentinel or invalid
                quant.slamSnpMask.pval = 0.001;
            }
            // Only log as user override if different from default
            if (quant.slamSnpMask.pval != 0.001) {
                userOverrides.push_back("pval=" + to_string(quant.slamSnpMask.pval));
            }
            
            // minTcRatio: STAR=0.3, GEDI=0.3 (same) - no sentinel needed
            if (quant.slamSnpMask.minTcRatio < 0.0) {  // sentinel or invalid
                quant.slamSnpMask.minTcRatio = 0.3;
            }
            // Only log as user override if different from default
            if (quant.slamSnpMask.minTcRatio != 0.3) {
                userOverrides.push_back("minTcRatio=" + to_string(quant.slamSnpMask.minTcRatio));
            }

            // kMode: STAR=conv, GEDI=any (GEDI counts any mismatch as "alt")
            {
                std::string km = quant.slamSnpMask.kMode;
                if (km.empty() || km == "-" || km == "None" || km == "none") {
                    km = isGediCompat ? "any" : "conv";
                    quant.slamSnpMask.kMode = km;
                    appliedSettings.push_back("kMode=" + km + (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
                } else {
                    for (auto& c : km) c = std::tolower(c);
                    if (km != "conv" && km != "any") {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL PARAMETER ERROR: "
                               << "--slamSnpMaskKMode must be 'conv' or 'any'\n"
                               << "Got: " << quant.slamSnpMask.kMode << "\n";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    }
                    quant.slamSnpMask.kMode = km;
                    userOverrides.push_back("kMode=" + km);
                }
            }
            
            // minCov: STAR=20, GEDI=6
            if (quant.slamSnpMask.minCov < 0) {  // sentinel (-1)
                quant.slamSnpMask.minCov = isGediCompat ? 6 : 20;
                appliedSettings.push_back("minCov=" + to_string(quant.slamSnpMask.minCov) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("minCov=" + to_string(quant.slamSnpMask.minCov));
            }
            
            // minAlt: STAR=3, GEDI=1
            if (quant.slamSnpMask.minAlt < 0) {  // sentinel (-1)
                quant.slamSnpMask.minAlt = isGediCompat ? 1 : 3;
                appliedSettings.push_back("minAlt=" + to_string(quant.slamSnpMask.minAlt) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("minAlt=" + to_string(quant.slamSnpMask.minAlt));
            }
            
            // minMapQ: STAR=20, GEDI=0
            if (quant.slamSnpMask.minMapQ < 0) {  // sentinel (-1)
                quant.slamSnpMask.minMapQ = isGediCompat ? 0 : 20;
                appliedSettings.push_back("minMapQ=" + to_string(quant.slamSnpMask.minMapQ) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("minMapQ=" + to_string(quant.slamSnpMask.minMapQ));
            }
            
            // minBaseQ: STAR=20, GEDI=0
            if (quant.slamSnpMask.minBaseQ < 0) {  // sentinel (-1)
                quant.slamSnpMask.minBaseQ = isGediCompat ? 0 : 20;
                appliedSettings.push_back("minBaseQ=" + to_string(quant.slamSnpMask.minBaseQ) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("minBaseQ=" + to_string(quant.slamSnpMask.minBaseQ));
            }
            
            // junctionFlank: STAR=6, GEDI=0
            if (quant.slamSnpMask.junctionFlank < 0) {  // sentinel (-1)
                quant.slamSnpMask.junctionFlank = isGediCompat ? 0 : 6;
                appliedSettings.push_back("junctionFlank=" + to_string(quant.slamSnpMask.junctionFlank) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("junctionFlank=" + to_string(quant.slamSnpMask.junctionFlank));
            }
            
            // indelFlank: STAR=3, GEDI=0
            if (quant.slamSnpMask.indelFlank < 0) {  // sentinel (-1)
                quant.slamSnpMask.indelFlank = isGediCompat ? 0 : 3;
                appliedSettings.push_back("indelFlank=" + to_string(quant.slamSnpMask.indelFlank) + 
                    (isGediCompat ? " (GEDI compat)" : " (STAR default)"));
            } else {
                userOverrides.push_back("indelFlank=" + to_string(quant.slamSnpMask.indelFlank));
            }
            
            // err: -1 means use computed snp_err_used (same for both modes)
            if (quant.slamSnpMask.err >= 0.0) {
                userOverrides.push_back("err=" + to_string(quant.slamSnpMask.err));
            } else {
                appliedSettings.push_back("err=<snp_err_used>");
            }
            
            // Log the settings
            if (isGediCompat) {
                inOut->logMain << "SNP mask GEDI compat mode enabled:\n";
            }
            for (const auto& s : appliedSettings) {
                inOut->logMain << "    " << s << "\n";
            }
            if (!userOverrides.empty()) {
                inOut->logMain << "  User overrides:\n";
                for (const auto& s : userOverrides) {
                    inOut->logMain << "    " << s << "\n";
                }
            }
        }
        
        // Validate SNP mask build parameters
        bool hasMaskIn = !quant.slamSnpMask.maskIn.empty() && quant.slamSnpMask.maskIn != "-" && quant.slamSnpMask.maskIn != "None";
        bool hasVcfIn = !quant.slamSnpMask.vcfIn.empty() && quant.slamSnpMask.vcfIn != "-" && quant.slamSnpMask.vcfIn != "None";
        bool hasBuildFastqs = !quant.slamSnpMask.buildFastqsFofn.empty() && quant.slamSnpMask.buildFastqsFofn != "-" && quant.slamSnpMask.buildFastqsFofn != "None";
        
        if (hasMaskIn && hasVcfIn) {
            inOut->logMain << "WARNING: --slamSnpMaskIn takes precedence over --slamSnpMaskVcfIn. "
                           << "Will load existing mask and skip VCF load.\n";
        }
        if (hasMaskIn && hasBuildFastqs) {
            inOut->logMain << "WARNING: --slamSnpMaskIn takes precedence over --slamSnpMaskBuildFastqs. "
                           << "Will load existing mask and skip build.\n";
        }
        if (hasVcfIn && hasBuildFastqs) {
            inOut->logMain << "WARNING: --slamSnpMaskVcfIn takes precedence over --slamSnpMaskBuildFastqs. "
                           << "Will load VCF mask and skip build.\n";
        }

        if (hasVcfIn) {
            if (quant.slamSnpMask.bedOut.empty()) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskVcfIn requires --slamSnpMaskBedOut\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            string modeLower = quant.slamSnpMask.vcfMode;
            for (auto& c : modeLower) c = std::tolower(c);
            quant.slamSnpMask.vcfMode = modeLower;
            if (quant.slamSnpMask.vcfMode != "gt" && quant.slamSnpMask.vcfMode != "any") {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskVcfMode must be 'gt' or 'any'\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            string filterLower = quant.slamSnpMask.vcfFilter;
            for (auto& c : filterLower) c = std::tolower(c);
            quant.slamSnpMask.vcfFilter = filterLower;
            if (quant.slamSnpMask.vcfFilter != "pass" && quant.slamSnpMask.vcfFilter != "all") {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskVcfFilter must be 'pass' or 'all'\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            if (quant.slamSnpMask.vcfSample == "-" || quant.slamSnpMask.vcfSample == "None" ||
                quant.slamSnpMask.vcfSample == "none") {
                quant.slamSnpMask.vcfSample.clear();
            }
        }
        
        if (hasBuildFastqs && !hasMaskIn && !hasVcfIn) {
            if (quant.slamSnpMask.bedOut.empty()) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskBuildFastqs requires --slamSnpMaskBedOut\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            if (quant.slamSnpMask.minCov == 0) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskMinCov must be > 0\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            // Validate model selection
            if (quant.slamSnpMask.model != "binom" && quant.slamSnpMask.model != "em") {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskModel must be 'binom' or 'em' (got: " << quant.slamSnpMask.model << ")\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            
            // Validate binomial parameters
            if (quant.slamSnpMask.model == "binom") {
                if (quant.slamSnpMask.pval <= 0.0 || quant.slamSnpMask.pval >= 1.0) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL PARAMETER ERROR: "
                           << "--slamSnpMaskPval must be in (0,1)\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
                if (quant.slamSnpMask.minTcRatio < 0.0 || quant.slamSnpMask.minTcRatio > 1.0) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL PARAMETER ERROR: "
                           << "--slamSnpMaskMinTcRatio must be in [0,1]\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }
            
            // Validate EM parameters
            if (quant.slamSnpMask.model == "em") {
                if (quant.slamSnpMask.posterior < 0.0 || quant.slamSnpMask.posterior > 1.0) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL PARAMETER ERROR: "
                           << "--slamSnpMaskPosterior must be in [0,1]\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
            }
            
            if (quant.slamSnpMask.maxIter == 0) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL PARAMETER ERROR: "
                       << "--slamSnpMaskMaxIter must be > 0\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
        }
        
        // Convert buildOnlyInt to bool
        quant.slamSnpMask.buildOnly = (quant.slamSnpMask.buildOnlyInt != 0);
        
        // Apply granular overrides (if explicitly set via int flags, they override mode defaults)
        // Sentinel is -1 (not set), so any value >= 0 is an explicit override (including 0 to disable)
        if (quant.slam.compatIntronicInt >= 0) {
            quant.slam.compatIntronic = (quant.slam.compatIntronicInt != 0);
        }
        if (quant.slam.compatLenientOverlapInt >= 0) {
            quant.slam.compatLenientOverlap = (quant.slam.compatLenientOverlapInt != 0);
        }
        if (quant.slam.compatOverlapWeightInt >= 0) {
            quant.slam.compatOverlapWeight = (quant.slam.compatOverlapWeightInt != 0);
        }
        if (quant.slam.compatIgnoreOverlapInt >= 0) {
            quant.slam.compatIgnoreOverlap = (quant.slam.compatIgnoreOverlapInt != 0);
        }
    }
    
    // Initialize transcriptVB defaults
    if (quant.transcriptVB.yes) {
        // Convert int flags to bool
        quant.transcriptVB.gcBias = (quant.transcriptVB.gcBiasInt != 0);
        // quantVBem flag: if set to 1, use EM (vb=false), otherwise use VB (vb=true, default)
        quant.transcriptVB.vb = (quant.transcriptVB.quantVBemInt == 0);
        quant.transcriptVB.geneOutput = (quant.transcriptVB.geneOutputInt != 0);
        if (quant.transcriptVB.outFile.empty()) {
            quant.transcriptVB.outFile = outFileNamePrefix + "quant.sf";
        }
        quant.transcriptVB.outFileGene = outFileNamePrefix + "quant.genes.sf";
        quant.transcriptVB.outFileGeneTximport = outFileNamePrefix + "quant.genes.tximport.sf";
        
        // Parse genesMode: Legacy or Tximport
        string genesModeLower = quant.transcriptVB.genesModeStr;
        for (auto& c : genesModeLower) c = std::tolower(c);
        if (genesModeLower == "legacy") {
            quant.transcriptVB.genesTximport = false;
        } else if (genesModeLower == "tximport") {
            quant.transcriptVB.genesTximport = true;
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--quantVBgenesMode must be 'Legacy' or 'Tximport'\n"
                   << "Got: " << quant.transcriptVB.genesModeStr << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        
        // Validate autoDetectWindow is positive
        if (quant.transcriptVB.autoDetectWindow <= 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--quantVBAutoDetectWindow must be > 0\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        
        // Validate libType is a known value and normalize to uppercase
        string upper = quant.transcriptVB.libType;
        for (auto& c : upper) c = std::toupper(c);
        if (upper != "A" && upper != "IU" && upper != "ISF" && upper != "ISR" && upper != "U") {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--quantVBLibType must be one of: A, IU, ISF, ISR, U\n"
                   << "Got: " << quant.transcriptVB.libType << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        // Normalize libType to uppercase so libType == "A" works in STAR.cpp
        quant.transcriptVB.libType = upper;
        
        // Normalize errorModelMode to lowercase
        string errorModelModeLower = quant.transcriptVB.errorModelMode;
        for (auto& c : errorModelModeLower) c = std::tolower(c);
        if (errorModelModeLower != "auto" && errorModelModeLower != "cigar" && 
            errorModelModeLower != "as" && errorModelModeLower != "off") {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: "
                   << "--quantVBErrorModel must be one of: auto, cigar, as, off\n"
                   << "Got: " << quant.transcriptVB.errorModelMode << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
        quant.transcriptVB.errorModelMode = errorModelModeLower;

        // Treat "None" and "-" as unset trace output
        if (quant.transcriptVB.traceFile == "-" || quant.transcriptVB.traceFile == "None") {
            quant.transcriptVB.traceFile.clear();
        }
        
        // If auto-detect enabled, validate window vs readMapNumber
        if (quant.transcriptVB.libType == "A") {
            if (readMapNumber > 0 && readMapNumber < (uint64_t)quant.transcriptVB.autoDetectWindow) {
                // User wants fewer reads than detection window - cap the window
                inOut->logMain << "Warning: --readMapNumber (" << readMapNumber 
                               << ") < autoDetectWindow (" << quant.transcriptVB.autoDetectWindow
                               << "). Capping detection window to " << readMapNumber << ".\n";
                quant.transcriptVB.autoDetectWindow = (int)readMapNumber;
                
                // If window is very small, warn about detection reliability
                if (quant.transcriptVB.autoDetectWindow < 100) {
                    inOut->logMain << "Warning: Detection window < 100 reads may be unreliable.\n"
                                   << "Consider using --quantVBLibType IU|ISF|ISR|U explicitly.\n";
                }
            }
        }
    }
    

    outSAMstrandField.type=0; //none
    if (outSAMstrandField.in=="None") {
        outSAMstrandField.type=0;
    } else if (outSAMstrandField.in=="intronMotif") {
        outSAMstrandField.type=1;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT error: unrecognized option in outSAMstrandField=" << outSAMstrandField.in << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outSAMstrandField : None or intronMotif \n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    //SAM attributes
    samAttributes();
    
    //solo
    pSolo.initialize(this);

    // Derive Y/noY FASTQ output paths (after Solo init so readNmates is final)
    if (emitYNoYFastqyes) {
        const bool hasYPrefix = !YFastqOutputPrefix.empty() && YFastqOutputPrefix != "-";
        const bool hasNoYPrefix = !noYFastqOutputPrefix.empty() && noYFastqOutputPrefix != "-";
        const string ext = (emitYNoYFastqCompression == "gz") ? ".fastq.gz" : ".fastq";
        const string outputDir = outputDirFromPrefix(outFileNamePrefix);
        bool useMateFallback = false;

        if (!hasYPrefix || !hasNoYPrefix) {
            if (readFilesN > 1) {
                warningMessage(" emitYNoYFastq: multiple input files detected; output FASTQ names are derived from the first file for each mate",
                               std::cerr, inOut->logMain, *this);
            }
        }

        if (!hasYPrefix || !hasNoYPrefix) {
            for (uint imate = 0; imate < readNmates; imate++) {
                if (readFilesNames.size() <= imate || readFilesNames[imate].empty()) {
                    useMateFallback = true;
                    break;
                }
                string base = pathBasename(readFilesNames[imate][0]);
                if (!hasReadToken(base)) {
                    useMateFallback = true;
                    break;
                }
            }
        }

        for (uint imate = 0; imate < readNmates; imate++) {
            if (hasYPrefix) {
                outYFastqFile[imate] = YFastqOutputPrefix + "mate" + to_string(imate + 1) + ext;
            } else if (!useMateFallback && readFilesNames.size() > imate && !readFilesNames[imate].empty()) {
                string base = pathBasename(readFilesNames[imate][0]);
                string tagged = insertTagBeforeReadToken(base, "_Y");
                if (!tagged.empty()) {
                    outYFastqFile[imate] = outputDir + adjustCompressionExt(tagged, emitYNoYFastqCompression);
                } else {
                    outYFastqFile[imate] = outFileNamePrefix + "Y_reads.mate" + to_string(imate + 1) + ext;
                }
            } else {
                outYFastqFile[imate] = outFileNamePrefix + "Y_reads.mate" + to_string(imate + 1) + ext;
            }

            if (hasNoYPrefix) {
                outNoYFastqFile[imate] = noYFastqOutputPrefix + "mate" + to_string(imate + 1) + ext;
            } else if (!useMateFallback && readFilesNames.size() > imate && !readFilesNames[imate].empty()) {
                string base = pathBasename(readFilesNames[imate][0]);
                string tagged = insertTagBeforeReadToken(base, "_noY");
                if (!tagged.empty()) {
                    outNoYFastqFile[imate] = outputDir + adjustCompressionExt(tagged, emitYNoYFastqCompression);
                } else {
                    outNoYFastqFile[imate] = outFileNamePrefix + "noY_reads.mate" + to_string(imate + 1) + ext;
                }
            } else {
                outNoYFastqFile[imate] = outFileNamePrefix + "noY_reads.mate" + to_string(imate + 1) + ext;
            }
        }
    }

    // Open Y/noY FASTQ output files (uncompressed only)
    if (runMode == "alignReads" && emitYNoYFastqyes && emitYNoYFastqCompression == "none") {
        for (uint imate = 0; imate < readNmates; imate++) {
            inOut->outYFastqStream[imate].open(outYFastqFile[imate].c_str());
            inOut->outNoYFastqStream[imate].open(outNoYFastqFile[imate].c_str());
            if (!inOut->outYFastqStream[imate].is_open() || !inOut->outNoYFastqStream[imate].is_open()) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: could not create Y/noY FASTQ output files\n";
                errOut << "Solution: check that you have permission to write and disk space\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
            }
        }
    }
    
    // Open BAM files after Solo initialization (so pSolo.addTagsToUnsorted is available)
    if (outBAMunsorted) {
        // Decide the write mode once, immediately after Solo init
        bool useSoloTmp = pSolo.addTagsToUnsorted && pSolo.samAttrYes;
        if (pSolo.skipProcessing) {
            useSoloTmp = false;
        } else if (pSolo.writeKeysBin) {
            useSoloTmp = true;
        }
        outBAMunsortedUseSoloTmp = useSoloTmp;
        // Safety: when two-pass unsorted injection is active, skip Solo post-map to avoid
        // a known multithreading crash in the post-map path. This still writes CB/UB sidecar
        // and performs pass-2 tag injection.
        if (outBAMunsortedUseSoloTmp && !pSolo.skipProcessing) {
            bool parityDebug = false;
#ifdef DEBUG_CB_UB_PARITY
            parityDebug = (std::getenv("STAR_DEBUG_CB_UB_PARITY") != nullptr);
#endif
            if (parityDebug) {
                inOut->logMain << "INFO: CB/UB parity validation forcing Solo counting despite two-pass unsorted injection" << std::endl;
            } else {
                inOut->logMain << "INFO: Enabling --soloSkipProcessing automatically for two-pass unsorted CB/UB injection" << std::endl;
                pSolo.skipProcessing = true;
            }
        }
        
        if (useSoloTmp) {
            // Two-pass mode: create temp file, skip opening final BAM for now
            outBAMfileUnsortedSoloTmpName = outFileNamePrefix + "Aligned.out.unsorted.solo.tmp";
            inOut->logMain << "DEBUG: Using solo tmp mode for unsorted BAM (CB/UB tags requested)" << std::endl;
            // Open shared tmp stream once
            inOut->outBAMfileUnsortedSoloTmp.open(outBAMfileUnsortedSoloTmpName.c_str(), ios::binary | ios::out | ios::trunc);
            if (!inOut->outBAMfileUnsortedSoloTmp.is_open()) {
                ostringstream errOut;
                errOut << "EXITING because of fatal OUTPUT ERROR: could not create solo tmp file: " << outBAMfileUnsortedSoloTmpName << "\n";
                errOut << "SOLUTION: check the path and permissions\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            // Final BAM will be opened in pass 2
            inOut->outBAMfileUnsorted = NULL;
        } else {
            // Direct mode: open final BAM directly, clear tmp filename
            outBAMfileUnsortedSoloTmpName = "";
            inOut->logMain << "DEBUG: Using direct mode for unsorted BAM (no CB/UB tags or not requested)" << std::endl;
            inOut->outBAMfileUnsorted = bgzf_open(outBAMfileUnsortedName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
            if (inOut->outBAMfileUnsorted == NULL) {
                ostringstream errOut;
                errOut << "EXITING because of fatal OUTPUT ERROR: could not create unsorted BAM file: " << outBAMfileUnsortedName << "\n";
                errOut << "SOLUTION: check the path and permissions\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            }
            
            // Open Y/noY BAM handles if Y-split is enabled
            if (emitNoYBAMyes) {
                inOut->outBAMfileY = bgzf_open(outBAMfileYName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
                if (inOut->outBAMfileY == NULL) {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal OUTPUT ERROR: could not create Y BAM file: " << outBAMfileYName << "\n";
                    errOut << "SOLUTION: check the path and permissions\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
                inOut->outBAMfileNoY = bgzf_open(outBAMfileNoYName.c_str(),("w"+to_string((long long) outBAMcompression)).c_str());
                if (inOut->outBAMfileNoY == NULL) {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal OUTPUT ERROR: could not create noY BAM file: " << outBAMfileNoYName << "\n";
                    errOut << "SOLUTION: check the path and permissions\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                }
                // Write headers to Y/noY BAMs (will be written by BAMoutput constructor)
            } else {
                inOut->outBAMfileY = nullptr;
                inOut->outBAMfileNoY = nullptr;
            }
        }
    }
    
    //clipping
    pClip.initialize(this);

        // Warn if trimCutadapt is enabled (ClipMate will be bypassed)
        if (trimCutadapt == "Yes") {
            inOut->logMain << "WARNING: --trimCutadapt is enabled. Existing clip* parameters will be IGNORED.\n";
            inOut->logMain << "         ClipMate clipping is bypassed; trimming uses trimCutadapt* parameters instead.\n";
        }
        
        // Warn if compatibility mode is enabled (non-default trimming)
        if (trimCutadapt == "Yes" && trimCutadaptCompat == "Cutadapt3") {
            inOut->logMain << "NOTICE: --trimCutadaptCompat Cutadapt3 is enabled. Using cutadapt 3.x compatibility mode.\n";
            inOut->logMain << "         This mode reproduces Trim Galore/cutadapt 3.2 behavior for adapter matching.\n";
        }

        trimQcEnabled = (!trimQcReport.empty() && trimQcReport != "-") ||
                        (!trimQcJson.empty() && trimQcJson != "-") ||
                        (!trimQcHtml.empty() && trimQcHtml != "-");

    //alignEnds
    alignEndsType.ext[0][0]=false;
    alignEndsType.ext[0][1]=false;
    alignEndsType.ext[1][0]=false;
    alignEndsType.ext[1][1]=false;

    if (alignEndsType.in=="EndToEnd") {
        alignEndsType.ext[0][0]=true;
        alignEndsType.ext[0][1]=true;
        alignEndsType.ext[1][0]=true;
        alignEndsType.ext[1][1]=true;
    } else if (alignEndsType.in=="Extend5pOfRead1" ) {
        alignEndsType.ext[0][0]=true;
    } else if (alignEndsType.in=="Extend5pOfReads12" ) {
        alignEndsType.ext[0][0]=true;
        alignEndsType.ext[1][0]=true;
    } else if (alignEndsType.in=="Extend3pOfRead1" ) {
        alignEndsType.ext[0][1]=true;
    } else if (alignEndsType.in=="Local") {
        //nothing to do for now
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --alignEndsType: "<<alignEndsType.in <<"\n";
        errOut <<"SOLUTION: re-run STAR with --alignEndsType Local OR EndToEnd OR Extend5pOfRead1 OR Extend3pOfRead1\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //open compilation-dependent streams
    #ifdef OUTPUT_localChains
            inOut->outLocalChains.open((outFileNamePrefix + "LocalChains.out.tab").c_str());
    #endif

    strcpy(genomeNumToNT,"ACGTN");

   //sjdb insert on the fly
    sjdbInsert.pass1=false;
    sjdbInsert.pass2=false;
    sjdbInsert.yes=false;
    if (pGe.sjdbFileChrStartEnd.at(0)!="-" || pGe.sjdbGTFfile!="-") {//will insert annotated sjdb on the fly
       sjdbInsert.pass1=true;
       sjdbInsert.yes=true;
    };
    if (twoPass.yes) {
       sjdbInsert.pass2=true;
       sjdbInsert.yes=true;
    };

    if (pGe.gLoad!="NoSharedMemory" && sjdbInsert.yes ) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: on the fly junction insertion and 2-pass mappng cannot be used with shared memory genome \n" ;
        errOut << "SOLUTION: run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (runMode=="alignReads" && sjdbInsert.yes )
    {//run-time genome directory, this is needed for genome files generated on the fly
        if (pGe.sjdbOverhang<=0) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: pGe.sjdbOverhang <=0 while junctions are inserted on the fly with --sjdbFileChrStartEnd or/and --sjdbGTFfile\n";
            errOut << "SOLUTION: specify pGe.sjdbOverhang>0, ideally readmateLength-1";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        sjdbInsert.outDir=outFileNamePrefix+"_STARgenome/";
        sysRemoveDir (sjdbInsert.outDir);
        if (mkdir (sjdbInsert.outDir.c_str(),runDirPerm)!=0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not make run-time genome directory directory: "<< sjdbInsert.outDir<<"\n";
            errOut <<"SOLUTION: please check the path and writing permissions \n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (outBAMcoord && limitBAMsortRAM==0) {//check limitBAMsortRAM
        if (pGe.gLoad!="NoSharedMemory") {
            ostringstream errOut;
            errOut <<"EXITING because of fatal PARAMETERS error: limitBAMsortRAM=0 (default) cannot be used with --genomeLoad="<<pGe.gLoad <<", or any other shared memory options\n";
            errOut <<"SOLUTION: please use default --genomeLoad NoSharedMemory, \n        OR specify --limitBAMsortRAM the amount of RAM (bytes) that can be allocated for BAM sorting in addition to shared memory allocated for the genome.\n        --limitBAMsortRAM typically has to be > 10000000000 (i.e 10GB).\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        inOut->logMain<<"WARNING: --limitBAMsortRAM=0, will use genome size as RAM limit for BAM sorting\n";
    };

    for (uint ii=0; ii<readNameSeparator.size(); ii++) {
        if (readNameSeparator.at(ii)=="space") {
            readNameSeparatorChar.push_back(' ');
        } else if (readNameSeparator.at(ii)=="none") {
            //nothing to do
        } else if (readNameSeparator.at(ii).size()==1) {
            readNameSeparatorChar.push_back(readNameSeparator.at(ii).at(0));
        } else{
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized value of --readNameSeparator="<<readNameSeparator.at(ii)<<"\n";
            errOut << "SOLUTION: use allowed values: space OR single characters";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    //outSAMunmapped
    outSAMunmapped.yes=false;
    outSAMunmapped.within=false;
    outSAMunmapped.keepPairs=false;
    if (outSAMunmapped.mode.at(0)=="None" && outSAMunmapped.mode.size()==1) {
        //nothing to do, all false
    } else if (outSAMunmapped.mode.at(0)=="Within" && outSAMunmapped.mode.size()==1) {
        outSAMunmapped.yes=true;
        outSAMunmapped.within=true;
    } else if (outSAMunmapped.mode.at(0)=="Within" && outSAMunmapped.mode.at(1)=="KeepPairs") {
        outSAMunmapped.yes=true;
        outSAMunmapped.within=true;
        if (readNmates==2) //not readNends, since this control output of alignments
            outSAMunmapped.keepPairs=true;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option for --outSAMunmapped=";
        for (uint ii=0; ii<outSAMunmapped.mode.size(); ii++) errOut <<" "<< outSAMunmapped.mode.at(ii);
        errOut << "\nSOLUTION: use allowed options: None OR Within OR Within KeepPairs";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    alignEndsProtrude.nBasesMax=stoi(alignEndsProtrude.in.at(0),nullptr);
    alignEndsProtrude.concordantPair=false;
    if (alignEndsProtrude.nBasesMax>0) {//allow ends protrusion
        if (alignEndsProtrude.in.at(1)=="ConcordantPair") {
            alignEndsProtrude.concordantPair=true;
        } else if (alignEndsProtrude.in.at(1)=="DiscordantPair") {
            alignEndsProtrude.concordantPair=false;
        } else  {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignEndsProtrude="<<alignEndsProtrude.in.at(1)<<"\n";
            errOut << "SOLUTION: use allowed option: ConcordantPair or DiscordantPair";
            exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if (alignInsertionFlush.in=="None") {
        alignInsertionFlush.flushRight=false;
    } else if (alignInsertionFlush.in=="Right") {
        alignInsertionFlush.flushRight=true;
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --alignInsertionFlush="<<alignInsertionFlush.in<<"\n";
        errOut << "SOLUTION: use allowed option: None or Right";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    //peOverlap
    if (peOverlap.NbasesMin>0) {
        peOverlap.yes=true;
    } else {
        peOverlap.yes=false;
    };

    //alignSoftClipAtReferenceEnds.in
    if (alignSoftClipAtReferenceEnds.in=="Yes") {
        alignSoftClipAtReferenceEnds.yes=true;
    } else if (alignSoftClipAtReferenceEnds.in=="No") {
        alignSoftClipAtReferenceEnds.yes=false;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in --alignSoftClipAtReferenceEnds   "<<alignSoftClipAtReferenceEnds.in<<"\n";
        errOut << "SOLUTION: use allowed option: Yes or No";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    outSAMreadIDnumber=false;
    if (outSAMreadID=="Number") {
        outSAMreadIDnumber=true;
    };


    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////// these parameters do not depend on other parameters
    /////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////// limitIObufferSize
    /* old before 2.7.9
    // in/out buffers
    #define BUFFER_InSizeFraction 0.5
    if (limitIObufferSize<limitOutSJcollapsed*Junction::dataSize+1000000) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: --limitIObufferSize="<<limitIObufferSize <<" is too small for ";
        errOut << "--limitOutSJcollapsed*"<<Junction::dataSize<<"="<< limitOutSJcollapsed<<"*"<<Junction::dataSize<<"="<<limitOutSJcollapsed*Junction::dataSize<<"\n";
        errOut <<"SOLUTION: re-run STAR with larger --limitIObufferSize or smaller --limitOutSJcollapsed\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    chunkInSizeBytesArray=(uint) int((limitIObufferSize-limitOutSJcollapsed*Junction::dataSize)*BUFFER_InSizeFraction)/2;
    chunkOutBAMsizeBytes= (uint) int((1.0/BUFFER_InSizeFraction-1.0)*chunkInSizeBytesArray*2.0);
    chunkInSizeBytes=chunkInSizeBytesArray-2*(DEF_readSeqLengthMax+1)-2*DEF_readNameLengthMax;//to prevent overflow
    */
    
    if (limitIObufferSize.size() != 2) 
        exitWithError("EXITING because of FATAL input ERROR: --limitIObufferSize requires 2 numbers since 2.7.9a.\n"
                      "SOLUTION: specify 2 numbers in --limitIObufferSize : size of input and output buffers in bytes.\n"
                        , std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    
    chunkInSizeBytesArray = limitIObufferSize[0]/readNends; //array size
    chunkInSizeBytes = chunkInSizeBytesArray-2*(DEF_readSeqLengthMax+1)-2*DEF_readNameLengthMax; //to prevent overflow - array is bigger to allow loading one read
    chunkOutBAMsizeBytes = limitIObufferSize[1];
    
    
    ///////////////////////////////////////////////////////// outSJ
    if (outSJ.type[0] == "None") {
        outSJ.yes = false;
    } else if (outSJ.type[0] == "Standard") {
        outSJ.yes = true;
    } else {
        exitWithError("EXITING because of FATAL input ERROR: unrecognized option in --outSJtype   " + outSJ.type[0] + '\n' +
                      "SOLUTION: use one of the allowed options: --outSJtype   Standard    OR    None\n"
                        , std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outFilterType=="Normal") {
        outFilterBySJoutStage=0;
    } else if (outFilterType=="BySJout") {
        if (!outSJ.yes)
            exitWithError("EXITING because of FATAL input ERROR: --outFilterType BySJout requires --outSJtype Standard\n"
                      "SOLUTION: --outFilterType Normal    OR   --outFilterType BySJout --outSJtype Standard\n"
                        , std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            
        outFilterBySJoutStage=1;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL input ERROR: unknown value of parameter outFilterType: " << outFilterType <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Normal | BySJout\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };    
    
    ////////////////////////////////////////////////
    inOut->logMain << "Finished loading and checking parameters\n" <<flush;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameters::scanAllLines (istream &streamIn, int inputLevel,  int inputLevelRequested) {//scan
//     istringstream stringInStream (stringIn);
    string lineIn;
    while (getline(streamIn,lineIn)) {
        scanOneLine(lineIn, inputLevel, inputLevelRequested);
    };
};

int Parameters::scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested) {//scan one line and load the parameters,
                                                             //0 if comment, 1 if OK
    if (lineIn=="") return 0; //empty line

    istringstream lineInStream (lineIn);

    if (inputLevel==0 && ( lineIn.substr(0,1)==" " || lineIn.substr(0,1)=="\t" ) ) return 0;//for Default input spaces also mark comments, for nice formatting

    string parIn("");
    lineInStream >> parIn;
    if (parIn=="" || parIn.substr(0,2)=="//" || parIn.substr(0,1)=="#") return 0; //this is a comment

    uint iPar;
    for (iPar=0; iPar<parArray.size(); iPar++) {
        if (parIn==parArray[iPar]->nameString) {//
            if (inputLevelRequested < 0 || inputLevelRequested == parArray[iPar]->inputLevelAllowed) {
                break;//will read this parameter values
            } else {
                return 1; //do not read inputs not requested at this level
            };
        };
    };

    string parV("");
    lineInStream >> parV;
    if (parV=="") {//parameter value cannot be empty
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: empty value for parameter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use non-empty value for this parameter\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    lineInStream.str(lineIn); lineInStream.clear(); lineInStream >> parIn; //get the correct state of stream, past reading parIn

    if (iPar==parArray.size()) {//string is not identified
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: unrecognized parameter name \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use correct parameter name (check the manual)\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    } else {//found the corresponding parameter
        if (inputLevel==0 && parArray[iPar]->inputLevel>0) {//this is one of the initial parameters, it was read from Command Line and should not be re-defined
            getline(lineInStream,parV);
            inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString <<parV<<" ... is RE-DEFINED on Command Line as: " << *(parArray[iPar]) <<"\n";
        } else if (parArray[iPar]->inputLevelAllowed>0 && parArray[iPar]->inputLevelAllowed < inputLevel) {//this is initial parameter and cannot be redefined
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: parameter \""<< parIn << "\" cannot be defined at the input level \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: define parameter \""<< parIn << "\" in \"" << parameterInputName.at(parArray[iPar]->inputLevelAllowed) <<"\"\n" <<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (parArray[iPar]->inputLevel==inputLevel) {//this parameter was already defined at this input level
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: duplicate parameter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: keep only one definition of input parameters in each input source\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else {//read values
            parArray[iPar]->inputValues(lineInStream);
            parArray[iPar]->inputLevel=inputLevel;
            if ( inOut->logMain.good() ) {
                inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString << *(parArray[iPar]);
                if ( parArray[iPar]->inputLevel > 0 ) inOut->logMain <<"     ~RE-DEFINED";
                inOut->logMain << endl;
            };
        };
    };
    return 0;
};
