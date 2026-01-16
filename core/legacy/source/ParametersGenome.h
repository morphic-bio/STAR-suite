#ifndef CODE_ParametersGenome
#define CODE_ParametersGenome

#include "IncludeDefine.h"
#include <unordered_set>

class Parameters;

class ParametersGenome {//"constant" genome parameters - user input
public:
    string gDir;
    string gLoad;
    
    uint32 gType;//type code
    string gTypeString;
    
    vector <string> gFastaFiles;
    vector <string> gChainFiles;
    //string gConsensusFile; DEPRECATED

    struct {
        int32 type;
        string typeString;
        string vcfFile;
        vector<string> output; //which output to transform
        bool outYes, outSAM, outSJ, outQuant;
    } transform;
    
    uint gSAindexNbases;//length of the SA pre-index strings
    uint gChrBinNbits;
    uint gSAsparseD;//SA sparsity
    uint gSuffixLengthMax;//maximum length of the suffixes, has to be longer than read length
    vector <uint> gFileSizes;//size of the genome files

    vector <string> sjdbFileChrStartEnd;
    string sjdbGTFfile;
    string sjdbGTFchrPrefix;
    string transcriptomeFasta;
    
    string sjdbGTFfeatureExon;
    string sjdbGTFtagExonParentTranscript;
    string sjdbGTFtagExonParentGene;
    vector<string> sjdbGTFtagExonParentGeneName;
    vector<string> sjdbGTFtagExonParentGeneType;
    
    string sjdbInsertSave;
    uint sjdbOverhang;
    int sjdbOverhang_par;
    int sjdbScore;

    struct {
        vector<string> mitoStrings;
        unordered_set<uint64> mito;
    } chrSet;

    // Flex gene probe parameters (50bp gene probes, distinct from 8bp sample tags)
    struct {
        string csvFile;           // Path to 50bp gene probe CSV
        uint32 enforceLength;     // Expected probe length (default: 50)
        string removeDeprecated;  // --removeDeprecated Yes/No (remove deprecated gene IDs from probe lists)
        bool enabled;             // Whether flex gene probe processing is enabled
        bool removeDeprecatedBool;// Converted from removeDeprecated string
    } flexGeneProbe;
    
    // CellRanger-style reference formatting parameters
    struct {
        string indexEnabled;   // --cellrangerStyleIndex Yes/No (gates GTF/FASTA formatting)
        string downloadOnly;    // --cellrangerStyleDownloadOnly Yes/No
        string cacheDir;       // --cellrangerStyleCacheDir
        string allUntrustedUrl;// --allUntrustedUrl Yes/No (required for URLs not in trusted table)
        string faUrl;          // --faUrl <url> (FASTA download URL)
        string gtfUrl;         // --gtfUrl <url> (GTF download URL)
        string autoCksumUpdate;// --autoCksumUpdate Yes/No (auto-fill missing cksum from CHECKSUMS files)
        string replaceUnverifiableFiles;// --replaceUnverifiableFiles Yes/No (replace unverifiable existing files)
        string refRelease;     // --cellrangerRefRelease (default: 2024-A)
        string refReleaseCanonical; // Normalized canonical release string (2024-A or 2020-A)
        bool indexEnabledBool; // Converted from indexEnabled string
        bool downloadOnlyBool; // Converted from downloadOnly string
        bool allUntrustedUrlBool;// Converted from allUntrustedUrl string
        bool autoCksumUpdateBool;// Converted from autoCksumUpdate string
        bool replaceUnverifiableFilesBool;// Converted from replaceUnverifiableFiles string
    } cellrangerStyle;
    
    // Auto-index workflow parameters
    struct {
        string autoIndex;      // --autoIndex Yes/No (automatically build index unless it exists)
        string forceIndex;     // --forceIndex Yes/No (rebuild index even if exists)
        string forceAllIndex;  // --forceAllIndex Yes/No (re-download and rebuild everything)
        bool autoIndexBool;    // Converted from autoIndex string
        bool forceIndexBool;   // Converted from forceIndex string
        bool forceAllIndexBool;// Converted from forceAllIndex string
    } autoIndexWorkflow;
    
    // Transcriptome FASTA generation parameters (index-time)
    struct {
        string generateTranscriptome;      // --genomeGenerateTranscriptome Yes/No (generate transcriptome.fa at index time)
        string transcriptomeFastaPath;     // --genomeGenerateTranscriptomeFasta <path> (output path, default: ${genomeDir}/transcriptome.fa)
        string overwrite;                  // --genomeGenerateTranscriptomeOverwrite Yes/No (overwrite existing file)
        bool generateTranscriptomeBool;   // Converted from generateTranscriptome string
        bool overwriteBool;                // Converted from overwrite string
    } transcriptomeGen;
    
    void initialize(Parameters *Pin);

private:
    Parameters *pP;
    
};

#endif
