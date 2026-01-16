#include <cmath>
#include <map>

#include "Genome.h"

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "GTF.h"
#include "SjdbClass.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "genomeParametersWrite.h"
#include "sjdbInsertJunctions.h"
#include "genomeScanFastaFiles.h"
#include "genomeSAindex.h"

#include "serviceFuns.cpp"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "FlexProbeIndex.h"
#include "CellRangerFormatter.h"
#include "TranscriptomeFasta.h"
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include <dirent.h>


char* globalG;
uint globalL;


inline int funCompareSuffixes ( const void *a, const void *b){

    uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
    uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));

    uint jj=0;
    int  ii=0;
    uint va=0,vb=0;
    uint8 *va1, *vb1;

    while (jj < globalL) {
        va=*(ga-jj);
        vb=*(gb-jj);

        #define has5(v) ((((v)^0x0505050505050505) - 0x0101010101010101) & ~((v)^0x0505050505050505) & 0x8080808080808080)

        if (has5(va) && has5(vb))
        {//there is 5 in the sequence - only compare bytes before 5
            va1=(uint8*) &va;
            vb1=(uint8*) &vb;
            for (ii=7;ii>=0;ii--)
            {
                if (va1[ii]>vb1[ii])
                {
                    return 1;
                } else if (va1[ii]<vb1[ii])
                {
                    return -1;
                } else if (va1[ii]==5)
                {//va=vb at the end of chr
                    if ( *((uint*)a) > *((uint*)b) )
                    {//anti-stable order,since indexes are sorted in the reverse order
                        return  -1;
                    } else
                    {//a cannot be equal to b
                        return 1;
                    };
                };
            };
        } else
        {//no 5, simple comparison
            if (va>vb)
            {
                return 1;
            } else if (va<vb)
            {
                return -1;
            };
        };
        jj++;
    };

    //suffixes are equal up to globalL, no simply compare the indexes
    if ( *((uint*)a) > *((uint*)b) )
    {//anti-stable order,since indexes are sorted in the reverse order
        return  -1;
    } else
    {//a cannot be equal to b
        return 1;
    };
};

inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
};

void Genome::genomeGenerate() {

    //check parameters
	createDirectory(pGe.gDir, P.runDirPerm, "--genomeDir", P);

	{//move Log.out file into genome directory
		string logfn=pGe.gDir+"Log.out";
		if ( rename( P.outLogFileName.c_str(), logfn.c_str() ) ) {
			warningMessage("Could not move Log.out file from " + P.outLogFileName + " into " + logfn + ". Will keep " + P.outLogFileName +"\n", \
						   std::cerr, P.inOut->logMain, P);
		} else {
			P.outLogFileName=logfn;
		};
	};
	
	// Auto-index workflow: check if index exists and handle force flags
	if (pGe.autoIndexWorkflow.autoIndexBool) {
	    // Check if index already exists (look for key sentinel files)
	    bool indexExists = false;
	    struct stat fileStat;
	    string genomeFile = pGe.gDir + "/Genome";
	    string saFile = pGe.gDir + "/SA";
	    string saIndexFile = pGe.gDir + "/SAindex";
	    string genomeParamsFile = pGe.gDir + "/genomeParameters.txt";
	    
	    // Index exists if all key files are present
	    if (stat(genomeFile.c_str(), &fileStat) == 0 &&
	        stat(saFile.c_str(), &fileStat) == 0 &&
	        stat(saIndexFile.c_str(), &fileStat) == 0 &&
	        stat(genomeParamsFile.c_str(), &fileStat) == 0) {
	        indexExists = true;
	    }
	    
	    // If index exists and no force flags, exit successfully
	    if (indexExists && !pGe.autoIndexWorkflow.forceIndexBool && !pGe.autoIndexWorkflow.forceAllIndexBool) {
	        P.inOut->logMain << "Index already exists in " << pGe.gDir << ". Skipping index generation.\n";
	        P.inOut->logMain << "Use --forceIndex Yes to rebuild the index, or --forceAllIndex Yes to re-download and rebuild everything.\n";
	        *(P.inOut->logStdOut) << "Index already exists in " << pGe.gDir << ". Skipping index generation.\n";
	        exit(0);
	    }
	    
	    // If forceAllIndex, delete cached downloads and formatted files
	    if (pGe.autoIndexWorkflow.forceAllIndexBool) {
	        P.inOut->logMain << "forceAllIndex enabled: cleaning up cached downloads and formatted files...\n";
	        
	        string cacheDir = pGe.cellrangerStyle.cacheDir.empty() 
	                          ? pGe.gDir + "cellranger_ref_cache" 
	                          : pGe.cellrangerStyle.cacheDir;
	        
	        // Delete cache directory contents (but keep the directory)
	        string rmCmd = "rm -rf \"" + cacheDir + "\"/* \"" + cacheDir + "\"/.* 2>/dev/null || true";
	        if (system(rmCmd.c_str()) != 0) {
	            // Ignore errors - cache dir might not exist
	        }
	        
	        // Delete formatted files in genome directory
	        string formattedFasta = pGe.gDir + "cellranger_ref/genome.fa";
	        string formattedGtf = pGe.gDir + "cellranger_ref/genes.gtf";
	        remove(formattedFasta.c_str());
	        remove(formattedGtf.c_str());
	        
	        P.inOut->logMain << "  Cache and formatted files cleaned up.\n";
	    }
	    
	    // Input resolution priority: files -> URLs -> default URLs
	    bool hasFiles = (pGe.gFastaFiles.size() > 0 && pGe.gFastaFiles.at(0) != "-") && 
	                    (pGe.sjdbGTFfile != "-");
	    bool hasUrls = !pGe.cellrangerStyle.faUrl.empty() && pGe.cellrangerStyle.faUrl != "-" &&
	                   !pGe.cellrangerStyle.gtfUrl.empty() && pGe.cellrangerStyle.gtfUrl != "-";
	    
	    // If no files and no URLs, use default URLs based on selected release
	    if (!hasFiles && !hasUrls) {
	        // Release mapping table (extensible design)
	        struct ReleaseUrls {
	            string release;
	            string ensemblRelease;
	            string gencodeRelease;
	            string faUrl;
	            string gtfUrl;
	        };
	        
	        static const map<string, ReleaseUrls> releaseMap = {
	            {"2024-A", {
	                "2024-A",
	                "110",
	                "44",
	                "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
	                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
	            }},
	            {"2020-A", {
	                "2020-A",
	                "98",
	                "32",
	                "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
	                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
	            }}
	        };
	        
	        // Get selected release (default to 2024-A if not set)
	        string selectedRelease = pGe.cellrangerStyle.refReleaseCanonical.empty() 
	                                ? "2024-A" 
	                                : pGe.cellrangerStyle.refReleaseCanonical;
	        
	        auto it = releaseMap.find(selectedRelease);
	        if (it == releaseMap.end()) {
	            ostringstream errOut;
	            errOut << "EXITING because of FATAL ERROR: Unknown CellRanger release: " << selectedRelease << "\n";
	            errOut << "SOLUTION: Use --cellrangerRefRelease with a supported release (2024-A or 2020-A)\n";
	            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
	        }
	        
	        const ReleaseUrls& urls = it->second;
	        
	        P.inOut->logMain << "No input files or URLs provided. Using default CellRanger reference release: " << urls.release << "\n";
	        P.inOut->logMain << "  Release mapping: Ensembl release-" << urls.ensemblRelease 
	                         << " / GENCODE release-" << urls.gencodeRelease << "\n";
	        
	        pGe.cellrangerStyle.faUrl = urls.faUrl;
	        pGe.cellrangerStyle.gtfUrl = urls.gtfUrl;
	        
	        P.inOut->logMain << "  Default FASTA URL: " << pGe.cellrangerStyle.faUrl << "\n";
	        P.inOut->logMain << "  Default GTF URL: " << pGe.cellrangerStyle.gtfUrl << "\n";
	        
	        // Default URLs are trusted, so no need for --allUntrustedUrl
	        hasUrls = true;
	    }
	    
	    if (hasUrls) {
	        P.inOut->logMain << "Input resolution: Using URLs for download.\n";
	    } else if (hasFiles) {
	        P.inOut->logMain << "Input resolution: Using provided files.\n";
	    }
	}
	
    // AutoIndex may set URLs without populating sjdbGTFfile yet; treat URL as annotation input.
    bool hasGtfInput = (pGe.sjdbGTFfile != "-");
    if (!hasGtfInput && pGe.autoIndexWorkflow.autoIndexBool) {
        if (!pGe.cellrangerStyle.gtfUrl.empty() && pGe.cellrangerStyle.gtfUrl != "-") {
            hasGtfInput = true;
        }
    }
    bool hasAnnoInput = (pGe.sjdbFileChrStartEnd.at(0)!="-" || hasGtfInput);

    if (sjdbOverhang<=0 && hasAnnoInput) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: for generating genome with annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
        errOut << "you need to specify >0 --sjdbOverhang\n";
        errOut << "SOLUTION: re-run genome generation specifying non-zero --sjdbOverhang, which ideally should be equal to OneMateLength-1, or could be chosen generically as ~100\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (pGe.sjdbFileChrStartEnd.at(0)=="-" && !hasGtfInput) {
        if (P.parArray.at(P.pGe.sjdbOverhang_par)->inputLevel>0 && sjdbOverhang>0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: when generating genome without annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
            errOut << "do not specify >0 --sjdbOverhang\n";
            errOut << "SOLUTION: re-run genome generation without --sjdbOverhang option\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        sjdbOverhang=0;
    };
    
    // Download if URLs are provided
    bool hasUrls = !pGe.cellrangerStyle.faUrl.empty() && pGe.cellrangerStyle.faUrl != "-" &&
                    !pGe.cellrangerStyle.gtfUrl.empty() && pGe.cellrangerStyle.gtfUrl != "-";
    
    if (hasUrls) {
        P.inOut->logMain << "Downloading reference files...\n";
        
        string cacheDir = pGe.cellrangerStyle.cacheDir.empty() 
                          ? pGe.gDir + "cellranger_ref_cache" 
                          : pGe.cellrangerStyle.cacheDir;
        
        // Create cache directory
        string mkdirCmd = "mkdir -p \"" + cacheDir + "\"";
        if (system(mkdirCmd.c_str()) != 0) {
            exitWithError("EXITING because of FATAL ERROR: Could not create cache directory: " + cacheDir + "\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        string fastaUrl = pGe.cellrangerStyle.faUrl;
        string gtfUrl = pGe.cellrangerStyle.gtfUrl;
        
        // Extract filenames for output paths
        string fastaFilename = fastaUrl.substr(fastaUrl.find_last_of('/') + 1);
        string gtfFilename = gtfUrl.substr(gtfUrl.find_last_of('/') + 1);
        // Remove .gz extension if present
        if (fastaFilename.length() > 3 && fastaFilename.substr(fastaFilename.length() - 3) == ".gz") {
            fastaFilename = fastaFilename.substr(0, fastaFilename.length() - 3);
        }
        if (gtfFilename.length() > 3 && gtfFilename.substr(gtfFilename.length() - 3) == ".gz") {
            gtfFilename = gtfFilename.substr(0, gtfFilename.length() - 3);
        }
        
        string downloadedFasta = cacheDir + "/" + fastaFilename;
        string downloadedGtf = cacheDir + "/" + gtfFilename;
        
        // Get expected cksums from EXPECTED_CKSUM map (if available)
        uint32_t expectedFastaCksum = 0;
        uint64_t expectedFastaSize = 0;
        uint32_t expectedGtfCksum = 0;
        uint64_t expectedGtfSize = 0;
        
        // Load cksum cache before downloads
        CellRangerFormatter::loadCksumCache(cacheDir);
        
        string errorMsg;
        P.inOut->logMain << "    Downloading FASTA from: " << fastaUrl << "\n";
        if (!CellRangerFormatter::downloadReference(fastaUrl, downloadedFasta, expectedFastaCksum, expectedFastaSize,
                                                    pGe.cellrangerStyle.allUntrustedUrlBool, cacheDir,
                                                    pGe.cellrangerStyle.autoCksumUpdateBool,
                                                    pGe.cellrangerStyle.replaceUnverifiableFilesBool, errorMsg)) {
            exitWithError("EXITING because of FATAL ERROR downloading FASTA: " + errorMsg + "\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        P.inOut->logMain << "    Downloading GTF from: " << gtfUrl << "\n";
        if (!CellRangerFormatter::downloadReference(gtfUrl, downloadedGtf, expectedGtfCksum, expectedGtfSize,
                                                    pGe.cellrangerStyle.allUntrustedUrlBool, cacheDir,
                                                    pGe.cellrangerStyle.autoCksumUpdateBool,
                                                    pGe.cellrangerStyle.replaceUnverifiableFilesBool, errorMsg)) {
            exitWithError("EXITING because of FATAL ERROR downloading GTF: " + errorMsg + "\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        // Update genome files to use downloaded versions
        pGe.gFastaFiles[0] = downloadedFasta;
        pGe.sjdbGTFfile = downloadedGtf;
        P.inOut->logMain << "    Download completed successfully\n";
        
        // Exit early if download-only mode
        if (pGe.cellrangerStyle.downloadOnlyBool) {
            P.inOut->logMain << "Download-only mode: skipping index generation\n";
            exit(0);
        }
    }
    
    // CellRanger-style reference formatting - format FASTA/GTF before indexing
    if (pGe.cellrangerStyle.indexEnabledBool) {
        P.inOut->logMain << "CellRanger-style reference formatting enabled...\n";
        
        string inputFasta = pGe.gFastaFiles.size() > 0 && pGe.gFastaFiles.at(0) != "-" 
                            ? pGe.gFastaFiles.at(0) : "";
        string inputGtf = pGe.sjdbGTFfile != "-" ? pGe.sjdbGTFfile : "";
        
        // Validate inputs
        if (inputFasta.empty()) {
            exitWithError("EXITING because of FATAL INPUT ERROR: --cellrangerStyleIndex requires either --genomeFastaFiles or --faUrl\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        if (inputGtf.empty()) {
            exitWithError("EXITING because of FATAL INPUT ERROR: --cellrangerStyleIndex requires either --sjdbGTFfile or --gtfUrl\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        // Format files (use cellranger_ref subdirectory to match expected output paths)
        string formattedFasta = pGe.gDir + "cellranger_ref/genome.fa";
        string formattedGtf = pGe.gDir + "cellranger_ref/genes.gtf";
        
        // Create cellranger_ref directory if it doesn't exist
        string mkdirCmd = "mkdir -p \"" + pGe.gDir + "cellranger_ref\"";
        if (system(mkdirCmd.c_str()) != 0) {
            exitWithError("EXITING because of FATAL ERROR: Could not create cellranger_ref directory: " + pGe.gDir + "cellranger_ref\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        CellRangerFormatter::Config formatConfig;
        formatConfig.inputFastaPath = inputFasta;
        formatConfig.outputFastaPath = formattedFasta;
        formatConfig.inputGtfPath = inputGtf;
        formatConfig.outputGtfPath = formattedGtf;
        
        P.inOut->logMain << "  Formatting FASTA and GTF files...\n";
        CellRangerFormatter::Result formatResult = CellRangerFormatter::format(formatConfig);
        if (!formatResult.success) {
            exitWithError("EXITING because of FATAL ERROR formatting reference: " + formatResult.errorMessage + "\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        // Update genome files to use formatted versions
        pGe.gFastaFiles[0] = formattedFasta;
        pGe.sjdbGTFfile = formattedGtf;
        P.inOut->logMain << "  Formatting completed successfully\n";
    }
    
    // Flex gene probe processing - generate hybrid reference if enabled
    if (pGe.flexGeneProbe.enabled) {
        P.inOut->logMain << "Flex gene probe processing enabled, generating hybrid reference...\n";
        
        // Validate inputs
        if (pGe.gFastaFiles.size() == 0 || pGe.gFastaFiles.at(0) == "-") {
            exitWithError("EXITING because of FATAL INPUT ERROR: --flexGeneProbeSet requires either --genomeFastaFiles or --faUrl\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        if (pGe.sjdbGTFfile == "-") {
            exitWithError("EXITING because of FATAL INPUT ERROR: --flexGeneProbeSet requires --sjdbGTFfile\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        // Configure FlexProbeIndex
        FlexProbeIndex::Config fpConfig;
        fpConfig.probeCSVPath = pGe.flexGeneProbe.csvFile;
        fpConfig.gtfPath = pGe.sjdbGTFfile;
        fpConfig.baseFastaPath = pGe.gFastaFiles.at(0);  // Use first FASTA file
        fpConfig.outputDir = pGe.gDir + "flex_probe_artifacts";
        fpConfig.enforceProbeLength = pGe.flexGeneProbe.enforceLength;
        fpConfig.removeDeprecated = pGe.flexGeneProbe.removeDeprecatedBool;
        
        // Run FlexProbeIndex
        P.inOut->logMain << "  Probe CSV: " << fpConfig.probeCSVPath << "\n";
        P.inOut->logMain << "  Target GTF: " << fpConfig.gtfPath << "\n";
        P.inOut->logMain << "  Base FASTA: " << fpConfig.baseFastaPath << "\n";
        P.inOut->logMain << "  Output dir: " << fpConfig.outputDir << "\n";
        
        FlexProbeIndex::Result fpResult = FlexProbeIndex::run(fpConfig);
        
        if (!fpResult.success) {
            exitWithError("EXITING because of FATAL ERROR in FlexProbeIndex: " + fpResult.errorMessage + "\n",
                         std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        
        P.inOut->logMain << "  FlexProbeIndex completed successfully:\n";
        P.inOut->logMain << "    Input probes: " << fpResult.stats.totalInput << "\n";
        P.inOut->logMain << "    Output probes: " << fpResult.stats.totalOutput << "\n";
        P.inOut->logMain << "    Unique genes: " << fpResult.stats.uniqueGenes << "\n";
        P.inOut->logMain << "    Dropped (DEPRECATED): " << fpResult.stats.droppedDeprecated << "\n";
        P.inOut->logMain << "    Dropped (included=FALSE): " << fpResult.stats.droppedIncludedFalse << "\n";
        P.inOut->logMain << "    Dropped (no GTF match): " << fpResult.stats.droppedNoMatch << "\n";
        if (pGe.flexGeneProbe.removeDeprecatedBool) {
            P.inOut->logMain << "    Note: Deprecated gene IDs were removed from probe_list.txt (--removeDeprecated enabled)\n";
        }
        
        // Update genome parameters to use hybrid files
        string hybridFasta = fpConfig.outputDir + "/genome.filtered.fa";
        string hybridGTF = fpConfig.outputDir + "/genes.filtered.gtf";
        
        P.inOut->logMain << "  Using hybrid FASTA: " << hybridFasta << "\n";
        P.inOut->logMain << "  Using hybrid GTF: " << hybridGTF << "\n";
        
        // Replace first FASTA file with hybrid
        pGe.gFastaFiles[0] = hybridFasta;
        // Replace GTF file with hybrid
        pGe.sjdbGTFfile = hybridGTF;
        
        // Copy probe_list.txt (generated by FlexProbeIndex) to genome index directory
        string srcProbeList = fpConfig.outputDir + "/probe_list.txt";
        string dstProbeList = pGe.gDir + "probe_gene_list.txt";
        
        std::ifstream srcFile(srcProbeList);
        std::ofstream dstFile(dstProbeList);
        if (srcFile.is_open() && dstFile.is_open()) {
            dstFile << srcFile.rdbuf();
            srcFile.close();
            dstFile.close();
            P.inOut->logMain << "  Copied probe gene list to: " << dstProbeList << "\n";
        } else {
            P.inOut->logMain << "  WARNING: Could not copy probe gene list to " << dstProbeList << "\n";
        }
    }
    
    //time
    time_t rawTime;
    string timeString;

    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to generate Genome files\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to generate Genome files\n" <<flush;

    //define some parameters from input parameters
    genomeChrBinNbases=1LLU << pGe.gChrBinNbits;

    nGenome = genomeScanFastaFiles(P,NULL,false,*this);//first scan the fasta file to find all the sizes
    genomeSequenceAllocate(nGenome, nG1alloc, G, G1);
    genomeScanFastaFiles(P,G,true,*this);    //load the genome sequence

    uint64 nGenomeTrue=0;
    for (auto &cl : chrLength)
    	nGenomeTrue += cl; //nGenomeTrue = sum of chr lengths

    P.inOut->logMain <<"Genome sequence total length = " << nGenomeTrue << "\n";
    P.inOut->logMain <<"Genome size with padding = "<< nGenome <<"\n";

    //consensusSequence(); //replace with consensus allele DEPRECATED
        
    SjdbClass sjdbLoci; //will be filled in transcriptGeneSJ below
    GTF mainGTF(*this, P, pGe.gDir, sjdbLoci); //this loads exonLoci and gene/transcript metadata only, sjdbLoci is not filled
    
    Genome::transformGenome(&mainGTF);
    
    mainGTF.superTranscript(); //this may change the genome into (Super)Transcriptome

    chrBinFill();//chrBin is first used in the transcriptGeneSJ below
    mainGTF.transcriptGeneSJ(pGe.gDir);
    
    // Generate transcriptome FASTA if enabled
    if (pGe.transcriptomeGen.generateTranscriptomeBool) {
        P.inOut->logMain << "Generating transcriptome FASTA...\n";
        
        // Determine output path
        string transcriptomePath = pGe.transcriptomeGen.transcriptomeFastaPath;
        
        // Check overwrite flag (also allow forceAllIndex to overwrite)
        bool overwrite = pGe.transcriptomeGen.overwriteBool || pGe.autoIndexWorkflow.forceAllIndexBool;
        
        // Generate standard transcriptome
        TranscriptomeFasta::writeTranscriptomeFasta(*this, mainGTF, transcriptomePath, P, 70, overwrite);
        
        // For CellRanger-style index, also write to cellranger_ref subdirectory
        // This is independent of the standard output - each has its own overwrite check
        if (pGe.cellrangerStyle.indexEnabledBool) {
            string cellrangerPath = pGe.gDir + "cellranger_ref/transcriptome.fa";
            P.inOut->logMain << "Also generating transcriptome FASTA for CellRanger-style reference...\n";
            
            // Create cellranger_ref directory if it doesn't exist
            string mkdirCmd = "mkdir -p \"" + pGe.gDir + "cellranger_ref\"";
            if (system(mkdirCmd.c_str()) != 0) {
                P.inOut->logMain << "WARNING: Could not create cellranger_ref directory: " << pGe.gDir << "cellranger_ref\n";
            } else {
                TranscriptomeFasta::writeTranscriptomeFasta(*this, mainGTF, cellrangerPath, P, 70, overwrite);
            }
        }
    }
    
    sjdbLoadFromFiles(P,sjdbLoci);//this will not be transformed. TODO prevent this parameter combination

    if (pGe.gSAindexNbases > log2(nGenomeTrue)/2-1) {
        ostringstream warnOut; 
        warnOut << "--genomeSAindexNbases " << pGe.gSAindexNbases << " is too large for the genome size=" << nGenomeTrue;
        warnOut << ", which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases " << int(log2(nGenomeTrue)/2-1);
        warningMessage(warnOut.str(),P.inOut->logMain,std::cerr,P);
    };    
    
    //output genome metadata
    writeChrInfo(pGe.gDir);

    //preparing to generate SA
    for (uint ii=0;ii<nGenome;ii++) {//- strand
        //if (G[ii]>5)
        //    cerr << ii <<" "<< G[ii]<<"\n";
        G[2*nGenome-1-ii]=G[ii]<4 ? 3-G[ii] : G[ii];
    };   
    nSA=0;
    for (uint ii=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {
        if (G[ii]<4) {
            nSA++;
        };
    };

    // GstrandBit
    GstrandBit = (char) (uint) floor(log(nGenome+P.limitSjdbInsertNsj*sjdbLength)/log(2))+1; //GstrandBit uses P.limitSjdbInsertNsj even if no insertion requested, in case it will be requested at the mapping stage
    if (GstrandBit<32) GstrandBit=32; //TODO: should not this be 31? Need to test for small genomes. TODO: use simple access function for SA
    P.inOut->logMain <<"Estimated genome size with padding and SJs: total=genome+SJ="<<nGenome+P.limitSjdbInsertNsj*sjdbLength<<" = "<<nGenome<<" + "<<P.limitSjdbInsertNsj*sjdbLength<<"\n";
    P.inOut->logMain << "GstrandBit=" << int(GstrandBit) <<"\n";
    GstrandMask = ~(1LLU<<GstrandBit);
    SA.defineBits(GstrandBit+1,nSA);

    if (P.sjdbInsert.yes) {//reserve space for junction insertion       
        SApass1.defineBits( GstrandBit+1, nSA+2*sjdbLength*min((uint64)sjdbLoci.chr.size(),P.limitSjdbInsertNsj) );//TODO: this allocation is wasteful, get a better estimate of the number of junctions
    } else {//same as SA
        SApass1.defineBits(GstrandBit+1,nSA);
    };

    P.inOut->logMain  << "Number of SA indices: "<< nSA << "\n"<<flush;

    //sort SA
    time ( &rawTime );
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to sort Suffix Array. This may take a long time...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to sort Suffix Array. This may take a long time...\n" <<flush;


//     if (false)
    {//sort SA chunks

        for (uint ii=0;ii<nGenome;ii++) {//re-fill the array backwards for sorting
            swap(G[2*nGenome-1-ii],G[ii]);
        };
        globalG=G;
        globalL=pGe.gSuffixLengthMax/sizeof(uint);
        //count the number of indices with 4nt prefix
        uint indPrefN=1LLU << 16;
        uint* indPrefCount = new uint [indPrefN];
        memset(indPrefCount,0,indPrefN*sizeof(indPrefCount[0]));
        nSA=0;
        for (uint ii=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {
            if (G[ii]<4) {
                uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                indPrefCount[p1]++;
                nSA++;
            };
        };

        uint saChunkSize=(P.limitGenomeGenerateRAM-nG1alloc)/8/P.runThreadN; //number of SA indexes per chunk
        saChunkSize=saChunkSize*6/10; //allow extra space for qsort
        //uint saChunkN=((nSA/saChunkSize+1)/P.runThreadN+1)*P.runThreadN;//ensure saChunkN is divisible by P.runThreadN
        //saChunkSize=nSA/saChunkN+100000;//final chunk size
        if (P.runThreadN>1) saChunkSize=min(saChunkSize,nSA/(P.runThreadN-1));
        uint saChunkN = nSA / saChunkSize + 1;//estimate
        uint* indPrefStart = new uint [saChunkN*2]; //start and stop, *2 just in case
        uint* indPrefChunkCount = new uint [saChunkN*2];
        indPrefStart[0]=0;
        saChunkN=0;//start counting chunks
        uint chunkSize1=indPrefCount[0];
        for (uint ii=1; ii<indPrefN; ii++) {
            chunkSize1 += indPrefCount[ii];
            if (chunkSize1 > saChunkSize) {
                saChunkN++;
                indPrefStart[saChunkN]=ii;
                indPrefChunkCount[saChunkN-1]=chunkSize1-indPrefCount[ii];
                chunkSize1=indPrefCount[ii];
            };
        };
        saChunkN++;
        indPrefStart[saChunkN]=indPrefN+1;
        indPrefChunkCount[saChunkN-1]=chunkSize1;

        P.inOut->logMain  << "Number of chunks: " << saChunkN <<";   chunks size limit: " << saChunkSize*8 <<" bytes\n" <<flush;

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;

        #pragma omp parallel for num_threads(P.runThreadN) ordered schedule(dynamic,1)
        for (int iChunk=0; iChunk < (int) saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
            uint* saChunk=new uint [indPrefChunkCount[iChunk]];//allocate local array for each chunk
            for (uint ii=0,jj=0;ii<2*nGenome;ii+=pGe.gSAsparseD) {//fill the chunk with SA indices
                if (G[ii]<4) {
                    uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                    if (p1>=indPrefStart[iChunk] && p1<indPrefStart[iChunk+1]) {
                        saChunk[jj]=ii;
                        jj++;
                    };
                    //TODO: if (jj==indPrefChunkCount[iChunk]) break;
                };
            };


            //sort the chunk
            qsort(saChunk,indPrefChunkCount[iChunk],sizeof(saChunk[0]),funCompareSuffixes);
            for (uint ii=0;ii<indPrefChunkCount[iChunk];ii++) {
                saChunk[ii]=2*nGenome-1-saChunk[ii];
            };
            //write files
            string chunkFileName=pGe.gDir+"/SA_"+to_string( (uint) iChunk);
            ofstream & saChunkFile = ofstrOpen(chunkFileName,ERROR_OUT, P);
            fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk],chunkFileName,ERROR_OUT,P);
            saChunkFile.close();
            delete [] saChunk;
            saChunk=NULL;
        };

        time ( &rawTime );
        P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;
        *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;

        //read chunks and pack into full SA
        SApass1.allocateArray();
        SA.pointArray(SApass1.charArray + SApass1.lengthByte-SA.lengthByte); //SA is shifted to have space for junction insertion
        uint N2bit= 1LLU << GstrandBit;
        uint packedInd=0;

        #define SA_CHUNK_BLOCK_SIZE 10000000
        uint* saIn=new uint[SA_CHUNK_BLOCK_SIZE]; //TODO make adjustable

        #ifdef genenomeGenerate_SA_textOutput
                ofstream SAtxtStream ((pGe.gDir + "/SAtxt").c_str());
        #endif

        for (uint iChunk=0;iChunk<saChunkN;iChunk++) {//load files one by one and convert to packed
            ostringstream saChunkFileNameStream("");
            saChunkFileNameStream<< pGe.gDir << "/SA_" << iChunk;
            ifstream saChunkFile(saChunkFileNameStream.str().c_str());
            while (! saChunkFile.eof()) {//read blocks from each file
                uint chunkBytesN=fstreamReadBig(saChunkFile,(char*) saIn,SA_CHUNK_BLOCK_SIZE*sizeof(saIn[0]));
                for (uint ii=0;ii<chunkBytesN/sizeof(saIn[0]);ii++) {
                    SA.writePacked( packedInd+ii, (saIn[ii]<nGenome) ? saIn[ii] : ( (saIn[ii]-nGenome) | N2bit ) );

                    #ifdef genenomeGenerate_SA_textOutput
                        SAtxtStream << saIn[ii] << "\n";
                    #endif
                };
                packedInd += chunkBytesN/sizeof(saIn[0]);
            };
            saChunkFile.close();
            remove(saChunkFileNameStream.str().c_str());//remove the chunk file
        };

        #ifdef genenomeGenerate_SA_textOutput
                SAtxtStream.close();
        #endif
        delete [] saIn;

        if (packedInd != nSA ) {//
            ostringstream errOut;
            errOut << "EXITING because of FATAL problem while generating the suffix array\n";
            errOut << "The number of indices read from chunks = "<<packedInd<<" is not equal to expected nSA="<<nSA<<"\n";
            errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        //DONE with suffix array generation

        for (uint ii=0;ii<nGenome;ii++) {//return to normal order for future use
            swap(G[2*nGenome-1-ii],G[ii]);
        };
        delete [] indPrefCount;
        delete [] indPrefStart;
        delete [] indPrefChunkCount;
    };

    time ( &rawTime );
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... finished generating suffix array\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... finished generating suffix array\n" <<flush;

    genomeSAindex(G, SA, P, SAi, *this);

    sjdbN=0;
    if (P.sjdbInsert.yes) {//insert junctions
        P.sjdbInsert.outDir=pGe.gDir;
        P.twoPass.pass2=false;

        Genome genome1(*this); //create copy here, *this will be changed by sjdbInsertJunctions
        sjdbInsertJunctions(P, *this, genome1, sjdbLoci);
    };

    pGe.gFileSizes.clear();
    pGe.gFileSizes.push_back(nGenome);
    pGe.gFileSizes.push_back(SA.lengthByte);

    //write genome parameters file
    genomeParametersWrite(pGe.gDir+("/genomeParameters.txt"), P, ERROR_OUT, *this);

    //write genome to disk
    time ( &rawTime );
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Genome to disk ...\n" <<flush;

    writeGenomeSequence(pGe.gDir);

    //write SA
    time ( &rawTime );
    P.inOut->logMain  << "SA size in bytes: "<<SA.lengthByte << "\n"<<flush;

    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;

    ofstream & SAout = ofstrOpen(pGe.gDir+"/SA",ERROR_OUT, P);
    fstreamWriteBig(SAout,(char*) SA.charArray, (streamsize) SA.lengthByte,pGe.gDir+"/SA",ERROR_OUT,P);
    SAout.close();

    //write SAi
    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;

    //write SAi to disk
    ofstream & SAiOut = ofstrOpen(pGe.gDir+"/SAindex",ERROR_OUT, P);

    fstreamWriteBig(SAiOut, (char*) &pGe.gSAindexNbases, sizeof(pGe.gSAindexNbases),pGe.gDir+"/SAindex",ERROR_OUT,P);
    fstreamWriteBig(SAiOut, (char*) genomeSAindexStart, sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1),pGe.gDir+"/SAindex",ERROR_OUT,P);
    fstreamWriteBig(SAiOut,  SAi.charArray, SAi.lengthByte,pGe.gDir+"/SAindex",ERROR_OUT,P);
    SAiOut.close();

    SApass1.deallocateArray();

    // Compute expected GC distribution for GC bias correction
    if (P.pGe.sjdbGTFfile != "-" || P.pGe.transcriptomeFasta != "-") {
        time(&rawTime);
        P.inOut->logMain    << timeMonthDayTime(rawTime) << " ... computing expected GC distribution ...\n" <<flush;
        *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ... computing expected GC distribution ...\n" <<flush;
        
        string compute_gc_cmd;
        // Try to find compute_expected_gc in common locations
        string tool_paths[] = {
            "compute_expected_gc",  // In PATH
            "../tools/compute_expected_gc/compute_expected_gc",  // Relative to STAR binary
            "../../tools/compute_expected_gc/compute_expected_gc",  // From build directory
            P.pGe.gDir + "/../tools/compute_expected_gc/compute_expected_gc"  // Relative to genome dir
        };
        string tool_path;
        bool tool_found = false;
        
        // Check which tool path exists
        for (int i = 0; i < 4; i++) {
            string test_cmd = "which " + tool_paths[i] + " > /dev/null 2>&1";
            if (i > 0) {
                test_cmd = "test -x " + tool_paths[i] + " > /dev/null 2>&1";
            }
            if (system(test_cmd.c_str()) == 0) {
                tool_path = tool_paths[i];
                tool_found = true;
                break;
            }
        }
        
        if (!tool_found) {
            P.inOut->logMain << "WARNING: compute_expected_gc tool not found. Skipping expected GC computation.\n" <<flush;
            P.inOut->logMain << "         GC bias correction will not be available. Install compute_expected_gc to enable.\n" <<flush;
        } else {
            string output_file = P.pGe.gDir + "/expected_gc.tsv";
            
            if (P.pGe.transcriptomeFasta != "-") {
                // Use provided transcriptome FASTA
                compute_gc_cmd = tool_path + " --transcriptome " + P.pGe.transcriptomeFasta 
                               + " --output " + output_file + " 2>&1";
            } else if (P.pGe.sjdbGTFfile != "-" && P.pGe.gFastaFiles.size() > 0) {
                // Extract from genome + GTF
                string genome_file = P.pGe.gFastaFiles[0];
                compute_gc_cmd = tool_path + " --genome " + genome_file 
                               + " --gtf " + P.pGe.sjdbGTFfile 
                               + " --output " + output_file + " 2>&1";
            }
            
            if (!compute_gc_cmd.empty()) {
                int ret = system(compute_gc_cmd.c_str());
                if (ret != 0) {
                    P.inOut->logMain << "WARNING: compute_expected_gc failed (exit code " << ret << "). GC bias correction may not be available.\n" <<flush;
                } else {
                    P.inOut->logMain << "Successfully computed expected GC distribution: " << output_file << "\n" <<flush;
                }
            }
        }
    }

    time(&rawTime);
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) << " ..... finished successfully\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... finished successfully\n" <<flush;
};

void Genome::writeChrInfo(const string dirOut) 
{//write chr information
    ofstream & chrN = ofstrOpen(dirOut+"/chrName.txt",ERROR_OUT, P);
    ofstream & chrS = ofstrOpen(dirOut+"/chrStart.txt",ERROR_OUT, P);
    ofstream & chrL = ofstrOpen(dirOut+"/chrLength.txt",ERROR_OUT, P);
    ofstream & chrNL = ofstrOpen(dirOut+"/chrNameLength.txt",ERROR_OUT, P);

    for (uint ii=0;ii<nChrReal;ii++) {//output names, starts, lengths
        chrN<<chrName[ii]<<"\n";
        chrS<<chrStart[ii]<<"\n";
        chrL<<chrLength.at(ii)<<"\n";
        chrNL<<chrName[ii]<<"\t"<<chrLength.at(ii)<<"\n";
    };
    chrS<<chrStart[nChrReal]<<"\n";//size of the genome
    chrN.close();chrL.close();chrS.close(); chrNL.close();
};
void Genome::writeGenomeSequence(const string dirOut) 
{//write genome sequence
    ofstream &genomeOut = ofstrOpen(dirOut+"/Genome",ERROR_OUT, P);
    fstreamWriteBig(genomeOut,G,nGenome,dirOut+"/Genome",ERROR_OUT,P);
    genomeOut.close();
};
