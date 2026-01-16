#include "ParametersGenome.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "CellRangerFormatter.h"
#include <algorithm>
#include <cctype>

void ParametersGenome::initialize(Parameters *pPin)
{
    pP=pPin;
    
    if (gDir.back()!='/') {
        gDir += '/';
    };
    
    //genome transformation
    if (transform.typeString=="None") {
        transform.type=0;
    } else if (transform.typeString=="Haploid") {
        transform.type=1;
    } else if (transform.typeString=="Diploid") {
        transform.type=2;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformType = " << transform.typeString << "\n";
        errOut << "SOLUTION: use one of the allowed values for --outTransformType : 'None' or 'Haploid' or 'Diploid' \n";
        exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    transform.outYes = transform.outSAM = transform.outSJ = false;
    if (transform.output.at(0) == "None") {
        //nothing to do
    } else {
        for (auto &ot: transform.output) {
            if (ot == "SAM") {
                transform.outYes = transform.outSAM = true;
            } else if (ot == "SJ") {
                transform.outYes = transform.outSJ = true;
            } else if (ot == "Quant") {
                transform.outYes = transform.outQuant = true;                
            } else {
                exitWithError("EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformOutput = " + ot + '\n'
                              + "SOLUTION: use allowed values for --outTransformOutput: None or SAM and/or SJ\n"
                              ,std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
        };
    };
    
    if (gTypeString!="Full" && gTypeString!="Transcriptome" && gTypeString!="SuperTranscriptome") {
        ostringstream errOut;
        errOut << "EXITING because of FATAL parameter error: --genomeType=" << gTypeString << "\n";
        errOut << "SOLUTION: use one of the allowed values for --genomeLoad : Full OR Transcriptome OR SuperTranscriptome\n" <<flush;
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    if (gLoad!="LoadAndKeep" && gLoad!="LoadAndRemove" && gLoad!="Remove" && gLoad!="LoadAndExit" && gLoad!="NoSharedMemory") {// find shared memory fragment
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: --genomeLoad=" << gLoad << "\n" <<flush;
        errOut << "SOLUTION: use one of the allowed values for --genomeLoad : NoSharedMemory,LoadAndKeep,LoadAndRemove,LoadAndExit,Remove.\n" <<flush;
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    // Initialize flex gene probe settings
    flexGeneProbe.enabled = (flexGeneProbe.csvFile != "-" && !flexGeneProbe.csvFile.empty());
    
    // Parse removeDeprecated flag (case-insensitive Yes/No)
    auto isYes = [](const string &value) -> bool {
        string normalized = value;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return (normalized == "yes");
    };
    flexGeneProbe.removeDeprecatedBool = isYes(flexGeneProbe.removeDeprecated);

    // Initialize cellranger-style settings (case-insensitive Yes/No)
    cellrangerStyle.indexEnabledBool = isYes(cellrangerStyle.indexEnabled);
    cellrangerStyle.downloadOnlyBool = isYes(cellrangerStyle.downloadOnly);
    cellrangerStyle.allUntrustedUrlBool = isYes(cellrangerStyle.allUntrustedUrl);
    cellrangerStyle.autoCksumUpdateBool = isYes(cellrangerStyle.autoCksumUpdate);
    cellrangerStyle.replaceUnverifiableFilesBool = isYes(cellrangerStyle.replaceUnverifiableFiles);
    
    if (cellrangerStyle.cacheDir.empty() || cellrangerStyle.cacheDir == "-") {
        cellrangerStyle.cacheDir = "";
    }
    
    // Parse and normalize cellrangerRefRelease
    {
        string release = cellrangerStyle.refRelease;
        if (release.empty()) {
            release = "2024-A"; // Default
        }
        
        // Normalize to lowercase for comparison
        string normalized = release;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        
        // Remove hyphens and spaces for alias matching
        string alias = normalized;
        alias.erase(std::remove(alias.begin(), alias.end(), '-'), alias.end());
        alias.erase(std::remove(alias.begin(), alias.end(), ' '), alias.end());
        
        // Map aliases to canonical forms
        if (normalized == "2024-a" || alias == "2024a" || alias == "2024") {
            cellrangerStyle.refReleaseCanonical = "2024-A";
        } else if (normalized == "2020-a" || alias == "2020a" || alias == "2020") {
            cellrangerStyle.refReleaseCanonical = "2020-A";
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: unrecognized option in --cellrangerRefRelease=" << release << "\n";
            errOut << "SOLUTION: use allowed option: 2024-A (or aliases: 2024A, 2024) OR 2020-A (or aliases: 2020A, 2020)\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        }
        
        // Store normalized input for logging
        cellrangerStyle.refRelease = cellrangerStyle.refReleaseCanonical;
    }
    
    // Initialize auto-index workflow settings (case-insensitive Yes/No)
    autoIndexWorkflow.autoIndexBool = isYes(autoIndexWorkflow.autoIndex);
    autoIndexWorkflow.forceIndexBool = isYes(autoIndexWorkflow.forceIndex);
    autoIndexWorkflow.forceAllIndexBool = isYes(autoIndexWorkflow.forceAllIndex);
    
    // Validate: forceAllIndex implies forceIndex
    if (autoIndexWorkflow.forceAllIndexBool && !autoIndexWorkflow.forceIndexBool) {
        autoIndexWorkflow.forceIndexBool = true; // forceAllIndex implies forceIndex
    }
    
    // Validate: if either faUrl or gtfUrl is set, both must be set
    bool hasFaUrl = !cellrangerStyle.faUrl.empty() && cellrangerStyle.faUrl != "-";
    bool hasGtfUrl = !cellrangerStyle.gtfUrl.empty() && cellrangerStyle.gtfUrl != "-";
    if (hasFaUrl != hasGtfUrl) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL PARAMETER ERROR: --faUrl and --gtfUrl must both be provided if either is set\n";
        exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    }
    
    // Validate: if URLs are provided and not trusted/have known cksum, require allUntrustedUrl flag
    if (hasFaUrl && hasGtfUrl) {
        // Check if URLs are trusted or have known cksum (will be checked again during download, but validate here too)
        bool faTrusted = CellRangerFormatter::isTrustedOrHasKnownCksum(cellrangerStyle.faUrl);
        bool gtfTrusted = CellRangerFormatter::isTrustedOrHasKnownCksum(cellrangerStyle.gtfUrl);
        
        if ((!faTrusted || !gtfTrusted) && !cellrangerStyle.allUntrustedUrlBool) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL PARAMETER ERROR: One or more URLs are not trusted and have no known cksum.\n"
                   << "  FASTA URL: " << cellrangerStyle.faUrl << (faTrusted ? " (trusted or has known cksum)" : " (NOT trusted, no known cksum)") << "\n"
                   << "  GTF URL: " << cellrangerStyle.gtfUrl << (gtfTrusted ? " (trusted or has known cksum)" : " (NOT trusted, no known cksum)") << "\n"
                   << "Untrusted URLs require --allUntrustedUrl Yes flag.\n"
                   << "WARNING: Using untrusted URLs disables integrity checking. Only use with trusted sources.\n";
            exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
        }
    }
    
    // Validate download-only requires URLs
    if (cellrangerStyle.downloadOnlyBool && (!hasFaUrl || !hasGtfUrl)) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL PARAMETER ERROR: --cellrangerStyleDownloadOnly requires --faUrl and --gtfUrl\n";
        exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    }
    
    // Initialize transcriptome generation settings (case-insensitive Yes/No)
    transcriptomeGen.generateTranscriptomeBool = isYes(transcriptomeGen.generateTranscriptome);
    transcriptomeGen.overwriteBool = isYes(transcriptomeGen.overwrite);
    
    // Set default transcriptome FASTA path if not provided
    if (transcriptomeGen.transcriptomeFastaPath.empty() || transcriptomeGen.transcriptomeFastaPath == "-") {
        transcriptomeGen.transcriptomeFastaPath = gDir + "transcriptome.fa";
    }
};
