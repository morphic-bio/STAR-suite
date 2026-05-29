#include "Parameters.h"
#include "ErrorWarning.h"
#include "input/CbqInputModule.h"
#include "input/FastxInputModule.h"
#include "streamFuns.h"
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
#include "serviceFuns.cpp"

namespace {
bool endsWithCaseInsensitive(const string& value, const string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    const size_t offset = value.size() - suffix.size();
    for (size_t ii = 0; ii < suffix.size(); ++ii) {
        const unsigned char c1 = static_cast<unsigned char>(value[offset + ii]);
        const unsigned char c2 = static_cast<unsigned char>(suffix[ii]);
        if (std::tolower(c1) != std::tolower(c2)) {
            return false;
        }
    }
    return true;
}

bool isFastqPath(const string& value) {
    return endsWithCaseInsensitive(value, ".fastq") ||
           endsWithCaseInsensitive(value, ".fq") ||
           endsWithCaseInsensitive(value, ".fastq.gz") ||
           endsWithCaseInsensitive(value, ".fq.gz");
}
}

void Parameters::readFilesInit() 
{//initialize read files - but do not open yet

    if (readFilesType.at(0) == "Fastx") {
        readFilesTypeN=1;
    } else if (readFilesType.at(0) == "Binseq") {
        readFilesTypeN=20;
    } else if (readFilesType.at(0) == "SAM"){
        readFilesTypeN=10;
        readFiles.samAttrKeepAll = false;
        readFiles.samAttrKeepNone = false;
        if (readFiles.samAttrKeepIn.at(0) == "All") {
            readFiles.samAttrKeepAll = true;
        } else if (readFiles.samAttrKeepIn.at(0) == "None") {
            readFiles.samAttrKeepNone = true;
        } else {
            for (auto &tag: readFiles.samAttrKeepIn) {
                if (tag.size()!=2) {
                    exitWithError("EXITING because of FATAL PARAMETER ERROR: each SAM tags in --readFilesSAMtagsKeep should contain two letters\n\
                                  SOLUTION: specify only two-letter tags in --readFilesSAMtagsKeep.",
                                  std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                };
                //array<char,2> taga = {tag[0], tag[1]};
                uint16_t tagn = * ( (uint16_t*) tag.c_str() );
                readFiles.samAttrKeep.insert(tagn);
            };
        };
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented value for --readFilesType: "<<readFilesType.at(0) <<"\n";
        errOut <<"SOLUTION: specify one of the allowed values: Fastx, Binseq, or SAM\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (emitYNoYFastqyes && readFilesTypeN != 1) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: --emitYNoYFastq currently requires FASTQ input through --readFilesType Fastx.\n";
        errOut << "SOLUTION: use FASTQ input files with --readFilesType Fastx or disable --emitYNoYFastq.\n";
        // TODO: Y-removal/Y-noY FASTQ emission is needed for non-FASTQ input formats.
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    }

    readFilesPrefixFinal=(readFilesPrefix=="-" ? "" : readFilesPrefix);
    
    if (readFilesManifest[0]=="-") {//no manifest, file names in readFilesIn
        readFilesNames.resize(readFilesIn.size());
        
        for (uint32 imate=0; imate<readFilesNames.size(); imate++) {
            splitString(readFilesIn[imate], ',', readFilesNames[imate]);
            if (readFilesNames[imate].back().empty()) {//extra comma at the end
                readFilesNames[imate].pop_back();
            };
        
            if (imate>0 && readFilesNames[imate].size() != readFilesNames[imate-1].size() ) {
                ostringstream errOut;
                errOut <<"EXITING: because of fatal INPUT ERROR: number of input files for mate" << imate+1 <<"="<< readFilesNames[imate].size()  <<" is not equal to that for mate"<< imate-1 <<"="<< readFilesNames[imate-1].size() <<"\n";
                errOut <<"Make sure that the number of files in --readFilesIn is the same for both mates\n";
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
            
            for ( auto &fn : readFilesNames[imate] )
                fn = readFilesPrefixFinal + fn; //add prefix
        };

        readFilesN = readFilesNames[0].size();

        //read groups
        if (outSAMattrRGline.at(0)!="-") {
            string linefull;
            for (uint ii=0;ii<outSAMattrRGline.size(); ii++) {//concatenate into one line
                if (ii==0 || outSAMattrRGline.at(ii)==",") {//start new entry
                    if (ii>0) ++ii;//skip comma
                    outSAMattrRGlineSplit.push_back(outSAMattrRGline.at(ii)); //start new RG line with the first field which must be ID:xxx
                    if (outSAMattrRGlineSplit.back().substr(0,3)!="ID:") {
                        ostringstream errOut;
                        errOut <<"EXITING because of FATAL INPUT ERROR: the first word of a line from --outSAMattrRGline="<<outSAMattrRGlineSplit.back()<<" does not start with ID:xxx read group identifier\n";
                        errOut <<"SOLUTION: re-run STAR with all lines in --outSAMattrRGline starting with ID:xxx\n";
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    };
                    outSAMattrRG.push_back(outSAMattrRGlineSplit.back().substr(3)); //this adds the ID field
                } else {//keep adding fields to this RG line, until the next comma
                    outSAMattrRGlineSplit.back()+="\t" + outSAMattrRGline.at(ii);
                };
            };
        };
        
        if (outSAMattrRG.size()>1 && outSAMattrRG.size()!=readFilesN) {
            ostringstream errOut;
            errOut <<"EXITING: because of fatal INPUT ERROR: number of input read files: "<< readFilesN << " does not agree with number of read group RG entries: "<< outSAMattrRG.size() <<"\n";
            errOut <<"Make sure that the number of RG lines in --outSAMattrRGline is equal to either 1, or the number of input read files in --readFilesIn\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (outSAMattrRG.size()==1) {//use the same read group for all files
            for (uint32 ifile=1; ifile<readFilesN; ifile++) {
                outSAMattrRG.push_back(outSAMattrRG.at(0));
            };
        };           
        
    } else {//read file names from manifest
        //TODO check that outSAMattrRGline and readFilesIn are not set, throw an error
        
        ifstream & rfM = ifstrOpen(readFilesManifest[0], ERROR_OUT, "SOLUTION: check the path and permissions for readFilesManifest = " + readFilesManifest[0], *this);
        inOut->logMain << "Reading input file names and read groups from readFileManifest " << readFilesManifest[0] << endl;

        readFilesNames.resize(2);
        string rfMline;
        while (getline(rfM, rfMline)) {
        	if (rfMline.find_first_not_of(" \t")>=rfMline.size())
        		continue; //skip blank lines

            uint32 itab1=0, itab2=0;
            for (uint32 imate=0; imate<2; imate++) {//SE manifest 2nd column contains "-"
                itab2=rfMline.find('\t',itab1);
                if (itab2>=rfMline.size()) {
                    ostringstream errOut;
                    errOut <<"EXITING because of FATAL INPUT FILE error: readFileManifest file " << readFilesManifest[0] <<  " has to contain at least 3 tab separated columns\n";
                    errOut <<"SOLUTION: fix the formatting of the readFileManifest file: Read1 <tab> Read2 <tab> ReadGroup. For single-end reads, use - in the 2nd column.\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
                };
                readFilesNames[imate].push_back( readFilesPrefixFinal + rfMline.substr(itab1,itab2-itab1) );
                itab1=itab2+1;
                
                inOut->logMain << readFilesNames[imate].back() <<'\t';
            };
            
            outSAMattrRGlineSplit.push_back(rfMline.substr(itab2+1));
            
            if (outSAMattrRGlineSplit.back().substr(0,3)!="ID:")
                outSAMattrRGlineSplit.back().insert(0,"ID:");
            
            itab2=outSAMattrRGlineSplit.back().find('\t');
            outSAMattrRG.push_back(outSAMattrRGlineSplit.back().substr(3,itab2-3));
            
            inOut->logMain <<  outSAMattrRGlineSplit.back() <<'\n';
            
        };
        rfM.close();
        
        if (readFilesTypeN == 20) {
            for (const auto& mate2Name : readFilesNames[1]) {
                if (mate2Name != "-") {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL INPUT FILE error: --readFilesType Binseq expects one CBQ file per lane.\n";
                    errOut << "In readFilesManifest, the second column must be '-' for Binseq input.\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
                }
            }
            readFilesNames.resize(1);
        } else {
            readNends = ( readFilesNames[1][0].back()=='-' ? 1 : 2);
            readFilesNames.resize(readNends);//resize if readFilesN=1
        }
        readFilesN = readFilesNames[0].size();
    };

    inOut->logMain << "Number of input files for each source = " << readFilesN << endl;
    
    // Parse legacy gzip override
    {
        string legacy = readFilesLegacyZcatStr;
        std::transform(legacy.begin(), legacy.end(), legacy.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (legacy == "yes") {
            readFilesLegacyZcat = true;
        } else if (legacy == "no" || legacy == "-" || legacy.empty()) {
            readFilesLegacyZcat = false;
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --readFilesLegacyZcat=" << readFilesLegacyZcatStr << "\n";
            errOut << "SOLUTION: use one of the allowed values: Yes or No\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    }

    const bool hasReadFilesCommand = (readFilesCommand.at(0) != "-");
    if (readFilesTypeN == 20 && hasReadFilesCommand) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: --readFilesCommand is not supported with --readFilesType Binseq.\n";
        errOut << "SOLUTION: pass CBQ files directly through --readFilesIn or --readFilesManifest.\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    }
    string explicitReadFilesCommandLower;
    if (hasReadFilesCommand) {
        for (const auto& token : readFilesCommand) {
            explicitReadFilesCommandLower += token;
            explicitReadFilesCommandLower.push_back(' ');
        }
        std::transform(explicitReadFilesCommandLower.begin(), explicitReadFilesCommandLower.end(),
                       explicitReadFilesCommandLower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    }
    const bool explicitLegacyZcatCmd =
        hasReadFilesCommand &&
        (explicitReadFilesCommandLower.find("zcat") != string::npos ||
         explicitReadFilesCommandLower.find("gzcat") != string::npos ||
         explicitReadFilesCommandLower.find("gunzip") != string::npos);
    bool allFastqGz = (readFilesTypeN == 1);
    if (allFastqGz) {
        for (const auto& mateFiles : readFilesNames) {
            for (const auto& fn : mateFiles) {
                if (!endsWithCaseInsensitive(fn, ".gz")) {
                    allFastqGz = false;
                    break;
                }
            }
            if (!allFastqGz) {
                break;
            }
        }
    }

    // Default behavior: when readFilesCommand is '-' and all Fastx inputs are .gz,
    // use internal zlib streaming path unless user forces legacy zcat helper mode.
    readFilesUseInternalGzip = (!hasReadFilesCommand && readFilesTypeN == 1 && allFastqGz && !readFilesLegacyZcat);

    readFilesCommandString = "";
    if (hasReadFilesCommand) {
        for (uint ii = 0; ii < readFilesCommand.size(); ii++) {
            readFilesCommandString += readFilesCommand.at(ii) + "   "; // concatenate into one string
        }
    } else if (readFilesUseInternalGzip) {
        // Non-empty sentinel to route openReadsFiles() through FIFO helper branch,
        // where internal zlib streaming is handled.
        readFilesCommandString = "INTERNAL_GZIP";
    } else if (readFilesTypeN == 1 && allFastqGz && readFilesLegacyZcat) {
        // Compatibility override: keep legacy FIFO + external zcat behavior.
        readFilesCommandString = "zcat   ";
    } else if (readFilesN > 1 && readFilesTypeN != 20) {
        readFilesCommandString = "cat   "; // concatenate multiple plain-text files
    }

    if (readFilesUseInternalGzip) {
        inOut->logMain << "NOTE: .gz Fastx input detected with --readFilesCommand - ; using internal zlib streaming path.\n"
                       << "      Set --readFilesLegacyZcat Yes to force legacy external zcat helper mode.\n";
    } else if (hasReadFilesCommand && readFilesLegacyZcat) {
        if (explicitLegacyZcatCmd) {
            inOut->logMain << "NOTE: --readFilesLegacyZcat Yes is redundant because --readFilesCommand already uses zcat/gunzip.\n"
                           << "      Continuing with the explicit --readFilesCommand path.\n";
        } else {
            inOut->logMain << "WARNING: --readFilesLegacyZcat Yes is ignored because an explicit --readFilesCommand is set.\n"
                           << "         Continuing with explicit --readFilesCommand: " << readFilesCommandString << "\n";
        }
    } else if (!hasReadFilesCommand && readFilesTypeN == 1 && allFastqGz && readFilesLegacyZcat) {
        inOut->logMain << "NOTE: --readFilesLegacyZcat Yes: forcing legacy external zcat helper for .gz Fastx input.\n";
    }
    
    if (readFilesTypeN==1) {
        readNends=readFilesNames.size(); //for now the number of mates is defined by the number of input files
    } else if (readFilesTypeN==20) {
        if (readFilesType.size()==2 && readFilesType.at(1)=="SE") {
            readNends=1;
        } else if (readFilesType.size()==2 && readFilesType.at(1)=="PE") {
            readNends=2;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --readFilesType Binseq requires specifying SE or PE reads"<<"\n";
            errOut <<"SOLUTION: specify --readFilesType Binseq SE for single-end CBQ or --readFilesType Binseq PE for paired-end CBQ\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        if (readFilesNames.size() != 1) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: --readFilesType Binseq expects one CBQ source list in --readFilesIn.\n";
            errOut << "SOLUTION: pass one comma-separated CBQ list, e.g. --readFilesType Binseq PE --readFilesIn lane1.cbq,lane2.cbq\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        }
    } else if (readFilesTypeN==10) {//find the number of mates from the SAM file
        if (readFilesType.size()==2 && readFilesType.at(1)=="SE") {
            readNends=1;
        } else if (readFilesType.size()==2 && readFilesType.at(1)=="PE") {
            readNends=2;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --readFilesType SAM requires specifying SE or PE reads"<<"\n";
            errOut <<"SOLUTION: specify --readFilesType SAM SE for single-end reads or --readFilesType SAM PE for paired-end reads\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };
    
    readNmates=readNends; //this may be changed later if one of the reads is barcode rea

    if (readFilesTypeN==1) {
        if (emitYNoYFastqyes) {
            for (uint32 imate = 0; imate < readFilesNames.size(); ++imate) {
                for (const auto& fastxName : readFilesNames[imate]) {
                    if (!isFastqPath(fastxName)) {
                        ostringstream errOut;
                        errOut << "EXITING because of FATAL INPUT ERROR: --emitYNoYFastq currently requires FASTQ input files.\n";
                        errOut << "Offending --readFilesIn entry: " << fastxName << "\n";
                        errOut << "SOLUTION: provide .fastq/.fq files, optionally gzip-compressed as .fastq.gz/.fq.gz, or disable --emitYNoYFastq.\n";
                        // TODO: Y-removal/Y-noY FASTQ emission is needed for other input formats.
                        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
                    }
                }
            }
        }

        star::input::InputSourcePlan fastxInputPlan =
            star::input::make_fastx_input_source_plan(
                readFilesNames,
                outSAMattrRG,
                readFilesCommandString,
                readFilesPrefixFinal,
                readFilesUseInternalGzip);
        fastxInputActive = false;
        fastxInputPendingRecordValid = false;
        fastxInputExhausted = false;
        fastxInputLastLoggedLane = -1;
        fastxInputPendingRecord.reset();
        fastxInputModule.reset(new star::input::FastxInputModule());
        string inputContractError;
        if (!fastxInputModule->configure(fastxInputPlan, &inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: invalid Fastx input source plan\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        fastxInputActive = true;
    } else if (readFilesTypeN==20) {
        vector<string> cbqReadGroups;
        if (!readFilesNames.empty() && outSAMattrRG.size() == readFilesNames.front().size()) {
            cbqReadGroups = outSAMattrRG;
        }
        star::input::InputSourcePlan cbqInputPlan =
            star::input::make_cbq_input_source_plan(
                readFilesNames,
                cbqReadGroups,
                readNends);
        cbqInputActive = false;
        cbqInputExhausted = false;
        cbqInputLastLoggedLane = -1;
        cbqInputPendingBatch.reset();
        cbqInputPendingBatchOffset = 0;
        cbqInputModule.reset(new star::input::CbqInputModule());
        string inputContractError;
        if (!cbqInputModule->configure(cbqInputPlan, &inputContractError)) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT ERROR: invalid Binseq/CBQ input source plan\n";
            errOut << inputContractError << "\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
        cbqInputActive = true;
    };
};
