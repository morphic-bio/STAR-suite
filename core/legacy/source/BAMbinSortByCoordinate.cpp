#include "BAMbinSortByCoordinate.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "BAMfunctions.h"
#include "SequenceFuns.h"

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &genome, Solo &solo) {

    if (binS==0) return; //nothing to do for empty bins
    //allocate arrays
    char *bamIn=new char[binS+1];
    // Use uint64 for iRead field to accommodate Y-bit in bit 63
    uint64 *startPos=new uint64[binN*3];

    uint bamInBytes=0;
    //load all aligns
    for (uint it=0; it<nThreads; it++) {
        string bamInFile=dirBAMsort+to_string(it)+"/"+to_string((uint) iBin);
        ifstream bamInStream;
        bamInStream.open(bamInFile.c_str(),std::ios::binary | std::ios::ate);//open at the end to get file size
        int64 s1=bamInStream.tellg();
        if (s1>0)         {
            bamInStream.seekg(std::ios::beg);
            bamInStream.read(bamIn+bamInBytes,s1);//read the whole file
        } else if (s1<0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: failed reading from temporary file: " << dirBAMsort+to_string(it)+"/"+to_string((uint) iBin);
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        bamInBytes += bamInStream.gcount();
        bamInStream.close();
        remove(bamInFile.c_str());
    };
    if (bamInBytes!=binS) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
        errOut << "Expected bin size=" <<binS <<" ; size on disk="<< bamInBytes <<" ; bin number="<< iBin <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
    };

    //extract coordinates

    for (uint ib=0,ia=0;ia<binN;ia++) {
        uint32 *bamIn32=(uint32*) (bamIn+ib);
        startPos[ia*3]  =( ((uint64) bamIn32[1]) << 32) | ( (uint64)bamIn32[2] );
        startPos[ia*3+2]=ib;
        ib+=bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        startPos[ia*3+1]=*( (uint64*) (bamIn+ib) ); //read order with Y-bit in bit 63
        ib+=sizeof(uint64);
    };

    //sort (using custom comparison that ignores Y-bit for sorting)
    qsort((void*) startPos, binN, sizeof(uint64)*3, funCompareArraysIgnoreYBit);

    // Open output handles
    BGZF *bgzfBin = nullptr;
    BGZF *bgzfBin_Y = nullptr;
    BGZF *bgzfBin_noY = nullptr;
    
    if (P.emitNoYBAMyes) {
        // Open Y and noY bin files
        bgzfBin_Y=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)+"_Y").c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
        if (bgzfBin_Y==NULL) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not open temporary Y bam file: " << dirBAMsort+"/b"+to_string((uint) iBin)+"_Y" << "\n";
            errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        };
        bgzfBin_noY=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)+"_noY").c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
        if (bgzfBin_noY==NULL) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not open temporary noY bam file: " << dirBAMsort+"/b"+to_string((uint) iBin)+"_noY" << "\n";
            errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        };
        outBAMwriteHeader(bgzfBin_Y,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
        outBAMwriteHeader(bgzfBin_noY,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
        
        // Open primary bin file if keepBAM is enabled
        if (P.keepBAMyes) {
            bgzfBin=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)).c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
            if (bgzfBin==NULL) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal ERROR: could not open temporary bam file: " << dirBAMsort+"/b"+to_string((uint) iBin) << "\n";
                errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            };
            outBAMwriteHeader(bgzfBin,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
        }
    } else {
        // Normal mode: open primary bin file
        bgzfBin=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)).c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
        if (bgzfBin==NULL) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal ERROR: could not open temporary bam file: " << dirBAMsort+"/b"+to_string((uint) iBin) << "\n";
            errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        };
        outBAMwriteHeader(bgzfBin,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
    }
    
    //send ordered aligns to bgzf one-by-one
    char bam1[BAM_ATTR_MaxSize];//temp array
    for (uint ia=0;ia<binN;ia++) {
        char* bam0=bamIn+startPos[ia*3+2];
        uint32 size0=*((uint32*) bam0)+sizeof(uint32);
        
        if (solo.pSolo.samAttrYes)
            solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->addBAMtags(bam0,size0,bam1);
        
        if (P.emitNoYBAMyes) {
            // Extract Y-bit from iRead
            uint64 iReadWithY = startPos[ia*3+1];
            bool hasY = (iReadWithY >> 63) & 1;
            
            // Route to appropriate handle
            BGZF *targetHandle = hasY ? bgzfBin_Y : bgzfBin_noY;
            bgzf_write(targetHandle, bam0, size0);
            
            // Also write to primary if enabled
            if (bgzfBin != nullptr) {
                bgzf_write(bgzfBin, bam0, size0);
            }
        } else {
            bgzf_write(bgzfBin, bam0, size0);
        }
    };

    if (bgzfBin != nullptr) {
        bgzf_flush(bgzfBin);
        bgzf_close(bgzfBin);
    }
    if (bgzfBin_Y != nullptr) {
        bgzf_flush(bgzfBin_Y);
        bgzf_close(bgzfBin_Y);
    }
    if (bgzfBin_noY != nullptr) {
        bgzf_flush(bgzfBin_noY);
        bgzf_close(bgzfBin_noY);
    }
    //release memory
    delete [] bamIn;
    delete [] startPos;
};
