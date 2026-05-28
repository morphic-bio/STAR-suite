#include "ReadAlignChunk.h"
#include "GlobalVariables.h"
#include "ThreadControl.h"
#include "ErrorWarning.h"
#include "input/CbqStarAdapter.h"
#include SAMTOOLS_BGZF_H

void ReadAlignChunk::mapChunk() {//map one chunk. Input reads stream has to be setup in RA->readInStream[ii]

    for (uint32 im=0; im<1; im++) {//hardcoded mate 1 5p onyl for now
        RA->clipMates[im][0].clipChunk(chunkIn[im], chunkInSizeBytesTotal[im]);
    };
    
    RA->statsRA.resetN();

    for (uint ii=0;ii<P.readNends;ii++) {//clear eof and rewind the input streams
        RA->readInStream[ii]->clear();
        RA->readInStream[ii]->seekg(0,ios::beg);
    };
    
    

    if ( P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {//open chunk file
        ostringstream name1("");
        name1 << P.outFileTmp + "/Aligned.tmp.sam.chunk"<<iChunkIn;
        chunkOutBAMfileName = name1.str();
        chunkOutBAMfile.open(chunkOutBAMfileName.c_str());
    };

    int readStatus=0;
    while (readStatus==0) {//main cycle over all reads

        readStatus=RA->oneRead(); //map one read

        if (readStatus==0) {//there was a read processed
            RA->iRead++;
//         chunkOutBAMtotal=(uint) RA->outSAMstream->tellp();
            chunkOutBAMtotal+=RA->outBAMbytes;
//             uint ddd=(uint) RA->outSAMstream->tellp();
        };

        //write SAM aligns to chunk buffer
        if (P.outSAMbool) {
            if ( chunkOutBAMtotal > P.chunkOutBAMsizeBytes ) {//this should not happen!
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SAM/BAM output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSAMoneReadBytes\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            } else if ( chunkOutBAMtotal + P.limitOutSAMoneReadBytes > P.chunkOutBAMsizeBytes || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
                if ( P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {//output chunks into separate files
                    chunkOutBAMfile.write(chunkOutBAM,chunkOutBAMtotal);
                    chunkOutBAMfile.clear(); //in case 0 bytes were written which could set fail bit
                    //chunkOutBAMfile.flush(); //not needed
                } else {//standard way, directly into Aligned.out.sam file
                    //SAM output
                    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
                    P.inOut->outSAM->write(chunkOutBAM,chunkOutBAMtotal);
                    P.inOut->outSAM->clear();//in case 0 bytes were written which could set fail bit
                    //P.inOut->outSAM->flush(); //not needed
                    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
                };
                RA->outSAMstream->seekp(0,ios::beg); //rewind the chunk storage
                chunkOutBAMtotal=0;
            };
        };

        //collapse SJ buffer if needed
        if ( !P.outSJ.yes ) {
            //do nothing
        } else if ( chunkOutSJ->N > chunkOutSJ->Nstore ) {//this means the number of collapsed junctions is larger than the chunks size
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if ( chunkOutSJ->N + P.limitOutSJoneRead > chunkOutSJ->Nstore || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
            chunkOutSJ->collapseSJ();
            if ( chunkOutSJ->N + 2*P.limitOutSJoneRead > chunkOutSJ->Nstore ) {
                /*
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                */
                chunkOutSJ->dataSizeIncrease();
                P.inOut->logMain << "Increased the size of chunkOutSJ to " << chunkOutSJ->Nstore <<'\n';
            };
        };

        //collapse SJ1 buffer if needed
        if ( P.outFilterBySJoutStage != 1 ) {//no outFilterBySJoutStage
            //do nothing
        } else if ( chunkOutSJ1->N > chunkOutSJ->Nstore ) {//this means the number of collapsed junctions is larger than the chunks size
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if ( chunkOutSJ1->N + P.limitOutSJoneRead > chunkOutSJ->Nstore || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
            chunkOutSJ1->collapseSJ();
            if ( chunkOutSJ1->N + 2*P.limitOutSJoneRead > chunkOutSJ->Nstore ) {
                /*
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
                */
                chunkOutSJ->dataSizeIncrease();
                P.inOut->logMain << "Increased the size of chunkOutSJ to " << chunkOutSJ->Nstore <<'\n';
            };
        };

    }; //reads cycle

    // Phase 2: Resolve accumulated ambiguous CBs using Bayesian inference
    RA->resolveAmbiguousCBs();

    if ( P.outSAMbool && P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {//write the remaining part of the buffer, close and rename chunk files
        chunkOutBAMfile.write(chunkOutBAM,chunkOutBAMtotal);
        chunkOutBAMfile.clear(); //in case 0 bytes were written which could set fail bit
        chunkOutBAMfile.close();
        RA->outSAMstream->seekp(0,ios::beg); //rewind the chunk storage
        chunkOutBAMtotal=0;
        ostringstream name2("");
        name2 << P.outFileTmp + "/Aligned.out.sam.chunk"<<iChunkIn;
        rename(chunkOutBAMfileName.c_str(),name2.str().c_str());//marks files as completedly written
    };

    //add stats, write progress if needed
    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexStats);
    g_statsAll.addStats(RA->statsRA);
    g_statsAll.progressReport(P.inOut->logProgress);
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexStats);
};

void ReadAlignChunk::mapCbqChunk() {//map one owned CBQ chunk through STAR read buffers

    RA->statsRA.resetN();

    if ( P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {//open chunk file
        ostringstream name1("");
        name1 << P.outFileTmp + "/Aligned.tmp.sam.chunk"<<iChunkIn;
        chunkOutBAMfileName = name1.str();
        chunkOutBAMfile.open(chunkOutBAMfileName.c_str());
    };

    auto flushOutputsIfNeeded = [&](int readStatus) {
        if (P.outSAMbool) {
            if ( chunkOutBAMtotal > P.chunkOutBAMsizeBytes ) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SAM/BAM output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSAMoneReadBytes\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            } else if ( chunkOutBAMtotal + P.limitOutSAMoneReadBytes > P.chunkOutBAMsizeBytes || (readStatus==-1 && noReadsLeft) ) {
                if ( P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {
                    chunkOutBAMfile.write(chunkOutBAM,chunkOutBAMtotal);
                    chunkOutBAMfile.clear();
                } else {
                    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
                    P.inOut->outSAM->write(chunkOutBAM,chunkOutBAMtotal);
                    P.inOut->outSAM->clear();
                    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
                };
                RA->outSAMstream->seekp(0,ios::beg);
                chunkOutBAMtotal=0;
            };
        };

        if ( !P.outSJ.yes ) {
        } else if ( chunkOutSJ->N > chunkOutSJ->Nstore ) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if ( chunkOutSJ->N + P.limitOutSJoneRead > chunkOutSJ->Nstore || (readStatus==-1 && noReadsLeft) ) {
            chunkOutSJ->collapseSJ();
            if ( chunkOutSJ->N + 2*P.limitOutSJoneRead > chunkOutSJ->Nstore ) {
                chunkOutSJ->dataSizeIncrease();
                P.inOut->logMain << "Increased the size of chunkOutSJ to " << chunkOutSJ->Nstore <<'\n';
            };
        };

        if ( P.outFilterBySJoutStage != 1 ) {
        } else if ( chunkOutSJ1->N > chunkOutSJ->Nstore ) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if ( chunkOutSJ1->N + P.limitOutSJoneRead > chunkOutSJ->Nstore || (readStatus==-1 && noReadsLeft) ) {
            chunkOutSJ1->collapseSJ();
            if ( chunkOutSJ1->N + 2*P.limitOutSJoneRead > chunkOutSJ->Nstore ) {
                chunkOutSJ->dataSizeIncrease();
                P.inOut->logMain << "Increased the size of chunkOutSJ to " << chunkOutSJ->Nstore <<'\n';
            };
        };
    };

    star::input::CbqStarAdapterOptions cbqOptions;
    cbqOptions.read_nends = P.readNends;
    cbqOptions.out_sam_read_id_number = P.outSAMreadIDnumber;
    cbqOptions.out_qs_conversion_add = P.outQSconversionAdd;
    cbqOptions.trim_cutadapt_enabled = (P.trimCutadapt == "Yes");
    cbqOptions.preserve_read_name_extra = false;

    star::input::CbqStarReadBuffers cbqBuffers;
    cbqBuffers.read_name_mates = RA->readNameMates;
    cbqBuffers.read0 = RA->Read0;
    cbqBuffers.read1 = RA->Read1;
    cbqBuffers.qual0 = RA->Qual0;
    cbqBuffers.read_name_extra = &RA->readNameExtra;
    cbqBuffers.read_length = RA->readLength;
    cbqBuffers.read_length_original = RA->readLengthOriginal;
    cbqBuffers.i_read_all = &RA->iReadAll;
    cbqBuffers.read_files_index = &RA->readFilesIndex;
    cbqBuffers.read_filter = &RA->readFilter;
    cbqBuffers.read_file_type = &RA->readFileType;
    string cbqLoadError;

    for (uint32 imate = 0; imate < P.readNends; ++imate) {
        cbqLoadError.clear();
        if (!star::input::prepare_cbq_star_chunk_clip_info(&cbqStarChunk,
                                                           imate,
                                                           &RA->clipMates[imate][0],
                                                           &cbqLoadError)) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in CBQ input adapter clip preparation\n";
            errOut << cbqLoadError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
    }

    for (size_t irecord = 0; irecord < cbqStarChunk.read_count(); ++irecord) {
        cbqLoadError.clear();
        if (!star::input::load_cbq_star_chunk_read_into_star_mates(cbqStarChunk,
                                                                   irecord,
                                                                   cbqOptions,
                                                                   &cbqBuffers,
                                                                   &RA->clipMates,
                                                                   &cbqLoadError)) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in CBQ input adapter\n";
            errOut << cbqLoadError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        }
        const int readStatus = RA->oneReadLoaded(RA->readFileType);

        if (readStatus==0) {
            RA->iRead++;
            chunkOutBAMtotal+=RA->outBAMbytes;
        };

        flushOutputsIfNeeded(readStatus);
    };

    flushOutputsIfNeeded(-1);

    if ( P.outSAMbool && P.outSAMorder == "PairedKeepInputOrder" && P.runThreadN>1 ) {
        chunkOutBAMfile.write(chunkOutBAM,chunkOutBAMtotal);
        chunkOutBAMfile.clear();
        chunkOutBAMfile.close();
        RA->outSAMstream->seekp(0,ios::beg);
        chunkOutBAMtotal=0;
        ostringstream name2("");
        name2 << P.outFileTmp + "/Aligned.out.sam.chunk"<<iChunkIn;
        rename(chunkOutBAMfileName.c_str(),name2.str().c_str());
    };

    if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexStats);
    g_statsAll.addStats(RA->statsRA);
    g_statsAll.progressReport(P.inOut->logProgress);
    if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexStats);
};
