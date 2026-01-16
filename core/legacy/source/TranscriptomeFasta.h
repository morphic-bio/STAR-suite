#ifndef H_TranscriptomeFasta
#define H_TranscriptomeFasta

#include "IncludeDefine.h"
#include "Genome.h"
#include "GTF.h"

class TranscriptomeFasta {
public:
    /**
     * Write transcriptome FASTA file from Genome and GTF data.
     * 
     * @param genome Reference to Genome object containing sequence data
     * @param gtf Reference to GTF object containing transcript annotations
     * @param outPath Output file path for transcriptome FASTA
     * @param P Reference to Parameters for logging
     * @param lineWidth Line width for FASTA output (default: 70)
     * @param overwrite If true, overwrite existing file; if false, skip if exists
     * @return true if successful, false otherwise
     */
    static bool writeTranscriptomeFasta(const Genome& genome, const GTF& gtf, 
                                        const std::string& outPath,
                                        Parameters& P,
                                        uint32_t lineWidth = 70,
                                        bool overwrite = false);
};

#endif
