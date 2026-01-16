/*
 * sample_fld - Extract fragment length distribution from paired-end BAM file
 * 
 * Usage:
 *   --bam FILE          Input BAM file (required)
 *   --output FILE       Output TSV file (default: stdout)
 *   --max-len INT       Maximum fragment length (default: 1000)
 *   --min-len INT       Minimum fragment length (default: 0)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include "htslib/sam.h"
#include "htslib/hts.h"

#define MAX_FRAG_LEN 2000

// Fragment length histogram
static uint64_t frag_hist[MAX_FRAG_LEN + 1] = {0};
static uint64_t total_fragments = 0;

// Compute fragment length from paired-end BAM records
// Matches Salmon's fragLengthPedantic logic
int32_t compute_fragment_length(bam1_t* read1, bam1_t* read2, int32_t ref_len) {
    // Check both reads are mapped
    if (read1->core.flag & BAM_FUNMAP || read2->core.flag & BAM_FUNMAP) {
        return 0;
    }
    
    // Check both reads map to same reference
    if (read1->core.tid != read2->core.tid) {
        return 0;
    }
    
    // Get strand orientations
    bool fwd1 = !(read1->core.flag & BAM_FREVERSE);
    bool fwd2 = !(read2->core.flag & BAM_FREVERSE);
    
    // Must be opposite strands (inward orientation)
    if (fwd1 == fwd2) {
        return 0;
    }
    
    // Get positions (0-based)
    int32_t pos1 = read1->core.pos;
    int32_t pos2 = read2->core.pos;
    int32_t len1 = bam_cigar2rlen(read1->core.n_cigar, bam_get_cigar(read1));
    int32_t len2 = bam_cigar2rlen(read2->core.n_cigar, bam_get_cigar(read2));
    
    // Clamp positions to valid range
    int32_t ref_len_signed = (int32_t)ref_len;
    pos1 = (pos1 < 0) ? 0 : (pos1 > ref_len_signed ? ref_len_signed : pos1);
    pos2 = (pos2 < 0) ? 0 : (pos2 > ref_len_signed ? ref_len_signed : pos2);
    
    // Compute fragment endpoints
    // Forward read's leftmost 5' end
    int32_t p1 = fwd1 ? pos1 : pos2;
    // Reverse read's rightmost 3' end
    int32_t p2 = fwd1 ? (pos2 + len2) : (pos1 + len1);
    
    // Clamp again
    p1 = (p1 < 0) ? 0 : (p1 > ref_len_signed ? ref_len_signed : p1);
    p2 = (p2 < 0) ? 0 : (p2 > ref_len_signed ? ref_len_signed : p2);
    
    // Fragment length = distance between endpoints
    int32_t frag_len = (p1 > p2) ? (p1 - p2) : (p2 - p1);
    
    return frag_len;
}

void print_usage(const char* progname) {
    fprintf(stderr, "Usage: %s [options]\n\n", progname);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  --bam FILE          Input BAM file\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --output FILE       Output TSV file (default: stdout)\n");
    fprintf(stderr, "  --max-len INT       Maximum fragment length (default: 1000)\n");
    fprintf(stderr, "  --min-len INT       Minimum fragment length (default: 0)\n");
    fprintf(stderr, "  -h, --help          Show this help\n");
}

int main(int argc, char** argv) {
    const char* bam_file = NULL;
    const char* output_file = NULL;
    int32_t max_len = 1000;
    int32_t min_len = 0;
    
    // Parse command line
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--bam") == 0 && i + 1 < argc) {
            bam_file = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--max-len") == 0 && i + 1 < argc) {
            max_len = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--min-len") == 0 && i + 1 < argc) {
            min_len = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    if (!bam_file) {
        fprintf(stderr, "Error: --bam is required\n");
        print_usage(argv[0]);
        return 1;
    }
    
    if (max_len > MAX_FRAG_LEN) {
        fprintf(stderr, "Error: --max-len cannot exceed %d\n", MAX_FRAG_LEN);
        return 1;
    }
    
    // Open BAM file
    samFile* fp = sam_open(bam_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open BAM file %s\n", bam_file);
        return 1;
    }
    
    bam_hdr_t* header = sam_hdr_read(fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header\n");
        sam_close(fp);
        return 1;
    }
    
    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();
    if (!read1 || !read2) {
        fprintf(stderr, "Error: Cannot allocate BAM records\n");
        bam_hdr_destroy(header);
        sam_close(fp);
        return 1;
    }
    
    // Read BAM records
    // Use TLEN field which is more reliable than trying to match mates
    while (sam_read1(fp, header, read1) >= 0) {
        // Skip unmapped reads
        if (read1->core.flag & BAM_FUNMAP) {
            continue;
        }
        
        // Check if properly paired
        if (!(read1->core.flag & BAM_FPROPER_PAIR)) {
            continue;
        }
        
        // Check if mate is unmapped
        if (read1->core.flag & BAM_FMUNMAP) {
            continue;
        }
        
        // Use TLEN (template length) field - BAM already computes this correctly
        // TLEN is signed: positive for forward orientation, negative for reverse
        int32_t tlen = read1->core.isize;
        int32_t frag_len = (tlen < 0) ? -tlen : tlen;
        
        // Only process if this is the "leftmost" read (TLEN is set on first read)
        // For paired-end, TLEN is typically set on read1
        if (frag_len > 0 && frag_len >= min_len && frag_len <= max_len && frag_len <= MAX_FRAG_LEN) {
            // Only count once per pair (use read1 flag to avoid double counting)
            if (read1->core.flag & BAM_FREAD1) {
                frag_hist[frag_len]++;
                total_fragments++;
            }
        }
    }
    bam_destroy1(read1);
    bam_hdr_destroy(header);
    sam_close(fp);
    
    // Output histogram
    FILE* out = output_file ? fopen(output_file, "w") : stdout;
    if (!out) {
        fprintf(stderr, "Error: Cannot open output file %s\n", output_file);
        return 1;
    }
    
    // Normalize to probabilities
    if (total_fragments > 0) {
        for (int i = 0; i <= max_len; i++) {
            double prob = (double)frag_hist[i] / total_fragments;
            fprintf(out, "%d\t%lu\t%.10f\n", i, frag_hist[i], prob);
        }
    } else {
        fprintf(stderr, "Warning: No properly paired fragments found\n");
        // Output zeros
        for (int i = 0; i <= max_len; i++) {
            fprintf(out, "%d\t0\t0.0\n", i);
        }
    }
    
    if (out != stdout) {
        fclose(out);
    }
    
    fprintf(stderr, "Processed %lu properly paired fragments\n", total_fragments);
    
    return 0;
}
