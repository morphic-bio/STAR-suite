/*
 * remove_y_reads - Split FASTQ files based on Y-only BAM
 *
 * Reads a Y-only BAM file, builds a hash set of read names, and splits
 * FASTQ files into Y/noY partitions while preserving read order.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <pthread.h>
#include <semaphore.h>
#include <zlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"

// Initialize kseq for gzFile
KSEQ_INIT(gzFile, gzread)

// Hash map: key=64-bit hash, value=struct with length + second hash (for collision detection)
typedef struct {
    uint32_t len;
    uint64_t hash2;  // Second hash for collision detection
} qname_info_t;

KHASH_MAP_INIT_INT64(qname_set, qname_info_t)

// FNV-1a 64-bit hash (primary)
static inline uint64_t hash_qname(const char *s, int len) {
    uint64_t h = 14695981039346656037ULL;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)(unsigned char)s[i];
        h *= 1099511628211ULL;
    }
    return h;
}

// djb2 64-bit hash (secondary, different algorithm for collision detection)
static inline uint64_t hash_qname2(const char *s, int len) {
    uint64_t h = 5381ULL;
    for (int i = 0; i < len; i++) {
        h = ((h << 5) + h) + (uint64_t)(unsigned char)s[i];
    }
    return h;
}

// Normalize qname:
// 1. Strip leading '@' (FASTQ format)
// 2. Strip trailing '/1' or '/2' (mate suffix)
// 3. Truncate at first whitespace or comment
// Returns normalized length
static int normalize_qname(const char *raw, int raw_len, char *out, int max_len) {
    const char *start = raw;
    int len = raw_len;
    
    // Skip leading '@' if present (FASTQ name lines)
    if (len > 0 && start[0] == '@') {
        start++;
        len--;
    }
    
    // Find end: stop at whitespace, tab, or comment
    int end = 0;
    while (end < len && start[end] != ' ' && start[end] != '\t' && 
           start[end] != '\n' && start[end] != '\r') {
        end++;
    }
    len = end;
    
    // Strip trailing /1 or /2 (mate suffix)
    if (len >= 2 && start[len-2] == '/' && 
        (start[len-1] == '1' || start[len-1] == '2')) {
        len -= 2;
    }
    
    // Copy to output
    if (len > max_len - 1) len = max_len - 1;
    memcpy(out, start, len);
    out[len] = '\0';
    return len;
}

// Load Y BAM into hash set
khash_t(qname_set) *load_y_bam(const char *bam_path, long *count) {
    khash_t(qname_set) *h = kh_init(qname_set);
    
    samFile *fp = sam_open(bam_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", bam_path);
        exit(1);
    }
    
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Error: Cannot read header from %s\n", bam_path);
        sam_close(fp);
        exit(1);
    }
    
    bam1_t *b = bam_init1();
    char norm[512];
    *count = 0;
    
    while (sam_read1(fp, hdr, b) >= 0) {
        const char *qname = bam_get_qname(b);
        int qlen = strlen(qname);
        int norm_len = normalize_qname(qname, qlen, norm, sizeof(norm));
        
        uint64_t hash = hash_qname(norm, norm_len);
        uint64_t hash2 = hash_qname2(norm, norm_len);
        
        int absent;
        khint_t k = kh_put(qname_set, h, hash, &absent);
        if (absent) {
            qname_info_t info;
            info.len = (uint32_t)norm_len;
            info.hash2 = hash2;
            kh_val(h, k) = info;  // Store length + second hash for collision check
            (*count)++;
        }
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    
    fprintf(stderr, "Loaded %ld unique Y read names\n", *count);
    return h;
}

// Parse FASTQ filename into stem and extension
// Input: "sample_R1.fastq.gz" -> stem="sample_R1", ext=".fastq.gz"
// Input: "sample_R1.fq" -> stem="sample_R1", ext=".fq"
static void parse_fastq_name(const char *filename, char *stem, char *ext) {
    const char *base = strrchr(filename, '/');
    base = base ? base + 1 : filename;
    
    // Check for compound extensions first
    const char *gz = strstr(base, ".fastq.gz");
    if (gz) {
        int stem_len = gz - base;
        strncpy(stem, base, stem_len);
        stem[stem_len] = '\0';
        strcpy(ext, ".fastq.gz");
        return;
    }
    
    gz = strstr(base, ".fq.gz");
    if (gz) {
        int stem_len = gz - base;
        strncpy(stem, base, stem_len);
        stem[stem_len] = '\0';
        strcpy(ext, ".fq.gz");
        return;
    }
    
    // Single extensions
    const char *dot = strrchr(base, '.');
    if (dot && (strcmp(dot, ".fastq") == 0 || strcmp(dot, ".fq") == 0 ||
                strcmp(dot, ".gz") == 0)) {
        int stem_len = dot - base;
        strncpy(stem, base, stem_len);
        stem[stem_len] = '\0';
        strcpy(ext, dot);
    } else {
        strcpy(stem, base);
        strcpy(ext, "");
    }
}

// Derive output paths
// Input: /path/to/sample_R1.fastq.gz, outdir=/out
// Y output: /out/sample_R1_Y.fastq.gz
// noY output: /out/sample_R1_noY.fastq.gz
static void derive_output_paths(const char *inpath, const char *outdir,
                                char *y_path, char *noy_path, int max_len) {
    char stem[256], ext[32];
    parse_fastq_name(inpath, stem, ext);
    
    const char *dir = outdir;
    char indir[256] = ".";
    if (!dir) {
        // Use input file's directory
        const char *last_slash = strrchr(inpath, '/');
        if (last_slash) {
            int dlen = last_slash - inpath;
            strncpy(indir, inpath, dlen);
            indir[dlen] = '\0';
        }
        dir = indir;
    }
    
    snprintf(y_path, max_len, "%s/%s_Y%s", dir, stem, ext);
    snprintf(noy_path, max_len, "%s/%s_noY%s", dir, stem, ext);
}

// Process single FASTQ file
static void process_fastq(const char *inpath, gzFile out_y, gzFile out_noy,
                          khash_t(qname_set) *y_set, long *total, long *y_cnt, long *noy_cnt) {
    gzFile fp = gzopen(inpath, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", inpath);
        return;
    }
    
    gzbuffer(fp, 128 * 1024);  // MANDATORY: large buffer
    
    kseq_t *seq = kseq_init(fp);
    char norm[512];
    
    while (kseq_read(seq) >= 0) {
        // seq->name.s is dynamically allocated, handles any length
        int norm_len = normalize_qname(seq->name.s, seq->name.l, norm, sizeof(norm));
        
        uint64_t h = hash_qname(norm, norm_len);
        uint64_t h2 = hash_qname2(norm, norm_len);
        khint_t k = kh_get(qname_set, y_set, h);
        
        // Collision check: verify length AND second hash match if found
        bool is_y = false;
        if (k != kh_end(y_set)) {
            qname_info_t *info = &kh_val(y_set, k);
            if (info->len == (uint32_t)norm_len && info->hash2 == h2) {
                is_y = true;
            }
        }
        
        gzFile out = is_y ? out_y : out_noy;
        
        // Write FASTQ record
        gzprintf(out, "@%s", seq->name.s);
        if (seq->comment.l) gzprintf(out, " %s", seq->comment.s);
        gzprintf(out, "\n%s\n+\n%s\n", seq->seq.s, seq->qual.s);
        
        (*total)++;
        if (is_y) (*y_cnt)++; else (*noy_cnt)++;
    }
    
    kseq_destroy(seq);
    gzclose(fp);
}

// Worker thread structure
typedef struct {
    char *inpath;
    char y_path[512];
    char noy_path[512];
    khash_t(qname_set) *y_set;
    int gzip_level;
    long total, y_count, noy_count;
    sem_t *thread_sem;  // Shared semaphore
} fastq_worker_t;

// Worker thread function
void *worker_thread(void *arg) {
    fastq_worker_t *w = (fastq_worker_t *)arg;
    
    // Acquire semaphore (blocks if at thread limit)
    sem_wait(w->thread_sem);
    
    // Open outputs with gzbuffer
    char mode[8];
    snprintf(mode, sizeof(mode), "wb%d", w->gzip_level);
    gzFile out_y = gzopen(w->y_path, mode);
    gzFile out_noy = gzopen(w->noy_path, mode);
    
    if (!out_y || !out_noy) {
        fprintf(stderr, "Error: Cannot create output files for %s\n", w->inpath);
        sem_post(w->thread_sem);
        return NULL;
    }
    
    gzbuffer(out_y, 128 * 1024);
    gzbuffer(out_noy, 128 * 1024);
    
    // Process FASTQ
    process_fastq(w->inpath, out_y, out_noy, w->y_set,
                  &w->total, &w->y_count, &w->noy_count);
    
    gzclose(out_y);
    gzclose(out_noy);
    
    // Release semaphore
    sem_post(w->thread_sem);
    
    return NULL;
}

// Split all FASTQ files
static void split_all_fastqs(int num_files, char **paths, const char *outdir,
                             khash_t(qname_set) *y_set, int num_threads, int gzip_level) {
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    fastq_worker_t *workers = malloc(num_files * sizeof(fastq_worker_t));
    int *file_indices = malloc(num_threads * sizeof(int));  // Map thread slot to file index
    
    if (!threads || !workers || !file_indices) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(1);
    }
    
    // Initialize semaphore to bound concurrency
    sem_t thread_sem;
    sem_init(&thread_sem, 0, num_threads);
    
    // Initialize all workers
    for (int i = 0; i < num_files; i++) {
        workers[i].inpath = paths[i];
        derive_output_paths(paths[i], outdir, 
                           workers[i].y_path, workers[i].noy_path, sizeof(workers[i].y_path));
        workers[i].y_set = y_set;
        workers[i].gzip_level = gzip_level;
        workers[i].thread_sem = &thread_sem;
        workers[i].total = workers[i].y_count = workers[i].noy_count = 0;
    }
    
    // Process files in batches to avoid creating too many parked threads
    int next_file = 0;
    int active_count = 0;
    
    while (next_file < num_files || active_count > 0) {
        // Start new threads up to the limit
        while (active_count < num_threads && next_file < num_files) {
            int slot = active_count;
            file_indices[slot] = next_file;
            
            if (pthread_create(&threads[slot], NULL, worker_thread, &workers[next_file]) != 0) {
                fprintf(stderr, "Error: Failed to create thread for %s\n", paths[next_file]);
                exit(1);
            }
            active_count++;
            next_file++;
        }
        
        // Wait for at least one thread to complete
        if (active_count > 0) {
            // Join the first active thread (blocks until it completes)
            pthread_join(threads[0], NULL);
            
            // Report results for the completed file
            int completed_file = file_indices[0];
            fprintf(stderr, "%s: total=%ld Y=%ld noY=%ld\n",
                    workers[completed_file].inpath, workers[completed_file].total,
                    workers[completed_file].y_count, workers[completed_file].noy_count);
            
            // Shift remaining threads and indices
            active_count--;
            for (int i = 0; i < active_count; i++) {
                threads[i] = threads[i + 1];
                file_indices[i] = file_indices[i + 1];
            }
        }
    }
    
    sem_destroy(&thread_sem);
    free(threads);
    free(workers);
    free(file_indices);
}

static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s -y <Y.bam> [-o <out_dir>] [--threads N] "
                    "[--gzip-level N] fastq1 [fastq2 ...]\n", prog);
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  -y, --ybam FILE     Y-only BAM file (required)\n");
    fprintf(stderr, "  -o, --outdir DIR    Output directory (default: alongside input)\n");
    fprintf(stderr, "  -t, --threads N     Number of threads (default: 1)\n");
    fprintf(stderr, "  -z, --gzip-level N  Gzip compression level 1-9 (default: 6)\n");
    fprintf(stderr, "  -h, --help          Show this help\n");
}

int main(int argc, char **argv) {
    char *ybam_path = NULL;
    char *outdir = NULL;
    int num_threads = 1;
    int gzip_level = 6;
    
    static struct option long_options[] = {
        {"ybam",       required_argument, 0, 'y'},
        {"outdir",     required_argument, 0, 'o'},
        {"threads",    required_argument, 0, 't'},
        {"gzip-level", required_argument, 0, 'z'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, "y:o:t:z:h", long_options, NULL)) != -1) {
        switch (c) {
            case 'y': ybam_path = optarg; break;
            case 'o': outdir = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            case 'z': gzip_level = atoi(optarg); break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }
    
    if (!ybam_path || optind >= argc) {
        fprintf(stderr, "Error: -y <Y.bam> and FASTQ files required\n");
        print_usage(argv[0]);
        return 1;
    }
    
    if (num_threads < 1) {
        fprintf(stderr, "Error: --threads must be >= 1\n");
        return 1;
    }
    
    if (gzip_level < 1 || gzip_level > 9) {
        fprintf(stderr, "Error: --gzip-level must be between 1 and 9\n");
        return 1;
    }
    
    // Load Y BAM
    long y_count;
    khash_t(qname_set) *y_set = load_y_bam(ybam_path, &y_count);
    
    if (y_count == 0) {
        fprintf(stderr, "Warning: No Y read names found in BAM\n");
    }
    
    // Split FASTQs
    int num_files = argc - optind;
    char **fastq_paths = &argv[optind];
    split_all_fastqs(num_files, fastq_paths, outdir, y_set, num_threads, gzip_level);
    
    kh_destroy(qname_set, y_set);
    return 0;
}

