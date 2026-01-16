/*
 * compute_expected_gc - Compute expected GC distribution from transcriptome
 * 
 * Usage:
 *   --transcriptome FILE    Transcriptome FASTA file
 *   --genome FILE           Genome FASTA file (requires --gtf)
 *   --gtf FILE              GTF annotation file (requires --genome)
 *   --fld FILE              Fragment length distribution (default: built-in)
 *   --output FILE           Output file (default: stdout)
 *   --min-frag INT          Minimum fragment length (default: 50)
 *   --max-frag INT          Maximum fragment length (default: 1000)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAX_LINE_LEN 1000000
#define MAX_SEQ_LEN 10000000
#define MAX_FRAG_LEN 2000
#define GC_BINS 101
#define MAX_TRANSCRIPTS 100000
#define MAX_EXONS_PER_TX 1000

// Fragment length distribution (default: normal-like, mean=200, sd=80)
static double default_fld[MAX_FRAG_LEN] = {0};

// Expected GC distribution
static double expected_gc[GC_BINS] = {0};

// Exon structure for GTF parsing
typedef struct {
    char* transcript_id;
    char* chrom;
    int start;  // 1-based, inclusive
    int end;    // 1-based, inclusive
    char strand;
} Exon;

// Transcript structure
typedef struct {
    char* id;
    char* seq;
    int len;
    Exon* exons;
    int num_exons;
} Transcript;

// Global arrays
static Transcript transcripts[MAX_TRANSCRIPTS];
static int num_transcripts = 0;
static char* genome_seqs[1000];
static int genome_lens[1000];
static char* genome_names[1000];
static int num_chromosomes = 0;

// Function prototypes
void init_default_fld(void);
int parse_fld_file(const char* filename, double* fld);
int parse_transcriptome_fasta(const char* filename);
int parse_genome_fasta(const char* filename);
int parse_gtf(const char* filename);
void extract_transcript_sequences(void);
void compute_expected_gc_distribution(double* fld, int min_frag, int max_frag);
void compute_quantile_bounds(double* fld, double quantile_low, double quantile_high, int* fld_low, int* fld_high);
void print_usage(const char* progname);
char* extract_transcript_id(const char* attributes);
int find_chromosome(const char* chrom_name);

// Initialize default fragment length distribution (normal-like)
void init_default_fld(void) {
    double mean = 200.0;
    double sd = 80.0;
    double sum = 0.0;
    
    for (int i = 0; i < MAX_FRAG_LEN; i++) {
        if (i < 50) {
            default_fld[i] = 0.0;
        } else {
            double x = (i - mean) / sd;
            default_fld[i] = exp(-0.5 * x * x);
            sum += default_fld[i];
        }
    }
    
    // Normalize
    for (int i = 0; i < MAX_FRAG_LEN; i++) {
        default_fld[i] /= sum;
    }
}

// Parse fragment length distribution file
// Supports two formats:
// 1. One value per line (legacy format)
// 2. TSV format: length <tab> count <tab> probability (from sample_fld)
int parse_fld_file(const char* filename, double* fld) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open FLD file %s\n", filename);
        return 0;
    }
    
    // Initialize to zero
    for (int i = 0; i < MAX_FRAG_LEN; i++) {
        fld[i] = 0.0;
    }
    
    char line[1024];
    int line_num = 0;
    int has_tab = 0;
    
    // Check first line to detect format
    if (fgets(line, sizeof(line), f)) {
        line_num++;
        if (strchr(line, '\t') != NULL) {
            has_tab = 1;  // TSV format
        }
        rewind(f);
    }
    
    double sum = 0.0;
    
    if (has_tab) {
        // TSV format: length <tab> count <tab> probability
        int len;
        unsigned long count;
        double prob;
        while (fgets(line, sizeof(line), f)) {
            if (sscanf(line, "%d %lu %lf", &len, &count, &prob) >= 2) {
                if (len >= 0 && len < MAX_FRAG_LEN) {
                    // Use probability if available, otherwise use count
                    if (prob > 0) {
                        fld[len] = prob;
                    } else if (count > 0) {
                        fld[len] = (double)count;
                    }
                    sum += fld[len];
                }
            }
        }
    } else {
        // Legacy format: one value per line
        for (int i = 0; i < MAX_FRAG_LEN; i++) {
            if (fscanf(f, "%lf", &fld[i]) != 1) {
                fld[i] = 0.0;
            }
            sum += fld[i];
        }
    }
    
    fclose(f);
    
    // Normalize if needed
    if (sum > 0) {
        for (int i = 0; i < MAX_FRAG_LEN; i++) {
            fld[i] /= sum;
        }
    }
    return 1;
}

// Parse transcriptome FASTA
int parse_transcriptome_fasta(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open transcriptome file %s\n", filename);
        return 0;
    }
    
    char line[MAX_LINE_LEN];
    char* current_seq = NULL;
    int current_len = 0;
    int current_alloc = 0;
    char* current_id = NULL;
    
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            // Save previous transcript
            if (current_seq && current_len > 0 && current_id) {
                if (num_transcripts >= MAX_TRANSCRIPTS) {
                    fprintf(stderr, "Error: Too many transcripts (max %d)\n", MAX_TRANSCRIPTS);
                    break;
                }
                transcripts[num_transcripts].id = strdup(current_id);
                transcripts[num_transcripts].seq = current_seq;
                transcripts[num_transcripts].len = current_len;
                transcripts[num_transcripts].exons = NULL;
                transcripts[num_transcripts].num_exons = 0;
                num_transcripts++;
            }
            
            // Start new transcript
            current_len = 0;
            current_alloc = 10000;
            current_seq = malloc(current_alloc);
            if (!current_seq) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                fclose(f);
                return 0;
            }
            
            // Extract transcript ID (first word after >)
            int id_start = 1;
            while (line[id_start] == ' ' || line[id_start] == '\t') id_start++;
            int id_end = id_start;
            while (line[id_end] != '\0' && line[id_end] != ' ' && 
                   line[id_end] != '\t' && line[id_end] != '\n') id_end++;
            int id_len = id_end - id_start;
            current_id = malloc(id_len + 1);
            strncpy(current_id, line + id_start, id_len);
            current_id[id_len] = '\0';
        } else {
            // Append sequence
            int line_len = strlen(line);
            if (line[line_len - 1] == '\n') line_len--;
            
            if (current_len + line_len + 1 > current_alloc) {
                current_alloc = (current_len + line_len) * 2;
                current_seq = realloc(current_seq, current_alloc);
                if (!current_seq) {
                    fprintf(stderr, "Error: Memory reallocation failed\n");
                    fclose(f);
                    return 0;
                }
            }
            
            memcpy(current_seq + current_len, line, line_len);
            current_len += line_len;
        }
    }
    
    // Save last transcript
    if (current_seq && current_len > 0 && current_id) {
        if (num_transcripts < MAX_TRANSCRIPTS) {
            transcripts[num_transcripts].id = strdup(current_id);
            transcripts[num_transcripts].seq = current_seq;
            transcripts[num_transcripts].len = current_len;
            transcripts[num_transcripts].exons = NULL;
            transcripts[num_transcripts].num_exons = 0;
            num_transcripts++;
        }
    } else {
        free(current_seq);
        free(current_id);
    }
    
    fclose(f);
    return 1;
}

// Parse genome FASTA
int parse_genome_fasta(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open genome file %s\n", filename);
        return 0;
    }
    
    char line[MAX_LINE_LEN];
    char* current_seq = NULL;
    int current_len = 0;
    int current_alloc = 0;
    char* current_name = NULL;
    
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            // Save previous chromosome
            if (current_seq && current_len > 0 && current_name) {
                if (num_chromosomes >= 1000) {
                    fprintf(stderr, "Error: Too many chromosomes (max 1000)\n");
                    break;
                }
                genome_seqs[num_chromosomes] = current_seq;
                genome_lens[num_chromosomes] = current_len;
                genome_names[num_chromosomes] = current_name;
                num_chromosomes++;
            }
            
            // Start new chromosome
            current_len = 0;
            current_alloc = 1000000;
            current_seq = malloc(current_alloc);
            if (!current_seq) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                fclose(f);
                return 0;
            }
            
            // Extract chromosome name
            int name_start = 1;
            while (line[name_start] == ' ' || line[name_start] == '\t') name_start++;
            int name_end = name_start;
            while (line[name_end] != '\0' && line[name_end] != ' ' && 
                   line[name_end] != '\t' && line[name_end] != '\n') name_end++;
            int name_len = name_end - name_start;
            current_name = malloc(name_len + 1);
            strncpy(current_name, line + name_start, name_len);
            current_name[name_len] = '\0';
        } else {
            // Append sequence
            int line_len = strlen(line);
            if (line[line_len - 1] == '\n') line_len--;
            
            if (current_len + line_len + 1 > current_alloc) {
                current_alloc = (current_len + line_len) * 2;
                current_seq = realloc(current_seq, current_alloc);
                if (!current_seq) {
                    fprintf(stderr, "Error: Memory reallocation failed\n");
                    fclose(f);
                    return 0;
                }
            }
            
            memcpy(current_seq + current_len, line, line_len);
            current_len += line_len;
        }
    }
    
    // Save last chromosome
    if (current_seq && current_len > 0 && current_name) {
        if (num_chromosomes < 1000) {
            genome_seqs[num_chromosomes] = current_seq;
            genome_lens[num_chromosomes] = current_len;
            genome_names[num_chromosomes] = current_name;
            num_chromosomes++;
        }
    } else {
        free(current_seq);
        free(current_name);
    }
    
    fclose(f);
    return 1;
}

// Extract transcript_id from GTF attributes column
char* extract_transcript_id(const char* attributes) {
    const char* key = "transcript_id";
    const char* p = strstr(attributes, key);
    if (!p) return NULL;
    
    p += strlen(key);
    while (*p == ' ' || *p == '\t') p++;
    if (*p != '"') return NULL;
    p++;
    
    int len = 0;
    while (p[len] != '"' && p[len] != '\0') len++;
    if (len == 0) return NULL;
    
    char* id = malloc(len + 1);
    strncpy(id, p, len);
    id[len] = '\0';
    return id;
}

// Find chromosome index by name
int find_chromosome(const char* chrom_name) {
    for (int i = 0; i < num_chromosomes; i++) {
        if (strcmp(genome_names[i], chrom_name) == 0) {
            return i;
        }
    }
    return -1;
}

// Parse GTF file
int parse_gtf(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open GTF file %s\n", filename);
        return 0;
    }
    
    char line[MAX_LINE_LEN];
    int line_num = 0;
    
    // First pass: count exons per transcript
    int* exon_counts = calloc(MAX_TRANSCRIPTS, sizeof(int));
    char** transcript_ids = calloc(MAX_TRANSCRIPTS, sizeof(char*));
    int num_unique_tx = 0;
    
    while (fgets(line, sizeof(line), f)) {
        line_num++;
        if (line[0] == '#') continue;
        
        char chrom[100], feature[100], strand;
        int start, end;
        char attributes[10000];
        
        if (sscanf(line, "%99s %*s %99s %d %d %*s %c %*s %9999s", 
                   chrom, feature, &start, &end, &strand, attributes) != 6) {
            continue;
        }
        
        if (strcmp(feature, "exon") != 0) continue;
        
        char* tx_id = extract_transcript_id(attributes);
        if (!tx_id) continue;
        
        // Find or create transcript index
        int tx_idx = -1;
        for (int i = 0; i < num_unique_tx; i++) {
            if (strcmp(transcript_ids[i], tx_id) == 0) {
                tx_idx = i;
                break;
            }
        }
        
        if (tx_idx == -1) {
            if (num_unique_tx >= MAX_TRANSCRIPTS) {
                free(tx_id);
                continue;
            }
            tx_idx = num_unique_tx;
            transcript_ids[tx_idx] = strdup(tx_id);
            num_unique_tx++;
        }
        
        exon_counts[tx_idx]++;
        free(tx_id);
    }
    
    // Allocate exon arrays
    for (int i = 0; i < num_unique_tx; i++) {
        transcripts[i].id = transcript_ids[i];
        transcripts[i].num_exons = exon_counts[i];
        transcripts[i].exons = malloc(exon_counts[i] * sizeof(Exon));
        transcripts[i].len = 0;
        transcripts[i].seq = NULL;
    }
    
    // Second pass: fill exon data
    rewind(f);
    int* exon_idx = calloc(num_unique_tx, sizeof(int));
    
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        
        char chrom[100], feature[100], strand;
        int start, end;
        char attributes[10000];
        
        if (sscanf(line, "%99s %*s %99s %d %d %*s %c %*s %9999s", 
                   chrom, feature, &start, &end, &strand, attributes) != 6) {
            continue;
        }
        
        if (strcmp(feature, "exon") != 0) continue;
        
        char* tx_id = extract_transcript_id(attributes);
        if (!tx_id) continue;
        
        // Find transcript index
        int tx_idx = -1;
        for (int i = 0; i < num_unique_tx; i++) {
            if (strcmp(transcript_ids[i], tx_id) == 0) {
                tx_idx = i;
                break;
            }
        }
        
        if (tx_idx >= 0 && exon_idx[tx_idx] < exon_counts[tx_idx]) {
            Exon* ex = &transcripts[tx_idx].exons[exon_idx[tx_idx]];
            ex->chrom = strdup(chrom);
            ex->start = start;
            ex->end = end;
            ex->strand = strand;
            ex->transcript_id = tx_id;
            exon_idx[tx_idx]++;
        } else {
            free(tx_id);
        }
    }
    
    num_transcripts = num_unique_tx;
    free(exon_counts);
    free(exon_idx);
    fclose(f);
    return 1;
}

// Extract transcript sequences from genome using exon coordinates
void extract_transcript_sequences(void) {
    for (int tx_idx = 0; tx_idx < num_transcripts; tx_idx++) {
        Transcript* tx = &transcripts[tx_idx];
        if (tx->num_exons == 0) continue;
        
        // Calculate total length
        int total_len = 0;
        for (int i = 0; i < tx->num_exons; i++) {
            total_len += tx->exons[i].end - tx->exons[i].start + 1;
        }
        
        tx->seq = malloc(total_len + 1);
        tx->len = total_len;
        int seq_pos = 0;
        
        // Extract exons in order (handle strand)
        if (tx->num_exons > 0) {
            // Sort exons by start position
            for (int i = 0; i < tx->num_exons - 1; i++) {
                for (int j = i + 1; j < tx->num_exons; j++) {
                    if (tx->exons[j].start < tx->exons[i].start) {
                        Exon tmp = tx->exons[i];
                        tx->exons[i] = tx->exons[j];
                        tx->exons[j] = tmp;
                    }
                }
            }
            
            // Extract sequence
            for (int i = 0; i < tx->num_exons; i++) {
                Exon* ex = &tx->exons[i];
                int chrom_idx = find_chromosome(ex->chrom);
                if (chrom_idx < 0) {
                    fprintf(stderr, "Warning: Chromosome %s not found in genome\n", ex->chrom);
                    continue;
                }
                
                // GTF coordinates are 1-based, inclusive
                int start = ex->start - 1;  // Convert to 0-based
                int end = ex->end;          // Still 1-based, but end is inclusive
                
                if (start < 0 || end > genome_lens[chrom_idx]) {
                    fprintf(stderr, "Warning: Exon coordinates out of range for %s\n", ex->chrom);
                    continue;
                }
                
                // Extract sequence (handle reverse strand)
                if (ex->strand == '+') {
                    for (int pos = start; pos < end; pos++) {
                        char c = toupper(genome_seqs[chrom_idx][pos]);
                        tx->seq[seq_pos++] = c;
                    }
                } else {
                    // Reverse complement
                    for (int pos = end - 1; pos >= start; pos--) {
                        char c = toupper(genome_seqs[chrom_idx][pos]);
                        // Reverse complement
                        if (c == 'A') c = 'T';
                        else if (c == 'T') c = 'A';
                        else if (c == 'G') c = 'C';
                        else if (c == 'C') c = 'G';
                        else if (c == 'N') c = 'N';
                        tx->seq[seq_pos++] = c;
                    }
                }
            }
        }
        
        tx->seq[seq_pos] = '\0';
    }
}

// Compute quantile bounds from FLD (matching Salmon's 0.5% and 99.5% quantiles)
void compute_quantile_bounds(double* fld, double quantile_low, double quantile_high, int* fld_low, int* fld_high) {
    // Build CDF
    double cdf[MAX_FRAG_LEN];
    cdf[0] = fld[0];
    for (int i = 1; i < MAX_FRAG_LEN; i++) {
        cdf[i] = cdf[i-1] + fld[i];
    }
    
    // Normalize CDF to [0, 1]
    double cdf_max = cdf[MAX_FRAG_LEN - 1];
    if (cdf_max <= 0) {
        *fld_low = 0;
        *fld_high = MAX_FRAG_LEN - 1;
        return;
    }
    
    for (int i = 0; i < MAX_FRAG_LEN; i++) {
        cdf[i] /= cdf_max;
    }
    
    // Find quantile bounds
    *fld_low = 0;
    *fld_high = MAX_FRAG_LEN - 1;
    
    for (int i = 0; i < MAX_FRAG_LEN; i++) {
        if (cdf[i] >= quantile_low && *fld_low == 0) {
            *fld_low = i;
        }
        if (cdf[i] >= quantile_high) {
            *fld_high = i;
            break;
        }
    }
}

// Compute context GC bin (matching Salmon's contextBin logic)
// Context is GC% of 5bp upstream + 5bp downstream of fragment
int compute_context_gc_bin(int context_gc_pct, int num_context_bins) {
    if (num_context_bins <= 1) {
        return 0;  // 1D model - no context
    }
    // Map context GC% (0-100) to bin (0 to num_context_bins-1)
    double bin_width = 100.0 / num_context_bins;
    int bin = (int)(context_gc_pct / bin_width);
    if (bin >= num_context_bins) bin = num_context_bins - 1;
    return bin;
}

// Compute expected GC distribution using CDF-based weighting with quantile bounds (matching Salmon)
void compute_expected_gc_distribution_with_quantiles(double* fld, int min_frag, int max_frag, double quantile_low, double quantile_high, int context_bins) {
    // Initialize expected_gc array
    int total_bins = context_bins * GC_BINS;
    for (int i = 0; i < total_bins && i < GC_BINS * 10; i++) {
        expected_gc[i] = 0.0;
    }
    // Compute quantile bounds (matching Salmon's 0.5% and 99.5%)
    int fld_low, fld_high;
    compute_quantile_bounds(fld, quantile_low, quantile_high, &fld_low, &fld_high);
    
    // Use quantile bounds to limit fragment length range
    if (min_frag < fld_low) min_frag = fld_low;
    if (max_frag > fld_high) max_frag = fld_high;
    
    // Build CDF from FLD (matching Salmon's approach)
    double cdf[MAX_FRAG_LEN];
    cdf[0] = fld[0];
    for (int i = 1; i < MAX_FRAG_LEN; i++) {
        cdf[i] = cdf[i-1] + fld[i];
    }
    
    // Normalize CDF to [0, 1]
    double cdf_max = cdf[MAX_FRAG_LEN - 1];
    if (cdf_max > 0) {
        for (int i = 0; i < MAX_FRAG_LEN; i++) {
            cdf[i] /= cdf_max;
        }
    }
    
    for (int tx_idx = 0; tx_idx < num_transcripts; tx_idx++) {
        Transcript* tx = &transcripts[tx_idx];
        if (!tx->seq || tx->len == 0) continue;
        
        // Build cumulative GC count array
        int* gc_count = malloc((tx->len + 1) * sizeof(int));
        gc_count[0] = 0;
        
        for (int i = 0; i < tx->len; i++) {
            char c = toupper(tx->seq[i]);
            int is_gc = (c == 'G' || c == 'C');
            gc_count[i + 1] = gc_count[i] + is_gc;
        }
        
        // For each start position
        for (int start = 0; start < tx->len; start++) {
            double prev_fl_mass = 0.0;  // Previous CDF value (matching Salmon)
            
            // For each fragment length (using CDF difference weighting)
            for (int fl = min_frag; fl <= max_frag && fl < MAX_FRAG_LEN && start + fl <= tx->len; fl++) {
                if (fl >= MAX_FRAG_LEN) break;
                
                // Use CDF difference instead of raw PDF (matching Salmon)
                double fl_mass = cdf[fl];
                double weight = fl_mass - prev_fl_mass;
                
                if (weight <= 0) {
                    prev_fl_mass = fl_mass;
                    continue;
                }
                
                int end = start + fl;
                int gc = gc_count[end] - gc_count[start];
                int gc_pct = (int)round(100.0 * gc / fl);
                
                if (gc_pct >= 0 && gc_pct < GC_BINS) {
                    expected_gc[gc_pct] += weight;
                }
                
                prev_fl_mass = fl_mass;
            }
        }
        
        free(gc_count);
    }
}

// Print usage
void print_usage(const char* progname) {
    fprintf(stderr, "Usage: %s [options]\n\n", progname);
    fprintf(stderr, "Input (one required):\n");
    fprintf(stderr, "  --transcriptome FILE    Transcriptome FASTA file\n");
    fprintf(stderr, "  --genome FILE           Genome FASTA file (requires --gtf)\n");
    fprintf(stderr, "  --gtf FILE              GTF annotation file (requires --genome)\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --fld FILE              Fragment length distribution (default: built-in)\n");
    fprintf(stderr, "  --output FILE           Output file (default: stdout)\n");
    fprintf(stderr, "  --min-frag INT          Minimum fragment length (default: 50)\n");
    fprintf(stderr, "  --max-frag INT          Maximum fragment length (default: 1000)\n");
    fprintf(stderr, "  --fld-quantile-low FLOAT  Lower quantile bound (default: 0.005)\n");
    fprintf(stderr, "  --fld-quantile-high FLOAT Upper quantile bound (default: 0.995)\n");
    fprintf(stderr, "  --context-bins INT      Number of context GC bins (default: 1 = 1D, 3 = 2D like Salmon)\n");
    fprintf(stderr, "  -h, --help              Show this help\n");
}

int main(int argc, char** argv) {
    const char* transcriptome_file = NULL;
    const char* genome_file = NULL;
    const char* gtf_file = NULL;
    const char* fld_file = NULL;
    const char* output_file = NULL;
    int min_frag = 50;
    int max_frag = 1000;
    double quantile_low = 0.005;
    double quantile_high = 0.995;
    int context_bins = 1;  // Default: 1D model
    
    // Parse command line
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--transcriptome") == 0 && i + 1 < argc) {
            transcriptome_file = argv[++i];
        } else if (strcmp(argv[i], "--genome") == 0 && i + 1 < argc) {
            genome_file = argv[++i];
        } else if (strcmp(argv[i], "--gtf") == 0 && i + 1 < argc) {
            gtf_file = argv[++i];
        } else if (strcmp(argv[i], "--fld") == 0 && i + 1 < argc) {
            fld_file = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "--min-frag") == 0 && i + 1 < argc) {
            min_frag = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-frag") == 0 && i + 1 < argc) {
            max_frag = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--fld-quantile-low") == 0 && i + 1 < argc) {
            quantile_low = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fld-quantile-high") == 0 && i + 1 < argc) {
            quantile_high = atof(argv[++i]);
        } else if (strcmp(argv[i], "--context-bins") == 0 && i + 1 < argc) {
            context_bins = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    // Validate input
    if (!transcriptome_file && (!genome_file || !gtf_file)) {
        fprintf(stderr, "Error: Must provide either --transcriptome or both --genome and --gtf\n");
        print_usage(argv[0]);
        return 1;
    }
    
    if (genome_file && !gtf_file) {
        fprintf(stderr, "Error: --genome requires --gtf\n");
        return 1;
    }
    
    // Initialize fragment length distribution
    double fld[MAX_FRAG_LEN] = {0};
    if (fld_file) {
        if (!parse_fld_file(fld_file, fld)) {
            return 1;
        }
    } else {
        init_default_fld();
        memcpy(fld, default_fld, sizeof(fld));
    }
    
    // Load input
    if (transcriptome_file) {
        if (!parse_transcriptome_fasta(transcriptome_file)) {
            return 1;
        }
    } else {
        if (!parse_genome_fasta(genome_file)) {
            return 1;
        }
        if (!parse_gtf(gtf_file)) {
            return 1;
        }
        extract_transcript_sequences();
    }
    
    fprintf(stderr, "Loaded %d transcripts\n", num_transcripts);
    fprintf(stderr, "Using %d context bin(s) for GC distribution\n", context_bins);
    
    // Compute expected GC distribution with quantile bounds
    compute_expected_gc_distribution_with_quantiles(fld, min_frag, max_frag, quantile_low, quantile_high, context_bins);
    
    // Normalize
    int total_bins = context_bins * GC_BINS;
    double total = 0.0;
    for (int i = 0; i < total_bins && i < GC_BINS * 10; i++) {
        total += expected_gc[i];
    }
    
    // Output
    FILE* out = output_file ? fopen(output_file, "w") : stdout;
    if (!out) {
        fprintf(stderr, "Error: Cannot open output file %s\n", output_file);
        return 1;
    }
    
    if (context_bins > 1) {
        // 2D output: context_bin <tab> gc_bin <tab> probability
        for (int ctx = 0; ctx < context_bins; ctx++) {
            for (int gc = 0; gc < GC_BINS; gc++) {
                int idx = ctx * GC_BINS + gc;
                double prob = total > 0 ? expected_gc[idx] / total : 0.0;
                fprintf(out, "%d\t%d\t%.10f\n", ctx, gc, prob);
            }
        }
    } else {
        // 1D output: gc_bin <tab> probability
        for (int i = 0; i < GC_BINS; i++) {
            double prob = total > 0 ? expected_gc[i] / total : 0.0;
            fprintf(out, "%d\t%.10f\n", i, prob);
        }
    }
    
    if (out != stdout) {
        fclose(out);
    }
    
    return 0;
}
