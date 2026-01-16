#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iomanip>

// htslib includes
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/tbx.h"

// Binomial math functions (defined in binomial_math.cpp)
double log_binom_pmf(uint32_t n, uint32_t k, double p);
double log_binom_tail_cdf(uint32_t n, uint32_t k, double p);
double logsumexp(double a, double b);

// Forward declarations
struct Locus {
    std::string chrom;
    uint32_t start;  // 0-based
    uint32_t end;     // 0-based, exclusive
};

struct PileupCounts {
    uint32_t n = 0;        // coverage
    uint32_t k_any = 0;    // any mismatch
    uint32_t k_conv = 0;   // T->C or A->G conversion
    char ref_base = 'N';
    // Track mismatch base counts (A/C/G/T) to report a reasonable ALT base.
    // Index mapping: 0=A,1=C,2=G,3=T
    uint32_t mis_acgt[4] = {0, 0, 0, 0};
};

struct Config {
    // Inputs
    std::string bam_path;
    std::string bed_path;
    std::string ref_fasta;
    
    // Counting modes
    std::string k_mode = "conv";  // "conv" or "any"
    bool include_secondary = false;
    bool include_supplementary = false;
    int min_mapq = 0;
    int min_baseq = 0;
    bool overlap_policy = false;  // bam_mplp_init_overlaps
    
    // Calling params (GEDI/STAR compat)
    double p_err = 0.001;
    double pval_threshold = 0.001;
    double min_tc_ratio = 0.3;
    uint32_t min_cov = 6;
    uint32_t min_alt = 1;
    
    // Outputs
    std::string output_bed;
    std::string output_debug_tsv;
};

// Parse BED file (0-based start, exclusive end)
std::vector<Locus> parseBed(const std::string& bed_path) {
    std::vector<Locus> loci;
    std::ifstream bed(bed_path);
    if (!bed.good()) {
        std::cerr << "ERROR: Cannot open BED file: " << bed_path << std::endl;
        return loci;
    }
    
    std::string line;
    while (std::getline(bed, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        Locus loc;
        if (iss >> loc.chrom >> loc.start >> loc.end) {
            loci.push_back(loc);
        }
    }
    return loci;
}

// Convert base char to integer (A=0, C=1, G=2, T=3, N=4)
int baseToInt(char b) {
    switch (b) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// Convert htslib BAM nt16 base encoding (bam_seqi) to A/C/G/T index (0..3) or 4 for non-ACGT.
// bam_seqi returns values in nt16 space (A=1,C=2,G=4,T=8,N=15,...).
// We use seq_nt16_str[] from htslib to map nt16->char, then baseToInt().
static inline int bamNt16ToBaseInt(uint8_t bam_nt16) {
    // seq_nt16_str is declared by htslib (sam.h/sam.c)
    char base = seq_nt16_str[bam_nt16];
    return baseToInt(base);
}

char intToBase(int i) {
    switch (i) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

// Check if base pair is a conversion (T->C or A->G)
bool isConversion(int ref_int, int read_int) {
    return (ref_int == 3 && read_int == 1) ||  // T->C
           (ref_int == 0 && read_int == 2);    // A->G
}

static inline char pickAltBase(const PileupCounts& counts) {
    // Pick the most frequent non-ref mismatch among A/C/G/T.
    int ref_i = baseToInt(counts.ref_base);
    uint32_t best = 0;
    int best_i = -1;
    for (int i = 0; i < 4; ++i) {
        if (i == ref_i) continue;
        if (counts.mis_acgt[i] > best) {
            best = counts.mis_acgt[i];
            best_i = i;
        }
    }
    return (best_i >= 0) ? intToBase(best_i) : 'N';
}

// Count bases at a single position using pileup
PileupCounts countPosition(htsFile* bam_file, bam_hdr_t* header, hts_idx_t* idx, faidx_t* fai,
                          const std::string& chrom, uint32_t pos0, const Config& cfg) {
    PileupCounts counts;
    
    // Get reference base (htslib faidx uses 0-based)
    int ref_len = 0;
    char* ref_seq = faidx_fetch_seq(fai, chrom.c_str(), static_cast<int>(pos0), static_cast<int>(pos0), &ref_len);
    if (!ref_seq || ref_len <= 0) {
        counts.ref_base = 'N';
        if (ref_seq) free(ref_seq);
        return counts;
    }
    counts.ref_base = ref_seq[0];
    int ref_int = baseToInt(counts.ref_base);
    free(ref_seq);
    
    // Build region string for hts_itr_query (htslib uses 1-based)
    std::ostringstream reg_oss;
    reg_oss << chrom << ":" << (pos0 + 1) << "-" << (pos0 + 1);
    std::string region = reg_oss.str();
    
    // Query iterator (sam_itr_querys takes idx first, then header)
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    if (!iter) {
        return counts;
    }
    
    bam1_t* b = bam_init1();
    
    // Iterate reads covering this position
    while (sam_itr_next(bam_file, iter, b) >= 0) {
        // Filter by alignment flags
        if (!cfg.include_secondary && (b->core.flag & BAM_FSECONDARY)) continue;
        if (!cfg.include_supplementary && (b->core.flag & BAM_FSUPPLEMENTARY)) continue;
        if (b->core.flag & BAM_FUNMAP) continue;
        
        // Filter by MAPQ
        if (b->core.qual < cfg.min_mapq) continue;
        
        // Check if read covers position (b->core.pos is 0-based, can be negative)
        if (static_cast<int32_t>(b->core.pos) > static_cast<int32_t>(pos0)) continue;
        
        // Compute reference position covered by this read
        uint32_t* cigar = bam_get_cigar(b);
        int32_t ref_pos = b->core.pos;
        int read_pos = 0;
        
        for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                // Check if our target position is in this match block
                if (ref_pos <= static_cast<int32_t>(pos0) && static_cast<int32_t>(pos0) < ref_pos + len) {
                    int offset = static_cast<int32_t>(pos0) - ref_pos;
                    uint8_t* seq = bam_get_seq(b);
                    uint8_t* qual = bam_get_qual(b);
                    
                    if (read_pos + offset < b->core.l_qseq) {
                        uint8_t bam_nt16 = bam_seqi(seq, read_pos + offset);
                        int read_int = bamNt16ToBaseInt(bam_nt16); // 0..3 for A/C/G/T, else 4
                        uint8_t base_qual = qual ? qual[read_pos + offset] : 0;
                        // Missing qualities are encoded as 255 in BAM
                        if (base_qual == 255) base_qual = 0;
                        
                        // Filter by base quality
                        if (base_qual >= cfg.min_baseq) {
                            // Only count A/C/G/T in coverage (exclude N and other codes)
                            if (read_int < 4) {
                                counts.n++;
                                if (read_int != ref_int) {
                                    counts.k_any++;
                                    counts.mis_acgt[read_int] += 1;
                                    if (isConversion(ref_int, read_int)) {
                                        counts.k_conv++;
                                    }
                                }
                            }
                        }
                    }
                    break;  // Found the position, done with this read
                }
                ref_pos += len;
                read_pos += len;
            } else if (op == BAM_CINS) {
                read_pos += len;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                if (ref_pos <= static_cast<int32_t>(pos0) && static_cast<int32_t>(pos0) < ref_pos + len) {
                    // Position is in a deletion - skip this read
                    break;
                }
                ref_pos += len;
            } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
                read_pos += (op == BAM_CSOFT_CLIP) ? len : 0;
            }
        }
    }
    
    bam_destroy1(b);
    hts_itr_destroy(iter);
    
    return counts;
}

// Binomial caller: determine if site passes filters
bool callSnp(const PileupCounts& counts, const Config& cfg) {
    uint32_t n = counts.n;
    uint32_t k = (cfg.k_mode == "any") ? counts.k_any : counts.k_conv;
    
    // Apply filters
    if (n < cfg.min_cov) return false;
    if (k < cfg.min_alt) return false;
    
    double f = (n > 0) ? (static_cast<double>(k) / static_cast<double>(n)) : 0.0;
    if (f < cfg.min_tc_ratio) return false;
    
    // Compute binomial p-value
    double log_pval = log_binom_tail_cdf(n, k, cfg.p_err);
    double pval = std::exp(log_pval);
    
    return pval < cfg.pval_threshold;
}

// Write BED output (bgzip + tabix)
bool writeBed(const std::vector<std::pair<Locus, PileupCounts>>& called_sites,
              const Config& cfg) {
    BGZF* bgzf = bgzf_open(cfg.output_bed.c_str(), "w");
    if (!bgzf) {
        std::cerr << "ERROR: Cannot open BED output: " << cfg.output_bed << std::endl;
        return false;
    }
    
    // Write header
    std::string header = "#chrom\tstart\tend\tref\talt\tn\tk\tf\tpval\n";
    bgzf_write(bgzf, header.c_str(), header.size());
    
    for (const auto& site : called_sites) {
        const Locus& loc = site.first;
        const PileupCounts& counts = site.second;
        
        uint32_t k = (cfg.k_mode == "any") ? counts.k_any : counts.k_conv;
        double f = (counts.n > 0) ? (static_cast<double>(k) / static_cast<double>(counts.n)) : 0.0;
        double log_pval = log_binom_tail_cdf(counts.n, k, cfg.p_err);
        double pval = std::exp(log_pval);
        
        char alt = pickAltBase(counts);
        
        std::ostringstream oss;
        oss << loc.chrom << "\t"
            << loc.start << "\t"
            << loc.end << "\t"
            << counts.ref_base << "\t"
            << alt << "\t"
            << counts.n << "\t"
            << k << "\t"
            << std::fixed << std::setprecision(6) << f << "\t"
            // pval in scientific to avoid printing 0.000000 for tiny values
            << std::scientific << std::setprecision(6) << pval << "\n";
        
        std::string line = oss.str();
        bgzf_write(bgzf, line.c_str(), line.size());
    }
    
    bgzf_close(bgzf);
    
    // Build tabix index
    // Use htslib's standard BED config (seq=0, beg=1, end=2, meta='#')
    int ret = tbx_index_build(cfg.output_bed.c_str(), 0, &tbx_conf_bed);
    if (ret != 0) {
        std::cerr << "WARNING: Failed to build tabix index for: " << cfg.output_bed << std::endl;
    }
    
    return true;
}

// Write debug TSV
bool writeDebugTsv(const std::vector<std::pair<Locus, PileupCounts>>& all_sites,
                   const Config& cfg) {
    if (cfg.output_debug_tsv.empty()) return true;
    
    std::ofstream out(cfg.output_debug_tsv);
    if (!out.good()) {
        std::cerr << "ERROR: Cannot open debug TSV: " << cfg.output_debug_tsv << std::endl;
        return false;
    }
    
    out << "chrom\tpos0\tref\talt\tn\tk_any\tk_conv\tf_any\tf_conv\tpval_any\tpval_conv\tcalled_any\tcalled_conv\n";
    
    for (const auto& site : all_sites) {
        const Locus& loc = site.first;
        const PileupCounts& counts = site.second;
        
        double f_any = (counts.n > 0) ? (static_cast<double>(counts.k_any) / static_cast<double>(counts.n)) : 0.0;
        double f_conv = (counts.n > 0) ? (static_cast<double>(counts.k_conv) / static_cast<double>(counts.n)) : 0.0;
        
        double log_pval_any = log_binom_tail_cdf(counts.n, counts.k_any, cfg.p_err);
        double pval_any = std::exp(log_pval_any);
        double log_pval_conv = log_binom_tail_cdf(counts.n, counts.k_conv, cfg.p_err);
        double pval_conv = std::exp(log_pval_conv);
        
        Config cfg_any = cfg;
        cfg_any.k_mode = "any";
        bool called_any = callSnp(counts, cfg_any);
        
        Config cfg_conv = cfg;
        cfg_conv.k_mode = "conv";
        bool called_conv = callSnp(counts, cfg_conv);
        
        char alt = pickAltBase(counts);
        out << loc.chrom << "\t"
            << loc.start << "\t"
            << counts.ref_base << "\t"
            << alt << "\t"
            << counts.n << "\t"
            << counts.k_any << "\t"
            << counts.k_conv << "\t"
            << std::fixed << std::setprecision(6) << f_any << "\t"
            << std::fixed << std::setprecision(6) << f_conv << "\t"
            << std::scientific << std::setprecision(6) << pval_any << "\t"
            << std::scientific << std::setprecision(6) << pval_conv << "\t"
            << (called_any ? "1" : "0") << "\t"
            << (called_conv ? "1" : "0") << "\n";
    }
    
    return true;
}

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " [OPTIONS]\n"
              << "Options:\n"
              << "  --bam PATH          Input BAM file (required)\n"
              << "  --bed PATH          Candidate BED file (required)\n"
              << "  --ref PATH          Reference FASTA (required)\n"
              << "  --output PATH        Output BED file (bgzip+tabix) (required)\n"
              << "  --debug-tsv PATH     Optional debug TSV with per-site counts\n"
              << "\n"
              << "Counting modes:\n"
              << "  --kMode MODE        'conv' (T->C/A->G only) or 'any' (any mismatch) [default: conv]\n"
              << "  --includeSecondary  0/1 Include secondary alignments [default: 0]\n"
              << "  --includeSupplementary 0/1 Include supplementary [default: 0]\n"
              << "  --minMapQ INT        Minimum MAPQ [default: 0]\n"
              << "  --minBaseQ INT       Minimum base quality [default: 0]\n"
              << "\n"
              << "Calling params (GEDI/STAR compat):\n"
              << "  --pErr FLOAT         Error rate for binomial test [default: 0.001]\n"
              << "  --pval FLOAT         P-value threshold [default: 0.001]\n"
              << "  --minTcRatio FLOAT   Minimum conversion ratio [default: 0.3]\n"
              << "  --minCov INT         Minimum coverage [default: 6]\n"
              << "  --minAlt INT         Minimum alternative count [default: 1]\n";
}

int main(int argc, char* argv[]) {
    Config cfg;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--bam" && i + 1 < argc) {
            cfg.bam_path = argv[++i];
        } else if (arg == "--bed" && i + 1 < argc) {
            cfg.bed_path = argv[++i];
        } else if (arg == "--ref" && i + 1 < argc) {
            cfg.ref_fasta = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            cfg.output_bed = argv[++i];
        } else if (arg == "--debug-tsv" && i + 1 < argc) {
            cfg.output_debug_tsv = argv[++i];
        } else if (arg == "--kMode" && i + 1 < argc) {
            cfg.k_mode = argv[++i];
            if (cfg.k_mode != "conv" && cfg.k_mode != "any") {
                std::cerr << "ERROR: --kMode must be 'conv' or 'any'\n";
                return 1;
            }
        } else if (arg == "--includeSecondary" && i + 1 < argc) {
            cfg.include_secondary = (std::atoi(argv[++i]) != 0);
        } else if (arg == "--includeSupplementary" && i + 1 < argc) {
            cfg.include_supplementary = (std::atoi(argv[++i]) != 0);
        } else if (arg == "--minMapQ" && i + 1 < argc) {
            cfg.min_mapq = std::atoi(argv[++i]);
        } else if (arg == "--minBaseQ" && i + 1 < argc) {
            cfg.min_baseq = std::atoi(argv[++i]);
        } else if (arg == "--pErr" && i + 1 < argc) {
            cfg.p_err = std::atof(argv[++i]);
        } else if (arg == "--pval" && i + 1 < argc) {
            cfg.pval_threshold = std::atof(argv[++i]);
        } else if (arg == "--minTcRatio" && i + 1 < argc) {
            cfg.min_tc_ratio = std::atof(argv[++i]);
        } else if (arg == "--minCov" && i + 1 < argc) {
            cfg.min_cov = std::atoi(argv[++i]);
        } else if (arg == "--minAlt" && i + 1 < argc) {
            cfg.min_alt = std::atoi(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "ERROR: Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate required args
    if (cfg.bam_path.empty() || cfg.bed_path.empty() || cfg.ref_fasta.empty() || cfg.output_bed.empty()) {
        std::cerr << "ERROR: Missing required arguments\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Load BED
    std::vector<Locus> loci = parseBed(cfg.bed_path);
    if (loci.empty()) {
        std::cerr << "ERROR: No loci found in BED file\n";
        return 1;
    }
    
    // Open BAM
    htsFile* bam_file = sam_open(cfg.bam_path.c_str(), "r");
    if (!bam_file) {
        std::cerr << "ERROR: Cannot open BAM: " << cfg.bam_path << std::endl;
        return 1;
    }
    
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (!header) {
        std::cerr << "ERROR: Cannot read BAM header\n";
        sam_close(bam_file);
        return 1;
    }
    
    // Load reference index
    faidx_t* fai = fai_load(cfg.ref_fasta.c_str());
    if (!fai) {
        std::cerr << "ERROR: Cannot load FASTA index. Run: samtools faidx " << cfg.ref_fasta << std::endl;
        bam_hdr_destroy(header);
        sam_close(bam_file);
        return 1;
    }
    
    // Load BAM index
    hts_idx_t* idx = sam_index_load(bam_file, cfg.bam_path.c_str());
    if (!idx) {
        std::cerr << "ERROR: Cannot load BAM index. Run: samtools index " << cfg.bam_path << std::endl;
        fai_destroy(fai);
        bam_hdr_destroy(header);
        sam_close(bam_file);
        return 1;
    }
    
    // Process each locus
    std::vector<std::pair<Locus, PileupCounts>> all_sites;
    std::vector<std::pair<Locus, PileupCounts>> called_sites;
    
    for (const Locus& loc : loci) {
        for (uint32_t pos = loc.start; pos < loc.end; pos++) {
            PileupCounts counts = countPosition(bam_file, header, idx, fai, loc.chrom, pos, cfg);
            Locus pos_locus = {loc.chrom, pos, pos + 1};
            all_sites.push_back({pos_locus, counts});
            
            if (callSnp(counts, cfg)) {
                called_sites.push_back({pos_locus, counts});
            }
        }
    }
    
    // Write outputs
    if (!writeBed(called_sites, cfg)) {
        std::cerr << "ERROR: Failed to write BED output\n";
        return 1;
    }
    
    if (!writeDebugTsv(all_sites, cfg)) {
        std::cerr << "ERROR: Failed to write debug TSV\n";
        return 1;
    }
    
    std::cerr << "Processed " << all_sites.size() << " positions, called " << called_sites.size() << " SNPs\n";
    
    // Cleanup
    hts_idx_destroy(idx);
    fai_destroy(fai);
    bam_hdr_destroy(header);
    sam_close(bam_file);
    
    return 0;
}
