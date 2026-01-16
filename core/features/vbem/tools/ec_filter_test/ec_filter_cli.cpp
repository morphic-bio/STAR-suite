// EC Filter CLI - BAM parsing and Salmon-format EC output
// Implements streaming BAM → alignment-mode EC building → eq_classes.txt (no pre-filtering)

#include "libem/alignment_filter.h"
#include "libem/ec_builder.h"
#include "libem/extended_pruning.h"
#include "libem/alignment_model.h"
#include "libem/forgetting_mass.h"
#include "libem/em_types.h"
#include "libem/gc_bias.h"
#include "htslib/htslib/sam.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <random>
#include <map>
#include <cstring>

using libem::AlignmentModel;
using libem::Transcriptome;
using libem::TranscriptSequence;

// Helper to get AS tag from BAM record
int32_t get_alignment_score(bam1_t* b) {
    uint8_t* as_tag = bam_aux_get(b, "AS");
    if (as_tag) {
        return bam_aux2i(as_tag);
    }
    return 0;  // Default to 0 if AS tag missing
}

// Aligner detection for error-model behavior
enum class AlignerType {
    UNKNOWN,
    STAR,
    BOWTIE2,
    RAPMAP,
    PUFFERFISH
};

AlignerType detect_aligner_type(const bam_hdr_t* header) {
    if (!header || !header->text) {
        return AlignerType::UNKNOWN;
    }
    std::string text(header->text);
    if (text.find("PN:STAR") != std::string::npos || text.find("ID:STAR") != std::string::npos) {
        return AlignerType::STAR;
    }
    if (text.find("PN:bowtie2") != std::string::npos || text.find("ID:bowtie2") != std::string::npos) {
        return AlignerType::BOWTIE2;
    }
    if (text.find("PN:RapMap") != std::string::npos || text.find("PN:rapmap") != std::string::npos) {
        return AlignerType::RAPMAP;
    }
    if (text.find("PN:pufferfish") != std::string::npos || text.find("PN:Pufferfish") != std::string::npos) {
        return AlignerType::PUFFERFISH;
    }
    return AlignerType::UNKNOWN;
}

const char* aligner_type_name(AlignerType type) {
    switch (type) {
        case AlignerType::STAR: return "STAR";
        case AlignerType::BOWTIE2: return "Bowtie2";
        case AlignerType::RAPMAP: return "RapMap";
        case AlignerType::PUFFERFISH: return "Pufferfish";
        default: return "Unknown";
    }
}

// Determine library format from paired-end read positions and orientations
// Ported from Salmon's hitType() function
LibraryFormat hitType(int32_t end1Start, bool end1Fwd, int32_t end2Start, bool end2Fwd) {
    // If reads come from opposite strands
    if (end1Fwd != end2Fwd) {
        // Read 1 from forward strand
        if (end1Fwd) {
            if (end1Start <= end2Start) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);  // ISF
            } else {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);     // OSF
            }
        }
        // Read 2 from forward strand
        if (end2Fwd) {
            if (end2Start <= end1Start) {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);  // ISR
            } else {
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);    // OSR
            }
        }
    } else {
        // Reads from same strand
        if (end1Fwd) {
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);  // MSF
        } else {
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);  // MSR
        }
    }
    // Default fallback
    return LibraryFormat::IU();
}

// Auto-detect library format by sampling BAM reads
// Note: This function reads from the file and should be called before main processing
LibraryFormat detectLibraryFormat(const char* bam_path, int sample_size = 1000) {
    std::map<uint8_t, int> type_counts;
    
    // Open BAM file for detection
    samFile* fp = sam_open(bam_path, "r");
    if (!fp) {
        return LibraryFormat::IU();  // Default fallback
    }
    
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        return LibraryFormat::IU();
    }
    
    // Map to store mates per qname: qname -> {read1, read2}
    std::unordered_map<std::string, std::pair<bam1_t*, bam1_t*>> mate_map;
    
    bam1_t* b = bam_init1();
    int count = 0;
    
    while (count < sample_size && sam_read1(fp, hdr, b) >= 0) {
        uint32_t flag = b->core.flag;
        
        // Skip unpaired, secondary, and supplementary alignments
        if (!(flag & BAM_FPAIRED) || (flag & BAM_FSECONDARY) || (flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }
        
        const char* qname = bam_get_qname(b);
        std::string qname_str(qname);
        
        // Determine if this is read1 or read2 using BAM_FREAD1/BAM_FREAD2 flags
        bool is_read1 = !(flag & BAM_FREAD2);  // BAM_FREAD2 is set for read2
        bool is_read2 = (flag & BAM_FREAD2) != 0;
        
        // Check if we already have a mate for this qname
        auto it = mate_map.find(qname_str);
        if (it == mate_map.end()) {
            // First mate seen - store it
            bam1_t* stored = bam_dup1(b);
            if (is_read1) {
                mate_map[qname_str] = std::make_pair(stored, nullptr);
            } else {
                mate_map[qname_str] = std::make_pair(nullptr, stored);
            }
        } else {
            // Second mate - check if we have the opposite read
            bam1_t*& stored_r1 = it->second.first;
            bam1_t*& stored_r2 = it->second.second;
            
            if (is_read1 && stored_r1 == nullptr && stored_r2 != nullptr) {
                // We have read2, now we have read1
                stored_r1 = bam_dup1(b);
            } else if (is_read2 && stored_r2 == nullptr && stored_r1 != nullptr) {
                // We have read1, now we have read2
                stored_r2 = bam_dup1(b);
            } else {
                // Skip r1+r1 or r2+r2 pairs, or already have both
                continue;
            }
            
            // Now we have both read1 and read2 - process the pair
            bam1_t* r1 = stored_r1;
            bam1_t* r2 = stored_r2;
            
            if (r1 && r2) {
                // Both mapped and proper pair
                bool mapped1 = !(r1->core.flag & BAM_FUNMAP);
                bool mapped2 = !(r2->core.flag & BAM_FUNMAP);
                bool proper1 = (r1->core.flag & BAM_FPROPER_PAIR) != 0;
                bool proper2 = (r2->core.flag & BAM_FPROPER_PAIR) != 0;
                
                if (mapped1 && mapped2 && proper1 && proper2 && r1->core.tid == r2->core.tid) {
                    bool fwd1 = !(r1->core.flag & BAM_FREVERSE);
                    bool fwd2 = !(r2->core.flag & BAM_FREVERSE);
                    int32_t pos1 = r1->core.pos;
                    int32_t pos2 = r2->core.pos;
                    
                    LibraryFormat observed = hitType(pos1, fwd1, pos2, fwd2);
                    type_counts[observed.typeId()]++;
                    count++;
                }
                
                // Clean up stored records
                bam_destroy1(stored_r1);
                bam_destroy1(stored_r2);
                mate_map.erase(it);
            }
        }
    }
    
    // Free any pending stored records
    for (auto& pair : mate_map) {
        if (pair.second.first) {
            bam_destroy1(pair.second.first);
        }
        if (pair.second.second) {
            bam_destroy1(pair.second.second);
        }
    }
    mate_map.clear();
    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    
    if (type_counts.empty()) {
        // Default to IU (inward unstranded) if no samples
        return LibraryFormat::IU();
    }
    
    // Find most common type
    uint8_t most_common_id = type_counts.begin()->first;
    int max_count = type_counts.begin()->second;
    for (const auto& pair : type_counts) {
        if (pair.second > max_count) {
            max_count = pair.second;
            most_common_id = pair.first;
        }
    }
    
    // Convert ID back to LibraryFormat
    ReadType rt = static_cast<ReadType>(most_common_id & 0x01);
    ReadOrientation ro = static_cast<ReadOrientation>((most_common_id >> 1) & 0x3);
    ReadStrandedness rs = static_cast<ReadStrandedness>((most_common_id >> 3) & 0x7);
    
    LibraryFormat detected(rt, ro, rs);
    
    // For unstranded detection, if we detect ISF/ISR but most reads are unstranded-compatible,
    // default to IU (common case)
    if (detected.orientation == ReadOrientation::TOWARD && 
        (detected.strandedness == ReadStrandedness::SA || detected.strandedness == ReadStrandedness::AS)) {
        // Check if IU would also be compatible with most reads
        int iu_compatible = 0;
        for (const auto& pair : type_counts) {
            LibraryFormat fmt = LibraryFormat::formatFromID(pair.first);
            if (fmt.orientation == ReadOrientation::TOWARD && fmt.strandedness == ReadStrandedness::U) {
                iu_compatible += pair.second;
            }
        }
        // If IU-compatible reads are significant, prefer IU
        if (iu_compatible > max_count * 0.5) {
            return LibraryFormat::IU();
        }
    }
    
    return detected;
}

// Helper to parse library format string (IU, ISF, ISR, U, etc.)
LibraryFormat parseLibraryFormatString(const std::string& fmt_str) {
    std::string fmt = fmt_str;
    for (auto& c : fmt) {
        c = std::toupper(c);
    }
    
    if (fmt == "IU") return LibraryFormat::IU();
    if (fmt == "ISF") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
    if (fmt == "ISR") return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
    if (fmt == "U") return LibraryFormat::U();
    
    // Default to IU
    return LibraryFormat::IU();
}

// Per-alignment BAM records for error model updates
struct AlignmentBamData {
    bam1_t* read1;
    bam1_t* read2;
    bool paired;
};

static bam1_t* clone_bam_record(bam1_t* src) {
    return src ? bam_dup1(src) : nullptr;
}

static void free_alignment_bams(std::vector<AlignmentBamData>& bam_data) {
    for (auto& entry : bam_data) {
        if (entry.read1) {
            bam_destroy1(entry.read1);
        }
        if (entry.read2) {
            bam_destroy1(entry.read2);
        }
    }
    bam_data.clear();
}

static double compute_err_like_single(
    bam1_t* b,
    const Transcriptome& transcriptome,
    AlignmentModel& model,
    bool is_left,
    bool use_model
) {
    if (!use_model || !b) {
        return 0.0;
    }
    int32_t tid = b->core.tid;
    if (tid < 0) {
        return 0.0;
    }
    const TranscriptSequence* ref = transcriptome.getTranscript(static_cast<uint32_t>(tid));
    if (!ref) {
        return 0.0;
    }

    uint32_t* cigar = bam_get_cigar(b);
    uint32_t cigar_len = b->core.n_cigar;
    uint8_t* seq = bam_get_seq(b);
    int32_t read_len = b->core.l_qseq;
    int32_t ref_pos = b->core.pos;

    return model.logLikelihood(cigar, cigar_len, seq, read_len, ref, ref_pos, is_left);
}

static double compute_err_like_paired(
    bam1_t* b1,
    bam1_t* b2,
    const Transcriptome& transcriptome,
    AlignmentModel& model,
    bool use_model
) {
    if (!use_model || !b1 || !b2) {
        return 0.0;
    }
    int32_t tid = b1->core.tid;
    if (tid < 0) {
        return 0.0;
    }
    const TranscriptSequence* ref = transcriptome.getTranscript(static_cast<uint32_t>(tid));
    if (!ref) {
        return 0.0;
    }

    uint32_t* cigar1 = bam_get_cigar(b1);
    uint32_t cigar_len1 = b1->core.n_cigar;
    uint8_t* seq1 = bam_get_seq(b1);
    int32_t read_len1 = b1->core.l_qseq;
    int32_t pos1 = b1->core.pos;

    uint32_t* cigar2 = bam_get_cigar(b2);
    uint32_t cigar_len2 = b2->core.n_cigar;
    uint8_t* seq2 = bam_get_seq(b2);
    int32_t read_len2 = b2->core.l_qseq;
    int32_t pos2 = b2->core.pos;

    return model.logLikelihoodPaired(
        cigar1, cigar_len1, seq1, read_len1,
        cigar2, cigar_len2, seq2, read_len2,
        ref, pos1, pos2
    );
}

static void update_error_model_single(
    bam1_t* b,
    const Transcriptome& transcriptome,
    AlignmentModel& model,
    bool is_left,
    double aligner_score,
    double log_forgetting_mass
) {
    if (!b) {
        return;
    }
    int32_t tid = b->core.tid;
    if (tid < 0) {
        return;
    }
    const TranscriptSequence* ref = transcriptome.getTranscript(static_cast<uint32_t>(tid));
    if (!ref) {
        return;
    }

    uint32_t* cigar = bam_get_cigar(b);
    uint32_t cigar_len = b->core.n_cigar;
    uint8_t* seq = bam_get_seq(b);
    int32_t read_len = b->core.l_qseq;
    int32_t ref_pos = b->core.pos;

    model.update(cigar, cigar_len, seq, read_len, ref, ref_pos, is_left,
                 aligner_score, log_forgetting_mass);
}

static void update_error_model_paired(
    bam1_t* b1,
    bam1_t* b2,
    const Transcriptome& transcriptome,
    AlignmentModel& model,
    double aligner_score,
    double log_forgetting_mass
) {
    if (!b1 || !b2) {
        return;
    }
    bool read1_is_left = (b1->core.pos <= b2->core.pos);
    bam1_t* left = read1_is_left ? b1 : b2;
    bam1_t* right = read1_is_left ? b2 : b1;

    update_error_model_single(left, transcriptome, model, true,
                              aligner_score, log_forgetting_mass);
    update_error_model_single(right, transcriptome, model, false,
                              aligner_score, log_forgetting_mass);
}

// Read decoy list from file (one transcript name per line)
std::unordered_set<std::string> load_decoy_list(const std::string& decoy_file) {
    std::unordered_set<std::string> decoys;
    if (decoy_file.empty()) {
        return decoys;
    }
    
    std::ifstream file(decoy_file);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open decoy file: " << decoy_file << std::endl;
        return decoys;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            decoys.insert(line);
        }
    }
    return decoys;
}

// Determine mate status from BAM flags
MateStatus determine_mate_status(uint32_t flag1, uint32_t flag2, bool both_mapped, bool proper_pair) {
    if (!both_mapped) {
        // Orphan: determine left/right from flag
        if (flag1 & BAM_FREAD1) {
            return MateStatus::PAIRED_END_LEFT;
        } else {
            return MateStatus::PAIRED_END_RIGHT;
        }
    }
    
    if (proper_pair && (flag1 & BAM_FPROPER_PAIR)) {
        return MateStatus::PAIRED_END_PAIRED;
    }
    
    // Discordant pair
    if (flag1 & BAM_FREAD1) {
        return MateStatus::PAIRED_END_LEFT;
    } else {
        return MateStatus::PAIRED_END_RIGHT;
    }
}

static void process_read_group_error_model(
    std::vector<RawAlignment>& alignments,
    std::vector<AlignmentBamData>& bam_data,
    const Transcriptome& transcriptome,
    AlignmentModel* alignment_model,
    const ECBuilderParams& ec_params,
    AlignerType aligner_type,
    double log_forgetting_mass,
    std::default_random_engine& rng,
    std::uniform_real_distribution<double>& uni,
    const std::string& qname,  // For tracing
    bool deterministic_updates  // If true, update ALL alignments
) {
    if (!alignment_model || alignments.empty()) {
        return;
    }

    const bool use_model = ec_params.use_error_model && alignment_model->useModel();
    
    // Compute errLike for each alignment (no updates here)
    // Only set has_err_like = true when the model is actually being used (after pre-burnin)
    for (size_t i = 0; i < alignments.size(); ++i) {
        RawAlignment& aln = alignments[i];
        const AlignmentBamData& bam_entry = bam_data[i];
        
        // Set trace context for this alignment
        if (alignment_model && !qname.empty()) {
            alignment_model->setTraceContext(qname, aln.transcript_id);
        }
        
        double err_like = 0.0;
        if (ec_params.use_error_model && use_model) {
            // Model is active - compute real errLike
            if (bam_entry.paired) {
                err_like = compute_err_like_paired(
                    bam_entry.read1, bam_entry.read2,
                    transcriptome, *alignment_model, use_model);
            } else {
                bool is_left = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                                aln.mate_status == MateStatus::SINGLE_END);
                err_like = compute_err_like_single(
                    bam_entry.read1, transcriptome, *alignment_model,
                    is_left, use_model);
            }
            aln.err_like = err_like;
            aln.has_err_like = true;
        } else {
            // Model not active yet (pre-burnin) or not enabled
            aln.has_err_like = false;
        }
    }

    // Compute auxProbs to get normalized weights for sampling
    ReadMapping mapping = computeAuxProbs(alignments, ec_params, false);
    if (!ec_params.use_error_model || !alignment_model->canUpdate()) {
        return;
    }

    for (size_t i = 0; i < mapping.alignment_indices.size(); ++i) {
        size_t aln_idx = mapping.alignment_indices[i];
        
        // Stochastic sampling: only update some alignments based on weight
        // Skip sampling if deterministic_updates is set
        if (!deterministic_updates) {
            double weight = mapping.aux_probs[i];
            double r = uni(rng);
            if (r >= weight) {
                continue;
            }
        }

        const RawAlignment& aln = alignments[aln_idx];
        const AlignmentBamData& bam_entry = bam_data[aln_idx];
        double aligner_score = (aligner_type == AlignerType::BOWTIE2)
            ? static_cast<double>(aln.score)
            : 0.0;

        if (bam_entry.paired) {
            update_error_model_paired(
                bam_entry.read1, bam_entry.read2,
                transcriptome, *alignment_model,
                aligner_score, log_forgetting_mass);
        } else {
            bool is_left = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                            aln.mate_status == MateStatus::SINGLE_END);
            update_error_model_single(
                bam_entry.read1, transcriptome, *alignment_model,
                is_left, aligner_score, log_forgetting_mass);
        }
    }
    
    // Trace read summary if enabled
    if (alignment_model && !qname.empty()) {
        // Sum errLike values (which are in log space)
        double errLikeSum = -std::numeric_limits<double>::infinity();  // LOG_0
        for (const auto& aln : alignments) {
            if (aln.has_err_like) {
                // Log-space addition
                if (errLikeSum == -std::numeric_limits<double>::infinity()) {
                    errLikeSum = aln.err_like;
                } else {
                    // logAdd: log(exp(a) + exp(b)) = max(a,b) + log(1 + exp(min(a,b) - max(a,b)))
                    if (aln.err_like == -std::numeric_limits<double>::infinity()) {
                        continue;  // Skip LOG_0 values
                    }
                    double max_val = std::max(errLikeSum, aln.err_like);
                    double min_val = std::min(errLikeSum, aln.err_like);
                    errLikeSum = max_val + std::log(1.0 + std::exp(min_val - max_val));
                }
            }
        }
        // If no valid errLike values, use LOG_1 (0.0)
        if (errLikeSum == -std::numeric_limits<double>::infinity()) {
            errLikeSum = 0.0;  // LOG_1
        }
        alignment_model->traceRead(qname, alignments.size(), errLikeSum, use_model, alignment_model->canUpdate());
        alignment_model->clearTraceContext();
    }
}

void print_usage(const char* prog_name) {
    std::cerr << "Usage: " << prog_name << " [options]\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  --input <file>              Input BAM file (required)\n";
    std::cerr << "  --transcripts <file>        Transcriptome FASTA file (optional, for validation)\n";
    std::cerr << "  --decoys <file>             Decoy transcript list (one name per line, optional)\n";
    std::cerr << "  --score-exp <float>         Score exponent (default: 1.0, used for AS-based errLike in computeAuxProbs)\n";
    std::cerr << "  --range-factorization-bins <int> Range factorization bins (default: 4)\n";
    std::cerr << "  --no-local-pruning          Disable local pruning (default: disabled)\n";
    std::cerr << "  --no-global-pruning         Disable global pruning (default: disabled)\n";
    std::cerr << "  --incompat-prior <float>    Incompatibility prior (default: 0.0 → skip incompatible)\n";
    std::cerr << "  --ignore-compat             Ignore compatibility check entirely (like --incompatPrior 1)\n";
    std::cerr << "  --lib-type <type>           Library type: A=auto-detect, IU, ISF, ISR, U (default: A)\n";
    std::cerr << "  --use-error-model           Use CIGAR-based error model (requires --transcripts)\n";
    std::cerr << "  --error-model-trace <file>  Trace error model computation (requires --use-error-model)\n";
    std::cerr << "  --trace-level <1|2|3>       Trace level: 1=read, 2=alignment, 3=transition (default: 1)\n";
    std::cerr << "  --dump-matrices <prefix>    Dump transition matrices at checkpoints (requires --use-error-model)\n";
    std::cerr << "  --read-order-file <file>    Process reads in order specified in file (for debugging)\n";
    std::cerr << "  --deterministic-updates     Update ALL alignments (disable stochastic sampling, for debugging)\n";
    std::cerr << "  --gc-bias                   Enable GC bias collection (requires --transcripts)\n";
    std::cerr << "  --expected-gc <file>        Expected GC distribution file (optional, for validation)\n";
    std::cerr << "  --observed-gc-out <file>   Output file for observed GC distribution\n";
    std::cerr << "  --discard-orphans           Discard orphan reads (default: false, keep orphans)\n";
    std::cerr << "  --paired                    Force paired-end mode\n";
    std::cerr << "  --single                    Force single-end mode\n";
    std::cerr << "  --trace-reads <file>        Trace read processing (alignment-mode format)\n";
    std::cerr << "  --trace-limit <N>           Stop tracing after N reads (default: no limit)\n";
    std::cerr << "  --threads <int>             Number of threads (default: 1)\n";
    std::cerr << "  --output-format <format>    Output format: salmon (default)\n";
    std::cerr << "  -o, --output <file>         Output EC file (required)\n";
    std::cerr << "  -h, --help                  Show this help\n";
    std::cerr << "\n";
    std::cerr << "Alignment-mode parity notes:\n";
    std::cerr << "  This CLI implements Salmon alignment-mode strictly (no pre-filtering).\n";
    std::cerr << "  Compatibility gating and auxProb computation happen in computeAuxProbs.\n";
    std::cerr << "  For Salmon parity, run Salmon with:\n";
    std::cerr << "    --noFragLengthDist --noLengthCorrection --noEffectiveLengthCorrection --incompatPrior 0\n";
    std::cerr << "  Use --lib-type A (auto-detect) unless you know the library type.\n";
    std::cerr << "  Default --incompat-prior 0 drops incompatible alignments (matches Salmon).\n";
    std::cerr << "  --ignore-compat keeps incompatible alignments (equivalent to --incompatPrior 1).\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse arguments
    std::string input_file;
    std::string transcripts_file;
    std::string decoy_file;
    std::string output_file;
    double score_exp = 1.0;
    uint32_t range_factorization_bins = 4;
    bool no_local_pruning = true;  // Default: disabled for Salmon parity
    bool no_global_pruning = true;  // Default: disabled for Salmon parity
    double incompat_prior = 0.0;    // Default: skip incompatible alignments (0.0 → ignore_incompat=true)
    bool ignore_compat = false;     // If true, bypass compatibility check entirely (like --incompatPrior 1)
    std::string lib_type = "A";     // Library type: A=auto, IU, ISF, ISR, U, etc.
    bool use_error_model = false;   // If true, use CIGAR-based error model
    std::string error_model_trace_file;  // Error model trace output file
    int trace_level = 1;            // Trace level: 1=read, 2=alignment, 3=transition
    std::string matrix_dump_prefix; // Prefix for matrix dump files
    size_t trace_limit = 0;         // 0 = no limit
    std::string output_format = "salmon";
    enum { AUTO, PAIRED, SINGLE } read_mode = AUTO;
    bool trace_reads = false;
    std::string trace_file;
    std::string read_order_file;    // File specifying read processing order (for debugging)
    bool deterministic_updates = false;  // If true, update ALL alignments (no stochastic sampling)
    bool gc_bias_enabled = false;   // If true, collect observed GC distribution
    std::string expected_gc_file;   // Expected GC distribution file (for loading)
    std::string observed_gc_out_file;  // Output file for observed GC distribution
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--input") == 0 && i + 1 < argc) {
            input_file = argv[++i];
        } else if (strcmp(argv[i], "--transcripts") == 0 && i + 1 < argc) {
            transcripts_file = argv[++i];
        } else if (strcmp(argv[i], "--decoys") == 0 && i + 1 < argc) {
            decoy_file = argv[++i];
        } else if (strcmp(argv[i], "--score-exp") == 0 && i + 1 < argc) {
            score_exp = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--range-factorization-bins") == 0 && i + 1 < argc) {
            range_factorization_bins = std::stoul(argv[++i]);
        } else if (strcmp(argv[i], "--no-local-pruning") == 0) {
            no_local_pruning = true;
        } else if (strcmp(argv[i], "--no-global-pruning") == 0) {
            no_global_pruning = true;
        } else if (strcmp(argv[i], "--incompat-prior") == 0 && i + 1 < argc) {
            incompat_prior = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--ignore-compat") == 0) {
            ignore_compat = true;
        } else if (strcmp(argv[i], "--lib-type") == 0 && i + 1 < argc) {
            lib_type = argv[++i];
        } else if (strcmp(argv[i], "--use-error-model") == 0) {
            use_error_model = true;
        } else if (strcmp(argv[i], "--error-model-trace") == 0 && i + 1 < argc) {
            error_model_trace_file = argv[++i];
        } else if (strcmp(argv[i], "--trace-level") == 0 && i + 1 < argc) {
            trace_level = std::stoi(argv[++i]);
            if (trace_level < 1 || trace_level > 3) {
                std::cerr << "Error: --trace-level must be 1, 2, or 3\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--dump-matrices") == 0 && i + 1 < argc) {
            matrix_dump_prefix = argv[++i];
        } else if (strcmp(argv[i], "--discard-orphans") == 0) {
            // TODO: implement discard_orphans flag (would set ec_params.discard_orphans = true)
        } else if (strcmp(argv[i], "--paired") == 0) {
            read_mode = PAIRED;
        } else if (strcmp(argv[i], "--single") == 0) {
            read_mode = SINGLE;
        } else if (strcmp(argv[i], "--trace-reads") == 0 && i + 1 < argc) {
            trace_reads = true;
            trace_file = argv[++i];
        } else if (strcmp(argv[i], "--read-order-file") == 0 && i + 1 < argc) {
            read_order_file = argv[++i];
        } else if (strcmp(argv[i], "--deterministic-updates") == 0) {
            deterministic_updates = true;
        } else if (strcmp(argv[i], "--gc-bias") == 0) {
            gc_bias_enabled = true;
        } else if (strcmp(argv[i], "--expected-gc") == 0 && i + 1 < argc) {
            expected_gc_file = argv[++i];
        } else if (strcmp(argv[i], "--observed-gc-out") == 0 && i + 1 < argc) {
            observed_gc_out_file = argv[++i];
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            // Threads parameter (not currently used, but accept for compatibility)
            std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--output-format") == 0 && i + 1 < argc) {
            output_format = argv[++i];
        } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    if (input_file.empty() || output_file.empty()) {
        std::cerr << "Error: --input and --output are required\n";
        print_usage(argv[0]);
        return 1;
    }
    
    if (output_format != "salmon") {
        std::cerr << "Error: Only 'salmon' output format is currently supported\n";
        return 1;
    }
    
    // Validate error model requirements
    if (use_error_model && transcripts_file.empty()) {
        std::cerr << "Error: --use-error-model requires --transcripts <fasta>\n";
        return 1;
    }
    
    // Validate GC bias requirements
    if (gc_bias_enabled && transcripts_file.empty()) {
        std::cerr << "Error: --gc-bias requires --transcripts <fasta>\n";
        return 1;
    }
    
    if (gc_bias_enabled && observed_gc_out_file.empty()) {
        std::cerr << "Error: --gc-bias requires --observed-gc-out <file>\n";
        return 1;
    }
    
    // Load transcriptome for error model and GC bias
    Transcriptome transcriptome;
    AlignmentModel* alignment_model = nullptr;
    GCFragModel* observed_gc = nullptr;
    
    if (use_error_model || gc_bias_enabled) {
        std::cerr << "Loading transcriptome...\n";
        if (!transcriptome.loadFromFasta(transcripts_file)) {
            std::cerr << "Error: Failed to load transcriptome from " << transcripts_file << std::endl;
            return 1;
        }
        
        if (use_error_model) {
            alignment_model = new AlignmentModel(1.0, 6);  // Salmon defaults: alpha=1.0, numErrorBins=6
            std::cerr << "Error model enabled (CIGAR-based)\n";
        }
        
        if (gc_bias_enabled) {
            observed_gc = new GCFragModel();
            std::cerr << "GC bias collection enabled\n";
        }
        
        // Set up error model tracing if requested
        if (!error_model_trace_file.empty()) {
            std::ofstream* trace_stream = new std::ofstream(error_model_trace_file);
            if (!trace_stream->is_open()) {
                std::cerr << "Error: Could not open error model trace file: " << error_model_trace_file << std::endl;
                delete trace_stream;
                delete alignment_model;
                return 1;
            }
            ErrorModelTraceLevel level = static_cast<ErrorModelTraceLevel>(trace_level);
            alignment_model->setTraceOutput(trace_stream, level);
            std::cerr << "Error model tracing enabled (level " << trace_level << ") -> " << error_model_trace_file << "\n";
        }
        
        // Set matrix dump prefix if requested
        if (!matrix_dump_prefix.empty()) {
            alignment_model->setMatrixDumpPrefix(matrix_dump_prefix);
            std::cerr << "Matrix dumping enabled -> " << matrix_dump_prefix << "_<checkpoint>_<left/right>_bin<N>.tsv\n";
        }
    }
    
    // Step 1: Open BAM file and read header
    samFile* bam_file = sam_open(input_file.c_str(), "r");
    if (!bam_file) {
        std::cerr << "Error: Failed to open BAM file: " << input_file << std::endl;
        return 1;
    }
    
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (!header) {
        std::cerr << "Error: Failed to read BAM header" << std::endl;
        sam_close(bam_file);
        return 1;
    }
    
    // Note: We assume BAM is unsorted or name-sorted for proper read grouping
    // Coordinate-sorted BAMs will break parity as read pairs won't be consecutive
    
    // Extract transcript names from BAM header
    size_t num_transcripts = header->n_targets;
    std::vector<std::string> transcript_names;
    transcript_names.reserve(num_transcripts);
    std::unordered_map<std::string, uint32_t> name_to_tid;
    
    for (size_t i = 0; i < num_transcripts; ++i) {
        std::string name(header->target_name[i]);
        transcript_names.push_back(name);
        name_to_tid[name] = static_cast<uint32_t>(i);
    }
    
    std::cerr << "Loaded " << num_transcripts << " transcripts from BAM header\n";

    // Detect aligner type from BAM header (affects error model / AS usage)
    AlignerType aligner_type = detect_aligner_type(header);
    std::cerr << "Detected aligner: " << aligner_type_name(aligner_type) << "\n";
    bool use_as_without_cigar = (aligner_type == AlignerType::RAPMAP ||
                                 aligner_type == AlignerType::PUFFERFISH);
    if (use_as_without_cigar && use_error_model) {
        std::cerr << "Warning: aligner does not provide CIGAR; disabling error model\n";
        use_error_model = false;
        if (alignment_model) {
            delete alignment_model;
            alignment_model = nullptr;
        }
    }
    
    // If error model or GC bias is enabled, reorder transcriptome to match BAM header order
    // This ensures transcript IDs from BAM match the transcriptome index
    if (alignment_model || gc_bias_enabled) {
        std::cerr << "Reordering transcriptome to match BAM header order...\n";
        if (!transcriptome.reorderByNames(transcript_names)) {
            std::cerr << "Warning: Some transcripts from BAM not found in FASTA\n";
        }
    }
    
    // Load decoy set
    std::unordered_set<std::string> decoy_names = load_decoy_list(decoy_file);
    std::vector<bool> is_decoy(num_transcripts, false);
    for (const auto& name : decoy_names) {
        auto it = name_to_tid.find(name);
        if (it != name_to_tid.end()) {
            is_decoy[it->second] = true;
        }
    }
    std::cerr << "Marked " << decoy_names.size() << " transcripts as decoys\n";
    
    // Step 2: Prepare filter and EC builder parameters
    // Alignment-mode: no pre-filtering (compatibility gating + auxProb only)
    
    ECBuilderParams ec_params;
    ec_params.use_range_factorization = (range_factorization_bins > 0);
    ec_params.range_factorization_bins = range_factorization_bins;
    ec_params.use_rank_eq_classes = false;  // Disabled for parity
    ec_params.use_frag_len_dist = false;   // Parity mode: no frag length dist (--noFragLengthDist)
    ec_params.use_error_model = use_error_model;   // Use CIGAR error model when enabled
    ec_params.use_as_without_cigar = use_as_without_cigar;
    ec_params.score_exp = score_exp;       // Score exponent for AS-based likelihood
    ec_params.incompat_prior = incompat_prior;  // Incompatibility prior (0.0 → skip incompatible)
    // If ignore_compat is set, don't skip any alignments (like --incompatPrior 1)
    // Otherwise, skip incompatible alignments when incompat_prior == 0.0
    ec_params.ignore_incompat = ignore_compat ? false : (incompat_prior == 0.0);
    // When ignore_compat is true, set incompat_prior to 1.0 to give no penalty
    if (ignore_compat) {
        ec_params.incompat_prior = 1.0;  // log(1.0) = 0.0, no penalty
    }
    // Detect or parse library format
    LibraryFormat detected_format = LibraryFormat::IU();  // Default
    if (lib_type == "A" || lib_type == "a" || lib_type == "auto") {
        std::cerr << "Auto-detecting library format from BAM...\n";
        detected_format = detectLibraryFormat(input_file.c_str(), 1000);
        std::cerr << "Detected library format: ";
        if (detected_format == LibraryFormat::IU()) std::cerr << "IU (Inward Unstranded)\n";
        else if (detected_format.strandedness == ReadStrandedness::SA && detected_format.orientation == ReadOrientation::TOWARD) std::cerr << "ISF (Inward Stranded Forward)\n";
        else if (detected_format.strandedness == ReadStrandedness::AS && detected_format.orientation == ReadOrientation::TOWARD) std::cerr << "ISR (Inward Stranded Reverse)\n";
        else std::cerr << "Other\n";
    } else {
        detected_format = parseLibraryFormatString(lib_type);
    }
    ec_params.lib_format = detected_format;
    ec_params.discard_orphans = false;    // Keep orphans (default false, matching Salmon)
    
    // Auto-detect compatibility gating: disable compat gating during auto-detect window
    // (matches Salmon's behavior: compat OFF during detection, ON after detection)
    // Only apply auto-detect gating when incompat_prior == 0.0 (default behavior)
    bool auto_detect_mode = (lib_type == "A" || lib_type == "a" || lib_type == "auto") && 
                            !ignore_compat && incompat_prior == 0.0;
    constexpr size_t AUTO_DETECT_SAMPLE_SIZE = 1000;  // Same as detectLibraryFormat sample size
    size_t reads_processed = 0;
    bool compat_gating_enabled = false;  // Will be set to true after auto-detect completes
    bool original_ignore_incompat = ec_params.ignore_incompat;  // Store original value
    double original_incompat_prior = ec_params.incompat_prior;  // Store original value
    
    if (auto_detect_mode) {
        // During auto-detect, disable compatibility gating
        // ignore_incompat = false means "don't skip incompatible" (gating OFF)
        // incompat_prior = 1.0 means no penalty for incompatible alignments
        ec_params.ignore_incompat = false;
        ec_params.incompat_prior = 1.0;
        std::cerr << "Auto-detect mode: compatibility gating disabled for first " 
                  << AUTO_DETECT_SAMPLE_SIZE << " reads\n";
    }
    
    ExtendedPruningParams pruning_params;
    pruning_params.enable_local_pruning = !no_local_pruning;
    pruning_params.enable_global_pruning = !no_global_pruning;
    
    // Step 3: Stream BAM records and group by read name
    std::vector<std::vector<RawAlignment>> read_alignments;
    std::vector<std::string> read_qnames;  // Track qnames for tracing
    bam1_t* b = bam_init1();
    
    std::string current_qname;
    std::vector<RawAlignment> current_read_alignments;
    std::vector<AlignmentBamData> current_read_bam_data;
    std::vector<bam1_t*> current_read_records;
    
    size_t total_reads = 0;
    size_t skipped_reads = 0;
    // Error model training schedule (Salmon-style minibatches)
    constexpr size_t mini_batch_size = 1000;  // Salmon default in alignment mode
    size_t batch_read_count = 0;
    double log_forgetting_mass = 0.0;
    uint64_t batch_timestep = 0;
    ForgettingMassCalculator fm_calc(0.65);  // Salmon default forgetting factor
    std::default_random_engine rng;
    std::uniform_real_distribution<double> uni(0.0, 1.0 + std::numeric_limits<double>::min());
    if (alignment_model) {
        // Use fixed seed (42) for deterministic updates - matches Salmon debug mode
        // This ensures the same alignments are selected for update in both implementations
        rng.seed(42);
        std::cerr << "Using fixed RNG seed (42) for deterministic error model updates\n";
    }

    auto build_alignments_from_records = [&](const std::vector<bam1_t*>& records) {
        current_read_alignments.clear();
        free_alignment_bams(current_read_bam_data);

        if (records.empty()) {
            return;
        }

        auto add_orphan = [&](bam1_t* rec, MateStatus ms) {
            int32_t tid = rec->core.tid;
            if (tid < 0) return;

            RawAlignment aln;
            aln.transcript_id = static_cast<uint32_t>(tid);
            aln.score = get_alignment_score(rec);
            aln.log_frag_prob = 0.0;
            aln.log_compat_prob = 0.0;
            aln.est_aln_prob = 1.0;
            aln.mate_status = ms;
            aln.is_forward = !(rec->core.flag & BAM_FREVERSE);
            aln.is_decoy = (tid >= 0 && tid < static_cast<int32_t>(is_decoy.size())) ? is_decoy[tid] : false;
            aln.fragment_len = -1;  // Orphan: no fragment length
            current_read_alignments.push_back(aln);
            if (alignment_model || gc_bias_enabled) {
                AlignmentBamData bam_entry{clone_bam_record(rec), nullptr, false};
                current_read_bam_data.push_back(bam_entry);
            }
        };

        if (read_mode == SINGLE) {
            for (bam1_t* rec : records) {
                if (!rec || (rec->core.flag & BAM_FUNMAP)) continue;
                add_orphan(rec, MateStatus::SINGLE_END);
            }
            return;
        }

        // Salmon-style pairing: treat each record independently
        // - If record is a proper pair on same transcript, create a paired alignment from this record
        //   using mate fields from BAM (mpos + mate strand)
        // - Otherwise treat as orphan based on BAM_FREAD1/BAM_FREAD2
        // - Keep secondary/supplementary alignments (do not filter them)

        for (bam1_t* rec : records) {
            if (!rec) continue;

            uint32_t flag = rec->core.flag;
            bool read_mapped = !(flag & BAM_FUNMAP);
            if (!read_mapped) continue;

            bool is_paired = (flag & BAM_FPAIRED) != 0;
            bool proper_pair = (flag & BAM_FPROPER_PAIR) != 0;
            int32_t tid = rec->core.tid;
            int32_t mtid = rec->core.mtid;

            if (is_paired && proper_pair && tid >= 0 && mtid >= 0 && tid == mtid) {
                // Only count one mate for proper pairs to avoid duplicate paired alignments
                if (flag & BAM_FREAD2) {
                    continue;
                }
                // Build a paired alignment from this single record using mate fields
                RawAlignment aln;
                aln.transcript_id = static_cast<uint32_t>(tid);
                aln.pos = rec->core.pos;
                aln.score = get_alignment_score(rec);
                aln.log_frag_prob = 0.0;
                aln.log_compat_prob = 0.0;
                aln.est_aln_prob = 1.0;
                aln.mate_status = MateStatus::PAIRED_END_PAIRED;
                aln.is_forward = !(flag & BAM_FREVERSE);
                aln.is_decoy = (tid >= 0 && tid < static_cast<int32_t>(is_decoy.size())) ? is_decoy[tid] : false;
                aln.mate_is_forward = !(flag & BAM_FMREVERSE);
                aln.mate_pos = rec->core.mpos;
                aln.mate_fields_set = true;
                
                // Compute fragment length: prefer TLEN if valid, else compute from mate positions + CIGAR
                int32_t frag_len = -1;
                if (rec->core.isize != 0) {
                    // TLEN is signed: positive for forward orientation, negative for reverse
                    int32_t tlen = rec->core.isize;
                    frag_len = (tlen < 0) ? -tlen : tlen;
                    // Validate TLEN is reasonable
                    if (frag_len < 20 || frag_len > 1000) {
                        frag_len = -1;
                    }
                }
                
                if (frag_len < 0) {
                    // Fallback: compute from mate positions + CIGAR reference length
                    int32_t pos1 = rec->core.pos;
                    int32_t pos2 = rec->core.mpos;
                    bool fwd1 = !(flag & BAM_FREVERSE);
                    bool fwd2 = !(flag & BAM_FMREVERSE);
                    
                    // Must be opposite strands (inward orientation) for valid fragment
                    if (fwd1 != fwd2 && pos1 >= 0 && pos2 >= 0) {
                        uint32_t* cigar = bam_get_cigar(rec);
                        if (cigar) {
                            int32_t len1 = bam_cigar2rlen(rec->core.n_cigar, cigar);
                            // Forward read's leftmost 5' end
                            int32_t p1 = fwd1 ? pos1 : pos2;
                            // Reverse read's rightmost 3' end
                            int32_t p2 = fwd1 ? (pos2 + len1) : (pos1 + len1);
                            frag_len = (p1 > p2) ? (p1 - p2) : (p2 - p1);
                            // Validate computed fragment length
                            if (frag_len < 20 || frag_len > 1000) {
                                frag_len = -1;
                            }
                        }
                    }
                }
                
                aln.fragment_len = frag_len;
                current_read_alignments.push_back(aln);
            } else {
                // Orphan path
                MateStatus ms;
                if (!is_paired) {
                    ms = MateStatus::SINGLE_END;
                } else {
                    bool is_read1 = !(flag & BAM_FREAD2);
                    ms = is_read1 ? MateStatus::PAIRED_END_LEFT : MateStatus::PAIRED_END_RIGHT;
                }
                add_orphan(rec, ms);
            }
        }
    };

    auto finalize_read_group = [&](void) {
        if (current_read_alignments.empty()) {
            return;
        }

        if (alignment_model) {
            if (batch_read_count == 0) {
                fm_calc.getLogMassAndTimestep(log_forgetting_mass, batch_timestep);
            }
            process_read_group_error_model(
                current_read_alignments, current_read_bam_data, transcriptome,
                alignment_model, ec_params, aligner_type, log_forgetting_mass,
                rng, uni, current_qname, deterministic_updates);
            batch_read_count++;
            if (batch_read_count >= mini_batch_size) {
                batch_read_count = 0;
            }
        }
        
        if (alignment_model) {
            alignment_model->incrementObserved();
            
            // Dump matrices at checkpoints (after incrementObserved so count is correct)
            if (!matrix_dump_prefix.empty() && alignment_model) {
                size_t num_obs = alignment_model->numObserved();
                if (num_obs == 1000 || num_obs == 2000 || num_obs == 3000 || num_obs == 4000 || num_obs == 5000) {
                    std::string suffix = "_pre" + std::to_string(num_obs);
                    alignment_model->dumpMatrices(matrix_dump_prefix, suffix);
                }
            }
        }
        
        // Collect observed GC from fragments (matching Salmon's approach)
        // Do this BEFORE freeing BAM records
        if (gc_bias_enabled && observed_gc && !current_read_alignments.empty() && 
            current_read_alignments.size() == current_read_bam_data.size()) {
            // Compute auxProbs to get alignment weights
            ReadMapping mapping = computeAuxProbs(current_read_alignments, ec_params, false);
            
            for (size_t i = 0; i < mapping.alignment_indices.size(); ++i) {
                size_t aln_idx = mapping.alignment_indices[i];
                if (aln_idx >= current_read_alignments.size() || aln_idx >= current_read_bam_data.size()) {
                    continue;
                }
                
                const RawAlignment& aln = current_read_alignments[aln_idx];
                const AlignmentBamData& bam_entry = current_read_bam_data[aln_idx];
                
                // Only collect from properly paired fragments
                if (!bam_entry.paired || aln.transcript_id >= transcriptome.size()) {
                    continue;
                }
                
                const TranscriptSequence* txp = transcriptome.getTranscript(aln.transcript_id);
                if (!txp) continue;
                
                // Get fragment boundaries from BAM records (must be valid here)
                bam1_t* r1 = bam_entry.read1;
                bam1_t* r2 = bam_entry.read2;
                if (!r1 || !r2 || r1->core.tid < 0 || r2->core.tid < 0) continue;
                
                // Compute fragment start and end positions using reference coordinates
                int32_t pos1 = r1->core.pos;
                int32_t pos2 = r2->core.pos;
                // Use reference length from CIGAR (matching Salmon's fragLengthPedantic)
                uint32_t* cigar1 = bam_get_cigar(r1);
                uint32_t* cigar2 = bam_get_cigar(r2);
                if (!cigar1 || !cigar2) continue;
                
                int32_t len1 = bam_cigar2rlen(r1->core.n_cigar, cigar1);
                int32_t len2 = bam_cigar2rlen(r2->core.n_cigar, cigar2);
                
                // Fragment boundaries: leftmost 5' to rightmost 3'
                int32_t frag_start = std::min(pos1, pos2);
                int32_t frag_end = std::max(pos1 + len1, pos2 + len2);
                
                // Clamp to transcript bounds
                int32_t txp_len = static_cast<int32_t>(txp->length());
                if (frag_start < 0) frag_start = 0;
                if (frag_end >= txp_len) frag_end = txp_len - 1;
                if (frag_start > frag_end || frag_start < 0 || frag_end < 0) continue;
                
                // Compute GC fraction
                int32_t gc_pct = txp->gcFrac(frag_start, frag_end);
                
                // Weight by alignment probability (in log space)
                if (i < mapping.aux_probs.size()) {
                    double log_weight = mapping.aux_probs[i];
                    if (log_weight > -std::numeric_limits<double>::infinity()) {
                        observed_gc->inc(gc_pct, log_weight);
                    }
                }
            }
        }
        
        read_alignments.push_back(std::move(current_read_alignments));
        read_qnames.push_back(current_qname);
        current_read_alignments.clear();
        free_alignment_bams(current_read_bam_data);
    };

    // =========================================================================
    // ORDERED PROCESSING MODE: Process reads in order specified by file
    // =========================================================================
    if (!read_order_file.empty()) {
        std::cerr << "Loading read order from: " << read_order_file << "\n";
        
        // Load the read order
        std::vector<std::string> read_order;
        {
            std::ifstream order_file(read_order_file);
            if (!order_file.is_open()) {
                std::cerr << "Error: Cannot open read order file: " << read_order_file << "\n";
                return 1;
            }
            std::string line;
            while (std::getline(order_file, line)) {
                if (!line.empty()) {
                    read_order.push_back(line);
                }
            }
        }
        std::cerr << "Loaded " << read_order.size() << " read names from order file\n";
        
        // First pass: Load all BAM records into a map keyed by qname
        // Each qname maps to a vector of BAM records (for multi-mapping reads)
        std::cerr << "Loading all BAM records into memory...\n";
        std::unordered_map<std::string, std::vector<bam1_t*>> bam_by_qname;
        
        while (true) {
            int ret = sam_read1(bam_file, header, b);
            if (ret < 0) break;
            
            // Skip unmapped
            if (b->core.flag & BAM_FUNMAP) continue;
            
            const char* qname = bam_get_qname(b);
            bam_by_qname[qname].push_back(bam_dup1(b));
        }
        std::cerr << "Loaded " << bam_by_qname.size() << " unique read names\n";
        
        // Auto-detect read mode from first read
        if (read_mode == AUTO && !bam_by_qname.empty()) {
            auto& first_records = bam_by_qname.begin()->second;
            if (!first_records.empty()) {
                read_mode = (first_records[0]->core.flag & BAM_FPAIRED) ? PAIRED : SINGLE;
            }
        }
        
        // Process reads in specified order
        std::cerr << "Processing reads in specified order...\n";
        for (const auto& target_qname : read_order) {
            auto it = bam_by_qname.find(target_qname);
            if (it == bam_by_qname.end()) {
                std::cerr << "Warning: Read " << target_qname << " not found in BAM\n";
                continue;
            }
            
            current_qname = target_qname;
            std::vector<bam1_t*>& records = it->second;
            build_alignments_from_records(records);
            
            finalize_read_group();
            total_reads++;
        }
        
        // Free all buffered BAM records
        for (auto& kv : bam_by_qname) {
            for (bam1_t* rec : kv.second) {
                bam_destroy1(rec);
            }
        }
        
        std::cerr << "Ordered processing complete. Processed " << total_reads << " reads.\n";
        
        // Skip the normal processing loop
        goto post_processing;
    }
    // =========================================================================
    // END ORDERED PROCESSING MODE
    // =========================================================================

    while (true) {
        int ret = sam_read1(bam_file, header, b);
        if (ret < 0) {
            break;  // EOF
        }
        
        if (b->core.flag & BAM_FUNMAP) {
            continue;
        }
        
        const char* qname = bam_get_qname(b);
        bool is_paired = (b->core.flag & BAM_FPAIRED) != 0;
        
        if (read_mode == AUTO) {
            read_mode = is_paired ? PAIRED : SINGLE;
        }
        
        if (current_qname.empty() || strcmp(qname, current_qname.c_str()) != 0) {
            if (!current_qname.empty()) {
                build_alignments_from_records(current_read_records);
                finalize_read_group();
                for (bam1_t* rec : current_read_records) {
                    bam_destroy1(rec);
                }
                current_read_records.clear();
            }
            current_qname = qname;
        }
        
        current_read_records.push_back(bam_dup1(b));
        total_reads++;
    }
    
    if (!current_read_records.empty()) {
        build_alignments_from_records(current_read_records);
        finalize_read_group();
        for (bam1_t* rec : current_read_records) {
            bam_destroy1(rec);
        }
        current_read_records.clear();
    }
    
    // Report error model status
    if (alignment_model) {
        std::cerr << "Error model: observed " << alignment_model->numObserved() << " fragments";
        std::cerr << " (useModel=" << (alignment_model->useModel() ? "yes" : "no");
        std::cerr << ", canUpdate=" << (alignment_model->canUpdate() ? "yes" : "no") << ")";
        std::cerr << "\n";
        
        // Dump final matrices if requested
        if (!matrix_dump_prefix.empty()) {
            alignment_model->dumpMatrices(matrix_dump_prefix, "_final");
            std::cerr << "Dumped final transition matrices to " << matrix_dump_prefix << "_final_*.tsv\n";
        }
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(bam_file);
    
post_processing:
    std::cerr << "Processed " << total_reads << " BAM records\n";
    std::cerr << "Grouped into " << read_alignments.size() << " reads\n";
    // Step 4: Build equivalence classes (alignment-mode: compatibility gating + auxProb only)
    std::cerr << "Building equivalence classes (alignment-mode, no pre-filtering)...\n";
    
    // Sort alignments by transcript_id per read (for consistent ordering)
    for (auto& alignments : read_alignments) {
        if (alignments.empty()) {
            skipped_reads++;
            continue;
        }
        
        // Sort by transcript_id to match Salmon order
        std::sort(alignments.begin(), alignments.end(),
                  [](const RawAlignment& a, const RawAlignment& b) {
                      return a.transcript_id < b.transcript_id;
                  });
    }
    
    // Remove empty read groups while keeping read_qnames aligned
    std::vector<std::vector<RawAlignment>> kept_read_alignments;
    std::vector<std::string> kept_read_qnames;
    kept_read_alignments.reserve(read_alignments.size());
    kept_read_qnames.reserve(read_qnames.size());
    
    for (size_t i = 0; i < read_alignments.size(); ++i) {
        if (!read_alignments[i].empty()) {
            kept_read_alignments.push_back(std::move(read_alignments[i]));
            if (i < read_qnames.size()) {
                kept_read_qnames.push_back(std::move(read_qnames[i]));
            }
        }
    }
    
    // Replace with kept vectors (maintains alignment between read_alignments and read_qnames)
    read_alignments = std::move(kept_read_alignments);
    read_qnames = std::move(kept_read_qnames);
    
    std::cerr << "Building ECs from " << read_alignments.size() << " reads with alignments\n";
    std::cerr << "Skipped " << skipped_reads << " empty read groups\n";
    
    // Sample size warning for error model stability
    if (alignment_model) {
        const size_t MIN_RECOMMENDED_READS = 200000;
        const size_t MIN_WARNING_READS = 50000;
        size_t total_reads = read_alignments.size();
        
        if (total_reads < MIN_WARNING_READS) {
            std::cerr << "\n⚠️  WARNING: Only " << total_reads << " reads for error model training.\n"
                      << "    Error model may be unstable at this sample size.\n"
                      << "    Recommend >= " << MIN_RECOMMENDED_READS << " reads for stable results.\n"
                      << "    Consider --no-error-model for small datasets.\n\n";
        } else if (total_reads < MIN_RECOMMENDED_READS) {
            std::cerr << "NOTE: " << total_reads << " reads is below recommended " 
                      << MIN_RECOMMENDED_READS << " for optimal error model stability.\n";
        }
    }
    
    // Optional: Trace read processing for debugging
    std::ofstream trace_out;
    size_t traced_reads = 0;
    size_t total_dropped_incompat = 0;
    size_t total_orphans = 0;
    size_t total_nonzero_errlike = 0;
    size_t zero_prob_frags = 0;  // Reads where all alignments dropped (sumOfAlignProbs == LOG_0)
    std::ofstream zero_prob_out;  // Output file for zero-prob reads
    
    if (trace_reads) {
        trace_out.open(trace_file);
        if (!trace_out.is_open()) {
            std::cerr << "Warning: Could not open trace file: " << trace_file << std::endl;
        } else {
            trace_out << std::setprecision(17);
            std::cerr << "Tracing reads to: " << trace_file;
            if (trace_limit > 0) {
                std::cerr << " (limit: " << trace_limit << " reads)";
            }
            std::cerr << std::endl;
        }
        
        // Also open file for zero-prob reads
        std::string zero_prob_file = trace_file + ".zero_prob";
        zero_prob_out.open(zero_prob_file);
        if (zero_prob_out.is_open()) {
            std::cerr << "Zero-prob reads will be written to: " << zero_prob_file << std::endl;
        }
    }
    
    // Step 5: Build equivalence classes (no filtering in alignment-mode)
    // Also trace reads if requested
    ECTable ec_table;
    ec_table.n_transcripts = num_transcripts;
    
    std::unordered_map<std::string, size_t> ec_map;
    size_t read_idx = 0;
    
    for (size_t i = 0; i < read_alignments.size(); ++i) {
        const auto& alignments = read_alignments[i];
        if (alignments.empty()) continue;
        
        // Check trace limit before processing
        if (trace_reads && trace_limit > 0 && traced_reads >= trace_limit) {
            break;
        }
        
        // Auto-detect compatibility gating: re-enable compat gating after detection window
        if (auto_detect_mode && !compat_gating_enabled && reads_processed >= AUTO_DETECT_SAMPLE_SIZE) {
            // Restore original values: ignore_incompat = true (gating ON), incompat_prior = 0.0
            ec_params.ignore_incompat = original_ignore_incompat;
            ec_params.incompat_prior = original_incompat_prior;
            compat_gating_enabled = true;
            std::cerr << "Auto-detect complete: compatibility gating re-enabled at read " 
                      << reads_processed << "\n";
        }
        reads_processed++;
        
        // Alignment-mode: no pre-filtering, pass all alignments to computeAuxProbs
        // Compatibility gating and auxProb computation happen inside computeAuxProbs
        const std::vector<RawAlignment>& alignments_to_process = alignments;
        
        // Compute auxProbs for this read (with tracing enabled if requested)
        ReadMapping mapping = computeAuxProbs(alignments_to_process, ec_params, trace_reads);
        
        // Check for zero-probability fragment (all alignments dropped by compatibility gating or auxProb computation)
        // This matches Salmon's check: if (sumOfAlignProbs == LOG_0) { continue; }
        if (mapping.transcript_ids.empty()) {
            zero_prob_frags++;
            if (trace_reads && zero_prob_out.is_open() && i < read_qnames.size()) {
                std::string qname = read_qnames[i];
                zero_prob_out << qname << "\t";
                zero_prob_out << "num_alignments=" << alignments.size() << "\t";
                zero_prob_out << "all_dropped=true";
                zero_prob_out << "\n";
            }
            continue;
        }
        
        // Update summary counters
        total_dropped_incompat += mapping.num_dropped_incompat;
        if (trace_reads && !mapping.trace_info.empty()) {
            for (const auto& trace : mapping.trace_info) {
                if (trace.is_orphan) total_orphans++;
                if (std::abs(trace.err_like) > 1e-10) total_nonzero_errlike++;
            }
        }
        
        // Trace: output structured per-read information
        if (trace_reads && trace_out.is_open() && i < read_qnames.size()) {
            std::string qname = read_qnames[i];
            
            // Build EC label (after range factorization/rank sorting)
            std::vector<uint32_t> final_txp_ids = mapping.transcript_ids;
            std::vector<double> final_weights = mapping.aux_probs;
            
            // Apply rank-based ECs if enabled
            if (ec_params.use_rank_eq_classes && final_txp_ids.size() > 1) {
                applyRankEqClasses(final_txp_ids, final_weights);
            }
            
            // Apply range factorization if enabled
            if (ec_params.use_range_factorization) {
                applyRangeFactorization(final_txp_ids, final_weights, ec_params.range_factorization_bins);
            }
            
            // Build EC label string
            std::string ec_label;
            for (size_t j = 0; j < final_txp_ids.size(); ++j) {
                if (j > 0) ec_label += ",";
                ec_label += std::to_string(final_txp_ids[j]);
            }
            
            // Output structured trace line
            trace_out << qname << "\t";
            trace_out << "txpIDs=";
            for (size_t j = 0; j < mapping.transcript_ids.size(); ++j) {
                if (j > 0) trace_out << ",";
                trace_out << mapping.transcript_ids[j];
            }
            trace_out << ";";
            
            trace_out << "as=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << mapping.trace_info[j].as_tag;
                }
            }
            trace_out << ";";
            
            trace_out << "bestAS=" << mapping.best_as << ";";
            
            trace_out << "logFragProb=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << mapping.trace_info[j].log_frag_prob;
                }
            }
            trace_out << ";";
            
            trace_out << "errLike=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << mapping.trace_info[j].err_like;
                }
            }
            trace_out << ";";
            
            trace_out << "logCompat=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << mapping.trace_info[j].log_compat_prob;
                }
            }
            trace_out << ";";
            
            trace_out << "droppedIncompat=" << mapping.num_dropped_incompat << ";";
            
            trace_out << "orphan=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << (mapping.trace_info[j].is_orphan ? "1" : "0");
                }
            }
            trace_out << ";";
            
            trace_out << "auxProb=";
            for (size_t j = 0; j < mapping.trace_info.size(); ++j) {
                if (j > 0) trace_out << ",";
                if (!mapping.trace_info[j].dropped_incompat && 
                    mapping.trace_info[j].aux_prob != LOG_0) {
                    trace_out << mapping.trace_info[j].aux_prob;
                }
            }
            trace_out << ";";
            
            trace_out << "weight=";
            for (size_t j = 0; j < mapping.aux_probs.size(); ++j) {
                if (j > 0) trace_out << ",";
                trace_out << mapping.aux_probs[j];
            }
            trace_out << ";";
            
            trace_out << "ec_label=" << ec_label;
            trace_out << "\n";
            
            traced_reads++;
        }
        
        // Check trace limit after processing this read
        // If limit reached, break out of main loop (stop processing reads)
        if (trace_reads && trace_limit > 0 && traced_reads >= trace_limit) {
            break;
        }
        
        // Note: rank-based ECs and range factorization are applied in trace output above
        // Apply rank-based ECs if enabled (for EC building)
        if (ec_params.use_rank_eq_classes && mapping.transcript_ids.size() > 1) {
            applyRankEqClasses(mapping.transcript_ids, mapping.aux_probs);
        }
        
        // Apply range factorization if enabled (for EC building)
        if (ec_params.use_range_factorization) {
            applyRangeFactorization(mapping.transcript_ids, mapping.aux_probs, ec_params.range_factorization_bins);
        }
        
        // Create canonical key for EC lookup (order-sensitive)
        std::string ec_key;
        for (size_t j = 0; j < mapping.transcript_ids.size(); ++j) {
            if (j > 0) ec_key += ",";
            ec_key += std::to_string(mapping.transcript_ids[j]);
        }
        
        // Find or create EC
        auto it = ec_map.find(ec_key);
        if (it == ec_map.end()) {
            // Create new EC
            EC ec;
            ec.transcript_ids = mapping.transcript_ids;
            ec.weights = mapping.aux_probs;
            ec.count = 1.0;
            ec_table.ecs.push_back(ec);
            ec_map[ec_key] = ec_table.ecs.size() - 1;
        } else {
            // Increment count and accumulate weights
            EC& ec = ec_table.ecs[it->second];
            ec.count += 1.0;
            if (ec.weights.size() == mapping.aux_probs.size()) {
                for (size_t j = 0; j < ec.weights.size(); ++j) {
                    ec.weights[j] += mapping.aux_probs[j];
                }
            }
        }
        
        read_idx++;
    }
    
    ec_table.n_ecs = ec_table.ecs.size();
    
    if (trace_reads && trace_out.is_open()) {
        trace_out.close();
        std::cerr << "Traced " << traced_reads << " reads to " << trace_file << std::endl;
        std::cerr << "Trace summary:" << std::endl;
        std::cerr << "  Dropped incompatible: " << total_dropped_incompat << std::endl;
        std::cerr << "  Orphan alignments: " << total_orphans << std::endl;
        std::cerr << "  Non-zero errLike: " << total_nonzero_errlike << std::endl;
        std::cerr << "  Zero-probability fragments: " << zero_prob_frags << std::endl;
    }
    if (zero_prob_out.is_open()) {
        zero_prob_out.close();
    }
    
    std::cerr << "Built " << ec_table.ecs.size() << " equivalence classes\n";
    std::cerr << "Zero-probability fragments (all alignments dropped): " << zero_prob_frags << "\n";
    
    // Step 6: Apply extended pruning if enabled
    if (!no_local_pruning || !no_global_pruning) {
        std::cerr << "Applying extended pruning...\n";
        applyExtendedPruning(ec_table, pruning_params);
        std::cerr << "After pruning: " << ec_table.ecs.size() << " equivalence classes\n";
    }
    
    // Step 7: Write Salmon format eq_classes.txt
    std::cerr << "Writing EC output to: " << output_file << "\n";
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Failed to open output file: " << output_file << std::endl;
        return 1;
    }
    
    // Set precision for weights
    out << std::setprecision(17);
    
    // Write header
    out << num_transcripts << "\n";
    out << ec_table.ecs.size() << "\n";
    
    // Write transcript names
    for (const auto& name : transcript_names) {
        out << name << "\n";
    }
    
    // Write ECs (preserve order, do NOT sort)
    for (const auto& ec : ec_table.ecs) {
        size_t k_labels = ec.transcript_ids.size();
        size_t k_weights = ec.weights.size();
        bool has_weights = (k_weights > 0);
        size_t group_size = has_weights ? k_weights : k_labels;
        if (group_size > k_labels) {
            group_size = k_labels;
        }
        out << group_size;
        
        // Write transcript IDs (Salmon writes only the first group_size IDs)
        for (size_t i = 0; i < group_size; ++i) {
            out << " " << ec.transcript_ids[i];
        }
        
        // Write weights if available (normalized by count to get average per-read weights)
        if (has_weights && ec.count > 0) {
            size_t w_count = std::min(group_size, k_weights);
            for (size_t i = 0; i < w_count; ++i) {
                // Normalize: accumulated_weight / count = average weight per read
                out << " " << (ec.weights[i] / ec.count);
            }
        }
        
        // Write count
        out << " " << ec.count << "\n";
    }
    
    out.close();
    std::cerr << "Done. Wrote " << ec_table.ecs.size() << " equivalence classes.\n";
    
    // Write observed GC distribution if collected
    if (gc_bias_enabled && observed_gc && !observed_gc_out_file.empty()) {
        std::cerr << "Writing observed GC distribution to: " << observed_gc_out_file << "\n";
        if (!observed_gc->writeToFile(observed_gc_out_file)) {
            std::cerr << "Warning: Failed to write observed GC distribution\n";
        } else {
            std::cerr << "Done. Wrote observed GC distribution.\n";
        }
    }
    
    // Cleanup
    if (alignment_model) {
        delete alignment_model;
    }
    
    if (observed_gc) {
        delete observed_gc;
    }
    
    return 0;
}
