#ifndef CODE_TranscriptQuantEC
#define CODE_TranscriptQuantEC

#include "IncludeDefine.h"
#include "Transcript.h"
#include "em_types.h"
#include "ec_builder.h"
#include "gc_bias.h"
#include "fld_accumulator.h"
#include "Parameters.h"
#include "LibFormatDetection.h"
#include "alignment_model.h"  // For AlignmentModel and Transcriptome
#include "forgetting_mass.h"  // For ForgettingMassCalculator
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <memory>
#include <random>

// Wrapper class for building equivalence classes during STAR alignment
// Maintains thread-local EC table and GC observations
class TranscriptQuantEC {
public:
    TranscriptQuantEC(uint num_transcripts, int threadId, const std::string& traceFile, int traceLimit, Parameters& P,
                      libem::AlignmentModel* alignmentModel = nullptr,
                      const libem::Transcriptome* transcriptome = nullptr);
    
    // Add alignments for a read (called after all alignments for a read are computed)
    // trMult: array of transcript pointers (alignments)
    // nTr: number of alignments
    // auxProbs: alignment probabilities (log-space, normalized)
    // qname: read name (for tracing)
    // packed_read1: packed read1 sequence (4-bit BAM encoding, nullptr if not available)
    // packed_read1_rc: packed read1 reverse-complement sequence (nullptr if not available)
    // read1_len: read1 length (0 if not available)
    // packed_read2: packed read2 sequence (4-bit BAM encoding, nullptr if not available)
    // packed_read2_rc: packed read2 reverse-complement sequence (nullptr if not available)
    // read2_len: read2 length (0 if not available)
    void addReadAlignments(Transcript* trMult, uint nTr, const std::vector<double>& auxProbs, const char* qname,
                          const uint8_t* packed_read1 = nullptr, const uint8_t* packed_read1_rc = nullptr, uint32_t read1_len = 0,
                          const uint8_t* packed_read2 = nullptr, const uint8_t* packed_read2_rc = nullptr, uint32_t read2_len = 0);
    
    // Get EC builder parameters (based on detection state)
    ECBuilderParams getParams() const;
    
    // Simplified method: add alignments using just transcript IDs and weights
    // transcriptIds: vector of transcript IDs this read maps to
    // weights: vector of alignment weights (uniform or from alignment scores)
    void addReadAlignmentsSimple(const std::vector<uint32_t>& transcriptIds, const std::vector<double>& weights);
    
    // Add GC observation from a properly-paired fragment
    // frag_start, frag_end: fragment boundaries in transcript coordinates
    // gc_pct: GC percentage (0-100)
    // weight: alignment probability weight (log-space)
    void addGCObservation(int32_t frag_start, int32_t frag_end, int32_t gc_pct, double weight);
    
    // Add FLD observation from a properly-paired fragment
    // frag_len: fragment length
    // weight: alignment probability weight (linear space, typically 1.0)
    // log_prob: log probability for stochastic acceptance (if <= LOG_0, always accept)
    // Only updates during burn-in phase (pre-burnin=0, burn-in=PMF, post-burnin=conditional PMF)
    void addFLDObservation(int32_t frag_len, double weight = 1.0, double log_prob = LOG_0);
    
    // Finalize EC table (called after all reads processed)
    void finalize();
    
    // Merge EC table from another thread
    void merge(const TranscriptQuantEC& other);
    
    // Get EC table (for quantification)
    const ECTable& getECTable() const { return ecTable_; }
    
    // Get observed GC model (for GC bias correction)
    const GCFragModel& getObservedGC() const { return observedGC_; }
    GCFragModel& getObservedGC() { return observedGC_; }
    
    // Get observed FLD (for fragment length distribution)
    const FLDAccumulator& getObservedFLD() const { return observedFLD_; }
    FLDAccumulator& getObservedFLD() { return observedFLD_; }
    
    // Set transcript lengths (called from ReadAlign.cpp constructor and STAR.cpp after merge)
    void setTranscriptLengths(const std::vector<int32_t>& lengths);
    
    // Get drop statistics (for logging)
    size_t getDroppedIncompat() const { return dropped_incompat_; }
    size_t getDroppedMissingMateFields() const { return dropped_missing_mate_fields_; }
    size_t getDroppedUnknownObsFmt() const { return dropped_unknown_obs_fmt_; }
    
private:
    ECTable ecTable_;
    GCFragModel observedGC_;
    FLDAccumulator observedFLD_;  // Fragment length distribution accumulator
    mutable std::vector<double> cached_fld_pmf_;  // Cached PMF (computed on demand)
    mutable bool fld_dirty_ = true;  // Flag to track if FLD PMF needs recomputation
    mutable std::vector<int32_t> cached_transcript_lengths_;  // Cached transcript lengths (populated on demand)
    uint64_t num_processed_fragments_ = 0;  // Track processed fragments for burn-in gating
    size_t local_batch_reads_ = 0;  // Thread-local batch counter (increments per read, flushed at batch size)
    uint64_t batch_start_count_ = 0;  // Batch start count (persisted across reads in a batch, only updated at batch boundaries)
    ForgettingMassCalculator forgetting_mass_;  // Forgetting mass calculator for FLD updates
    ForgettingMassCalculator error_model_forgetting_mass_;  // Error model forgetting mass (minibatches)
    double cached_log_forgetting_mass_ = 0.0;  // Cached log forgetting mass for current mini-batch
    uint64_t cached_forgetting_batch_ = UINT64_MAX;  // Mini-batch number for cached forgetting mass (sentinel = never cached)
    size_t error_model_batch_count_ = 0;
    double error_model_log_mass_ = 0.0;
    std::mt19937 rng_;  // Random number generator for stochastic acceptance
    uint num_transcripts_;
    
    // Current read sequences (packed, 4-bit BAM encoding) - set per read in addReadAlignments()
    const uint8_t* current_packed_read1_ = nullptr;
    const uint8_t* current_packed_read1_rc_ = nullptr;
    uint32_t current_read1_len_ = 0;
    const uint8_t* current_packed_read2_ = nullptr;
    const uint8_t* current_packed_read2_rc_ = nullptr;
    uint32_t current_read2_len_ = 0;
    
    // Helper: Extract CIGAR from Transcript and convert to BAM-packed format
    std::vector<uint32_t> extractCIGAR(Transcript* tr, bool is_read1) const;
    
    // Helper: Parse CIGAR string into BAM-packed format
    std::vector<uint32_t> parseCIGARString(const std::string& cigar_str) const;
    int threadId_;
    int trace_limit_;
    int traced_count_;
    std::ofstream trace_out_;
    Parameters& P_;  // Reference to Parameters for accessing detection state
    libem::AlignmentModel* alignment_model_;  // Error model (optional, null if not available)
    const libem::Transcriptome* transcriptome_;  // Transcript sequences cache (optional, null if not available)
    
    // Temporary storage for current read's alignments
    std::vector<RawAlignment> current_read_alignments_;
    
    // Drop statistics (aggregated across all reads)
    size_t dropped_incompat_;
    size_t dropped_missing_mate_fields_;
    size_t dropped_unknown_obs_fmt_;
    
    // Convert STAR Transcript to RawAlignment for EC building
    RawAlignment transcriptToRawAlignment(Transcript* tr, uint32_t transcript_id, uint32_t read1_len, uint32_t read2_len);
    
    // Write trace line for a read
    void writeTraceLine(const char* qname, const ReadMapping& mapping, const ECBuilderParams& params, 
                        const std::vector<RawAlignment>& alignments);
    
    // Hash function for EC signature (sorted transcript IDs only - matches Salmon)
    // NOTE: Salmon keys ECs by transcript IDs only, not by weights
    // This ensures ECs with the same transcripts are merged even if weights differ slightly
    struct ECSignature {
        std::vector<uint32_t> transcript_ids;
        
        bool operator==(const ECSignature& other) const {
            if (transcript_ids.size() != other.transcript_ids.size()) return false;
            for (size_t i = 0; i < transcript_ids.size(); ++i) {
                if (transcript_ids[i] != other.transcript_ids[i]) return false;
            }
            return true;
        }
    };
    
    struct ECSignatureHash {
        size_t operator()(const ECSignature& sig) const {
            size_t h = 0;
            for (uint32_t id : sig.transcript_ids) {
                h ^= std::hash<uint32_t>{}(id) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
    
    // Map from EC signature to EC index
    std::unordered_map<ECSignature, size_t, ECSignatureHash> ec_signature_map_;
};

#endif // CODE_TranscriptQuantEC
