#include "TranscriptQuantEC.h"
#include "Transcriptome.h"
#include "ec_builder.h"
#include "Parameters.h"
#include "LibFormatDetection.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <random>
#include <atomic>

// Global atomic counter for processed read groups (matches Salmon's processedReads)
// Used for pre-burn-in gating: aux params are enabled when this count >= numPreBurninFrags (5000)
extern std::atomic<uint64_t> global_processed_fragments;

TranscriptQuantEC::TranscriptQuantEC(uint num_transcripts, int threadId, const std::string& traceFile, int traceLimit, Parameters& P,
                                      libem::AlignmentModel* alignmentModel, const libem::Transcriptome* transcriptome)
    : num_transcripts_(num_transcripts), threadId_(threadId), trace_limit_(traceLimit), traced_count_(0), P_(P),
      alignment_model_(alignmentModel), transcriptome_(transcriptome),
      forgetting_mass_(0.65), error_model_forgetting_mass_(0.65),
      rng_(static_cast<unsigned>(threadId + 1)),  // Seed RNG with thread ID
      dropped_incompat_(0), dropped_missing_mate_fields_(0), dropped_unknown_obs_fmt_(0),
      local_batch_reads_(0), batch_start_count_(0) {  // Initialize batch counter and batch start count
    ecTable_.n_transcripts = static_cast<size_t>(num_transcripts);
    ecTable_.n_ecs = 0;
    
    // Open trace file if specified
    if (!traceFile.empty()) {
        std::string fname = traceFile + "." + std::to_string(threadId_) + ".tsv";
        trace_out_.open(fname);
        if (trace_out_.is_open()) {
            // Write header (trace format matches Salmon/CLI alignment-mode trace)
            // Added instrumentation: fragLen, useAuxParams, globalFragCount, batchReads, batchStartCount, inDetectionMode
            trace_out_ << "#qname\ttxpIDs=...;as=...;bestAS=...;logFragProb=...;errLike=...;"
                       << "logCompat=...;expFmt=...;obsFmt=...;isCompat=...;mateStatus=...;mateFields=...;"
                       << "fwd=...;mateFwd=...;pos=...;matePos=...;primary=...;dropped=...;orphan=...;auxProb=...;"
                       << "fragLen=...;useAuxParams=...;globalFragCount=...;batchReads=...;batchStartCount=...;inDetectionMode=...;\n";
        }
    }
}

RawAlignment TranscriptQuantEC::transcriptToRawAlignment(Transcript* tr, uint32_t transcript_id, uint32_t read1_len, uint32_t read2_len) {
    RawAlignment aln;
    aln.transcript_id = transcript_id;
    aln.score = tr->maxScore;
    aln.est_aln_prob = 1.0;  // Will be normalized later
    aln.log_frag_prob = 0.0;  // Not using fragment length dist for now
    aln.log_compat_prob = 0.0;  // Compatibility already checked
    aln.err_like = 0.0;
    aln.has_err_like = false;
    aln.is_decoy = false;  // TODO: check if transcript is decoy

    // Identify read1/read2 exon indices (if present) for proper orientation
    uint iExRead1 = UINT_MAX;
    uint iExRead2 = UINT_MAX;
    for (uint i = 0; i < tr->nExons; i++) {
        if (tr->exons[i][EX_iFrag] == 0 && iExRead1 == UINT_MAX) {
            iExRead1 = i;
        }
        if (tr->exons[i][EX_iFrag] == 1 && iExRead2 == UINT_MAX) {
            iExRead2 = i;
        }
    }

    // Default to leftmost exon if read1/read2 not found
    uint iExPos = (iExRead1 != UINT_MAX) ? iExRead1 : 0;
    aln.pos = tr->exons[iExPos][EX_G];  // TRANSCRIPTOME coordinate

    // Determine read1 forwardness based on which mate is leftmost
    // Match detection logic: tr->Str is the strand of the leftmost mate
    // (see LibFormatDetection.cpp observeFormatFromTranscript)
    bool read1_fwd;
    bool read2_fwd;
    
    if (tr->readNmates == 2 && iExRead1 != UINT_MAX && iExRead2 != UINT_MAX) {
        // Determine which mate is leftmost based on positions
        int32_t pos1 = tr->exons[iExRead1][EX_G];
        int32_t pos2 = tr->exons[iExRead2][EX_G];
        bool read1_leftmost = (pos1 <= pos2);
        
        if (read1_leftmost) {
            // Read1 is leftmost: read1 strand = Str
            read1_fwd = (tr->Str == 0);
            // Assume opposite strands (as Str does for proper pairs)
            read2_fwd = !read1_fwd;
        } else {
            // Read2 is leftmost: read2 strand = Str
            read2_fwd = (tr->Str == 0);
            // Assume opposite strands
            read1_fwd = !read2_fwd;
        }
    } else {
        // Single-end or missing mate info: use Str directly
        read1_fwd = (tr->Str == 0);
        read2_fwd = !read1_fwd;  // Not used for single-end
    }
    
    aln.is_forward = read1_fwd;
    
    if (tr->readNmates == 2) {
        // Find mate boundary using canonSJ == -3
        uint iExMate = 0;
        bool foundMate = false;
        for (uint i = 0; i < tr->nExons - 1; i++) {
            if (tr->canonSJ[i] == -3) {
                iExMate = i;
                foundMate = true;
                break;
            }
        }
        
        bool proper = foundMate && 
            (tr->exons[0][EX_iFrag] != tr->exons[tr->nExons-1][EX_iFrag]);
        
        if (proper) {
            aln.mate_status = MateStatus::PAIRED_END_PAIRED;
            if (iExRead2 != UINT_MAX) {
                aln.mate_pos = tr->exons[iExRead2][EX_G];  // TRANSCRIPTOME coord
            } else {
                aln.mate_pos = tr->exons[iExMate + 1][EX_G];  // fallback
            }
            aln.mate_is_forward = read2_fwd;  // Use computed read2_fwd
            // Set mate fields for proper pairs (regardless of primaryFlag)
            // Only skip compat gating for orphans or missing mates
            aln.mate_fields_set = true;
            
            // Compute fragment length using Salmon's pedantic definition (ReadPair::fragLengthPedantic)
            // Only for proper pairs and opposite orientation (inward/outward)
            // Formula: rightmost mate 3' end - leftmost mate 5' start + 1
            // For same-orientation mates: set fragment_len = -1 (invalid, skip FLD/logFragProb)
            int32_t fragLen = -1;
            if (read1_fwd != read2_fwd && read1_len > 0 && read2_len > 0) {
                // Opposite strands: compute pedantic fragment length
                // Get transcript length from cached lengths (for clamping)
                int32_t transcript_len = 0;
                if (!cached_transcript_lengths_.empty() && transcript_id < cached_transcript_lengths_.size()) {
                    transcript_len = cached_transcript_lengths_[transcript_id];
                } else {
                    // Fallback: compute from exons
                    for (uint32_t iex = 0; iex < tr->nExons; iex++) {
                        int32_t ex_end = tr->exons[iex][EX_G] + tr->exons[iex][EX_L];
                        if (ex_end > transcript_len) {
                            transcript_len = ex_end;
                        }
                    }
                }
                
                // Get positions (transcriptomic coordinates)
                int32_t pos1 = aln.pos;  // Read1 5' start (transcriptomic)
                int32_t pos2 = aln.mate_pos;  // Read2 5' start (transcriptomic)
                
                // Determine 5' start and 3' end for each mate
                // Forward mate: 5' start = pos, 3' end = pos + len - 1
                // Reverse mate: 5' start = pos + len - 1, 3' end = pos
                int32_t read1_5prime, read1_3prime, read2_5prime, read2_3prime;
                if (read1_fwd) {
                    read1_5prime = pos1;
                    read1_3prime = pos1 + static_cast<int32_t>(read1_len) - 1;
                } else {
                    read1_5prime = pos1 + static_cast<int32_t>(read1_len) - 1;
                    read1_3prime = pos1;
                }
                
                if (read2_fwd) {
                    read2_5prime = pos2;
                    read2_3prime = pos2 + static_cast<int32_t>(read2_len) - 1;
                } else {
                    read2_5prime = pos2 + static_cast<int32_t>(read2_len) - 1;
                    read2_3prime = pos2;
                }
                
                // Find leftmost 5' start and rightmost 3' end
                int32_t leftmost_5prime = std::min(read1_5prime, read2_5prime);
                int32_t rightmost_3prime = std::max(read1_3prime, read2_3prime);
                
                // Clamp to transcript length if needed
                if (transcript_len > 0) {
                    leftmost_5prime = std::max(0, std::min(leftmost_5prime, transcript_len - 1));
                    rightmost_3prime = std::max(0, std::min(rightmost_3prime, transcript_len - 1));
                }
                
                // Pedantic fragment length: rightmost 3' end - leftmost 5' start + 1
                fragLen = rightmost_3prime - leftmost_5prime + 1;
                
                // Validate fragment length (must be positive and within reasonable bounds)
                if (fragLen <= 0 || fragLen > FLDAccumulator::MAX_FRAG_LEN) {
                    fragLen = -1;  // Invalid
                }
            } else {
                // Same orientation or missing read lengths: Salmon's fragLengthPedantic returns -1 (invalid)
                fragLen = -1;
            }
            aln.fragment_len = fragLen;
        } else {
            // Orphan - DO NOT set mate fields
            aln.mate_status = (tr->exons[0][EX_iFrag] == 0)
                ? MateStatus::PAIRED_END_LEFT : MateStatus::PAIRED_END_RIGHT;
            aln.mate_fields_set = false;
        }
    } else {
        aln.mate_status = MateStatus::SINGLE_END;
        aln.mate_fields_set = false;
        aln.fragment_len = -1;  // Single-end: no fragment length
    }
    
    return aln;
}

ECBuilderParams TranscriptQuantEC::getParams() const {
    ECBuilderParams params;
    // Range factorization: Enabled by default to match Salmon (default bins=4)
    // This splits multi-mapper ECs into bins based on alignment positions
    // For strict parity with Salmon, both tools now use range factorization
    params.use_range_factorization = true;
    params.range_factorization_bins = 4;  // Salmon default
    params.use_rank_eq_classes = false;
    params.use_frag_len_dist = false;  // Will be set based on Parameters
    params.use_error_model = false;  // Will be set based on errorModelMode below
    params.use_as_without_cigar = false;
    
    // Set error model mode based on Parameters
    // auto: CIGAR-only (no AS fallback), cigar: CIGAR only, as: AS only, off: disabled
    // NOTE: Error model requires:
    //   1. CIGAR strings in alignments (extracted via extractCIGAR())
    //   2. Burn-in threshold: Error model becomes active after PRE_BURNIN_FRAGS=5000 observations
    //   3. AlignmentModel must be created (alignment_model_ != nullptr)
    //   4. Transcriptome sequences must be available (transcriptome_ != nullptr)
    // After burn-in, err_like is computed from CIGAR-based error model if has_err_like=true
    // If prerequisites are missing in auto mode, errLike is disabled (no AS fallback)
    string errorModelMode = P_.quant.transcriptVB.errorModelMode;
    if (errorModelMode == "auto") {
        // Auto mode: CIGAR-only (no AS fallback)
        // If prerequisites are missing, errLike will be disabled and a warning logged
        params.use_error_model = true;
        params.use_as_without_cigar = false;
    } else if (errorModelMode == "cigar") {
        params.use_error_model = true;
        params.use_as_without_cigar = false;
    } else if (errorModelMode == "as") {
        params.use_error_model = false;
        params.use_as_without_cigar = true;
    } else {  // "off"
        params.use_error_model = false;
        params.use_as_without_cigar = false;
    }
    
    // Wire FLD into params if observations exist
    // Prefer O(1) log-space accumulator over PMF to avoid recomputation
    bool has_fld = observedFLD_.isValid();
    if (has_fld) {
        // Prefer accumulator for O(1) logFragProb computation
        params.fld_accumulator = &observedFLD_;
        // Only recompute PMF when needed (for conditional PMF in post-burnin phase)
        // Since fld_accumulator is preferred, we don't need to recompute PMF on every observation
        // Only recompute when dirty AND conditional PMF is needed (post-burnin)
        bool need_pmf = params.use_conditional_pmf && 
                        params.num_processed_fragments >= params.num_burnin_frags;
        if (need_pmf && (fld_dirty_ || cached_fld_pmf_.empty())) {
            cached_fld_pmf_ = observedFLD_.getPMF();
            fld_dirty_ = false;
        }
        // Provide PMF only if computed (may be empty if not needed)
        if (!cached_fld_pmf_.empty()) {
            params.fld_pmf = &cached_fld_pmf_;
        } else {
            params.fld_pmf = nullptr;
        }
        params.use_frag_len_dist = true;  // Enable FLD when observations exist
    } else {
        params.fld_accumulator = nullptr;
        params.fld_pmf = nullptr;
        params.use_frag_len_dist = false;
    }
    
    // Wire transcript lengths and start position probability
    bool has_transcript_lengths = !cached_transcript_lengths_.empty();
    if (has_transcript_lengths) {
        params.transcript_lengths = &cached_transcript_lengths_;
        // Enable start position probability when FLD + fragment lengths + transcript lengths are available
        params.use_start_pos_prob = has_fld;  // Only if FLD is also enabled
    } else {
        params.transcript_lengths = nullptr;
        params.use_start_pos_prob = false;
    }
    
    // Wire alignment filtering parameters.
    // Salmon uses AS-based score filtering only when using AS (RapMap/Pufferfish).
    if (params.use_as_without_cigar) {
        params.min_score_fraction = 0.65;  // Salmon default when AS is used
        params.apply_filtering = true;
    } else {
        params.min_score_fraction = 0.0;
        params.apply_filtering = false;
    }
    
    // Wire pre-burn-in fragment count threshold from Parameters
    params.num_pre_burnin_frags = static_cast<uint64_t>(P_.quant.transcriptVB.preBurninFrags);
    
    // Check detection mode flag (set in Parameters)
    if (P_.quant.transcriptVB.inDetectionMode) {
        // DETECTION MODE: gating OFF
        params.lib_format = LibraryFormat::IU();
        params.ignore_incompat = false;  // Don't skip ANY alignments
        params.incompat_prior = 1.0;     // No penalty
    } else if (P_.quant.transcriptVB.detectionComplete) {
        // POST-DETECTION: gating ON with detected format
        // Use detected format when libType=A (auto-detect), otherwise use explicit format
        params.lib_format = LibraryFormat::formatFromID(
            P_.quant.transcriptVB.detectedLibFormatId);
        params.incompat_prior = 0.0;
        // Match Salmon: drop incompatible when ignoreIncompat is enabled
        params.ignore_incompat = true;
    } else {
        // Explicit format specified (not auto-detect)
        params.lib_format = parseLibFormat(P_.quant.transcriptVB.libType);
        params.incompat_prior = 0.0;
        // Match Salmon: drop incompatible when ignoreIncompat is enabled
        params.ignore_incompat = true;
    }
    
    return params;
}

void TranscriptQuantEC::addReadAlignments(Transcript* trMult, uint nTr, const std::vector<double>& auxProbs, const char* qname,
                                          const uint8_t* packed_read1, const uint8_t* packed_read1_rc, uint32_t read1_len,
                                          const uint8_t* packed_read2, const uint8_t* packed_read2_rc, uint32_t read2_len) {
    if (nTr == 0) return;
    
    // Store packed read sequences for error model computation (per-read, not per-alignment)
    current_packed_read1_ = packed_read1;
    current_packed_read1_rc_ = packed_read1_rc;
    current_read1_len_ = read1_len;
    current_packed_read2_ = packed_read2;
    current_packed_read2_rc_ = packed_read2_rc;
    current_read2_len_ = read2_len;
    
    // >>> VOTE ON FORMAT DURING DETECTION MODE <<<
    if (P_.quant.transcriptVB.inDetectionMode && 
        P_.quant.transcriptVB.libFormatDetector != nullptr) {
        
        // Prefer primary alignments, but do not require them for detection.
        // primaryFlag is only set when TranscriptomeSAM output is enabled.
        Transcript* vote_tr = nullptr;
        for (uint i = 0; i < nTr; ++i) {
            Transcript* tr = &trMult[i];
            if (tr->readNmates != 2) {
                continue;
            }
            // Find mate boundary
            bool proper = false;
            for (uint j = 0; j < tr->nExons - 1; j++) {
                if (tr->canonSJ[j] == -3) {
                    proper = (tr->exons[0][EX_iFrag] != tr->exons[tr->nExons-1][EX_iFrag]);
                    break;
                }
            }
            if (!proper) {
                continue;
            }
            if (tr->primaryFlag) {
                vote_tr = tr;
                break;  // Use first primary proper pair if present
            }
            if (vote_tr == nullptr) {
                vote_tr = tr;  // Fallback to first proper pair
            }
        }
        if (vote_tr != nullptr) {
            P_.quant.transcriptVB.libFormatDetector->vote(vote_tr);
        }
    }
    
    // Get params (gating OFF during detection, ON after)
    ECBuilderParams params = getParams();
    
    const bool error_model_available = params.use_error_model &&
        alignment_model_ != nullptr &&
        transcriptome_ != nullptr &&
        current_packed_read1_ != nullptr &&
        current_read1_len_ > 0;
    
    // Warn if auto mode is requested but prerequisites are missing (log once per thread)
    static thread_local bool warned_auto_prerequisites = false;
    if (params.use_error_model && !error_model_available && !warned_auto_prerequisites) {
        string errorModelMode = P_.quant.transcriptVB.errorModelMode;
        if (errorModelMode == "auto") {
            P_.inOut->logMain << "WARNING: Auto error model requested but transcriptome/CIGAR unavailable; "
                              << "errLike disabled (no AS fallback)." << "\n";
            warned_auto_prerequisites = true;
        }
    }
    
    const bool use_error_model = error_model_available;
    const bool use_model = use_error_model && alignment_model_->useModel();
    
    struct ErrModelData {
        const libem::TranscriptSequence* ref = nullptr;
        std::vector<uint32_t> cigar1;
        std::vector<uint32_t> cigar2;
        int32_t pos1 = -1;
        int32_t pos2 = -1;
        bool paired = false;
    };
    std::vector<ErrModelData> err_data;
    err_data.reserve(nTr);
    
    current_read_alignments_.clear();
    current_read_alignments_.reserve(nTr);
    
    // Convert Transcript objects to RawAlignment
    for (uint i = 0; i < nTr; ++i) {
        Transcript* tr = &trMult[i];
        if (tr->nExons == 0) {
            continue;
        }
        // Mirror TranscriptomeSAM filtering for EC input.
        if (!P_.quant.trSAM.indel && (tr->nDel > 0 || tr->nIns > 0)) {
            continue;
        }
        if (!P_.quant.trSAM.singleEnd && tr->readNmates == 2) {
            if (tr->exons[0][EX_iFrag] == tr->exons[tr->nExons - 1][EX_iFrag]) {
                continue;
            }
        }
        uint transcript_id = tr->Chr;  // Chr field contains transcript index
        
        if (transcript_id >= num_transcripts_) {
            continue;  // Skip invalid transcript IDs
        }
        
        RawAlignment aln = transcriptToRawAlignment(tr, static_cast<uint32_t>(transcript_id), 
                                                     current_read1_len_, current_read2_len_);
        aln.is_primary = tr->primaryFlag;
        
        // Set alignment probability if provided
        if (i < auxProbs.size()) {
            aln.est_aln_prob = auxProbs[i];
        }
        
        ErrModelData err;
        if (use_error_model) {
            err.ref = transcriptome_->getTranscript(transcript_id);
            err.pos1 = aln.pos;
            err.pos2 = aln.mate_pos;
            err.paired = (aln.mate_status == MateStatus::PAIRED_END_PAIRED &&
                          current_packed_read2_ != nullptr &&
                          current_read2_len_ > 0);
            err.cigar1 = extractCIGAR(tr, true);
            if (err.paired) {
                err.cigar2 = extractCIGAR(tr, false);
            }
            
            if (use_model && err.ref != nullptr && !err.cigar1.empty()) {
                const uint8_t* read1_seq = current_packed_read1_;
                if (!aln.is_forward && current_packed_read1_rc_ != nullptr) {
                    read1_seq = current_packed_read1_rc_;
                }
                const uint8_t* read2_seq = current_packed_read2_;
                if (err.paired && !aln.mate_is_forward && current_packed_read2_rc_ != nullptr) {
                    read2_seq = current_packed_read2_rc_;
                }
                if (err.paired && !err.cigar2.empty() &&
                    current_packed_read2_ != nullptr && current_read2_len_ > 0) {
                    aln.err_like = alignment_model_->logLikelihoodPaired(
                        err.cigar1.data(), static_cast<uint32_t>(err.cigar1.size()),
                        read1_seq, static_cast<int32_t>(current_read1_len_),
                        err.cigar2.data(), static_cast<uint32_t>(err.cigar2.size()),
                        read2_seq, static_cast<int32_t>(current_read2_len_),
                        err.ref, err.pos1, err.pos2);
                } else {
                    bool is_left = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                                    aln.mate_status == MateStatus::SINGLE_END);
                    aln.err_like = alignment_model_->logLikelihood(
                        err.cigar1.data(), static_cast<uint32_t>(err.cigar1.size()),
                        read1_seq, static_cast<int32_t>(current_read1_len_),
                        err.ref, err.pos1, is_left);
                }
                aln.has_err_like = true;
            } else {
                // CIGAR missing or other prerequisites not met
                aln.has_err_like = false;
                aln.err_like = 0.0;
                
                // Warn once per thread if CIGAR is missing in auto mode (only after burn-in check)
                // This warning is for when error model is available but CIGAR is missing per alignment
                static thread_local bool warned_cigar_missing = false;
                if (!warned_cigar_missing && use_model && err.ref != nullptr && err.cigar1.empty()) {
                    string errorModelMode = P_.quant.transcriptVB.errorModelMode;
                    if (errorModelMode == "auto") {
                        P_.inOut->logMain << "WARNING: CIGAR missing for alignment; "
                                          << "errLike disabled for this read (auto mode)." << "\n";
                        warned_cigar_missing = true;
                    }
                }
            }
        } else {
            aln.has_err_like = false;
            aln.err_like = 0.0;
        }
        
        current_read_alignments_.push_back(aln);
        err_data.push_back(std::move(err));
    }
    
    if (current_read_alignments_.empty()) return;
    
    // Match Salmon/CLI ordering: sort alignments by transcript_id
    std::vector<size_t> sort_idx(current_read_alignments_.size());
    for (size_t i = 0; i < sort_idx.size(); ++i) {
        sort_idx[i] = i;
    }
    std::sort(sort_idx.begin(), sort_idx.end(),
              [&](size_t a, size_t b) {
                  return current_read_alignments_[a].transcript_id <
                         current_read_alignments_[b].transcript_id;
              });
    std::vector<RawAlignment> sorted_alignments;
    std::vector<ErrModelData> sorted_err_data;
    sorted_alignments.reserve(current_read_alignments_.size());
    sorted_err_data.reserve(err_data.size());
    for (size_t idx : sort_idx) {
        sorted_alignments.push_back(current_read_alignments_[idx]);
        if (idx < err_data.size()) {
            sorted_err_data.push_back(std::move(err_data[idx]));
        } else {
            sorted_err_data.emplace_back();
        }
    }
    current_read_alignments_.swap(sorted_alignments);
    err_data.swap(sorted_err_data);
    
    // Apply Salmon's filterAndCollectAlignments if enabled (min_score_fraction > 0)
    std::vector<RawAlignment> filtered_alignments;
    std::vector<ErrModelData> filtered_err_data;
    filtered_alignments.reserve(current_read_alignments_.size());
    filtered_err_data.reserve(current_read_alignments_.size());
    if (params.apply_filtering && params.min_score_fraction > 0.0) {
        // Compute best score for filtering
        int32_t best_score = std::numeric_limits<int32_t>::min();
        for (const auto& aln : current_read_alignments_) {
            if (aln.score > best_score) {
                best_score = aln.score;
            }
        }
        
        // Apply min_score_fraction filtering (Salmon's score fraction gating)
        int32_t min_score_threshold = static_cast<int32_t>(static_cast<double>(best_score) * params.min_score_fraction);
        for (size_t i = 0; i < current_read_alignments_.size(); ++i) {
            const auto& aln = current_read_alignments_[i];
            if (aln.score < min_score_threshold) {
                continue;
            }
            filtered_alignments.push_back(aln);
            if (i < err_data.size()) {
                filtered_err_data.push_back(std::move(err_data[i]));
            } else {
                filtered_err_data.emplace_back();
            }
        }
    } else {
        filtered_alignments = current_read_alignments_;
        filtered_err_data = std::move(err_data);
    }
    
    if (filtered_alignments.empty()) return;
    
    // Mini-batch gating (matches Salmon's processedReads += batchReads behavior)
    // Salmon gates on processed reads per mini-batch (size 1000), not per-read
    // This ensures useAuxParams stays constant for the whole batch
    // Key: batch_start_count_ is persisted across reads and only updated at batch boundaries
    
    if (!P_.quant.transcriptVB.inDetectionMode) {
        // If starting a new batch (local_batch_reads_ == 0), capture batch start count
        // This count persists for all reads in this batch
        if (local_batch_reads_ == 0) {
            batch_start_count_ = global_processed_fragments.load(std::memory_order_relaxed);
        }
    }
    
    // Set burn-in tracking in params using persisted batch start count
    // This matches Salmon's behavior: gates on processed read groups per mini-batch
    // useAuxParams stays constant for the whole batch (matches Salmon)
    params.num_processed_fragments = batch_start_count_;
    
    // Store batch info in params for trace instrumentation (before incrementing)
    params.batch_reads = local_batch_reads_;
    params.batch_start_count = batch_start_count_;
    
    // Compute auxProbs for this read (enable trace info if tracing)
    bool do_trace = trace_out_.is_open() && qname != nullptr &&
                    (trace_limit_ == 0 || traced_count_ < trace_limit_);
    ReadMapping mapping = computeAuxProbs(filtered_alignments, params, do_trace);
    
    // AFTER computeAuxProbs: increment batch counter and flush if needed
    if (!P_.quant.transcriptVB.inDetectionMode) {
        // Increment thread-local batch counter AFTER computing auxProbs
        local_batch_reads_++;
        
        // Check if we've reached mini-batch size
        int mini_batch_size = P_.quant.transcriptVB.miniBatchSize;
        if (mini_batch_size > 0 && local_batch_reads_ >= static_cast<size_t>(mini_batch_size)) {
            // Flush batch: add batch count to global counter AFTER processing all reads in batch
            global_processed_fragments.fetch_add(local_batch_reads_, std::memory_order_relaxed);
            num_processed_fragments_ += local_batch_reads_;  // Keep thread-local for merge
            local_batch_reads_ = 0;  // Reset batch counter (next read will start new batch)
        }
    }
    
    // Increment global FLD observation counter AFTER computing auxProbs (for FLD stats only)
    // Only increment if this read has valid FLD observations (valid fragment length + compatible)
    // SKIP during detection mode (Salmon doesn't count detection reads)
    if (mapping.num_valid_fld_obs > 0 && !P_.quant.transcriptVB.inDetectionMode) {
        Parameters::global_fld_obs_count.fetch_add(1, std::memory_order_relaxed);
    }
    
    // Salmon-style FLD and error model updates with stochastic acceptance
    // Salmon uses: if (!burnedIn && r < exp(aln->logProb)) { update FLD and error model }
    // Where aln->logProb is NORMALIZED (sum to 1 across alignments for a read)
    // Updates use logForgettingMass as weight for FLD
    // Only update during burn-in (processedReads < 5e6), not post-burn-in
    static constexpr uint64_t kNumBurninFrags = 5000000;  // Salmon default
    uint64_t processed_reads = batch_start_count_;
    bool burned_in = (processed_reads >= kNumBurninFrags);
    
    // Skip updates during detection mode and post burn-in
    if (!P_.quant.transcriptVB.inDetectionMode && !burned_in && 
        !mapping.alignment_indices.empty() && !mapping.log_probs.empty() &&
        mapping.log_prob_denom != LOG_0) {
        
        // Get forgetting mass for updates (Salmon uses same mass for FLD and error model)
        // IMPORTANT: Only advance the forgetting mass calculator once per mini-batch, not per read!
        // Salmon calls getLogMassAndTimestep once per mini-batch, not once per read.
        // The mini-batch number is (batch_start_count_ / miniBatchSize)
        int mini_batch_size = P_.quant.transcriptVB.miniBatchSize;
        uint64_t current_batch_num = (mini_batch_size > 0) ? (batch_start_count_ / mini_batch_size) : 0;
        
        double log_forgetting_mass = 0.0;
        if (current_batch_num != cached_forgetting_batch_) {
            // New mini-batch: advance the forgetting mass calculator
            uint64_t timestep = 0;
            forgetting_mass_.getLogMassAndTimestep(log_forgetting_mass, timestep);
            cached_log_forgetting_mass_ = log_forgetting_mass;
            cached_forgetting_batch_ = current_batch_num;
        } else {
            // Same mini-batch: reuse cached forgetting mass
            log_forgetting_mass = cached_log_forgetting_mass_;
        }
        
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        
        for (size_t i = 0; i < mapping.alignment_indices.size(); ++i) {
            size_t aln_idx = mapping.alignment_indices[i];
            if (aln_idx >= filtered_alignments.size()) {
                continue;
            }
            
            // Get logProb for this alignment and NORMALIZE it
            // Salmon normalizes: aln->logProb -= denom, so exp(logProb) sums to 1
            double log_prob = (i < mapping.log_probs.size()) ? mapping.log_probs[i] : LOG_0;
            double normalized_log_prob = log_prob - mapping.log_prob_denom;
            
            // Stochastic acceptance: r < exp(normalized_logProb)
            // Since logProb is normalized, exp(normalized_logProb) is a valid probability [0,1]
            double r = uni(rng_);
            double acceptance_prob = std::exp(normalized_log_prob);
            if (r >= acceptance_prob) {
                continue;  // Reject this update
            }
            
            // FLD update (Salmon: fragLengthDist.addVal(fragLength, logForgettingMass))
            // Only for paired-end with valid fragment length
            int32_t frag_len = (i < mapping.frag_lens.size()) ? mapping.frag_lens[i] : -1;
            const auto& aln = filtered_alignments[aln_idx];
            bool is_paired = (aln.mate_status == MateStatus::PAIRED_END_PAIRED);
            if (is_paired && frag_len > 0 && frag_len <= static_cast<int32_t>(FLDAccumulator::MAX_FRAG_LEN)) {
                // Apply forgetting mass as weight (LOG space - Salmon style)
                // Salmon: fragLengthDist.addVal(fragLength, logForgettingMass);
                observedFLD_.add(frag_len, log_forgetting_mass);
                fld_dirty_ = true;
            }
            
            // Error model update (same stochastic acceptance pattern)
            if (use_error_model && alignment_model_ != nullptr && alignment_model_->canUpdate()) {
                if (aln_idx < filtered_err_data.size()) {
                    const ErrModelData& err = filtered_err_data[aln_idx];
                    if (err.ref != nullptr && !err.cigar1.empty()) {
                        const uint8_t* read1_seq = current_packed_read1_;
                        if (!aln.is_forward && current_packed_read1_rc_ != nullptr) {
                            read1_seq = current_packed_read1_rc_;
                        }
                        const uint8_t* read2_seq = current_packed_read2_;
                        if (err.paired && !aln.mate_is_forward && current_packed_read2_rc_ != nullptr) {
                            read2_seq = current_packed_read2_rc_;
                        }
                        if (err.paired && !err.cigar2.empty() &&
                            current_packed_read2_ != nullptr && current_read2_len_ > 0) {
                            alignment_model_->logLikelihoodPairedAndUpdate(
                                err.cigar1.data(), static_cast<uint32_t>(err.cigar1.size()),
                                read1_seq, static_cast<int32_t>(current_read1_len_),
                                err.cigar2.data(), static_cast<uint32_t>(err.cigar2.size()),
                                read2_seq, static_cast<int32_t>(current_read2_len_),
                                err.ref, err.pos1, err.pos2,
                                0.0, log_forgetting_mass);
                        } else {
                            bool is_left = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                                            aln.mate_status == MateStatus::SINGLE_END);
                            alignment_model_->logLikelihoodAndUpdate(
                                err.cigar1.data(), static_cast<uint32_t>(err.cigar1.size()),
                                read1_seq, static_cast<int32_t>(current_read1_len_),
                                err.ref, err.pos1, is_left,
                                0.0, log_forgetting_mass);
                        }
                    }
                }
            }
        }
    }
    
    if (use_error_model && alignment_model_ != nullptr) {
        alignment_model_->incrementObserved();
    }
    
    // Accumulate drop statistics
    dropped_incompat_ += mapping.num_dropped_incompat;
    dropped_missing_mate_fields_ += mapping.num_dropped_missing_mate_fields;
    dropped_unknown_obs_fmt_ += mapping.num_dropped_unknown_obs_fmt;
    
    // Set detection mode flag in trace info (for debugging)
    if (do_trace) {
        for (auto& tr : mapping.trace_info) {
            tr.in_detection_mode = P_.quant.transcriptVB.inDetectionMode;
        }
        writeTraceLine(qname, mapping, params, current_read_alignments_);
        traced_count_++;
    }
    
    if (mapping.transcript_ids.empty()) return;
    
    // mapping.aux_probs is already normalized to linear weights by computeAuxProbs
    // (see ec_builder.cpp line 465: mapping.aux_probs[i] = normalized_weight)
    // Use directly without exponentiation
    const std::vector<double>& linear_weights = mapping.aux_probs;
    
    // Sort transcript IDs for canonical EC representation
    std::vector<size_t> sort_indices(mapping.transcript_ids.size());
    for (size_t i = 0; i < sort_indices.size(); ++i) {
        sort_indices[i] = i;
    }
    std::sort(sort_indices.begin(), sort_indices.end(), 
        [&](size_t a, size_t b) {
            return mapping.transcript_ids[a] < mapping.transcript_ids[b];
        });
    
    // Reorder transcript IDs and weights
    std::vector<uint32_t> sorted_ids;
    std::vector<double> sorted_weights;
    sorted_ids.reserve(mapping.transcript_ids.size());
    sorted_weights.reserve(linear_weights.size());
    for (size_t idx : sort_indices) {
        sorted_ids.push_back(mapping.transcript_ids[idx]);
        sorted_weights.push_back(linear_weights[idx]);
    }
    
    // Create EC signature (sorted transcript IDs only - matches Salmon)
    // NOTE: Salmon keys ECs by transcript IDs only, weights are accumulated separately
    ECSignature sig;
    sig.transcript_ids = sorted_ids;
    
    // Find or create EC
    auto it = ec_signature_map_.find(sig);
    if (it != ec_signature_map_.end()) {
        // Existing EC: increment count and accumulate weights
        size_t ec_idx = it->second;
        EC& ec = ecTable_.ecs[ec_idx];
        ec.count += 1.0;
        // Accumulate weights (Salmon sums combinedWeights per EC)
        if (ec.weights.size() == sorted_weights.size()) {
            for (size_t i = 0; i < ec.weights.size(); ++i) {
                ec.weights[i] += sorted_weights[i];
            }
        }
    } else {
        // New EC: create entry
        EC ec;
        ec.transcript_ids = sorted_ids;
        ec.weights = sorted_weights;
        ec.count = 1.0;
        
        size_t ec_idx = ecTable_.ecs.size();
        ecTable_.ecs.push_back(ec);
        ecTable_.n_ecs++;
        ec_signature_map_[sig] = ec_idx;
    }
}

void TranscriptQuantEC::addReadAlignmentsSimple(const std::vector<uint32_t>& transcriptIds, const std::vector<double>& weights) {
    if (transcriptIds.empty()) return;
    
    // Create sorted indices
    std::vector<size_t> sort_indices(transcriptIds.size());
    for (size_t i = 0; i < sort_indices.size(); ++i) {
        sort_indices[i] = i;
    }
    std::sort(sort_indices.begin(), sort_indices.end(), 
        [&](size_t a, size_t b) {
            return transcriptIds[a] < transcriptIds[b];
        });
    
    // Build sorted transcript IDs and weights
    std::vector<uint32_t> sorted_ids;
    std::vector<double> sorted_weights;
    sorted_ids.reserve(transcriptIds.size());
    sorted_weights.reserve(weights.size());
    for (size_t idx : sort_indices) {
        sorted_ids.push_back(transcriptIds[idx]);
        if (idx < weights.size()) {
            sorted_weights.push_back(weights[idx]);
        } else {
            sorted_weights.push_back(1.0 / transcriptIds.size());  // Default uniform
        }
    }
    
    // Create EC signature (sorted transcript IDs only - matches Salmon)
    ECSignature sig;
    sig.transcript_ids = sorted_ids;
    
    // Find or create EC
    auto it = ec_signature_map_.find(sig);
    if (it != ec_signature_map_.end()) {
        // Existing EC: increment count and accumulate weights
        size_t ec_idx = it->second;
        EC& ec = ecTable_.ecs[ec_idx];
        ec.count += 1.0;
        // Accumulate weights
        if (ec.weights.size() == sorted_weights.size()) {
            for (size_t i = 0; i < ec.weights.size(); ++i) {
                ec.weights[i] += sorted_weights[i];
            }
        }
    } else {
        // New EC: create entry
        EC ec;
        ec.transcript_ids = sorted_ids;
        ec.weights = sorted_weights;
        ec.count = 1.0;
        
        size_t ec_idx = ecTable_.ecs.size();
        ecTable_.ecs.push_back(ec);
        ecTable_.n_ecs++;
        ec_signature_map_[sig] = ec_idx;
    }
}

void TranscriptQuantEC::addGCObservation(int32_t frag_start, int32_t frag_end, int32_t gc_pct, double weight) {
    observedGC_.inc(gc_pct, weight);
}

void TranscriptQuantEC::addFLDObservation(int32_t frag_len, double weight, double log_prob) {
    // Salmon updates FLD starting immediately (no pre-burn-in gating for FLD updates)
    // Only stop updating after post-burn-in (>= 5e6 processed read groups)
    static constexpr uint64_t num_burnin_frags = 5000000;
    
    // Get current global processed read groups count (matches Salmon's processedReads)
    uint64_t processed_reads = global_processed_fragments.load(std::memory_order_relaxed);
    
    if (processed_reads >= num_burnin_frags) {
        return;  // Post-burnin: no FLD updates (use conditional PMF)
    }
    
    // Pre-burnin: continue updating FLD (Salmon behavior - FLD updates start immediately)
    
    // Burn-in: apply stochastic acceptance if log_prob is provided
    if (log_prob > LOG_0 && log_prob < 0.0) {
        // Stochastic acceptance: accept if log(u) < logProb
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double u = uniform(rng_);
        double log_u = std::log(u);
        if (log_u >= log_prob) {
            return;  // Reject update
        }
    }
    
    // Apply forgetting mass (LOG space - Salmon style)
    double log_forgetting_mass = 0.0;
    uint64_t timestep = 0;
    forgetting_mass_.getLogMassAndTimestep(log_forgetting_mass, timestep);
    
    // Adjust weight by forgetting mass in LOG space
    // log(weight * forgetting_mass) = log(weight) + log_forgetting_mass
    double log_weight = (weight > 0.0) ? std::log(weight) : LOG_0;
    double log_adjusted_weight = log_weight + log_forgetting_mass;
    
    observedFLD_.add(frag_len, log_adjusted_weight);
    fld_dirty_ = true;  // Mark FLD as dirty after adding observation
}

void TranscriptQuantEC::setTranscriptLengths(const std::vector<int32_t>& lengths) {
    cached_transcript_lengths_ = lengths;
}

void TranscriptQuantEC::finalize() {
    // Normalize observed GC if needed
    // (Will be normalized later when computing bias ratios)
}

void TranscriptQuantEC::writeTraceLine(const char* qname, const ReadMapping& mapping, 
                                        const ECBuilderParams& params,
                                        const std::vector<RawAlignment>& alignments) {
    if (!trace_out_.is_open() || qname == nullptr) return;
    (void)alignments;
    
    auto keep_trace = [](const AlignmentTrace& tr) {
        return true;
    };
    
    // Build EC label (after rank/range factorization)
    std::vector<uint32_t> final_txp_ids = mapping.transcript_ids;
    std::vector<double> final_weights = mapping.aux_probs;
    if (params.use_rank_eq_classes && final_txp_ids.size() > 1) {
        applyRankEqClasses(final_txp_ids, final_weights);
    }
    if (params.use_range_factorization) {
        applyRangeFactorization(final_txp_ids, final_weights, params.range_factorization_bins);
    }
    
    // Output trace line in Salmon/CLI alignment-mode format
    trace_out_ << qname << "\t";
    
    trace_out_ << "txpIDs=";
    for (size_t i = 0; i < mapping.transcript_ids.size(); ++i) {
        if (i > 0) trace_out_ << ",";
        trace_out_ << mapping.transcript_ids[i];
    }
    trace_out_ << ";";
    
    trace_out_ << "as=";
    bool first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.as_tag;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "bestAS=" << mapping.best_as << ";";
    
    trace_out_ << "logFragProb=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.log_frag_prob;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "errLike=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.err_like;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "logCompat=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.log_compat_prob;
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "expFmt=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << static_cast<uint32_t>(tr.expected_format_id);
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "obsFmt=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << static_cast<uint32_t>(tr.observed_format_id);
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "isCompat=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.is_compat ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "mateStatus=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << static_cast<uint32_t>(tr.mate_status);
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "mateFields=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.mate_fields_set ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "fwd=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.is_forward ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "mateFwd=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.mate_is_forward ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "pos=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.pos;
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "matePos=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.mate_pos;
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "primary=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.is_primary ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";

    trace_out_ << "dropped=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.dropped ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "droppedIncompat=" << mapping.num_dropped_incompat << ";";
    
    trace_out_ << "orphan=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.is_orphan ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "auxProb=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.aux_prob;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "weight=";
    for (size_t i = 0; i < mapping.aux_probs.size(); ++i) {
        if (i > 0) trace_out_ << ",";
        trace_out_ << mapping.aux_probs[i];
    }
    trace_out_ << ";";
    
    trace_out_ << "ec_label=";
    for (size_t i = 0; i < final_txp_ids.size(); ++i) {
        if (i > 0) trace_out_ << ",";
        trace_out_ << final_txp_ids[i];
    }
    trace_out_ << ";";
    
    // Trace instrumentation fields
    trace_out_ << "fragLen=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.frag_len;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "useAuxParams=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.use_aux_params ? "1" : "0");
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "globalFragCount=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.global_frag_count;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "batchReads=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.batch_reads;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "batchStartCount=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << tr.batch_start_count;
        first = false;
    }
    trace_out_ << ";";
    
    trace_out_ << "inDetectionMode=";
    first = true;
    for (const auto& tr : mapping.trace_info) {
        if (!keep_trace(tr)) continue;
        if (!first) trace_out_ << ",";
        trace_out_ << (tr.in_detection_mode ? "1" : "0");
        first = false;
    }
    trace_out_ << "\n";
    
    trace_out_.flush();
}

void TranscriptQuantEC::merge(const TranscriptQuantEC& other) {
    // Merge drop statistics
    dropped_incompat_ += other.dropped_incompat_;
    dropped_missing_mate_fields_ += other.dropped_missing_mate_fields_;
    dropped_unknown_obs_fmt_ += other.dropped_unknown_obs_fmt_;
    
    // Merge EC tables
    for (const EC& other_ec : other.ecTable_.ecs) {
        // Sort transcript IDs for canonical representation
        std::vector<size_t> sort_indices(other_ec.transcript_ids.size());
        for (size_t i = 0; i < sort_indices.size(); ++i) {
            sort_indices[i] = i;
        }
        std::sort(sort_indices.begin(), sort_indices.end(),
            [&](size_t a, size_t b) {
                return other_ec.transcript_ids[a] < other_ec.transcript_ids[b];
            });
        
        std::vector<uint32_t> sorted_ids;
        std::vector<double> sorted_weights;
        sorted_ids.reserve(other_ec.transcript_ids.size());
        sorted_weights.reserve(other_ec.weights.size());
        for (size_t idx : sort_indices) {
            sorted_ids.push_back(other_ec.transcript_ids[idx]);
            if (idx < other_ec.weights.size()) {
                sorted_weights.push_back(other_ec.weights[idx]);
            }
        }
        
        // Create signature (transcript IDs only - matches Salmon)
        ECSignature sig;
        sig.transcript_ids = sorted_ids;
        
        // Find or create EC
        auto it = ec_signature_map_.find(sig);
        if (it != ec_signature_map_.end()) {
            // Existing EC: add counts and accumulate weights
            size_t ec_idx = it->second;
            EC& ec = ecTable_.ecs[ec_idx];
            ec.count += other_ec.count;
            // Accumulate weights
            if (ec.weights.size() == sorted_weights.size()) {
                for (size_t i = 0; i < ec.weights.size(); ++i) {
                    ec.weights[i] += sorted_weights[i];
                }
            }
        } else {
            // New EC: add entry
            EC ec;
            ec.transcript_ids = sorted_ids;
            ec.weights = sorted_weights;
            ec.count = other_ec.count;
            size_t ec_idx = ecTable_.ecs.size();
            ecTable_.ecs.push_back(ec);
            ecTable_.n_ecs++;
            ec_signature_map_[sig] = ec_idx;
        }
    }
    
    // Merge GC observations
    observedGC_.combineCounts(other.observedGC_);
    
    // Merge FLD observations
    observedFLD_.combine(other.observedFLD_);
    fld_dirty_ = true;  // Mark FLD as dirty after merging observations
    
    // Merge processed fragment counts (for burn-in tracking)
    num_processed_fragments_ += other.num_processed_fragments_;
}

std::vector<uint32_t> TranscriptQuantEC::parseCIGARString(const std::string& cigar_str) const {
    std::vector<uint32_t> cigar;
    int64_t num = 0;
    int sign = 1;
    bool in_num = false;

    auto push_op = [&](int32_t len, char op) {
        if (len <= 0) {
            return;
        }
        uint8_t code = 0;
        switch (op) {
            case 'M': code = BAM_CMATCH; break;
            case 'I': code = BAM_CINS; break;
            case 'D': code = BAM_CDEL; break;
            case 'N': code = BAM_CREF_SKIP; break;
            case 'S': code = BAM_CSOFT_CLIP; break;
            case 'H': code = BAM_CHARD_CLIP; break;
            case 'P': code = BAM_CPAD; break;
            case '=': code = BAM_CEQUAL; break;
            case 'X': code = BAM_CDIFF; break;
            default: return;
        }
        cigar.push_back((static_cast<uint32_t>(len) << BAM_CIGAR_SHIFT) | code);
    };

    for (char c : cigar_str) {
        if (c == '-') {
            sign = -1;
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c))) {
            num = num * 10 + (c - '0');
            in_num = true;
            continue;
        }
        if (!in_num) {
            continue;
        }
        int32_t len = static_cast<int32_t>(num * sign);
        num = 0;
        sign = 1;
        in_num = false;
        if (c == 'p' || c == 'P') {
            continue;  // gap between mates, not part of per-mate CIGAR
        }
        if (len < 0) {
            len = -len;
        }
        push_op(len, c);
    }

    return cigar;
}

std::vector<uint32_t> TranscriptQuantEC::extractCIGAR(Transcript* tr, bool is_read1) const {
    std::string cigar_str = tr->cigarString;
    if (cigar_str.empty()) {
        cigar_str = tr->generateCigarP();
    }

    // Single-end or no mate separator: parse entire CIGAR
    if (tr->readNmates < 2 || cigar_str.find('p') == std::string::npos) {
        return parseCIGARString(cigar_str);
    }

    std::vector<uint32_t> left_cigar;
    std::vector<uint32_t> right_cigar;
    std::vector<uint32_t>* cur = &left_cigar;
    int64_t num = 0;
    int sign = 1;
    bool in_num = false;

    auto push_op = [&](int32_t len, char op) {
        if (len <= 0) {
            return;
        }
        uint8_t code = 0;
        switch (op) {
            case 'M': code = BAM_CMATCH; break;
            case 'I': code = BAM_CINS; break;
            case 'D': code = BAM_CDEL; break;
            case 'N': code = BAM_CREF_SKIP; break;
            case 'S': code = BAM_CSOFT_CLIP; break;
            case 'H': code = BAM_CHARD_CLIP; break;
            case 'P': code = BAM_CPAD; break;
            case '=': code = BAM_CEQUAL; break;
            case 'X': code = BAM_CDIFF; break;
            default: return;
        }
        cur->push_back((static_cast<uint32_t>(len) << BAM_CIGAR_SHIFT) | code);
    };

    for (char c : cigar_str) {
        if (c == '-') {
            sign = -1;
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c))) {
            num = num * 10 + (c - '0');
            in_num = true;
            continue;
        }
        if (!in_num) {
            continue;
        }
        int32_t len = static_cast<int32_t>(num * sign);
        num = 0;
        sign = 1;
        in_num = false;
        if (c == 'p' || c == 'P') {
            // Switch to right mate; ignore gap length
            cur = &right_cigar;
            continue;
        }
        if (len < 0) {
            len = -len;
        }
        push_op(len, c);
    }

    if (right_cigar.empty()) {
        return left_cigar;
    }

    // Determine whether read1 is leftmost based on transcript exons
    bool read1_leftmost = true;
    uint iExRead1 = UINT_MAX;
    uint iExRead2 = UINT_MAX;
    for (uint i = 0; i < tr->nExons; i++) {
        if (tr->exons[i][EX_iFrag] == 0 && iExRead1 == UINT_MAX) {
            iExRead1 = i;
        }
        if (tr->exons[i][EX_iFrag] == 1 && iExRead2 == UINT_MAX) {
            iExRead2 = i;
        }
    }
    if (iExRead1 != UINT_MAX && iExRead2 != UINT_MAX) {
        int32_t pos1 = tr->exons[iExRead1][EX_G];
        int32_t pos2 = tr->exons[iExRead2][EX_G];
        read1_leftmost = (pos1 <= pos2);
    }

    if (is_read1) {
        return read1_leftmost ? left_cigar : right_cigar;
    }
    return read1_leftmost ? right_cigar : left_cigar;
}
