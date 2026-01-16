#include "ec_builder.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

// Check if alignment is compatible with library format (for single-end or orphaned reads)
// Ported from Salmon's compatibleHit() function
bool compatibleHit(const LibraryFormat& expected, bool isForward, MateStatus ms) {
    auto expectedStrand = expected.strandedness;
    auto expectedType = expected.type;
    
    switch (ms) {
    case MateStatus::SINGLE_END:
        if (isForward) { // U, SF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::S);
        } else { // U, SR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::A);
        }
        break;
        
    case MateStatus::PAIRED_END_LEFT:
        if (expected.orientation == ReadOrientation::SAME) {
            return (expectedStrand == ReadStrandedness::U ||
                    (expectedStrand == ReadStrandedness::S && isForward) ||
                    (expectedStrand == ReadStrandedness::A && !isForward));
        } else if (isForward) { // IU, ISF, OU, OSF, MU, MSF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::SA);
        } else { // IU, ISR, OU, OSR, MU, MSR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::AS);
        }
        break;
        
    case MateStatus::PAIRED_END_RIGHT:
        if (expected.orientation == ReadOrientation::SAME) {
            return (expectedStrand == ReadStrandedness::U ||
                    (expectedStrand == ReadStrandedness::S && isForward) ||
                    (expectedStrand == ReadStrandedness::A && !isForward));
        } else if (isForward) { // IU, ISR, OU, OSR, MU, MSR
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::AS);
        } else { // IU, ISF, OU, OSF, MU, MSF
            return (expectedStrand == ReadStrandedness::U || expectedStrand == ReadStrandedness::SA);
        }
        break;
        
    default:
        return false;
    }
}

// Check if paired-end alignment is compatible (check observed vs expected format)
bool compatibleHit(const LibraryFormat& expected, const LibraryFormat& observed) {
    if (observed.type != ReadType::PAIRED_END) {
        return false;
    }
    
    auto es = expected.strandedness;
    auto eo = expected.orientation;
    auto os = observed.strandedness;
    auto oo = observed.orientation;
    
    // If orientations are different, incompatible
    if (eo != oo) {
        return false;
    } else {
        // Orientations match - check strandedness
        return (es == ReadStrandedness::U || es == os);
    }
}

// Helper function to compute observed format from paired-end alignment
// Matches hitType() logic but uses RawAlignment fields
static LibraryFormat computeObservedFormat(const RawAlignment& aln) {
    if (aln.mate_status != MateStatus::PAIRED_END_PAIRED) {
        // Not a proper pair - return default
        return LibraryFormat::IU();
    }
    
    bool end1Fwd = aln.is_forward;
    bool end2Fwd = aln.mate_is_forward;
    int32_t end1Start = aln.pos;
    int32_t end2Start = aln.mate_pos;
    
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

// Main compatibility check function (matching Salmon's isCompatible)
bool isCompatible(const RawAlignment& aln, const LibraryFormat& expected_format, bool isForward) {
    // For unstranded libraries, all alignments are compatible
    if (expected_format.strandedness == ReadStrandedness::U) {
        return true;
    }
    
    // For paired-end paired reads, compute observed format from both mates and check compatibility
    if (aln.mate_status == MateStatus::PAIRED_END_PAIRED) {
        // If mate fields are missing for paired-end alignments, treat as incompatible.
        // Salmon drops these from paired-end ECs, so do not auto-accept.
        if (!aln.mate_fields_set) {
            return false;
        }
        
        // Compute observed format from both mates' orientations and positions
        LibraryFormat observed = computeObservedFormat(aln);
        // Check compatibility between expected and observed formats
        return compatibleHit(expected_format, observed);
    } else {
        // Single-end or orphan - use the compatibility check
        return compatibleHit(expected_format, isForward, aln.mate_status);
    }
}

// Compute auxProb for each transcript in a read's alignments
// Matches Salmon's SalmonQuantifyAlignments.cpp alignment-mode calculation
ReadMapping computeAuxProbs(
    const std::vector<RawAlignment>& alignments,
    const ECBuilderParams& params,
    bool enable_trace
) {
    ReadMapping mapping;
    mapping.aux_denom = LOG_0;  // -infinity in log space
    
    if (alignments.empty()) {
        return mapping;
    }
    
    // Step 1: Find bestAS across all alignments for this read
    // For paired-end: AS is sum of both mates; for single-end: AS tag value
    int32_t best_as = 0;
    if (params.use_as_without_cigar) {
        best_as = std::numeric_limits<int32_t>::min();
        for (const auto& aln : alignments) {
            int32_t aln_as = aln.score;  // Already contains AS (sum for paired, single AS for single)
            if (aln_as > best_as) {
                best_as = aln_as;
            }
        }
    }
    mapping.best_as = best_as;
    
    // For post-burnin conditional PMF, use cached log-CMF from FLD accumulator (O(1) lookup)
    // No need to compute CMF per read batch - FLD accumulator caches it
    
    // Step 2: Compute auxProb for each alignment
    for (size_t aln_idx = 0; aln_idx < alignments.size(); ++aln_idx) {
        const auto& aln = alignments[aln_idx];
        AlignmentTrace trace;
        trace.transcript_id = aln.transcript_id;
        trace.as_tag = params.use_as_without_cigar ? aln.score : 0;
        trace.best_as = best_as;
        // Use Salmon's formatID encoding for trace output parity
        trace.expected_format_id = params.lib_format.salmonFormatID();
        trace.mate_status = static_cast<uint8_t>(aln.mate_status);
        trace.mate_fields_set = aln.mate_fields_set;
        trace.is_forward = aln.is_forward;
        trace.mate_is_forward = aln.mate_is_forward;
        trace.pos = aln.pos;
        trace.mate_pos = aln.mate_pos;
        trace.is_primary = aln.is_primary;
        trace.dropped = false;
        trace.observed_format_id = 255;
        if (aln.mate_status == MateStatus::PAIRED_END_PAIRED && aln.mate_fields_set) {
            trace.observed_format_id = computeObservedFormat(aln).salmonFormatID();
        }
        
        // Initialize trace instrumentation fields (will be populated below)
        trace.frag_len = aln.fragment_len;
        trace.use_aux_params = false;  // Will be set below
        trace.global_frag_count = params.num_processed_fragments;  // Batch start count (for gating)
        trace.batch_reads = params.batch_reads;
        trace.batch_start_count = params.batch_start_count;
        trace.in_detection_mode = false;  // Will be set by caller if needed
        
        // Skip decoys (if needed - for now, process all)
        if (aln.is_decoy) {
            // Could skip decoys here if desired
        }
        
        // Step 2a: Check for missing mate fields in paired-end mode (per trace plan Step 2)
        const bool missing_mate_fields =
            (params.lib_format.type == ReadType::PAIRED_END) &&
            (aln.mate_status == MateStatus::PAIRED_END_PAIRED) &&
            (!aln.mate_fields_set);
        
        // Step 2b: Check for unknown observed format in paired-end mode (per trace plan Step 2)
        const bool unknown_obs_fmt =
            (params.lib_format.type == ReadType::PAIRED_END) &&
            (aln.mate_status == MateStatus::PAIRED_END_PAIRED) &&
            (trace.observed_format_id == 255);
        
        // Drop alignments with missing mate fields or unknown observed format in paired-end mode
        // (Match Salmon behavior: these are treated as incompatible and dropped)
        if (missing_mate_fields || unknown_obs_fmt) {
            trace.dropped_incompat = true;
            trace.dropped = true;
            trace.is_compat = false;
            trace.log_compat_prob = LOG_1;  // logCompat = 0 (LOG_1) when incompatPrior=0
            if (missing_mate_fields) {
                mapping.num_dropped_missing_mate_fields++;
            }
            if (unknown_obs_fmt) {
                mapping.num_dropped_unknown_obs_fmt++;
            }
            // Count as incompatible drop (for backward compatibility)
            mapping.num_dropped_incompat++;
            if (enable_trace) {
                mapping.trace_info.push_back(trace);
            }
            continue;
        }
        
        // Check compatibility for alignments with valid mate fields
        double log_compat_prob;
        bool is_compat = isCompatible(aln, params.lib_format, aln.is_forward);
        trace.dropped_incompat = false;

        // 3. logAlignCompatProb: Compatibility probability (compute before skipping)
        if (!is_compat) {
            if (params.incompat_prior > 0.0) {
                // Use log(incompatPrior) as penalty
                log_compat_prob = std::log(params.incompat_prior);
            } else {
                // Match Salmon: incompatPrior == 0 => logCompat = 0
                log_compat_prob = LOG_1;
            }
        } else {
            log_compat_prob = LOG_1;
        }
        trace.log_compat_prob = log_compat_prob;
        trace.is_compat = is_compat;
        
        // If ignore_incompat=true and not compatible, skip this alignment
        if (params.ignore_incompat && !is_compat) {
            trace.dropped_incompat = true;
            trace.dropped = true;
            mapping.num_dropped_incompat++;
            if (enable_trace) {
                mapping.trace_info.push_back(trace);
            }
            continue;
        }
        
        // 1. logFragProb: Fragment length probability
        double log_frag_prob = 0.0;
        bool is_orphan = (aln.mate_status == MateStatus::PAIRED_END_LEFT ||
                         aln.mate_status == MateStatus::PAIRED_END_RIGHT);
        trace.is_orphan = is_orphan;
        
        // Salmon FLD phase-based logFragProb (using GLOBAL FLD observation count):
        // - Pre-burnin (<5000 FLD obs): logFragProb = 0, errLike = 0 (useAuxParams = false)
        // - Burn-in (5000-5e6 FLD obs): logFragProb = log PMF(len)
        // - Post-burnin (>5e6 FLD obs): logFragProb = log PMF(len) - log CMF(refLen) [conditional]
        // Note: Use FLD observation count (not read count) for gating - matches Salmon's behavior
        // FLD observations only increment when fragment length is valid and alignment is compatible
        // params.num_processed_fragments is now the GLOBAL FLD observation count (atomic counter)
        // Get current FLD count BEFORE processing this alignment (for gating decision)
        uint64_t current_fld_count = params.num_processed_fragments;
        bool use_aux_params = current_fld_count >= params.num_pre_burnin_frags;
        trace.use_aux_params = use_aux_params;  // Set trace instrumentation
        bool in_pre_burnin = !use_aux_params;
        
        // Track if this alignment should count toward FLD observations
        // Increment FLD counter when: valid fragment length, compatible, not orphaned, not missing fields
        // Note: unknown_obs_fmt and missing_mate_fields are already computed above
        bool should_count_fld = !is_orphan && 
                                aln.fragment_len >= 0 && 
                                aln.fragment_len <= FLDAccumulator::MAX_FRAG_LEN &&
                                is_compat &&
                                !missing_mate_fields &&
                                !unknown_obs_fmt;
        bool in_burnin = use_aux_params && params.num_processed_fragments < params.num_burnin_frags;
        bool post_burnin = params.num_processed_fragments >= params.num_burnin_frags;
        
        if (is_orphan) {
            // For orphans: use LOG_EPSILON (orphanProb = LOG_EPSILON when discardOrphansAln=false)
            log_frag_prob = LOG_EPSILON;
        } else if (in_pre_burnin) {
            // Pre-burnin: logFragProb = 0 (Salmon behavior)
            log_frag_prob = 0.0;
        } else if (params.use_frag_len_dist && aln.fragment_len >= 0) {
            // Burn-in or post-burnin: compute logFragProb from FLD
            double log_pmf = 0.0;
            
            // Prefer O(1) log-space computation from FLD accumulator if available
            if (params.fld_accumulator != nullptr) {
                log_pmf = params.fld_accumulator->getLogProb(aln.fragment_len);
            } else if (params.fld_pmf != nullptr) {
                // Fallback to PMF lookup (requires PMF recomputation)
                int32_t frag_len = aln.fragment_len;
                if (frag_len < static_cast<int32_t>(params.fld_pmf->size())) {
                    double prob = (*params.fld_pmf)[frag_len];
                    if (prob > 0.0) {
                        log_pmf = std::log(prob);
                    } else {
                        log_pmf = LOG_0;  // Zero probability = -infinity
                    }
                } else {
                    log_pmf = 0.0;  // Fragment length out of range
                }
            } else {
                log_pmf = 0.0;  // No FLD available
            }
            
            if (post_burnin && params.use_conditional_pmf && 
                params.transcript_lengths != nullptr &&
                aln.transcript_id < params.transcript_lengths->size()) {
                // Post-burnin: use conditional PMF (PMF - CMF)
                // logFragProb = log PMF(len) - log CMF(refLen)
                // CMF(refLen) = sum(PMF[i] for i <= refLen - fragLen)
                int32_t ref_len = (*params.transcript_lengths)[aln.transcript_id];
                int32_t frag_len = aln.fragment_len;
                int32_t max_valid_len = ref_len - frag_len + 1;
                
                if (max_valid_len > 0 && params.fld_accumulator != nullptr) {
                    // Use cached log-CMF from FLD accumulator (O(1) lookup)
                    double log_cmf = params.fld_accumulator->getLogCMF(max_valid_len);
                    if (log_cmf != LOG_0) {
                        log_frag_prob = log_pmf - log_cmf;
                    } else {
                        log_frag_prob = LOG_0;  // Invalid CMF
                    }
                } else {
                    log_frag_prob = log_pmf;  // Fallback to unconditional PMF
                }
            } else {
                // Burn-in: use unconditional PMF
                log_frag_prob = log_pmf;
            }
        } else {
            // For parity with --noFragLengthDist: set to 0.0
            log_frag_prob = 0.0;
        }
        trace.log_frag_prob = log_frag_prob;
        
        // 2. errLike: Error model likelihood
        // Salmon gates aux params (errLike + logFragProb) during pre-burn-in using GLOBAL counter:
        // - Pre-burnin (<5000 frags): errLike = 0, logFragProb = 0 (useAuxParams = false)
        // - Post-burnin (>=5000 frags): use actual error model/FLD values (useAuxParams = true)
        // The error model still updates during pre-burnin, but auxProb doesn't use it until post-burn-in
        int32_t aln_as = params.use_as_without_cigar ? aln.score : 0;
        double err_like = 0.0;
        if (use_aux_params) {
            // Only compute errLike after pre-burnin threshold (GLOBAL gating)
            // Auto mode: prefer CIGAR if available (has_err_like), fallback to AS
            // CIGAR mode: only use CIGAR-based (has_err_like)
            // AS mode: only use AS-based
            // Off mode: err_like = 0.0
            if (params.use_error_model && aln.has_err_like) {
                // Use pre-computed CIGAR-based error likelihood
                err_like = aln.err_like;
            } else if (params.use_as_without_cigar && (!params.use_error_model || !aln.has_err_like)) {
                // AS-based likelihood (RapMap/Pufferfish style)
                // Use when: AS mode enabled, OR auto mode but CIGAR not available
                double score_diff = static_cast<double>(best_as - aln_as);
                err_like = -params.score_exp * score_diff;
            }
        }
        // Pre-burnin: err_like remains 0.0 (Salmon behavior, GLOBAL gating)
        trace.err_like = err_like;
        
        // 3. startPosProb: Start position probability (for paired alignments)
        // startPosProb = -log(refLen - fragLen + 1) when paired and length correction enabled
        // Only applied when frag_len is valid and ref_len - frag_len + 1 > 0
        // If invalid, set to LOG_0 and drop alignment (consistent with Salmon)
        double start_pos_prob = 0.0;
        if (params.use_start_pos_prob && 
            aln.mate_status == MateStatus::PAIRED_END_PAIRED &&
            aln.fragment_len >= 0 &&
            params.transcript_lengths != nullptr &&
            aln.transcript_id < params.transcript_lengths->size()) {
            int32_t ref_len = (*params.transcript_lengths)[aln.transcript_id];  // Transcript length (not genomic)
            int32_t frag_len = aln.fragment_len;
            int32_t valid_start_positions = ref_len - frag_len + 1;
            if (valid_start_positions > 0) {
                start_pos_prob = -std::log(static_cast<double>(valid_start_positions));
            } else {
                start_pos_prob = LOG_0;  // Invalid: fragment longer than transcript or invalid positions
            }
        }
        
        // Combine: auxProb = logFragProb + errLike + logAlignCompatProb
        // Note: Salmon does NOT include startPosProb in auxProb (it's only in aln->logProb)
        // Remove startPosProb from auxProb for Salmon parity
        double aux_prob = log_frag_prob + err_like + log_compat_prob;
        trace.aux_prob = aux_prob;
        trace.start_pos_prob = start_pos_prob;
        
        // Compute full alignment log probability (Salmon's aln->logProb)
        // logProb = transcriptLogCount + auxProb + startPosProb
        // For now, transcriptLogCount = 0 (uniform prior mass)
        double transcript_log_count = 0.0;  // TODO: use actual transcript mass when available
        double log_prob = transcript_log_count + aux_prob + start_pos_prob;
        trace.log_prob = log_prob;
        
        // Skip if aux_prob is LOG_0 (incompatible and incompatPrior == 0.0)
        if (aux_prob == LOG_0) {
            if (enable_trace) {
                mapping.trace_info.push_back(trace);
            }
            continue;
        }
        
        // Match Salmon's filtering checks:
        // 1. transcriptLogCount != LOG_0 (transcript mass check)
        //    In initial round, all transcripts should have prior mass, so this always passes.
        //    We don't track transcript mass, so we skip this check.
        //    TODO: If we add transcript mass tracking, add: if (transcriptLogCount == LOG_0) continue;
        
        // 2. startPosProb != LOG_0 (start position probability)
        //    With --noLengthCorrection: startPosProb = -logRefLength = -1.0 (not LOG_0)
        //    So this check always passes. We skip it since we don't use length correction.
        //    TODO: If we add length correction, add: if (startPosProb == LOG_0) continue;
        
        mapping.transcript_ids.push_back(aln.transcript_id);
        mapping.aux_probs.push_back(aux_prob);
        mapping.log_probs.push_back(log_prob);
        mapping.frag_lens.push_back(aln.fragment_len);
        mapping.alignment_indices.push_back(aln_idx);
        mapping.aux_denom = logAdd(mapping.aux_denom, aux_prob);  // log-sum-exp
        mapping.log_prob_denom = logAdd(mapping.log_prob_denom, log_prob);  // for logProb normalization
        
        // Count valid FLD observations (for pre-burn-in gating)
        // Only count once per read (not per alignment), so count only for the first valid alignment
        if (should_count_fld && mapping.num_valid_fld_obs == 0) {
            mapping.num_valid_fld_obs = 1;  // Count once per read
        }
        
        if (enable_trace) {
            mapping.trace_info.push_back(trace);
        }
    }
    
    // Step 3: Normalize weights: weight = exp(auxProb - auxDenom)
    // Note: Salmon keeps duplicate transcript IDs for multi-mapping reads
    // (each alignment location is a separate entry in the EC)
    if (mapping.aux_denom != LOG_0) {
        for (size_t i = 0; i < mapping.aux_probs.size(); ++i) {
            double normalized_weight = std::exp(mapping.aux_probs[i] - mapping.aux_denom);
            mapping.aux_probs[i] = normalized_weight;
        }
    }
    
    return mapping;
}

// Apply range factorization by appending bin IDs to transcript list
// Matches Salmon's SalmonQuantify.cpp lines 578-586
void applyRangeFactorization(
    std::vector<uint32_t>& txp_ids,
    const std::vector<double>& aux_probs,
    uint32_t range_factorization_bins
) {
    if (range_factorization_bins == 0) return;
    
    int32_t txps_size = static_cast<int32_t>(txp_ids.size());
    // rangeCount = sqrt(n) + bins (Salmon's formula)
    int32_t range_count = static_cast<int32_t>(std::sqrt(static_cast<double>(txps_size))) + static_cast<int32_t>(range_factorization_bins);
    
    // Append range bin numbers as pseudo-transcript IDs
    for (int32_t i = 0; i < txps_size; i++) {
        int32_t range_number = static_cast<int32_t>(aux_probs[i] * range_count);
        txp_ids.push_back(static_cast<uint32_t>(range_number));
    }
}

// Sort transcripts by conditional probability for rank-based ECs
// Matches Salmon's SalmonQuantify.cpp lines 557-575
void applyRankEqClasses(
    std::vector<uint32_t>& txp_ids,
    std::vector<double>& aux_probs
) {
    if (txp_ids.size() <= 1) return;
    
    std::vector<size_t> inds(txp_ids.size());
    std::iota(inds.begin(), inds.end(), 0);
    
    // Sort indices by aux_probs (ascending order)
    std::sort(inds.begin(), inds.end(), [&aux_probs](size_t i, size_t j) {
        return aux_probs[i] < aux_probs[j];
    });
    
    std::vector<uint32_t> txp_ids_new(txp_ids.size());
    std::vector<double> aux_probs_new(aux_probs.size());
    for (size_t r = 0; r < inds.size(); ++r) {
        txp_ids_new[r] = txp_ids[inds[r]];
        aux_probs_new[r] = aux_probs[inds[r]];
    }
    
    std::swap(txp_ids, txp_ids_new);
    std::swap(aux_probs, aux_probs_new);
}

// Build equivalence classes from filtered alignments
// This is the main entry point that combines all the steps
ECTable buildEquivalenceClasses(
    const std::vector<std::vector<RawAlignment>>& read_alignments,
    const ECBuilderParams& params,
    size_t num_transcripts
) {
    ECTable ec_table;
    ec_table.n_transcripts = num_transcripts;
    
    // Map from EC label (order-sensitive) to EC index
    // Salmon treats TranscriptGroup ordering as significant.
    std::unordered_map<std::string, size_t> ec_map;
    
    for (const auto& alignments : read_alignments) {
        if (alignments.empty()) continue;
        
        // Compute auxProbs for this read (tracing disabled in EC building)
        ReadMapping mapping = computeAuxProbs(alignments, params, false);
        
        if (mapping.transcript_ids.empty()) continue;
        
        // Apply rank-based ECs if enabled
        if (params.use_rank_eq_classes && mapping.transcript_ids.size() > 1) {
            applyRankEqClasses(mapping.transcript_ids, mapping.aux_probs);
        }
        
        // Apply range factorization if enabled
        if (params.use_range_factorization) {
            applyRangeFactorization(mapping.transcript_ids, mapping.aux_probs, params.range_factorization_bins);
        }
        
        std::string ec_key;
        for (size_t i = 0; i < mapping.transcript_ids.size(); ++i) {
            if (i > 0) ec_key += ",";
            ec_key += std::to_string(mapping.transcript_ids[i]);
        }
        
        // Find or create EC
        auto it = ec_map.find(ec_key);
        if (it == ec_map.end()) {
            // Create new EC
            EC ec;
            ec.transcript_ids = mapping.transcript_ids;
            // mapping.aux_probs is already normalized to linear weights by computeAuxProbs
            // (see line 465: mapping.aux_probs[i] = normalized_weight)
            // Store directly without exponentiation
            ec.weights = mapping.aux_probs;
            ec.count = 1.0;
            ec_table.ecs.push_back(ec);
            ec_map[ec_key] = ec_table.ecs.size() - 1;
        } else {
            // Increment count and accumulate weights (Salmon sums linear weights per EC)
            EC& ec = ec_table.ecs[it->second];
            ec.count += 1.0;
            if (ec.weights.size() == mapping.aux_probs.size()) {
                for (size_t i = 0; i < ec.weights.size(); ++i) {
                    // Add linear weight (mapping.aux_probs is already linear)
                    ec.weights[i] += mapping.aux_probs[i];
                }
            }
        }
    }
    
    ec_table.n_ecs = ec_table.ecs.size();
    return ec_table;
}
