// Sparse-gate parity harness for the VB/EM finalizer.
//
// Proves that restricting VB/EM to transcripts with non-zero equivalence-class
// support is EXACT, not an approximation. It runs the same VB twice over the
// same merged sidecar evidence:
//
//   dense  - the full transcript index space, as today
//   gated  - only transcripts appearing in >=1 EC, with the index space
//            compacted and EC transcript ids remapped, then the result
//            expanded back to the full index space
//
// The two count vectors must agree exactly. Any disagreement means the gate is
// not exact and must not be used, since it would perturb quantification for
// reasons unrelated to sharding.
//
// The gate is only a compaction of the transcript index space, so it leaves the
// deterministic fixed-thread-order reduction in vb_engine untouched.
//
// usage: sparse_gate_parity QUANT_SF SIDECAR...
//   QUANT_SF supplies transcript names, lengths and effective lengths, and its
//   NumReads column is reported alongside for context.

#include "TranscriptVBSidecar.h"
#include "em_types.h"
#include "vb_engine.h"
#include "fld_accumulator.h"
#include "effective_length_wrapper.h"
#include "alignment_model.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <map>
#include <unordered_map>
#include <functional>

namespace {

double nowSec() {
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}


struct QuantRow {
    std::string name;
    double length = 0.0, eff_length = 0.0, num_reads = 0.0;
};

bool loadQuant(const std::string& path, std::vector<QuantRow>& rows) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    if (!std::getline(in, line)) return false;  // header
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        QuantRow r;
        std::string tpm;
        if (!std::getline(ss, r.name, '\t')) continue;
        std::string len, eff, nr;
        std::getline(ss, len, '\t');
        std::getline(ss, eff, '\t');
        std::getline(ss, tpm, '\t');
        std::getline(ss, nr, '\t');
        r.length = std::atof(len.c_str());
        r.eff_length = std::atof(eff.c_str());
        r.num_reads = std::atof(nr.c_str());
        rows.push_back(std::move(r));
    }
    return true;
}

EMParams starLikeParams() {
    // STAR's in-process VB converged at 1863 iterations, far past em_quant's
    // 200 default, so the limit must be generous or the comparison is invalid.
    EMParams p;
    p.use_vb = true;
    p.max_iters = 10000;
    p.tolerance = 0.01;
    p.vb_prior = 0.01;
    p.zero_threshold = 1e-8;
    return p;
}


// Verbatim replica of TranscriptQuantEC::finalize (TranscriptQuantEC.cpp:991).
// The sidecar is written from getECTable() BEFORE finalize() runs, so the
// evidence it carries is unnormalised and every consumer must normalise before
// VB. No Parameters dependency; the class itself is not needed.
void finalizeECWeights(ECTable& t) {
    for (EC& ec : t.ecs) {
        if (ec.weights.empty()) continue;
        double sum = 0.0;
        for (double w : ec.weights)
            if (std::isfinite(w) && w > 0.0) sum += w;
        if (sum > 0.0) {
            const double inv = 1.0 / sum;
            for (double& w : ec.weights)
                w = (std::isfinite(w) && w > 0.0) ? w * inv : 0.0;
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: sparse_gate_parity QUANT_SF SIDECAR...\n";
        return 2;
    }
    std::vector<QuantRow> quant;
    if (!loadQuant(argv[1], quant)) {
        std::cerr << "FAIL cannot read " << argv[1] << "\n";
        return 1;
    }
    std::cout << "reference rows: " << quant.size() << "\n";

    // Optional transcriptome FASTA as argv[2], enabling the GC axis.
    std::string tx_fasta;
    int first_sidecar = 2;
    if (argc > 2) {
        std::string a2 = argv[2];
        if (a2.size() > 3 && a2.rfind(".fa") == a2.size() - 3) { tx_fasta = a2; first_sidecar = 3; }
    }
    std::vector<transcript_vb_sidecar::Sidecar> shards;
    std::string err;
    for (int i = first_sidecar; i < argc; ++i) {
        transcript_vb_sidecar::Sidecar s;
        if (!transcript_vb_sidecar::read(argv[i], s, err)) {
            std::cerr << "FAIL read " << argv[i] << ": " << err << "\n";
            return 1;
        }
        shards.push_back(std::move(s));
    }
    transcript_vb_sidecar::Sidecar gathered;
    if (!transcript_vb_sidecar::merge(shards, gathered, err)) {
        std::cerr << "FAIL merge: " << err << "\n";
        return 1;
    }
    ECTable& ecs = gathered.evidence.ec_table;
    finalizeECWeights(ecs);   // STAR calls mergedEC.finalize() before VB
    const std::size_t n_tx = ecs.n_transcripts;
    std::cout << "merged ECs: " << ecs.ecs.size()
              << "  transcripts: " << n_tx << "\n";
    if (quant.size() != n_tx) {
        std::cerr << "FAIL reference row count != sidecar transcript count\n";
        return 1;
    }

    auto fillState = [&](TranscriptState& st, const std::vector<std::size_t>& idx) {
        st.resize(idx.size());
        for (std::size_t k = 0; k < idx.size(); ++k) {
            const QuantRow& r = quant[idx[k]];
            st.names[k] = r.name;
            st.lengths[k] = r.length;
            st.eff_lengths[k] = r.eff_length > 0.0 ? r.eff_length : r.length;
            st.abundances[k] = 0.0;
            st.counts[k] = 0.0;
        }
    };

    // ---------------------------------------------------------------- dense --
    std::vector<std::size_t> all(n_tx);
    for (std::size_t i = 0; i < n_tx; ++i) all[i] = i;
    TranscriptState dense_state;
    fillState(dense_state, all);
    ECTable dense_ecs = ecs;
    const EMParams params = starLikeParams();
    std::cout << "\nrunning dense VB over " << n_tx << " transcripts...\n";
    EMResult dense = run_vb(dense_ecs, dense_state, params);
    std::cout << "  converged=" << (dense.converged ? "yes" : "no")
              << " iterations=" << dense.iterations << "\n";

    // ------------------------------------------------- dense repeat control --
    // The engine's E-step uses schedule(dynamic), so thread-to-EC assignment
    // can vary between runs. Repeating dense establishes the engine's own
    // reproducibility floor: the gate can only be judged exact relative to it.
    // Characterise the engine's own reproducibility floor over several repeats.
    // schedule(dynamic) in the E-step means dense does not reproduce itself
    // bit-for-bit, so a single pair proves nothing; the gate must be judged
    // against the distribution.
    double repro_floor = 0.0;
    {
        const int repeats = 4;
        std::vector<double> floors;
        std::vector<std::vector<double>> runs;
        runs.push_back(dense.counts);
        for (int r = 1; r < repeats; ++r) {
            TranscriptState st; fillState(st, all);
            ECTable e = ecs;
            EMResult res = run_vb(e, st, params);
            runs.push_back(res.counts);
        }
        for (std::size_t a = 0; a < runs.size(); ++a)
            for (std::size_t b = a + 1; b < runs.size(); ++b) {
                double mr = 0.0;
                for (std::size_t i = 0; i < n_tx; ++i) {
                    const double x = runs[a][i], y = runs[b][i];
                    const double m = std::max(std::fabs(x), std::fabs(y));
                    if (m > 0.0) mr = std::max(mr, std::fabs(x - y) / m);
                }
                floors.push_back(mr);
            }
        std::sort(floors.begin(), floors.end());
        repro_floor = floors.back();
        std::printf("\nengine reproducibility floor over %d dense runs (%zu pairs):\n",
                    repeats, floors.size());
        for (double f : floors) std::printf("  max rel delta %.6g\n", f);
        std::printf("  worst observed floor = %.6g\n", repro_floor);
    }

    // ---------------------------------------------------------------- gated --
    std::vector<char> supported(n_tx, 0);
    for (const auto& ec : ecs.ecs)
        for (std::uint32_t tid : ec.transcript_ids)
            if (tid < n_tx) supported[tid] = 1;
    std::vector<std::size_t> keep;
    std::vector<std::uint32_t> remap(n_tx, UINT32_MAX);
    for (std::size_t i = 0; i < n_tx; ++i)
        if (supported[i]) { remap[i] = static_cast<std::uint32_t>(keep.size()); keep.push_back(i); }
    std::cout << "gate: " << keep.size() << " of " << n_tx
              << " transcripts have >=1 EC ("
              << (100.0 * keep.size() / n_tx) << "%)\n";

    ECTable gated_ecs = ecs;
    gated_ecs.n_transcripts = keep.size();
    for (auto& ec : gated_ecs.ecs)
        for (auto& tid : ec.transcript_ids) tid = remap[tid];
    TranscriptState gated_state;
    fillState(gated_state, keep);
    // Salmon derives uniformPrior from ALL transcripts, so the compacted run
    // must divide by the original count, not the compacted one.
    EMParams gated_params = params;
    gated_params.num_active_override = n_tx;
    std::cout << "running gated VB over " << keep.size() << " transcripts...\n";
    EMResult gated = run_vb(gated_ecs, gated_state, gated_params);
    std::cout << "  converged=" << (gated.converged ? "yes" : "no")
              << " iterations=" << gated.iterations << "\n";

    // Expand gated back into the full index space.
    std::vector<double> expanded(n_tx, 0.0);
    for (std::size_t k = 0; k < keep.size(); ++k) expanded[keep[k]] = gated.counts[k];

    // ------------------------------------------------------------ comparison --
    std::size_t differing = 0, dropped_nonzero = 0;
    double max_abs = 0.0, max_rel = 0.0, sum_dense = 0.0, sum_gated = 0.0;
    for (std::size_t i = 0; i < n_tx; ++i) {
        const double a = dense.counts[i], b = expanded[i];
        sum_dense += a; sum_gated += b;
        if (a != b) {
            ++differing;
            const double d = std::fabs(a - b);
            max_abs = std::max(max_abs, d);
            const double m = std::max(std::fabs(a), std::fabs(b));
            if (m > 0.0) max_rel = std::max(max_rel, d / m);
        }
        if (!supported[i] && a != 0.0) ++dropped_nonzero;
    }
    std::printf("\ndense vs gated:\n"
                "  iterations           %u vs %u\n"
                "  transcripts differing %zu / %zu\n"
                "  max abs delta         %.17g\n"
                "  max rel delta         %.17g\n"
                "  total counts          %.6f vs %.6f (delta %.6g)\n"
                "  unsupported transcripts with non-zero dense count: %zu\n",
                dense.iterations, gated.iterations, differing, n_tx,
                max_abs, max_rel, sum_dense, sum_gated, sum_dense - sum_gated,
                dropped_nonzero);

    // Bit-equality is the wrong bar: the engine's E-step uses schedule(dynamic),
    // so thread-to-EC assignment varies and dense does not reproduce itself
    // bit-for-bit either. The meaningful test is whether gating perturbs results
    // by more than the engine's own run-to-run floor, measured above.
    int rc = 0;
    if (dropped_nonzero != 0) {
        std::cout << "FAIL gate dropped transcripts that carried counts\n";
        rc = 1;
    } else if (dense.iterations != gated.iterations) {
        std::cout << "FAIL gate changed the iteration count\n";
        rc = 1;
    } else if (max_rel <= repro_floor) {
        std::printf("PASS gate is exact within the engine's reproducibility floor "
                    "(gate %.3g <= floor %.3g)\n", max_rel, repro_floor);
    } else {
        std::printf("FAIL gate perturbs results beyond the engine floor "
                    "(gate %.3g > floor %.3g)\n", max_rel, repro_floor);
        rc = 1;
    }

    // Axis 3: harness dense output against STAR's own quant.sf NumReads.
    // Not expected to be exact: STAR updates effective lengths dynamically
    // during VB (its log records a GC update at iteration 11) whereas this
    // harness holds STAR's FINAL effective lengths fixed throughout. The
    // established acceptance band is the bar, not equality.
    {
        double sx = 0, sy = 0; std::size_t n = quant.size();
        for (std::size_t i = 0; i < n; ++i) { sx += dense.counts[i]; sy += quant[i].num_reads; }
        const double mx = sx / n, my = sy / n;
        double cov = 0, vx = 0, vy = 0, maxrel = 0; std::size_t diff = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const double x = dense.counts[i], y = quant[i].num_reads;
            cov += (x - mx) * (y - my); vx += (x - mx) * (x - mx); vy += (y - my) * (y - my);
            if (x != y) {
                ++diff;
                const double m = std::max(std::fabs(x), std::fabs(y));
                if (m > 0.0) maxrel = std::max(maxrel, std::fabs(x - y) / m);
            }
        }
        const double r = cov / (std::sqrt(vx) * std::sqrt(vy));
        std::printf("\naxis 3 - harness dense vs STAR quant.sf NumReads:\n"
                    "  transcripts differing %zu / %zu\n"
                    "  Pearson r             %.15f\n"
                    "  max rel delta         %.6g\n"
                    "  totals                %.3f vs %.3f\n",
                    diff, n, r, maxrel, sx, sy);
        const double band = 0.99998;   // chr22 TranscriptVB smoke NumReads band
        if (r >= band) std::printf("  PASS Pearson >= %.5f (chr22 smoke band)\n", band);
        else { std::printf("  FAIL Pearson < %.5f\n", band); rc = 1; }
    }
    // Axis 4: reconstruct effective lengths from the sidecar's own FLD state and
    // compare against STAR's EffectiveLength column. This is the piece axis 3
    // borrows from STAR, so quantifying it directly says how much of STAR's
    // effective length the sidecar alone can reproduce.
    //
    // Only the FLD term is reconstructed here. STAR additionally applies a GC
    // correction, which needs transcript sequences this harness does not load,
    // so a residual gap is expected and its size is the point of the measurement.
    {
        FLDAccumulator fld;
        std::string ferr;
        if (!fld.restore(gathered.evidence.fld_state, ferr)) {
            std::printf("\naxis 4 - SKIP: cannot restore FLD state: %s\n", ferr.c_str());
        } else {
            std::vector<double> pmf = fld.getPMF();
            std::vector<std::int32_t> raw_len(quant.size());
            for (std::size_t i = 0; i < quant.size(); ++i)
                raw_len[i] = static_cast<std::int32_t>(quant[i].length);
            std::vector<double> eff = computeEffectiveLengthsFromPMFWrapper(pmf, raw_len);
            if (eff.size() != quant.size()) {
                std::printf("\naxis 4 - SKIP: size mismatch %zu vs %zu\n",
                            eff.size(), quant.size());
            } else {
                double sx=0, sy=0; std::size_t n=quant.size(), cmp=0;
                for (std::size_t i=0;i<n;++i){ sx+=eff[i]; sy+=quant[i].eff_length; }
                const double mx=sx/n, my=sy/n;
                double cov=0, vx=0, vy=0, maxrel=0, sumabs=0;
                for (std::size_t i=0;i<n;++i){
                    const double x=eff[i], y=quant[i].eff_length;
                    cov+=(x-mx)*(y-my); vx+=(x-mx)*(x-mx); vy+=(y-my)*(y-my);
                    const double m=std::max(std::fabs(x),std::fabs(y));
                    if (m>0.0){ const double rel=std::fabs(x-y)/m;
                                maxrel=std::max(maxrel,rel); sumabs+=std::fabs(x-y); ++cmp; }
                }
                std::printf("\naxis 4 - FLD-only effective lengths vs STAR EffectiveLength:\n"
                            "  fld pmf support       %zu\n"
                            "  transcripts compared  %zu / %zu\n"
                            "  Pearson r             %.15f\n"
                            "  mean abs delta        %.6g nt\n"
                            "  max rel delta         %.6g\n"
                            "  mean eff length       %.3f (harness) vs %.3f (STAR)\n",
                            pmf.size(), cmp, n,
                            cov/(std::sqrt(vx)*std::sqrt(vy)),
                            cmp? sumabs/cmp : 0.0, maxrel, mx, my);
            }
        }
    }
    // Axis 4b: add the GC term, using the sidecar's own observed GC counts and
    // the transcript sequences. This distinguishes the two candidates for the
    // axis-4 gap: a missing GC correction, or an FLD reconstruction that
    // diverges from STAR's internal use of the same state.
    if (!tx_fasta.empty()) {
        FLDAccumulator fld; std::string ferr;
        if (fld.restore(gathered.evidence.fld_state, ferr)) {
            libem::Transcriptome txome;
            if (!txome.loadFromFasta(tx_fasta)) {
                std::printf("\naxis 4b - SKIP: cannot load %s\n", tx_fasta.c_str());
            } else {
                std::vector<double> pmf = fld.getPMF();
                std::vector<std::int32_t> raw_len(quant.size());
                for (std::size_t i = 0; i < quant.size(); ++i)
                    raw_len[i] = static_cast<std::int32_t>(quant[i].length);
                std::vector<double> eff_in = computeEffectiveLengthsFromPMFWrapper(pmf, raw_len);
                std::vector<double> obs_gc(gathered.evidence.gc_counts.begin(),
                                           gathered.evidence.gc_counts.end());
                DynamicGCEffectiveLengthResult res =
                    computeDynamicGCBiasedEffectiveLengthsWrapper(
                        txome, pmf, raw_len, dense.counts, eff_in, obs_gc);
                const std::vector<double>& eff = res.effective_lengths;
                if (!res.applied || eff.size() != quant.size()) {
                    std::printf("\naxis 4b - GC update applied=%d, size %zu vs %zu\n",
                                (int)res.applied, eff.size(), quant.size());
                } else {
                    double sx=0, sy=0; std::size_t n=quant.size(), cmp=0;
                    for (std::size_t i=0;i<n;++i){ sx+=eff[i]; sy+=quant[i].eff_length; }
                    const double mx=sx/n, my=sy/n;
                    double cov=0, vx=0, vy=0, maxrel=0, sumabs=0;
                    for (std::size_t i=0;i<n;++i){
                        const double x=eff[i], y=quant[i].eff_length;
                        cov+=(x-mx)*(y-my); vx+=(x-mx)*(x-mx); vy+=(y-my)*(y-my);
                        const double m=std::max(std::fabs(x),std::fabs(y));
                        if (m>0.0){ maxrel=std::max(maxrel,std::fabs(x-y)/m);
                                    sumabs+=std::fabs(x-y); ++cmp; }
                    }
                    std::printf("\naxis 4b - FLD+GC effective lengths vs STAR EffectiveLength:\n"
                                "  Pearson r             %.15f\n"
                                "  mean abs delta        %.6g nt\n"
                                "  max rel delta         %.6g\n"
                                "  mean eff length       %.3f (harness) vs %.3f (STAR)\n",
                                cov/(std::sqrt(vx)*std::sqrt(vy)),
                                cmp? sumabs/cmp : 0.0, maxrel, mx, my);
                }
            }
        }
    }
    // Axis 5: the fully traced reproduction of STAR's TranscriptVB path.
    //
    // Traced from STAR.cpp:2551-2915. Two defects in earlier attempts, both
    // found only by reading the source rather than approximating it:
    //   1. the sidecar is written from getECTable() BEFORE finalize(), so its
    //      EC weights are unnormalised and must be normalised before VB;
    //   2. the transcript index space is transcriptInfo.tab / quant.sf order,
    //      NOT transcriptome FASTA order -- same 199,138 names, different
    //      permutation -- so sequences must be reordered by name.
    // TranscriptQuantEC is not needed: finalize() is pure EC-weight
    // normalisation and importEvidence only restores FLD/GC state.
    if (!tx_fasta.empty()) {
        std::vector<std::string> index_names(quant.size());
        std::vector<std::int32_t> raw_len(quant.size());
        for (std::size_t i = 0; i < quant.size(); ++i) {
            index_names[i] = quant[i].name;
            raw_len[i] = static_cast<std::int32_t>(quant[i].length);
        }
        libem::Transcriptome txome;
        if (!txome.loadFromFasta(tx_fasta)) {
            std::printf("\naxis 5 - SKIP: cannot load %s\n", tx_fasta.c_str());
        } else if (!txome.reorderByNames(index_names)) {
            std::printf("\naxis 5 - FAIL: cannot reorder FASTA to index order\n");
            rc = 1;
        } else {
            std::printf("\naxis 5 - transcriptome reordered to index order (%zu)\n",
                        txome.size());
            FLDAccumulator fld; std::string ferr;
            if (!fld.restore(gathered.evidence.fld_state, ferr)) {
                std::printf("  SKIP: fld restore: %s\n", ferr.c_str());
            } else {
                std::vector<double> pmf = fld.getPMF();
                std::vector<double> obs_gc(gathered.evidence.gc_counts.begin(),
                                           gathered.evidence.gc_counts.end());
                std::vector<double> eff0 =
                    computeEffectiveLengthsFromPMFWrapper(pmf, raw_len);
                TranscriptState st; fillState(st, all);
                for (std::size_t i = 0; i < st.n; ++i) st.eff_lengths[i] = eff0[i];
                EMParams gcp = params;
                int applied_at = -1;
                gcp.effective_length_update =
                    [&](std::uint32_t it, TranscriptState& cb,
                        const std::vector<double>& a) -> bool {
                        DynamicGCEffectiveLengthResult r =
                            computeDynamicGCBiasedEffectiveLengthsWrapper(
                                txome, pmf, raw_len, a, cb.eff_lengths, obs_gc);
                        if (!r.applied) return false;
                        cb.eff_lengths.swap(r.effective_lengths);
                        applied_at = static_cast<int>(it);
                        return true;
                    };
                ECTable e = ecs;   // already normalised by finalizeECWeights
                EMResult r = run_vb(e, st, gcp);
                std::printf("  GC applied at iteration %d, converged=%s, iterations=%u\n",
                            applied_at, r.converged ? "yes" : "no", r.iterations);
                auto report = [&](const char* label, const std::vector<double>& x,
                                  bool use_eff) {
                    std::size_t n = quant.size(); double sx=0, sy=0;
                    for (std::size_t i=0;i<n;++i){
                        sx += x[i]; sy += use_eff ? quant[i].eff_length : quant[i].num_reads; }
                    const double mx=sx/n, my=sy/n;
                    double cov=0,vx=0,vy=0,sab=0;
                    for (std::size_t i=0;i<n;++i){
                        const double a=x[i];
                        const double b= use_eff ? quant[i].eff_length : quant[i].num_reads;
                        cov+=(a-mx)*(b-my); vx+=(a-mx)*(a-mx); vy+=(b-my)*(b-my);
                        sab+=std::fabs(a-b);
                    }
                    std::printf("  %-22s Pearson %.15f  mean abs %.6g  totals %.3f vs %.3f\n",
                                label, cov/(std::sqrt(vx)*std::sqrt(vy)), sab/n, sx, sy);
                };
                report("effective lengths:", st.eff_lengths, true);
                report("NumReads:", r.counts, false);
            }
        }
    }
    // Per-component convergence experiment.
    //
    // The EC-transcript graph decomposes into ~24k independent components, so
    // each is a separate VB problem. The global convergence check makes every
    // component iterate until the slowest one settles; per component they stop
    // on their own schedule. This measures how much iteration work that saves.
    //
    // NOTE: initialization is NOT matched. Salmon derives uniformPrior from
    // global totalWeight/numActive; running components in isolation changes it.
    // So these results are not parity-comparable to quant.sf -- the purpose here
    // is the iteration histogram, not the counts.
    {
        const std::size_t n_tx = ecs.n_transcripts;
        const double t_uf = nowSec();
        std::vector<std::size_t> parent(n_tx);
        for (std::size_t i = 0; i < n_tx; ++i) parent[i] = i;
        std::function<std::size_t(std::size_t)> find =
            [&](std::size_t x){ while(parent[x]!=x){ parent[x]=parent[parent[x]]; x=parent[x]; } return x; };
        for (const auto& ec : ecs.ecs)
            for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k) {
                std::size_t a=find(ec.transcript_ids[0]), b=find(ec.transcript_ids[k]);
                if (a!=b) parent[a]=b;
            }
        std::unordered_map<std::size_t, std::vector<std::size_t>> comp_ecs;
        for (std::size_t e = 0; e < ecs.ecs.size(); ++e)
            if (!ecs.ecs[e].transcript_ids.empty())
                comp_ecs[find(ecs.ecs[e].transcript_ids[0])].push_back(e);
        const double t_comp = nowSec() - t_uf;

        std::printf("\nper-component VB over %zu components\n", comp_ecs.size());
        std::printf("  component decomposition (one-time): %.1f ms\n", 1000.0 * t_comp);
        std::map<int,std::size_t> hist;
        std::uint64_t work_pc = 0, work_global = 0;
        std::size_t max_iter = 0;
        const double t0 = nowSec();
        for (const auto& kv : comp_ecs) {
            // local transcript space for this component
            std::vector<std::uint32_t> local;
            for (std::size_t e : kv.second)
                for (std::uint32_t t : ecs.ecs[e].transcript_ids) local.push_back(t);
            std::sort(local.begin(), local.end());
            local.erase(std::unique(local.begin(), local.end()), local.end());
            std::unordered_map<std::uint32_t,std::uint32_t> remap;
            for (std::uint32_t i = 0; i < local.size(); ++i) remap[local[i]] = i;

            ECTable sub; sub.n_transcripts = local.size();
            sub.ecs.reserve(kv.second.size());
            for (std::size_t e : kv.second) {
                EC c = ecs.ecs[e];
                for (auto& t : c.transcript_ids) t = remap[t];
                sub.ecs.push_back(std::move(c));
            }
            sub.n_ecs = sub.ecs.size();
            TranscriptState st; st.resize(local.size());
            for (std::size_t i = 0; i < local.size(); ++i) {
                const QuantRow& r = quant[local[i]];
                st.names[i]=r.name; st.lengths[i]=r.length;
                st.eff_lengths[i]= r.eff_length>0.0 ? r.eff_length : r.length;
                st.abundances[i]=0.0; st.counts[i]=0.0;
            }
            EMParams pp = params; pp.threads = 1;
            EMResult r = run_vb(sub, st, pp);
            int bucket = r.iterations <= 10 ? 10 : r.iterations <= 100 ? 100
                       : r.iterations <= 500 ? 500 : r.iterations <= 1000 ? 1000 : 2000;
            ++hist[bucket];
            work_pc     += (std::uint64_t)r.iterations * sub.ecs.size();
            work_global += (std::uint64_t)dense.iterations * sub.ecs.size();
            max_iter = std::max(max_iter, (std::size_t)r.iterations);
        }
        const double dt = nowSec() - t0;
        std::printf("  iteration histogram (components by iterations to converge):\n");
        for (const auto& kv : hist)
            std::printf("    <=%-5d %8zu components\n", kv.first, kv.second);
        std::printf("  slowest component            : %zu iterations\n", max_iter);
        std::printf("  EC-iterations, per-component : %llu\n", (unsigned long long)work_pc);
        std::printf("  EC-iterations, global (%u it): %llu\n",
                    dense.iterations, (unsigned long long)work_global);
        std::printf("  reduction in iteration work  : %.2fx\n",
                    work_global ? double(work_global)/double(work_pc) : 0.0);
        std::printf("  wall for all components, 1 thread each: %.1f s\n", dt);
    }
    // Is the static-range architecture worth building?
    //
    // Simulates it without the file layout: partition components into W equal
    // cost ranges, give each range the global effective lengths (so no GC
    // rendezvous is needed), and solve each range over ONLY its own transcripts.
    // That is the point of the design -- a range's per-iteration O(n) passes
    // become O(n_range), which is what the Lustre benchmark showed the global
    // run cannot escape.
    //
    // Iteration count is pinned identical for every arm so this measures the
    // architecture, not a different stopping rule.
    {
        const std::size_t n_tx = ecs.n_transcripts;
        const uint32_t FIXED = 300;   // same for global and ranges
        std::vector<std::size_t> par(n_tx);
        for (std::size_t i = 0; i < n_tx; ++i) par[i] = i;
        std::function<std::size_t(std::size_t)> fnd =
            [&](std::size_t x){ while(par[x]!=x){par[x]=par[par[x]];x=par[x];} return x; };
        for (const auto& ec : ecs.ecs)
            for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k) {
                std::size_t a=fnd(ec.transcript_ids[0]), b=fnd(ec.transcript_ids[k]);
                if (a!=b) par[a]=b;
            }
        std::unordered_map<std::size_t,std::vector<std::size_t>> comp_ecs;
        for (std::size_t e = 0; e < ecs.ecs.size(); ++e)
            if (!ecs.ecs[e].transcript_ids.empty())
                comp_ecs[fnd(ecs.ecs[e].transcript_ids[0])].push_back(e);

        EMParams fx = params;
        fx.min_iters = FIXED; fx.max_iters = FIXED; fx.num_active_override = n_tx;

        // global baseline, single thread, fixed iterations
        {
            TranscriptState st; fillState(st, all);
            ECTable e = ecs;
            EMParams g = fx; g.threads = 1;
            const double t0 = nowSec();
            run_vb(e, st, g);
            std::printf("\nstatic-range architecture, %u fixed iterations\n", FIXED);
            std::printf("  global, 1 thread, all %zu transcripts : %.2f s\n",
                        n_tx, nowSec()-t0);
        }

        for (int W : {8, 16}) {
            // greedy equal-cost assignment of components to ranges
            std::vector<std::pair<std::size_t,std::size_t>> cs;  // (cost, root)
            std::size_t total = 0;
            for (auto& kv : comp_ecs) {
                std::size_t c = 0;
                for (std::size_t e : kv.second) c += ecs.ecs[e].transcript_ids.size();
                cs.emplace_back(c, kv.first); total += c;
            }
            std::sort(cs.rbegin(), cs.rend());
            std::vector<std::vector<std::size_t>> range_roots(W);
            std::vector<std::size_t> load(W,0);
            for (auto& pr : cs) {
                const int w = std::min_element(load.begin(), load.end()) - load.begin();
                range_roots[w].push_back(pr.second); load[w] += pr.first;
            }
            double worst = 0.0, sum = 0.0;
            std::size_t worst_tx = 0;
            for (int w = 0; w < W; ++w) {
                std::vector<std::size_t> ecidx;
                for (std::size_t r : range_roots[w])
                    for (std::size_t e : comp_ecs[r]) ecidx.push_back(e);
                std::vector<std::uint32_t> local;
                for (std::size_t e : ecidx)
                    for (std::uint32_t t : ecs.ecs[e].transcript_ids) local.push_back(t);
                std::sort(local.begin(), local.end());
                local.erase(std::unique(local.begin(), local.end()), local.end());
                std::unordered_map<std::uint32_t,std::uint32_t> rm;
                for (std::uint32_t i = 0; i < local.size(); ++i) rm[local[i]] = i;
                ECTable sub; sub.n_transcripts = local.size();
                for (std::size_t e : ecidx) {
                    EC c = ecs.ecs[e];
                    for (auto& t : c.transcript_ids) t = rm[t];
                    sub.ecs.push_back(std::move(c));
                }
                sub.n_ecs = sub.ecs.size();
                TranscriptState st; st.resize(local.size());
                for (std::size_t i = 0; i < local.size(); ++i) {
                    const QuantRow& r = quant[local[i]];
                    st.names[i]=r.name; st.lengths[i]=r.length;
                    st.eff_lengths[i]= r.eff_length>0.0 ? r.eff_length : r.length;
                    st.abundances[i]=0.0; st.counts[i]=0.0;
                }
                EMParams rp = fx; rp.threads = 1;
                const double t0 = nowSec();
                run_vb(sub, st, rp);
                const double dt = nowSec()-t0;
                sum += dt;
                if (dt > worst) { worst = dt; worst_tx = local.size(); }
            }
            std::printf("  %2d ranges, 1 thread each: slowest %.2f s (%zu transcripts),"
                        " total work %.2f s\n", W, worst, worst_tx, sum);
        }
    }
    // Cost of the global setup pass (projected counts, cluster bounds, alpha
    // initialization) versus one VB iteration.
    //
    // If setup is cheap, the gather can compute alpha globally and write each
    // range's slice into its file, so workers start already correctly seeded and
    // never need global state. That removes the initialization coupling that
    // made isolated components disagree (the num_active_override failure).
    {
        auto timed = [&](uint32_t iters) {
            TranscriptState st; fillState(st, all);
            ECTable e = ecs;
            EMParams q = params; q.threads = 1;
            q.min_iters = iters; q.max_iters = iters;
            const double t0 = nowSec();
            run_vb(e, st, q);
            return nowSec() - t0;
        };
        const double t1  = timed(1);
        const double t11 = timed(11);
        const double per_iter = (t11 - t1) / 10.0;
        const double setup = t1 - per_iter;
        std::printf("\nglobal setup pass vs iteration cost (1 thread)\n");
        std::printf("  1 iteration total     : %.3f s\n", t1);
        std::printf("  11 iterations total   : %.3f s\n", t11);
        std::printf("  => per-iteration       : %.4f s\n", per_iter);
        std::printf("  => setup (projected counts, cluster bounds, alpha init): %.3f s\n", setup);
        std::printf("  setup as %% of a 1863-iteration run: %.2f%%\n",
                    100.0 * setup / (setup + 1863.0 * per_iter));
    }
    return rc;
}
