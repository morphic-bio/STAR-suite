// Standalone sidecar gather parity harness.
//
// Reads N TranscriptVB shard sidecars, merges them, and reports the merged
// evidence, using only TranscriptVBSidecar + libem. It links no Genome, no
// suffix array and no STAR index, which is the point: it demonstrates that the
// gather half of the finalizer runs without the 29.7 GB index that the current
// in-STAR finalizer loads and then frees.
//
// This is the first parity milestone for a standalone VB/EM finalizer
// (docs/STANDALONE_VBEM_FINALIZER_20260807.md). It checks the merge, not the
// quantification: compare the reported figures against what STAR logged for the
// same sidecars, e.g.
//
//     TranscriptVB gathered 8 complete shard sidecars (pairs=50917353)
//     Merged 232633 equivalence classes from 32 threads
//
// usage: sidecar_gather_parity [--expect-pairs N] [--expect-ecs N] SIDECAR...

#include "TranscriptVBSidecar.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <set>
#include <numeric>
#include <fstream>

namespace {

std::uint64_t sumGC(const transcript_vb_sidecar::Evidence& e) {
    double total = 0.0;
    for (double v : e.gc_counts) total += v;
    return static_cast<std::uint64_t>(total);
}

int fail(const std::string& what) {
    std::cerr << "FAIL " << what << "\n";
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> paths;
    std::string tx_info;
    std::uint64_t expect_pairs = 0, expect_ecs = 0;
    bool have_pairs = false, have_ecs = false;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--expect-pairs") == 0 && i + 1 < argc) {
            expect_pairs = std::strtoull(argv[++i], nullptr, 10);
            have_pairs = true;
        } else if (std::strcmp(argv[i], "--expect-ecs") == 0 && i + 1 < argc) {
            expect_ecs = std::strtoull(argv[++i], nullptr, 10);
            have_ecs = true;
        } else if (std::strcmp(argv[i], "--tx-info") == 0 && i + 1 < argc) {
            tx_info = argv[++i];
        } else {
            paths.emplace_back(argv[i]);
        }
    }
    if (paths.empty()) {
        std::cerr << "usage: sidecar_gather_parity [--expect-pairs N] "
                     "[--expect-ecs N] SIDECAR...\n";
        return 2;
    }

    std::vector<transcript_vb_sidecar::Sidecar> shards;
    shards.reserve(paths.size());
    std::string err;
    for (const std::string& p : paths) {
        transcript_vb_sidecar::Sidecar s;
        if (!transcript_vb_sidecar::read(p, s, err))
            return fail("read " + p + ": " + err);
        std::cout << "shard ordinal=" << s.metadata.shard_ordinal
                  << " of " << s.metadata.shard_count
                  << "  first_pair=" << s.metadata.first_pair
                  << "  pairs=" << s.metadata.pair_count
                  << "  ecs=" << s.evidence.ec_table.ecs.size()
                  << "  transcripts=" << s.metadata.transcript_count << "\n";
        shards.push_back(std::move(s));
    }

    // Every ordinal exactly once, which is the invariant the gather relies on.
    const std::uint32_t n = shards.front().metadata.shard_count;
    std::vector<int> seen(n, 0);
    for (const auto& s : shards) {
        if (s.metadata.shard_count != n)
            return fail("shard_count disagrees across sidecars");
        if (s.metadata.shard_ordinal >= n)
            return fail("shard ordinal out of range");
        if (seen[s.metadata.shard_ordinal]++)
            return fail("duplicate shard ordinal");
    }
    for (std::uint32_t i = 0; i < n; ++i)
        if (!seen[i]) return fail("missing shard ordinal");
    if (shards.size() != n)
        return fail("sidecar count does not equal shard_count");
    std::cout << "PASS ordinals 0.." << (n - 1) << " each present exactly once\n";

    transcript_vb_sidecar::Sidecar gathered;
    if (!transcript_vb_sidecar::merge(shards, gathered, err))
        return fail("merge: " + err);

    const auto& ev = gathered.evidence;
    std::cout << "\nmerged:\n"
              << "  pairs                = " << gathered.metadata.pair_count << "\n"
              << "  equivalence classes  = " << ev.ec_table.ecs.size() << "\n"
              << "  ec_table.n_ecs       = " << ev.ec_table.n_ecs << "\n"
              << "  transcripts          = " << ev.ec_table.n_transcripts << "\n"
              << "  processed_fragments  = " << ev.processed_fragments << "\n"
              << "  gc observations      = " << sumGC(ev) << "\n"
              << "  dropped_incompat     = " << ev.dropped_incompat << "\n"
              << "  dropped_missing_mate = " << ev.dropped_missing_mate_fields << "\n"
              << "  dropped_unknown_fmt  = " << ev.dropped_unknown_obs_fmt << "\n"
              << "  library_format_id    = " << gathered.metadata.library_format_id << "\n";

    // Non-zero support census.
    //
    // A transcript appearing in no equivalence class contributes nothing to the
    // E-step and receives nothing from it, so its expected count stays 0 and it
    // is truncated to 0 on output regardless. Excluding such transcripts from
    // the dense per-thread accumulation buffer is therefore exact, not an
    // approximation: no mass moves and no EC is emptied, so the downstream
    // dropping mechanics are untouched.
    //
    // The dense buffer is num_threads * n_transcripts doubles, zeroed and
    // summed every iteration, so the saving is proportional to the fraction of
    // transcripts with no support.
    {
        const std::size_t n_tx = ev.ec_table.n_transcripts;
        std::vector<std::uint32_t> ec_occurrences(n_tx, 0);
        std::vector<double> global_weight(n_tx, 0.0);
        std::uint64_t incidence = 0;
        for (const auto& ec : ev.ec_table.ecs) {
            const std::size_t g = ec.transcript_ids.size();
            for (std::size_t i = 0; i < g; ++i) {
                const std::uint32_t tid = ec.transcript_ids[i];
                if (tid >= n_tx) continue;
                ++ec_occurrences[tid];
                const double w = ec.weights.empty() ? 1.0 / static_cast<double>(g)
                                                    : ec.weights[i];
                global_weight[tid] += ec.count * w;
            }
            incidence += g;
        }
        std::size_t supported = 0, positive_weight = 0;
        for (std::size_t i = 0; i < n_tx; ++i) {
            if (ec_occurrences[i] > 0) ++supported;
            if (global_weight[i] > 0.0) ++positive_weight;
        }
        const double frac = n_tx ? 100.0 * static_cast<double>(supported) /
                                       static_cast<double>(n_tx) : 0.0;
        std::cout << "\nnon-zero support census (gate = appears in >=1 EC):\n"
                  << "  transcripts total            = " << n_tx << "\n"
                  << "  with >=1 EC                  = " << supported
                  << "  (" << frac << "%)\n"
                  << "  with global weight > 0       = " << positive_weight << "\n"
                  << "  no support (exactly zero)    = " << (n_tx - supported) << "\n"
                  << "  total EC-transcript incidence= " << incidence << "\n"
                  << "  mean transcripts per EC      = "
                  << (ev.ec_table.ecs.empty() ? 0.0
                        : static_cast<double>(incidence) /
                          static_cast<double>(ev.ec_table.ecs.size())) << "\n";
        for (int threads : {8, 32, 256}) {
            const double dense = static_cast<double>(threads) * n_tx * 8.0 / (1024*1024);
            const double gated = static_cast<double>(threads) * supported * 8.0 / (1024*1024);
            std::cout << "  reduction buffer @" << threads << " threads: "
                      << dense << " MiB dense -> " << gated << " MiB gated\n";
        }
    }

    // Dependency structure: connected components of the EC-transcript graph.
    //
    // Two equivalence classes interact only if they share a transcript, so the
    // connected components of the bipartite EC/transcript graph are fully
    // independent VB problems. If the graph decomposes, a distributed VB needs
    // no per-iteration exchange at all -- workers take whole components and run
    // to convergence locally. If one giant component dominates, it does.
    {
        const std::size_t n_tx = ev.ec_table.n_transcripts;
        const std::size_t n_ec = ev.ec_table.ecs.size();
        std::vector<std::size_t> parent(n_tx);
        for (std::size_t i = 0; i < n_tx; ++i) parent[i] = i;
        std::function<std::size_t(std::size_t)> find =
            [&](std::size_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
        auto unite = [&](std::size_t a, std::size_t b) {
            a = find(a); b = find(b); if (a != b) parent[a] = b;
        };
        for (const auto& ec : ev.ec_table.ecs)
            for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k)
                unite(ec.transcript_ids[0], ec.transcript_ids[k]);

        // component sizes over SUPPORTED transcripts, and EC counts per component
        std::vector<char> sup(n_tx, 0);
        for (const auto& ec : ev.ec_table.ecs)
            for (std::uint32_t t : ec.transcript_ids) if (t < n_tx) sup[t] = 1;
        std::unordered_map<std::size_t, std::size_t> tx_per, ec_per;
        std::unordered_map<std::size_t, double> mass_per;
        for (std::size_t i = 0; i < n_tx; ++i) if (sup[i]) ++tx_per[find(i)];
        double total_mass = 0.0;
        for (const auto& ec : ev.ec_table.ecs) {
            if (ec.transcript_ids.empty()) continue;
            const std::size_t r = find(ec.transcript_ids[0]);
            ++ec_per[r]; mass_per[r] += ec.count; total_mass += ec.count;
        }
        std::vector<std::pair<std::size_t,std::size_t>> comps;  // (ecs, root)
        for (const auto& kv : ec_per) comps.emplace_back(kv.second, kv.first);
        std::sort(comps.rbegin(), comps.rend());

        std::printf("\ndependency structure: connected components of the EC-transcript graph\n");
        std::printf("  components (with >=1 EC) = %zu\n", comps.size());
        std::printf("  largest components:\n");
        std::printf("    %-10s %-12s %-14s %s\n", "ECs", "transcripts", "EC count mass", "%% of mass");
        for (std::size_t i = 0; i < comps.size() && i < 6; ++i) {
            const std::size_t r = comps[i].second;
            std::printf("    %-10zu %-12zu %-14.0f %.2f%%\n",
                        comps[i].first, tx_per[r], mass_per[r],
                        100.0 * mass_per[r] / total_mass);
        }
        // how much sits outside the single largest component
        if (!comps.empty()) {
            const std::size_t big = comps[0].second;
            std::size_t rest_ec = 0; double rest_mass = 0.0;
            for (const auto& kv : ec_per) if (kv.first != big) rest_ec += kv.second;
            for (const auto& kv : mass_per) if (kv.first != big) rest_mass += kv.second;
            std::printf("  outside the largest component: %zu ECs (%.1f%%), %.2f%% of mass\n",
                        rest_ec, 100.0 * rest_ec / n_ec, 100.0 * rest_mass / total_mass);
        }
    }

    // Are the observed components a refinement of gene boundaries?
    //
    // An EC edge requires a read compatible with both transcripts, hence shared
    // sequence, so the observed graph is a subgraph of a static sequence
    // compatibility graph. Gene assignment is a static partition already stored
    // in transcriptInfo.tab. If components sit inside genes, genes are a ready
    // made precomputed partition; if components span genes, they do not.
    if (!tx_info.empty()) {
        std::ifstream ti(tx_info);
        std::size_t nTr = 0;
        if (ti >> nTr) {
            std::vector<std::uint32_t> gene(nTr);
            bool ok = true;
            for (std::size_t i = 0; i < nTr && ok; ++i) {
                std::string nm; std::uint64_t a,b,c; unsigned d,e,f,g;
                if (!(ti >> nm >> a >> b >> c >> d >> e >> f >> g)) ok = false;
                else gene[i] = g;
            }
            const std::size_t n_tx = ev.ec_table.n_transcripts;
            if (ok && nTr == n_tx) {
                std::vector<std::size_t> par(n_tx);
                for (std::size_t i = 0; i < n_tx; ++i) par[i] = i;
                std::function<std::size_t(std::size_t)> fnd =
                    [&](std::size_t x){ while(par[x]!=x){par[x]=par[par[x]];x=par[x];} return x; };
                for (const auto& ec : ev.ec_table.ecs)
                    for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k) {
                        std::size_t A=fnd(ec.transcript_ids[0]), B=fnd(ec.transcript_ids[k]);
                        if (A!=B) par[A]=B;
                    }
                std::unordered_map<std::size_t, std::set<std::uint32_t>> genes_in_comp;
                std::unordered_map<std::uint32_t, std::set<std::size_t>> comps_in_gene;
                std::vector<char> sup(n_tx, 0);
                for (const auto& ec : ev.ec_table.ecs)
                    for (std::uint32_t t : ec.transcript_ids) if (t < n_tx) sup[t]=1;
                for (std::size_t i = 0; i < n_tx; ++i) if (sup[i]) {
                    genes_in_comp[fnd(i)].insert(gene[i]);
                    comps_in_gene[gene[i]].insert(fnd(i));
                }
                std::size_t single_gene = 0, multi_gene = 0, max_genes = 0;
                for (const auto& kv : genes_in_comp) {
                    if (kv.second.size() == 1) ++single_gene; else ++multi_gene;
                    max_genes = std::max(max_genes, kv.second.size());
                }
                std::size_t split_genes = 0;
                for (const auto& kv : comps_in_gene) if (kv.second.size() > 1) ++split_genes;
                std::printf("\ncomponents vs gene boundaries (static partition candidate)\n");
                std::printf("  components confined to one gene : %zu\n", single_gene);
                std::printf("  components spanning >1 gene     : %zu\n", multi_gene);
                std::printf("  widest component                : %zu genes\n", max_genes);
                std::printf("  genes split across >1 component : %zu of %zu\n",
                            split_genes, comps_in_gene.size());
            } else {
                std::printf("\n(transcriptInfo.tab has %zu rows, evidence has %zu; skipped)\n",
                            nTr, n_tx);
            }
        }
    }

    // Static range partition feasibility.
    //
    // A Lustre-friendly layout writes components contiguously with an index, so
    // each worker pulls one byte range, solves it, and writes one byte range --
    // no scheduler, no crosstalk. That only works if fixed ranges balance
    // without dynamic stealing, so measure the imbalance directly.
    {
        const auto& t = ev.ec_table;
        const std::size_t n_tx = t.n_transcripts;
        std::vector<std::size_t> par(n_tx);
        for (std::size_t i = 0; i < n_tx; ++i) par[i] = i;
        std::function<std::size_t(std::size_t)> fnd =
            [&](std::size_t x){ while(par[x]!=x){par[x]=par[par[x]];x=par[x];} return x; };
        for (const auto& ec : t.ecs)
            for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k) {
                std::size_t a=fnd(ec.transcript_ids[0]), b=fnd(ec.transcript_ids[k]);
                if (a!=b) par[a]=b;
            }
        // cost proxy: EC-transcript incidences, which is what the E-step touches
        std::unordered_map<std::size_t,std::size_t> cost, ntx_of;
        for (const auto& ec : t.ecs)
            if (!ec.transcript_ids.empty())
                cost[fnd(ec.transcript_ids[0])] += ec.transcript_ids.size();
        std::vector<char> sup(n_tx,0);
        for (const auto& ec : t.ecs) for (auto x : ec.transcript_ids) if (x<n_tx) sup[x]=1;
        for (std::size_t i=0;i<n_tx;++i) if (sup[i]) ++ntx_of[fnd(i)];
        std::vector<std::size_t> costs;
        for (auto& kv : cost) costs.push_back(kv.second);
        std::sort(costs.rbegin(), costs.rend());
        const std::size_t total = std::accumulate(costs.begin(), costs.end(), (std::size_t)0);
        std::printf("\nstatic range partition (no scheduler)\n");
        std::printf("  components=%zu  total cost=%zu  largest single component=%zu (%.3f%%)\n",
                    costs.size(), total, costs.front(), 100.0*costs.front()/total);
        for (int W : {8, 16, 64}) {
            // greedy contiguous split of the size-sorted order into W equal-cost runs
            const double target = double(total)/W;
            std::vector<std::size_t> load(W,0);
            std::size_t w=0, acc=0;
            for (std::size_t c : costs) {
                if (w+1<(std::size_t)W && acc+c > target*(w+1)) ++w;
                load[w]+=c; acc+=c;
            }
            const std::size_t mx=*std::max_element(load.begin(),load.end());
            const std::size_t mn=*std::min_element(load.begin(),load.end());
            std::printf("    %2d ranges: max/mean = %.3f   min/mean = %.3f   (perfect = 1.000)\n",
                        W, mx/target, mn/target);
        }
        std::printf("  => a fixed range assignment is viable when max/mean is near 1\n");
    }

    // Representation cost of a shared EC signature dictionary.
    //
    // The distinct EC signature set is a property of the transcriptome and the
    // reads, not of any cell partition, so the bulk sidecar measures it
    // directly. Under a per-cell design the signatures are stored once globally
    // and each cell holds only (signature id, count) pairs.
    {
        const auto& t = ev.ec_table;
        std::size_t inc = 0, with_w = 0;
        for (const auto& ec : t.ecs) {
            inc += ec.transcript_ids.size();
            if (!ec.weights.empty()) with_w += ec.weights.size();
        }
        const double MB = 1024.0 * 1024.0;
        // current: struct-of-vectors, one EC object per entry
        const double cur_struct = t.ecs.size() * (double)sizeof(EC) / MB;
        const double cur_payload = (inc * 4.0 + with_w * 8.0) / MB;
        const double cur_allocs  = t.ecs.size() * 2.0 * 16.0 / MB;   // ~2 heap blocks each
        // compact CSR dictionary: offsets + ids + weights, no per-EC objects
        const double csr = ((t.ecs.size() + 1) * 4.0 + inc * 4.0 + with_w * 8.0) / MB;
        std::printf("\nEC representation cost (this dataset)\n");
        std::printf("  distinct signatures      = %zu\n", t.ecs.size());
        std::printf("  transcript incidences    = %zu\n", inc);
        std::printf("  current struct-of-vectors: %.1f MB struct + %.1f MB payload"
                    " + ~%.1f MB alloc overhead = %.1f MB\n",
                    cur_struct, cur_payload, cur_allocs,
                    cur_struct + cur_payload + cur_allocs);
        std::printf("  compact CSR dictionary   : %.1f MB (shared across all cells)\n", csr);
        std::printf("  per-cell cost would be   : 8 B per (signature id, count) entry\n");
        for (int cells : {1000, 10000, 100000}) {
            for (int per : {2000, 10000}) {
                const double tot = csr + (double)cells * per * 8.0 / MB;
                std::printf("    %6d cells x %5d ECs/cell -> %.2f GB\n",
                            cells, per, tot / 1024.0);
            }
        }
    }

    // Partition analysis for a distributed VB.
    //
    // A worker owning a subset of equivalence classes only ever writes the
    // transcripts those ECs touch, so its published partial can be sparse. How
    // small depends entirely on how the ECs are partitioned: a split that keeps
    // co-occurring transcripts together cuts the exchange, a split that
    // scatters them does not.
    {
        const std::size_t n_tx = ev.ec_table.n_transcripts;
        const std::size_t n_ec = ev.ec_table.ecs.size();
        std::printf("\ndistributed exchange: distinct transcripts touched per worker\n");
        std::printf("  %-8s %-12s %8s %8s %8s   %10s %10s\n",
                    "workers", "partition", "min", "mean", "max",
                    "dense MB", "sparse MB");
        for (int nw : {2, 4, 8, 16}) {
            for (int mode = 0; mode < 2; ++mode) {   // 0 = contiguous, 1 = strided
                std::vector<std::size_t> distinct(nw, 0);
                std::vector<char> touched(n_tx, 0);
                for (int w = 0; w < nw; ++w) {
                    std::fill(touched.begin(), touched.end(), 0);
                    for (std::size_t e = 0; e < n_ec; ++e) {
                        const int owner = (mode == 0)
                            ? static_cast<int>(e * nw / n_ec)   // contiguous blocks
                            : static_cast<int>(e % nw);          // round robin
                        if (owner != w) continue;
                        for (std::uint32_t t : ev.ec_table.ecs[e].transcript_ids)
                            if (t < n_tx) touched[t] = 1;
                    }
                    std::size_t c = 0;
                    for (std::size_t i = 0; i < n_tx; ++i) c += touched[i];
                    distinct[w] = c;
                }
                std::size_t mn = distinct[0], mx = distinct[0], sum = 0;
                for (std::size_t d : distinct) { mn = std::min(mn,d); mx = std::max(mx,d); sum += d; }
                const double mean = double(sum) / nw;
                // dense: every worker publishes the whole vector.
                // sparse: (index,value) pairs for touched transcripts only.
                const double dense_mb  = nw * n_tx * 8.0 / 1048576.0;
                const double sparse_mb = sum * 12.0 / 1048576.0;
                std::printf("  %-8d %-12s %8zu %8.0f %8zu   %10.2f %10.2f\n",
                            nw, mode == 0 ? "contiguous" : "round-robin",
                            mn, mean, mx, dense_mb, sparse_mb);
            }
        }
        std::printf("  (per iteration; x1863 iterations for the run total)\n");
    }

    int rc = 0;
    if (have_pairs) {
        if (gathered.metadata.pair_count == expect_pairs)
            std::cout << "PASS pairs == " << expect_pairs << "\n";
        else { std::cout << "FAIL pairs got " << gathered.metadata.pair_count
                         << " want " << expect_pairs << "\n"; rc = 1; }
    }
    if (have_ecs) {
        if (ev.ec_table.ecs.size() == expect_ecs)
            std::cout << "PASS equivalence classes == " << expect_ecs << "\n";
        else { std::cout << "FAIL equivalence classes got " << ev.ec_table.ecs.size()
                         << " want " << expect_ecs << "\n"; rc = 1; }
    }
    return rc;
}
