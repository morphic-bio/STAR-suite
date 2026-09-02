// transcriptvb_finalize - standalone TranscriptVB shard finalizer.
//
// Gathers N TranscriptVB shard sidecars and produces quant.sf without loading
// the STAR genome index. It reproduces the path traced from STAR.cpp:2551-2915;
// see docs/STANDALONE_VBEM_FINALIZER_20260807.md.
//
// Two things are non-obvious and were the source of earlier wrong results:
//
//   1. The sidecar is written from getECTable() BEFORE TranscriptQuantEC::
//      finalize(), so its EC weights are UNNORMALISED. Every consumer must
//      normalise them before VB or mass is conserved but its assignment is
//      scrambled.
//   2. The transcript index space is transcriptInfo.tab order, NOT
//      transcriptome.fa order. The FASTA holds the same transcripts in a
//      different permutation, so sequences must be reordered by name before the
//      GC correction is applied, or each correction lands on the wrong
//      transcript.
//
// Inputs are the transcriptome and the small annotation tables; Genome, SA and
// SAindex are never read.

#include "TranscriptVBSidecar.h"
#include "em_types.h"
#include "vb_engine.h"
#include "em_engine.h"
#include "fld_accumulator.h"
#include "effective_length_wrapper.h"
#include "alignment_model.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <functional>
#include <iomanip>
#include <sstream>

namespace {

double nowSec() {
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
}

// Verbatim replica of TranscriptQuantEC::finalize (TranscriptQuantEC.cpp:991).
void normaliseECWeights(ECTable& t) {
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

// Mirrors Transcriptome.cpp:135-176: names and gene ids from transcriptInfo.tab,
// lengths derived from exonInfo.tab.
struct TxIndex {
    std::vector<std::string> names;
    std::vector<std::int32_t> lengths;
    std::vector<std::uint32_t> gene;
    std::size_t n = 0;
};

bool loadTxIndex(const std::string& genome_dir, TxIndex& out, std::string& err) {
    std::ifstream tri(genome_dir + "/transcriptInfo.tab");
    if (!tri) { err = "cannot open transcriptInfo.tab"; return false; }
    std::size_t nTr = 0;
    if (!(tri >> nTr) || nTr == 0) { err = "bad transcriptInfo.tab header"; return false; }
    out.n = nTr;
    out.names.resize(nTr); out.lengths.assign(nTr, 0); out.gene.resize(nTr);
    std::vector<std::uint32_t> trExI(nTr), trExN(nTr);
    for (std::size_t i = 0; i < nTr; ++i) {
        std::uint64_t s, e, emax; unsigned str, exN, exI, gene;
        if (!(tri >> out.names[i] >> s >> e >> emax >> str >> exN >> exI >> gene)) {
            err = "transcriptInfo.tab truncated"; return false;
        }
        trExN[i] = exN; trExI[i] = exI; out.gene[i] = gene;
    }
    std::ifstream exi(genome_dir + "/exonInfo.tab");
    if (!exi) { err = "cannot open exonInfo.tab"; return false; }
    std::size_t nEx = 0;
    if (!(exi >> nEx)) { err = "bad exonInfo.tab header"; return false; }
    std::vector<std::uint32_t> exS(nEx), exE(nEx), exCum(nEx);
    for (std::size_t i = 0; i < nEx; ++i)
        if (!(exi >> exS[i] >> exE[i] >> exCum[i])) { err = "exonInfo.tab truncated"; return false; }
    for (std::size_t i = 0; i < nTr; ++i) {
        const std::size_t last = trExI[i] + trExN[i] - 1;
        if (last >= nEx) { err = "exon index out of range"; return false; }
        out.lengths[i] = static_cast<std::int32_t>(exCum[last] + exE[last] - exS[last] + 1);
    }
    return true;
}

// Verbatim replica of TranscriptQuantOutput.cpp's formatNum: fixed precision
// then trailing-zero trimming, matching Salmon's quant.genes.sf ("510" not
// "510.000000"). quant.sf uses plain fixed precision instead; the two files
// genuinely format differently.
std::string formatNum(double val, int maxPrecision = 6) {
    if (val == 0.0) return "0";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(maxPrecision) << val;
    std::string t = oss.str();
    const std::size_t dot = t.find('.');
    if (dot != std::string::npos) {
        const std::size_t last = t.find_last_not_of('0');
        if (last != std::string::npos && last > dot) t = t.substr(0, last + 1);
        else if (last == dot) t = t.substr(0, dot);
    }
    return t;
}

// geneInfo.tab: count, then "geneID geneName biotype" per line.
bool loadGeneNames(const std::string& genome_dir, std::vector<std::string>& names,
                   std::string& err) {
    std::ifstream gi(genome_dir + "/geneInfo.tab");
    if (!gi) { err = "cannot open geneInfo.tab"; return false; }
    std::size_t nGe = 0;
    if (!(gi >> nGe)) { err = "bad geneInfo.tab header"; return false; }
    names.resize(nGe);
    for (std::size_t i = 0; i < nGe; ++i) {
        std::string rest;
        if (!(gi >> names[i])) { err = "geneInfo.tab truncated"; return false; }
        std::getline(gi, rest);
    }
    return true;
}

// Mirrors writeQuantGeneSF (TranscriptQuantOutput.cpp:69). Length and effective
// length are TPM-weighted averages where a gene has TPM > 0, and uniform
// averages otherwise, which is what Salmon's aggregateEstimatesToGeneLevel does.
bool writeGeneQuant(const std::string& path,
                    const std::vector<std::string>& gene_names,
                    const std::vector<std::uint32_t>& tr_gene,
                    const std::vector<std::int32_t>& tr_len,
                    const std::vector<double>& eff_len,
                    const std::vector<double>& tpm,
                    const std::vector<double>& counts) {
    const std::size_t nGe = gene_names.size();
    std::vector<double> gTPM(nGe,0.0), gCnt(nGe,0.0),
                        gLenW(nGe,0.0), gEffW(nGe,0.0),
                        gLenU(nGe,0.0), gEffU(nGe,0.0);
    std::vector<std::uint32_t> gN(nGe, 0);
    for (std::size_t t = 0; t < tr_gene.size(); ++t) {
        const std::uint32_t g = tr_gene[t];
        if (g >= nGe) continue;
        gTPM[g] += tpm[t];
        gCnt[g] += counts[t];
        gLenW[g] += tpm[t] * tr_len[t];
        gEffW[g] += tpm[t] * eff_len[t];
        gLenU[g] += tr_len[t];
        gEffU[g] += eff_len[t];
        ++gN[g];
    }
    std::ofstream out(path);
    if (!out) return false;
    out << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
    for (std::size_t g = 0; g < nGe; ++g) {
        double len, efflen;
        if (gTPM[g] > 0)       { len = gLenW[g]/gTPM[g];  efflen = gEffW[g]/gTPM[g]; }
        else if (gN[g] > 0)    { len = gLenU[g]/gN[g];    efflen = gEffU[g]/gN[g];   }
        else                   { len = 0.0;               efflen = 0.0;              }
        out << gene_names[g] << '\t' << formatNum(len,3) << '\t'
            << formatNum(efflen,4) << '\t' << formatNum(gTPM[g],6) << '\t'
            << formatNum(gCnt[g],3) << '\n';
    }
    return true;
}

// Shard outputs in sorted name order. Sorting matters only for reproducible
// diagnostics -- the scatter itself is order-independent because every entry
// carries its global transcript id.
static std::vector<std::string> globShardOutputs(const std::string& dir) {
    std::vector<std::string> out;
    DIR* d = ::opendir(dir.c_str());
    if (d == nullptr) return out;
    while (struct dirent* e = ::readdir(d)) {
        const std::string n = e->d_name;
        if (n.size() > 4 && n.compare(n.size()-4, 4, ".bin") == 0)
            out.push_back(dir + "/" + n);
    }
    ::closedir(d);
    std::sort(out.begin(), out.end());
    return out;
}

void usage() {
    std::cerr <<
      "usage: transcriptvb_finalize --genome-dir DIR --transcriptome FA \\\n"
      "                            --out-prefix PREFIX [--threads N] [--sparse]\\\n"
      "                            SIDECAR...\n"
      "   split:    transcriptvb_finalize --emit-ranges W --range-dir DIR \\\n"
      "                            --genome-dir DIR --transcriptome FA \\\n"
      "                            --out-prefix PREFIX SIDECAR...\n"
      "   worker:   transcriptvb_finalize --range-file F --range-out O\n"
      "   assemble: transcriptvb_finalize --assemble --genome-dir DIR \\\n"
      "                            --range-dir DIR --shard-out DIR \\\n"
      "                            --out-prefix PREFIX\n"
      "  --genome-dir     directory holding transcriptInfo.tab and exonInfo.tab\n"
      "                   (Genome, SA and SAindex are never read)\n"
      "  --transcriptome  transcriptome FASTA, for GC sequence context\n"
      "\n"
      "Non-zero-support gating is applied inside libem automatically and is\n"
      "exact; there is no caller-side flag.\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string genome_dir, tx_fasta, out_prefix;
    int threads = 0;
    double tolerance = 0.0;   // 0 = engine default
    bool per_component = false;
    bool comp_partition = false;
    int emit_ranges = 0;
    std::string range_dir;
    std::string range_file;
    bool parallel_gc = false;        // opt in: trades byte-identity for ~12x
                                     // on the GC background pass
    bool assemble = false;           // gather shard outputs -> quant.sf
    std::string shard_out;
    bool no_gc = false;              // isolates the shard architecture from the
                                     // not-yet-sharded GC effective-length pass
    int range_iters = 0;
    std::string range_out;
    std::vector<std::string> sidecars;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--genome-dir")    genome_dir = next();
        else if (a == "--transcriptome") tx_fasta   = next();
        else if (a == "--out-prefix")    out_prefix = next();
        else if (a == "--threads")       threads    = std::atoi(next().c_str());
        else if (a == "--tolerance")     tolerance  = std::atof(next().c_str());
        else if (a == "--per-component") per_component = true;
        else if (a == "--component-parallel") comp_partition = true;
        else if (a == "--emit-ranges")   emit_ranges = std::atoi(next().c_str());
        else if (a == "--range-dir")     range_dir  = next();
        else if (a == "--range-file")    range_file = next();
        else if (a == "--range-iters")   range_iters = std::atoi(next().c_str());
        else if (a == "--range-out")     range_out  = next();
        else if (a == "--no-gc")         no_gc = true;
        else if (a == "--parallel-gc")   parallel_gc = true;
        else if (a == "--assemble")      assemble = true;
        else if (a == "--shard-out")     shard_out = next();
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else sidecars.push_back(a);
    }
    // ---- range worker ----------------------------------------------------
    // Consume one self-contained range file: read it contiguously, rebuild the
    // local EC table and transcript state from it, and iterate. No global state,
    // no other input.
    if (!range_file.empty()) {
        const double t0 = nowSec();
        std::ifstream in(range_file, std::ios::binary);
        if (!in) { std::cerr << "EXITING: cannot read " << range_file << "\n"; return 1; }
        auto R32=[&]{ std::uint32_t v; in.read((char*)&v,4); return v; };
        auto R64=[&]{ std::uint64_t v; in.read((char*)&v,8); return v; };
        auto RD =[&]{ double v; in.read((char*)&v,8); return v; };
        if (R32() != 0x52425654u) { std::cerr << "EXITING: bad range magic\n"; return 1; }
        R32();
        const std::uint32_t nR = R32(), idx = R32();
        const std::uint64_t nTxL = R64(), nEcL = R64(), nTxG = R64();
        EMParams rp; rp.use_vb = true; rp.threads = 1;
        rp.vb_prior = RD(); rp.tolerance = RD();
        rp.min_iters = R32(); rp.max_iters = R32();
        rp.num_active_override = (std::size_t)nTxG;
        std::vector<std::uint32_t> l2g(nTxL);
        for (auto& v : l2g) v = R32();
        std::vector<double> a0(nTxL), efl(nTxL);
        for (auto& v : a0)  v = RD();
        for (auto& v : efl) v = RD();
        std::vector<std::uint32_t> rawl(nTxL);
        for (auto& v : rawl) v = R32();
        const std::uint64_t nnz = R64();
        std::vector<std::uint32_t> off(nEcL + 1);
        for (auto& v : off) v = R32();
        std::vector<std::uint32_t> ids(nnz);
        for (auto& v : ids) v = R32();
        std::vector<double> wts(nnz), cnts(nEcL);
        for (auto& v : wts)  v = RD();
        for (auto& v : cnts) v = RD();
        if (!in) { std::cerr << "EXITING: range file truncated\n"; return 1; }
        const double t_read = nowSec();

        ECTable e; e.n_transcripts = nTxL; e.ecs.resize(nEcL); e.n_ecs = nEcL;
        for (std::size_t k = 0; k < nEcL; ++k) {
            EC& c = e.ecs[k];
            c.transcript_ids.assign(ids.begin()+off[k], ids.begin()+off[k+1]);
            c.weights.assign(wts.begin()+off[k], wts.begin()+off[k+1]);
            c.count = cnts[k];
        }
        TranscriptState st; st.resize(nTxL);
        for (std::size_t i = 0; i < nTxL; ++i) {
            st.lengths[i] = rawl[i]; st.eff_lengths[i] = efl[i];
            st.abundances[i] = 0.0; st.counts[i] = 0.0;
        }
        // Use the globally-computed seed rather than letting the engine derive
        // one from this slice; a shard's local totalWeight/numActive are not
        // Salmon's, and the difference is silent.
        rp.initial_alpha = a0;
        if (tolerance > 0.0) rp.tolerance = tolerance;
        if (range_iters > 0) { rp.min_iters = range_iters; rp.max_iters = range_iters; }
        const double t1 = nowSec();
        EMResult r = run_vb(e, st, rp);
        const double t2 = nowSec();
        // Write this shard's results keyed by GLOBAL transcript id, so the
        // gather is a plain concatenation with no coordination.
        double t_write = 0.0;
        if (!range_out.empty()) {
            const double tw = nowSec();
            std::ofstream o(range_out, std::ios::binary);
            if (!o) { std::cerr << "EXITING: cannot write " << range_out << "\n"; return 1; }
            const std::uint64_t n = nTxL;
            o.write((const char*)&n, 8);
            o.write((const char*)l2g.data(), (std::streamsize)(nTxL * 4));
            o.write((const char*)r.counts.data(), (std::streamsize)(nTxL * 8));
            o.write((const char*)st.eff_lengths.data(), (std::streamsize)(nTxL * 8));
            o.close();
            t_write = nowSec() - tw;
        }
        std::printf("range %u/%u: %llu tx, %llu ECs | read %.3f s | vb %.3f s (%u iters) | write %.3f s | total %.3f s\n",
                    idx, nR, (unsigned long long)nTxL, (unsigned long long)nEcL,
                    t_read - t0, t2 - t1, r.iterations, t_write, nowSec() - t0);
        return 0;
    }

    if (assemble) {
        if (genome_dir.empty() || out_prefix.empty() || range_dir.empty() || shard_out.empty()) {
            std::cerr << "EXITING: --assemble needs --genome-dir --range-dir --shard-out --out-prefix\n";
            return 2;
        }
        const double ta0 = nowSec();
        std::string aerr;
        TxIndex atx;
        if (!loadTxIndex(genome_dir, atx, aerr)) {
            std::cerr << "EXITING: " << aerr << "\n"; return 1;
        }
        // Effective lengths come from the setup manifest, counts from the
        // workers. No sidecars and no transcriptome are read here: the gather
        // is a pure scatter-back, which is what keeps it cheap.
        std::vector<double> aeff(atx.n, 0.0);
        {
            const std::string mpath = range_dir + "/manifest.tvbm";
            std::ifstream m(mpath, std::ios::binary);
            if (!m) { std::cerr << "EXITING: cannot read " << mpath << "\n"; return 1; }
            std::uint32_t magic=0, ver=0; std::uint64_t n64=0;
            m.read((char*)&magic,4); m.read((char*)&ver,4); m.read((char*)&n64,8);
            if (magic != 0x4D425654u || n64 != atx.n) {
                std::cerr << "EXITING: manifest does not match this index ("
                          << n64 << " vs " << atx.n << ")\n"; return 1;
            }
            m.read((char*)aeff.data(), (std::streamsize)(atx.n * 8));
        }
        std::vector<double> acnt(atx.n, 0.0);
        std::vector<char> filled(atx.n, 0);
        std::size_t nfiles = 0, ndup = 0;
        for (const std::string& f : globShardOutputs(shard_out)) {
            std::ifstream in(f, std::ios::binary);
            if (!in) { std::cerr << "EXITING: cannot read " << f << "\n"; return 1; }
            std::uint64_t n = 0;
            in.read((char*)&n, 8);
            std::vector<std::uint32_t> l2g(n);
            std::vector<double> c(n), e(n);
            in.read((char*)l2g.data(), (std::streamsize)(n*4));
            in.read((char*)c.data(),   (std::streamsize)(n*8));
            in.read((char*)e.data(),   (std::streamsize)(n*8));
            if (!in) { std::cerr << "EXITING: truncated " << f << "\n"; return 1; }
            for (std::uint64_t k = 0; k < n; ++k) {
                const std::uint32_t g = l2g[k];
                if (g >= atx.n) { std::cerr << "EXITING: bad global id in " << f << "\n"; return 1; }
                if (filled[g]++) ++ndup;      // ranges partition components, so
                acnt[g] = c[k];               // an overlap means the split is wrong
            }
            ++nfiles;
        }
        if (ndup) {
            std::cerr << "EXITING: " << ndup << " transcripts claimed by more than one shard\n";
            return 1;
        }
        // Salmon TPM, identical to the in-process formula.
        double total_norm = 0.0;
        for (std::size_t i = 0; i < atx.n; ++i)
            if (aeff[i] > 0 && acnt[i] > 0) total_norm += acnt[i] / aeff[i];
        std::vector<double> atpm(atx.n, 0.0);
        if (total_norm > 0)
            for (std::size_t i = 0; i < atx.n; ++i)
                if (aeff[i] > 0 && acnt[i] > 0)
                    atpm[i] = (acnt[i] / aeff[i]) / total_norm * 1e6;

        const std::string qpath = out_prefix + "quant.sf";
        std::ofstream q(qpath);
        if (!q) { std::cerr << "EXITING: cannot write " << qpath << "\n"; return 1; }
        q << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
        q.setf(std::ios::fixed);
        for (std::size_t i = 0; i < atx.n; ++i) {
            q << atx.names[i] << '\t' << atx.lengths[i] << '\t';
            q.precision(3); q << aeff[i] << '\t';
            q.precision(6); q << atpm[i] << '\t';
            q.precision(3); q << acnt[i] << '\n';
        }
        q.close();
        std::vector<std::string> agenes;
        if (!loadGeneNames(genome_dir, agenes, aerr)) {
            std::cerr << "EXITING: " << aerr << "\n"; return 1;
        }
        const std::string gpath = out_prefix + "quant.genes.sf";
        if (!writeGeneQuant(gpath, agenes, atx.gene, atx.lengths, aeff, atpm, acnt)) {
            std::cerr << "EXITING: cannot write " << gpath << "\n"; return 1;
        }
        std::size_t covered = 0;
        for (std::size_t i = 0; i < atx.n; ++i) covered += filled[i] ? 1 : 0;
        std::printf("assembled %zu shard outputs: %zu of %zu transcripts covered, "
                    "%zu zero-support (%.2f s)\n",
                    nfiles, covered, atx.n, atx.n - covered, nowSec() - ta0);
        return 0;
    }

    if (genome_dir.empty() || tx_fasta.empty() || out_prefix.empty() || sidecars.empty()) {
        usage(); return 2;
    }

    const double t_start = nowSec();
    std::string err;

    TxIndex tx;
    if (!loadTxIndex(genome_dir, tx, err)) {
        std::cerr << "EXITING: " << err << "\n"; return 1;
    }
    const double t_index = nowSec();
    std::printf("transcript index: %zu transcripts (%.2f s)\n", tx.n, t_index - t_start);

    std::vector<transcript_vb_sidecar::Sidecar> shards;
    for (const std::string& p : sidecars) {
        transcript_vb_sidecar::Sidecar s;
        if (!transcript_vb_sidecar::read(p, s, err)) {
            std::cerr << "EXITING: read " << p << ": " << err << "\n"; return 1;
        }
        shards.push_back(std::move(s));
    }
    // Every ordinal exactly once.
    const std::uint32_t nshard = shards.front().metadata.shard_count;
    std::vector<int> seen(nshard, 0);
    for (const auto& s : shards) {
        if (s.metadata.shard_count != nshard || s.metadata.shard_ordinal >= nshard ||
            seen[s.metadata.shard_ordinal]++) {
            std::cerr << "EXITING: shard ordinals are not a complete set 0.." << (nshard-1) << "\n";
            return 1;
        }
    }
    if (shards.size() != nshard) {
        std::cerr << "EXITING: expected " << nshard << " sidecars, got " << shards.size() << "\n";
        return 1;
    }
    transcript_vb_sidecar::Sidecar gathered;
    if (!transcript_vb_sidecar::merge(shards, gathered, err)) {
        std::cerr << "EXITING: merge: " << err << "\n"; return 1;
    }
    const double t_gather = nowSec();
    ECTable& ecs = gathered.evidence.ec_table;
    if (ecs.n_transcripts != tx.n) {
        std::cerr << "EXITING: sidecar transcript count " << ecs.n_transcripts
                  << " != index " << tx.n << "\n"; return 1;
    }
    std::printf("gathered %zu sidecars: %zu pairs, %zu ECs (%.2f s)\n",
                shards.size(), (std::size_t)gathered.metadata.pair_count,
                ecs.ecs.size(), t_gather - t_index);

    normaliseECWeights(ecs);   // the sidecar is written pre-finalize

    libem::Transcriptome txome;
    if (!txome.loadFromFasta(tx_fasta)) {
        std::cerr << "EXITING: cannot load " << tx_fasta << "\n"; return 1;
    }
    if (!txome.reorderByNames(tx.names)) {
        std::cerr << "EXITING: transcriptome does not match the index order\n"; return 1;
    }
    // NOTE: under --sparse the transcript index space is compacted below, and
    // the GC correction indexes sequences by that same space. The transcriptome
    // must therefore be compacted to match, or every correction lands on the
    // wrong transcript.
    const double t_txome = nowSec();
    std::printf("transcriptome loaded and reordered (%.2f s)\n", t_txome - t_gather);

    FLDAccumulator fld;
    if (!fld.restore(gathered.evidence.fld_state, err)) {
        std::cerr << "EXITING: fld restore: " << err << "\n"; return 1;
    }
    std::vector<double> pmf = fld.getPMF();
    std::vector<double> obs_gc(gathered.evidence.gc_counts.begin(),
                                gathered.evidence.gc_counts.end());
    std::vector<double> eff0 = computeEffectiveLengthsFromPMFWrapper(pmf, tx.lengths);

    // Non-zero-support gating now lives inside vb_engine, where the full
    // transcript space remains available for initialization, the convergence
    // check and the GC background. Compacting here instead shifted the GC bias,
    // because background selection admits unsupported transcripts carrying only
    // the uniform prior.
    std::vector<std::size_t> keep(tx.n);
    for (std::size_t i = 0; i < tx.n; ++i) keep[i] = i;

    // ---- Lustre range emission -------------------------------------------
    // Split the merged evidence into W self-contained range files, one per
    // worker. Each carries its own local transcript space, its globally seeded
    // alpha, its effective lengths and its ECs, so a worker needs no global
    // state and no scheduler: it reads one contiguous file and iterates.
    if (emit_ranges > 0) {
        if (range_dir.empty()) { std::cerr << "EXITING: --emit-ranges needs --range-dir\n"; return 2; }
        const double t_r0 = nowSec();
        // Iteration parameters are recorded in each range file so a worker
        // needs nothing but its own file.
        EMParams rp;
        rp.use_vb = true;
        if (tolerance > 0.0) rp.tolerance = tolerance;
        const std::size_t n_tx = ecs.n_transcripts;
        // components
        std::vector<std::size_t> par(n_tx);
        for (std::size_t i = 0; i < n_tx; ++i) par[i] = i;
        std::function<std::size_t(std::size_t)> fnd =
            [&](std::size_t x){ while(par[x]!=x){par[x]=par[par[x]];x=par[x];} return x; };
        for (const auto& ec : ecs.ecs)
            for (std::size_t k = 1; k < ec.transcript_ids.size(); ++k) {
                std::size_t a=fnd(ec.transcript_ids[0]), b=fnd(ec.transcript_ids[k]);
                if (a!=b) par[a]=b;
            }
        std::unordered_map<std::size_t,std::vector<std::size_t>> comp;
        for (std::size_t e = 0; e < ecs.ecs.size(); ++e)
            if (!ecs.ecs[e].transcript_ids.empty())
                comp[fnd(ecs.ecs[e].transcript_ids[0])].push_back(e);
        // greedy: assign each component to the currently smallest range
        std::vector<std::pair<std::size_t,std::size_t>> byCost;
        for (auto& kv : comp) {
            std::size_t c = 0;
            for (std::size_t e : kv.second) c += ecs.ecs[e].transcript_ids.size();
            byCost.emplace_back(c, kv.first);
        }
        std::sort(byCost.rbegin(), byCost.rend());
        const int W = emit_ranges;
        std::vector<std::vector<std::size_t>> rng(W);
        std::vector<std::size_t> load(W, 0);
        for (auto& pr : byCost) {
            const int w = std::min_element(load.begin(), load.end()) - load.begin();
            for (std::size_t e : comp[pr.second]) rng[w].push_back(e);
            load[w] += pr.first;
        }
        // Setup burn-in. The engine applies the dynamic GC effective-length
        // correction exactly once, at iter > effective_length_update_target_iter,
        // from the alpha of that moment. Its only global quantity is a 101-bin
        // histogram -- tiny, but collecting it mid-run would force every worker
        // through a barrier. So we pay those ~11 iterations of ~1270 once, here,
        // and hand each worker the corrected lengths plus the alpha the global
        // run held at the fire point. A worker's iteration 0 is then the global
        // run's iteration 11, and workers need no rendezvous at all.
        //
        // The firing rule is not restated here: we let the engine call us, so
        // setup cannot drift from what an in-process run would have done.
        std::vector<double> alpha0, effS;
        const double t_burn = nowSec();
        {
            TranscriptState st0; st0.resize(n_tx);
            for (std::size_t i = 0; i < n_tx; ++i) {
                st0.names[i]=tx.names[i]; st0.lengths[i]=tx.lengths[i];
                st0.eff_lengths[i]=eff0[i]; st0.abundances[i]=0.0; st0.counts[i]=0.0;
            }
            EMParams pb = rp;
            pb.threads = threads;
            // One iteration past the trigger, with early exit disabled so the
            // callback is guaranteed to fire.
            pb.min_iters = pb.max_iters = pb.effective_length_update_target_iter + 2;
            pb.tolerance = 0.0;
            pb.effective_length_update =
                [&](std::uint32_t, TranscriptState& cb, const std::vector<double>& alpha) -> bool {
                    alpha0 = alpha;
                    DynamicGCEffectiveLengthResult r =
                        computeDynamicGCBiasedEffectiveLengthsWrapper(
                            txome, pmf, tx.lengths, alpha, cb.eff_lengths, obs_gc,
                            threads);   // setup is a one-shot stage; take the speedup
                    if (!r.applied) { effS = cb.eff_lengths; return false; }
                    cb.eff_lengths.swap(r.effective_lengths);
                    effS = cb.eff_lengths;
                    return true;
                };
            run_vb(ecs, st0, pb);
        }
        if (effS.empty())   effS = eff0;                 // GC correction not applicable
        if (alpha0.empty()) {                            // callback never fired
            TranscriptState st0; st0.resize(n_tx);
            for (std::size_t i = 0; i < n_tx; ++i) {
                st0.names[i]=tx.names[i]; st0.lengths[i]=tx.lengths[i];
                st0.eff_lengths[i]=effS[i]; st0.abundances[i]=0.0; st0.counts[i]=0.0;
            }
            alpha0 = compute_initial_alpha(ecs, st0, rp);
        }
        std::printf("  setup burn-in (%u iters + GC length correction): %.2f s\n",
                    rp.effective_length_update_target_iter + 2, nowSec() - t_burn);
        std::size_t bytes_total = 0;
        for (int w = 0; w < W; ++w) {
            std::vector<std::uint32_t> loc;
            for (std::size_t e : rng[w])
                for (std::uint32_t t : ecs.ecs[e].transcript_ids) loc.push_back(t);
            std::sort(loc.begin(), loc.end());
            loc.erase(std::unique(loc.begin(), loc.end()), loc.end());
            std::unordered_map<std::uint32_t,std::uint32_t> rm;
            for (std::uint32_t i = 0; i < loc.size(); ++i) rm[loc[i]] = i;
            char path[1024];
            std::snprintf(path, sizeof path, "%s/range_%04d.tvbr", range_dir.c_str(), w);
            std::ofstream o(path, std::ios::binary);
            if (!o) { std::cerr << "EXITING: cannot write " << path << "\n"; return 1; }
            auto W32=[&](std::uint32_t v){ o.write((const char*)&v,4); };
            auto W64=[&](std::uint64_t v){ o.write((const char*)&v,8); };
            auto WD =[&](double v){ o.write((const char*)&v,8); };
            const std::uint32_t magic = 0x52425654u; // "TVBR"
            W32(magic); W32(1); W32((std::uint32_t)W); W32((std::uint32_t)w);
            W64(loc.size()); W64(rng[w].size()); W64(n_tx);
            WD(rp.vb_prior); WD(rp.tolerance);
            W32(rp.min_iters); W32(rp.max_iters);
            for (std::uint32_t g : loc) W32(g);                     // local -> global
            for (std::uint32_t g : loc) WD(alpha0[g]);              // seeded alpha
            for (std::uint32_t g : loc) WD(effS[g]);                // GC-corrected eff lengths
            for (std::uint32_t g : loc) W32((std::uint32_t)tx.lengths[g]);
            std::uint64_t nnz = 0;
            for (std::size_t e : rng[w]) nnz += ecs.ecs[e].transcript_ids.size();
            W64(nnz);
            std::uint32_t off = 0;
            for (std::size_t e : rng[w]) { W32(off); off += ecs.ecs[e].transcript_ids.size(); }
            W32(off);
            for (std::size_t e : rng[w]) for (std::uint32_t t : ecs.ecs[e].transcript_ids) W32(rm[t]);
            for (std::size_t e : rng[w])
                for (std::size_t k = 0; k < ecs.ecs[e].transcript_ids.size(); ++k)
                    WD(ecs.ecs[e].weights.empty() ? 1.0/ecs.ecs[e].transcript_ids.size()
                                                  : ecs.ecs[e].weights[k]);
            for (std::size_t e : rng[w]) WD(ecs.ecs[e].count);
            bytes_total += (std::size_t)o.tellp();
            o.close();
        }
        std::printf("emitted %d range files to %s (%.1f MB total, %.2f s)\n",
                    W, range_dir.c_str(), bytes_total/1048576.0, nowSec()-t_r0);
        // Effective lengths for every transcript, including the zero-support
        // ones no shard owns. Workers never touch eff_lengths (they run without
        // the GC callback), so this is the whole picture and the assembler
        // needs nothing else to reconstruct quant.sf.
        {
            const std::string mpath = range_dir + "/manifest.tvbm";
            std::ofstream m(mpath, std::ios::binary);
            if (!m) { std::cerr << "EXITING: cannot write " << mpath << "\n"; return 1; }
            const std::uint32_t magic = 0x4D425654u; // "TVBM"
            const std::uint32_t ver = 1;
            const std::uint64_t n64 = n_tx;
            m.write((const char*)&magic,4); m.write((const char*)&ver,4);
            m.write((const char*)&n64,8);
            m.write((const char*)effS.data(), (std::streamsize)(n_tx * 8));
        }
        std::size_t mn=load[0], mx=load[0];
        for (std::size_t l : load) { mn=std::min(mn,l); mx=std::max(mx,l); }
        std::printf("  balance: max/mean %.4f  min/mean %.4f\n",
                    (double)mx*W/(double)std::accumulate(load.begin(),load.end(),(std::size_t)0),
                    (double)mn*W/(double)std::accumulate(load.begin(),load.end(),(std::size_t)0));
        return 0;
    }

    TranscriptState st;
    st.resize(keep.size());
    for (std::size_t k = 0; k < keep.size(); ++k) {
        const std::size_t i = keep[k];
        st.names[k] = tx.names[i];
        st.lengths[k] = tx.lengths[i];
        st.eff_lengths[k] = eff0[i];
        st.abundances[k] = 0.0;
        st.counts[k] = 0.0;
    }


    EMParams params;
    params.use_vb = true;
    params.threads = threads;
    if (tolerance > 0.0) params.tolerance = tolerance;
    params.per_component_convergence = per_component;
    params.component_partition = comp_partition;
    // Salmon derives uniformPrior from ALL transcripts, so a compacted run must
    // still divide by the original count.
    params.num_active_override = 0;
    int gc_iter = -1;
    if (!no_gc)
    params.effective_length_update =
        [&](std::uint32_t it, TranscriptState& cb, const std::vector<double>& alpha) -> bool {
            DynamicGCEffectiveLengthResult r =
                computeDynamicGCBiasedEffectiveLengthsWrapper(
                    txome, pmf, tx.lengths, alpha, cb.eff_lengths, obs_gc,
                    parallel_gc ? threads : 1);
            if (!r.applied) return false;
            cb.eff_lengths.swap(r.effective_lengths);
            gc_iter = (int)it;
            return true;
        };

    const double t_pre = nowSec();
    EMResult res = run_vb(ecs, st, params);
    const double t_vb = nowSec();
    std::printf("VB: GC update at iteration %d, converged=%s, iterations=%u (%.2f s)\n",
                gc_iter, res.converged ? "yes" : "no", res.iterations, t_vb - t_pre);

    // Expand back to the full index space if gated.
    std::vector<double> counts(tx.n, 0.0), efflen(tx.n, 0.0), tpm(tx.n, 0.0);
    for (std::size_t k = 0; k < keep.size(); ++k) {
        counts[keep[k]] = res.counts[k];
        efflen[keep[k]] = st.eff_lengths[k];
        if (k < res.tpm.size()) tpm[keep[k]] = res.tpm[k];
    }
    for (std::size_t i = 0; i < tx.n; ++i)
        if (efflen[i] == 0.0) efflen[i] = eff0[i];

    const std::string qpath = out_prefix + "quant.sf";
    std::ofstream q(qpath);
    if (!q) { std::cerr << "EXITING: cannot write " << qpath << "\n"; return 1; }
    q << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
    q.setf(std::ios::fixed);
    for (std::size_t i = 0; i < tx.n; ++i) {
        q << tx.names[i] << '\t' << tx.lengths[i] << '\t';
        q.precision(3); q << efflen[i] << '\t';
        q.precision(6); q << tpm[i] << '\t';
        q.precision(3); q << counts[i] << '\n';
    }
    q.close();

    std::vector<std::string> gene_names;
    if (!loadGeneNames(genome_dir, gene_names, err)) {
        std::cerr << "EXITING: " << err << "\n"; return 1;
    }
    const std::string gpath = out_prefix + "quant.genes.sf";
    if (!writeGeneQuant(gpath, gene_names, tx.gene, tx.lengths, efflen, tpm, counts)) {
        std::cerr << "EXITING: cannot write " << gpath << "\n"; return 1;
    }
    const double t_end = nowSec();
    std::printf("wrote %s and %s (%.2f s)\ntotal %.2f s\n",
            qpath.c_str(), gpath.c_str(), t_end - t_vb, t_end - t_start);
    return 0;
}
