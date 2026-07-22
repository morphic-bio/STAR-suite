#!/usr/bin/env python3
"""Compare STAR Suite vs the nf-core star_salmon reference chain and write a report.

Reads the run tree produced by run_compare.sh and emits compare/report.md plus
TSVs. Baseline = the reference chain's Salmon quant (what nf-core/rnaseq
--aligner star_salmon produces). We compare it to STAR Suite's integrated quant
and to STAR Suite's BAM -> Salmon (exact) quant, and we compare wall time, peak
RSS, read counts, and BAM stats. numpy-only (no scipy).
"""
import argparse, os, re, sys
import numpy as np

# ----------------------------- small helpers --------------------------------
def parse_time_v(path):
    """Return (wall_seconds, max_rss_gb) from a /usr/bin/time -v log, or (None,None)."""
    if not os.path.exists(path):
        return None, None
    wall = rss = None
    for line in open(path):
        if "Elapsed (wall clock) time" in line:
            t = line.strip().split()[-1]            # last token, e.g. "0:02.34" or "1:23:45"
            parts = [float(x) for x in t.split(":")]
            wall = parts[0]*3600+parts[1]*60+parts[2] if len(parts) == 3 else \
                   parts[0]*60+parts[1] if len(parts) == 2 else parts[0]
        elif "Maximum resident set size" in line:
            rss = float(line.strip().split()[-1])/1e6  # KB -> GB
    return wall, rss

def fmt_secs(s):
    if s is None: return "n/a"
    m, sec = divmod(int(round(s)), 60); h, m = divmod(m, 60)
    return f"{h:d}:{m:02d}:{sec:02d}" if h else f"{m:d}:{sec:02d}"

def load_quant(path):
    """Load a salmon-format quant.sf -> dict name->(TPM, NumReads). None if missing."""
    if not os.path.exists(path): return None
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[ix["Name"]]] = (float(f[ix["TPM"]]), float(f[ix["NumReads"]]))
    return out

def pearson(x, y):
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0: return float("nan")
    return float(np.corrcoef(x, y)[0, 1])

def spearman(x, y):
    if len(x) < 2: return float("nan")
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return pearson(rx.astype(float), ry.astype(float))

def corr_block(ref, other, tx2gene):
    """Return a dict of correlations (transcript + gene, all + expressed) for NumReads & TPM."""
    keys = sorted(set(ref) & set(other))
    res = {"n_shared": len(keys), "n_ref_only": len(set(ref)-set(other)),
           "n_other_only": len(set(other)-set(ref))}
    all_keys = sorted(set(ref) | set(other))
    if all_keys:
        nr_r_all = np.array([ref.get(k, (0.0, 0.0))[1] for k in all_keys])
        nr_o_all = np.array([other.get(k, (0.0, 0.0))[1] for k in all_keys])
        abs_delta = np.abs(nr_o_all - nr_r_all)
        total_ref = float(nr_r_all.sum())
        total_other = float(nr_o_all.sum())
        total_delta = total_other - total_ref
        l1_delta = float(abs_delta.sum())
        res["numreads_total_ref"] = total_ref
        res["numreads_total_other"] = total_other
        res["numreads_total_delta"] = total_delta
        res["numreads_abs_delta_sum"] = l1_delta
        # If totals are equal, half-L1 is the number of read-equivalents moved
        # between transcripts. With unequal totals, include the assigned/unassigned
        # delta as an extra bin.
        res["numreads_reassigned_est"] = 0.5 * (l1_delta + abs(total_delta))
        res["tx_numreads_absdiff_gt_0p01"] = int((abs_delta > 0.01).sum())
        res["tx_numreads_absdiff_gt_0p1"] = int((abs_delta > 0.1).sum())
        res["tx_numreads_absdiff_gt_1"] = int((abs_delta > 1.0).sum())
        res["tx_numreads_max_abs_delta"] = float(abs_delta.max())
    if not keys: return res
    tpm_r = np.array([ref[k][0] for k in keys]); tpm_o = np.array([other[k][0] for k in keys])
    nr_r  = np.array([ref[k][1] for k in keys]); nr_o  = np.array([other[k][1] for k in keys])
    expr = (nr_r > 0) | (nr_o > 0)
    res["tx_pearson_numreads"]   = pearson(nr_r, nr_o)
    res["tx_spearman_numreads"]  = spearman(nr_r, nr_o)
    res["tx_pearson_tpm_expr"]   = pearson(tpm_r[expr], tpm_o[expr])
    res["tx_spearman_tpm_expr"]  = spearman(tpm_r[expr], tpm_o[expr])
    res["n_expressed"] = int(expr.sum())
    # gene-level aggregation of NumReads
    if tx2gene:
        gr, go = {}, {}
        for k in keys:
            g = tx2gene.get(k)
            if g is None: continue
            gr[g] = gr.get(g, 0.0) + ref[k][1]; go[g] = go.get(g, 0.0) + other[k][1]
        gk = sorted(gr)
        if gk:
            a = np.array([gr[g] for g in gk]); b = np.array([go[g] for g in gk])
            res["gene_pearson_numreads"]  = pearson(a, b)
            res["gene_spearman_numreads"] = spearman(a, b)
            res["n_genes"] = len(gk)
    return res

def flagstat_mapped(path):
    if not os.path.exists(path): return None
    for line in open(path):
        m = re.match(r"(\d+) \+ \d+ mapped \(", line)
        if m: return int(m.group(1))
    return None

def star_input_reads(log_final):
    if not os.path.exists(log_final): return None
    for line in open(log_final):
        if "Number of input reads" in line:
            return int(line.rsplit("|", 1)[1].strip())
    return None

# ----------------------------------- main -----------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True); ap.add_argument("--mode", required=True)
ap.add_argument("--threads", required=True); ap.add_argument("--salmon-threads", default=None)
ap.add_argument("--tx2gene", required=True)
a = ap.parse_args()
O = a.outdir

tx2gene = {}
if os.path.exists(a.tx2gene):
    for line in open(a.tx2gene):
        t, g = line.rstrip("\n").split("\t"); tx2gene[t] = g

# timings
t_trim   = parse_time_v(f"{O}/A_reference_chain/trim/time.log")
t_star   = parse_time_v(f"{O}/A_reference_chain/star/time.log")
t_salmon = parse_time_v(f"{O}/A_reference_chain/salmon/time.log")
t_suite  = parse_time_v(f"{O}/B_starsuite/time.log")
t_b2     = parse_time_v(f"{O}/B2_starsuite_bam_salmon/time.log")
ref_wall = sum(x[0] for x in (t_trim, t_star, t_salmon) if x[0] is not None)
ref_rss  = max([x[1] for x in (t_trim, t_star, t_salmon) if x[1] is not None] or [0])
suite_wall, suite_rss = t_suite
speedup = ref_wall/suite_wall if suite_wall else float("nan")

# quant
ref_q   = load_quant(f"{O}/A_reference_chain/salmon/quant.sf")
suite_q = load_quant(f"{O}/B_starsuite/quant.sf")
b2_q    = load_quant(f"{O}/B2_starsuite_bam_salmon/quant.sf")
cmp_int = corr_block(ref_q, suite_q, tx2gene) if ref_q and suite_q else {}
cmp_exact = corr_block(ref_q, b2_q, tx2gene) if ref_q and b2_q else {}

# bam + reads
ref_mapped   = flagstat_mapped(f"{O}/A_reference_chain/star/flagstat.txt")
suite_mapped = flagstat_mapped(f"{O}/B_starsuite/flagstat.txt")
ref_in   = star_input_reads(f"{O}/A_reference_chain/star/Log.final.out")
suite_in = star_input_reads(f"{O}/B_starsuite/Log.final.out")

def g(d, k):
    v = d.get(k); return f"{v:.4f}" if isinstance(v, float) and v == v else ("n/a" if v is None else str(v))

R = []
R.append(f"# STAR Suite vs nf-core star_salmon — output parity check ({a.mode})\n")
R.append(f"Reference chain = Trim Galore -> STAR -> Salmon (the tools "
         f"`nf-core/rnaseq --aligner star_salmon` runs); its Salmon quant is the "
         f"baseline. The question here is **parity** — does swapping in STAR Suite "
         f"change the *answer*. STAR threads: {a.threads}; Salmon quant threads: "
         f"{a.salmon_threads or a.threads}. All intermediates under `{O}`.\n")

R.append("## Quant parity vs the reference chain's Salmon (the headline)\n")
R.append("Two STAR Suite quant modes: *integrated* (STAR Suite's own quant) and "
         "*exact* (STAR Suite's BAM handed to Salmon).\n")
R.append("| Comparison | shared tx | tx Pearson (NumReads) | tx Spearman (NumReads) | tx Pearson TPM (expr) | gene Pearson (NumReads) |")
R.append("|---|---|---|---|---|---|")
R.append(f"| integrated vs Salmon | {g(cmp_int,'n_shared')} | {g(cmp_int,'tx_pearson_numreads')} | "
         f"{g(cmp_int,'tx_spearman_numreads')} | {g(cmp_int,'tx_pearson_tpm_expr')} | {g(cmp_int,'gene_pearson_numreads')} |")
R.append(f"| exact (BAM->Salmon) vs Salmon | {g(cmp_exact,'n_shared')} | {g(cmp_exact,'tx_pearson_numreads')} | "
         f"{g(cmp_exact,'tx_spearman_numreads')} | {g(cmp_exact,'tx_pearson_tpm_expr')} | {g(cmp_exact,'gene_pearson_numreads')} |")
R.append("")
R.append(f"- integrated: {g(cmp_int,'n_expressed')} expressed transcripts, {g(cmp_int,'n_genes')} genes.")
R.append("- exact mode is the *same* Salmon on STAR Suite's BAM, so ~1.0 confirms the alignment matches.")
R.append("- parity is depth-robust (both pipelines see the same reads), so the chr22 fixture is a valid parity test.\n")

R.append("## Read-equivalent differences\n")
R.append("NumReads are fractional after multimapper allocation. `Half-L1 moved` is the "
         "approximate number of read-equivalents assigned to a different transcript "
         "(including assigned/unassigned imbalance when totals differ).\n")
R.append("| Comparison | Salmon total | Other total | total delta | abs delta sum | half-L1 moved | tx >0.1 reads | tx >1 read | max tx delta |")
R.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
R.append(f"| integrated vs Salmon | {g(cmp_int,'numreads_total_ref')} | {g(cmp_int,'numreads_total_other')} | "
         f"{g(cmp_int,'numreads_total_delta')} | {g(cmp_int,'numreads_abs_delta_sum')} | "
         f"{g(cmp_int,'numreads_reassigned_est')} | {g(cmp_int,'tx_numreads_absdiff_gt_0p1')} | "
         f"{g(cmp_int,'tx_numreads_absdiff_gt_1')} | {g(cmp_int,'tx_numreads_max_abs_delta')} |")
R.append(f"| exact (BAM->Salmon) vs Salmon | {g(cmp_exact,'numreads_total_ref')} | {g(cmp_exact,'numreads_total_other')} | "
         f"{g(cmp_exact,'numreads_total_delta')} | {g(cmp_exact,'numreads_abs_delta_sum')} | "
         f"{g(cmp_exact,'numreads_reassigned_est')} | {g(cmp_exact,'tx_numreads_absdiff_gt_0p1')} | "
         f"{g(cmp_exact,'tx_numreads_absdiff_gt_1')} | {g(cmp_exact,'tx_numreads_max_abs_delta')} |")
R.append("")

R.append("## Alignment (genome BAM) and reads\n")
R.append("| Metric | Reference chain (STAR) | STAR Suite |")
R.append("|---|---|---|")
R.append(f"| Input reads into alignment | {ref_in if ref_in is not None else 'n/a'} | {suite_in if suite_in is not None else 'n/a'} |")
R.append(f"| Mapped alignments (flagstat) | {ref_mapped if ref_mapped is not None else 'n/a'} | {suite_mapped if suite_mapped is not None else 'n/a'} |")
R.append("")
R.append("Trimming: the reference chain trims with Trim Galore; STAR Suite trims "
         "internally (cutadapt-compatible) in the same pass. The 'input reads into "
         "alignment' row is the post-trim count for each — close agreement indicates "
         "equivalent trimming; trimmed FASTQs are kept under `A_reference_chain/trim/`.\n")

R.append("## Timing (indicative only — NOT the point of this run)\n")
R.append("chr22 is a parity slice, far too small to represent the integrated pass's "
         "speed advantage. Speed is a full-dataset question — see STAR Suite's "
         "published full-genome benchmarks, which anyone can re-run. Shown only for "
         "completeness.\n")
R.append("| Stage | Wall | Peak RSS (GB) |")
R.append("|---|---|---|")
R.append(f"| A. reference chain (Trim Galore+STAR+Salmon) | {fmt_secs(ref_wall)} | {ref_rss:.2f} |")
R.append(f"| B. STAR Suite (integrated, one binary) | {fmt_secs(suite_wall)} | {suite_rss or 0:.2f} |")
R.append("")

R.append("## What is on disk (every intermediate is kept)\n")
R.append("```")
R.append("fixtures/chr22/            cached chr22-filtered FASTQs (the parity fixture)")
R.append("A_reference_chain/trim/    trimmed FASTQs + Trim Galore reports")
R.append("A_reference_chain/star/    genome BAM, Aligned.toTranscriptome.out.bam, Log.final.out")
R.append("A_reference_chain/salmon/  quant.sf (BASELINE)")
R.append("B_starsuite/               genome BAM, transcriptome BAM, quant.sf + quant.genes.sf, Log.final.out")
R.append("B2_starsuite_bam_salmon/   quant.sf (exact star_salmon-equivalent)")
R.append("compare/                   this report + TSVs")
R.append("```")

os.makedirs(f"{O}/compare", exist_ok=True)
open(f"{O}/compare/report.md", "w").write("\n".join(R) + "\n")

# machine-readable summary
with open(f"{O}/compare/summary.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    fh.write(f"mode\t{a.mode}\nthreads\t{a.threads}\n")
    fh.write(f"salmon_threads\t{a.salmon_threads or a.threads}\n")
    fh.write(f"ref_chain_wall_s\t{ref_wall:.1f}\nstarsuite_wall_s\t{suite_wall or 0:.1f}\n")
    fh.write(f"speedup\t{speedup:.3f}\n")
    for k, v in cmp_int.items():   fh.write(f"integrated.{k}\t{v}\n")
    for k, v in cmp_exact.items(): fh.write(f"exact.{k}\t{v}\n")
print(f"wrote {O}/compare/report.md")
