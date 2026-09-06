#!/usr/bin/env python3
"""Concordance vs Cell Ranger under the manuscript's stated Methods protocol:
Pearson and Spearman on log1p counts, restricted to genes with >=20 total
counts in BOTH outputs and detected in >=1% of common cells (in both), over
common called cells. Reports mean and median per-cell correlations and a
gene-level Pearson on log1p per-gene totals. Genes matched by Ensembl ID.
Usage: paper_protocol_concordance.py OUTPUT_ROOT LABEL [--input-kind star|cyto]
"""
import sys, os, argparse
import numpy as np, scipy.stats, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import concordance_vs_cr as cc
from collections import Counter
from pathlib import Path

ap = argparse.ArgumentParser(); ap.add_argument("output_root", type=Path); ap.add_argument("label")
ap.add_argument("--input-kind", choices=("star", "cyto"), default="star"); ap.add_argument("--cr-root", type=Path, default=cc.DEFAULT_CRROOT)
a = ap.parse_args()
groups = cc.DEFAULT_GROUPS
first_cr = a.cr_root / groups[0][0] / "count" / "sample_filtered_feature_bc_matrix"
with cc.open_maybe_gz(first_cr, "features.tsv") as h: cr_rows = [l.rstrip("\n").split("\t") for l in h]
sc = Counter(r[1] for r in cr_rows if len(r) > 1); sym2id = {r[1]: r[0] for r in cr_rows if len(r) > 1 and sc[r[1]] == 1}
print(f"# paper-protocol concordance vs Cell Ranger 9.0.1 - {a.label}  (log1p; genes >=20 counts in both & detected in >=1% of cells; Ensembl-ID matched)")
print(f"{'sample':20s} {'cells':>11s} {'genes':>6s} {'Jaccard':>8s} {'cellPear mean':>14s} {'median':>8s} {'cellSpear mean':>15s} {'median':>8s} {'genePear':>9s}")
J=[];PM=[];PMd=[];SM=[];SMd=[];G=[]
for sample, bcs in groups:
    if a.input_kind == "star": qx, qb, qs = cc.read_star_group(a.output_root, bcs)
    else:
        qx, qb, qs = cc.read_cyto_group(a.output_root, bcs); qs = [sym2id.get(s, s) for s in qs]
    cx, cb, cs = cc.read_mex(a.cr_root / sample / "count" / "sample_filtered_feature_bc_matrix"); cx, cb = cc.collapse_duplicate_cells(cx, cb)
    qi = {g: i for i, g in enumerate(qs) if Counter(qs)[g] == 1}; ci = {g: i for i, g in enumerate(cs) if Counter(cs)[g] == 1}
    cg = sorted(set(qi) & set(ci)); common = sorted(set(qb) & set(cb)); jac = len(common) / len(set(qb) | set(cb))
    iq = {b: i for i, b in enumerate(qb)}; ic = {b: i for i, b in enumerate(cb)}
    Q = qx[[iq[b] for b in common]][:, [qi[g] for g in cg]].tocsc(); C = cx[[ic[b] for b in common]][:, [ci[g] for g in cg]].tocsc()
    tq = np.asarray(Q.sum(0)).ravel(); tc = np.asarray(C.sum(0)).ravel(); dq = np.asarray((Q > 0).sum(0)).ravel(); dc = np.asarray((C > 0).sum(0)).ravel()
    keep = (tq >= 20) & (tc >= 20) & (dq >= 0.01 * len(common)) & (dc >= 0.01 * len(common))
    Ql = np.log1p(Q[:, keep].toarray()); Cl = np.log1p(C[:, keep].toarray())
    pr = []; sr = []
    for i in range(len(common)):
        x, y = Ql[i], Cl[i]
        if x.std() > 0 and y.std() > 0: pr.append(scipy.stats.pearsonr(x, y)[0]); sr.append(scipy.stats.spearmanr(x, y)[0])
    pr = np.array(pr); sr = np.array(sr); gp = scipy.stats.pearsonr(np.log1p(tq[keep]), np.log1p(tc[keep]))[0]
    print(f"{sample:20s} {len(common):>11d} {int(keep.sum()):>6d} {jac:8.4f} {pr.mean():14.6f} {np.median(pr):8.6f} {sr.mean():15.6f} {np.median(sr):8.6f} {gp:9.6f}")
    J.append(jac); PM.append(pr.mean()); PMd.append(np.median(pr)); SM.append(sr.mean()); SMd.append(np.median(sr)); G.append(gp)
print(f"{'mean over samples':20s} {'':>11s} {'':>6s} {np.mean(J):8.4f} {np.mean(PM):14.6f} {np.mean(PMd):8.6f} {np.mean(SM):15.6f} {np.mean(SMd):8.6f} {np.mean(G):9.6f}")
