#!/usr/bin/env python3
"""Concordance of STAR Suite with Cell Ranger, per sample and with all samples pooled.

For a dataset made of one or more samples (pairs of Cell Ranger / STAR Suite filtered MEX directories),
every metric is reported at three levels:
  per-sample  each sample compared on its own, then averaged over samples (sample-aware: per-sample
              totals also carry any disagreement in assigning reads to samples)
  pooled      all samples concatenated and each metric computed once (sample-ignorant; the values the
              manuscript quotes). A cell is still its 16-base barcode plus sample tag, because in
              multiplexed Flex one barcode can hold cells of several samples.
  barcode     diagnostic only: tags dropped and counts summed per 16-base barcode. This compares
              partitions rather than cells in multiplexed Flex (the 320k dataset has 325,410 Cell Ranger
              cells on 128,394 barcodes), so it is not quoted.
For a single-sample dataset the levels coincide and only the per-sample line is printed.

Metrics, each on the cells both tools called ("shared"):
  jaccard            |A & B| / |A | B| of the called-cell sets
  cell_total_r       cell Pearson: Pearson, across shared cells, of each cell's total count
                     (all Gene Expression features)
  cell_profile_r     mean per-cell Pearson: Pearson, across genes, of one shared cell's log1p counts in the
                     two outputs, averaged over shared cells (median also printed); genes with >=20 counts
                     in both outputs and detected in >=1% of shared cells in both; cells with no variation
                     skipped
  cell_profile_rho   the same with Spearman correlation (ties averaged)
  gene_spearman, gene_pearson   raw per-gene totals over shared cells, every gene ID in both outputs,
                     zero-count genes included

Barcodes: a trailing "-1" is stripped. Cell Ranger and STAR Suite Flex outputs both write CB16 + tag8
barcodes, so the pooled level keys cells correctly without further mapping.

SPEC.json is a list of datasets:
  [{"name": "JAX",
    "translation": null,                          # optional two-column barcode translation file
    "samples": [{"name": "WT-Day-7",
                 "cr":   ".../per_sample_outs/WT-Day-7/count/sample_filtered_feature_bc_matrix",
                 "star": ".../per_sample/BC004/Gene/filtered"}, ...]}]
A sample may give "star_h5ad" (cyto's filtered h5ad: cells CB16-BCnnn, genes by symbol) instead of "star";
its cells are keyed CB16 + the sample's tag 8-mer (taken from the Cell Ranger barcodes) and its gene
symbols are mapped to Cell Ranger's gene IDs where the symbol is unique (requires anndata).
For scRNA-seq and Perturb-seq, "cr" is outs/filtered_feature_bc_matrix and "star" is
Solo.out/GeneFull/filtered.

Requires numpy, scipy and pyarrow. Output is a tab-separated table on stdout.
Usage: concordance_levels.py SPEC.json
"""
import gzip, json, os, sys
import numpy as np, scipy.sparse as sp, scipy.stats as ss
import pyarrow as pa, pyarrow.csv as pc

pa.set_cpu_count(8); pa.set_io_thread_count(8)

def fpath(d, stem):
    for n in (stem, stem + ".gz"):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    raise FileNotFoundError(os.path.join(d, stem))

def lines(p):
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as h:
        return [l.rstrip("\n") for l in h]

def read_mex(d, tr=None):
    """-> (cells x genes CSR float64, barcodes without '-1', gene ids) restricted to Gene Expression rows."""
    feats = [l.split("\t") for l in lines(fpath(d, "features.tsv"))]
    bcs = [b[:-2] if b.endswith("-1") else b for b in lines(fpath(d, "barcodes.tsv"))]
    if tr: bcs = [tr.get(b, b) for b in bcs]
    mp = fpath(d, "matrix.mtx"); skip = 0
    with (gzip.open(mp, "rt") if mp.endswith(".gz") else open(mp)) as h:
        for l in h:
            skip += 1
            if not l.startswith("%"): break   # the dimension line is skipped too
    t = pc.read_csv(mp, read_options=pc.ReadOptions(skip_rows=skip, column_names=["g", "c", "v"]),
                    parse_options=pc.ParseOptions(delimiter=" "),
                    convert_options=pc.ConvertOptions(column_types={"g": pa.int32(), "c": pa.int32(), "v": pa.float64()}))
    g = t.column("g").to_numpy() - 1; c = t.column("c").to_numpy() - 1; v = t.column("v").to_numpy()
    m = sp.csr_matrix((v, (c, g)), shape=(len(bcs), len(feats)))
    ge = np.array([len(f) < 3 or f[2] == "Gene Expression" for f in feats])
    return m[:, np.flatnonzero(ge)].tocsr(), bcs, [f[0] for f, k in zip(feats, ge) if k]

def read_h5ad(path, cr_bcs, cr_dir):
    import anndata as ad
    tags = {b[16:] for b in cr_bcs}; assert len(tags) == 1, "h5ad input needs a single-tag sample"
    feats = [l.split("\t") for l in lines(fpath(cr_dir, "features.tsv"))]
    n = {}
    for f in feats: n[f[1]] = n.get(f[1], 0) + 1
    sym2id = {f[1]: f[0] for f in feats if n[f[1]] == 1}
    a = ad.read_h5ad(path)
    return (sp.csr_matrix(a.X, dtype=np.float64), [str(x)[:16] + next(iter(tags)) for x in a.obs_names],
            [sym2id.get(str(g), str(g)) for g in a.var_names])

def collapse(m, keys):
    """Sum rows sharing a key; returns (CSR, sorted unique keys)."""
    uk, inv = np.unique(np.array(keys, dtype=object), return_inverse=True)
    if len(uk) == len(keys):
        o = np.argsort(inv); return m[o].tocsr(), list(uk)
    agg = sp.csr_matrix((np.ones(len(keys)), (inv, np.arange(len(keys)))), shape=(len(uk), len(keys)))
    return (agg @ m).tocsr(), list(uk)

def row_pearson(a, b):
    a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
    na = np.sqrt((a * a).sum(1)); nb = np.sqrt((b * b).sum(1)); ok = (na > 0) & (nb > 0)
    return (a * b).sum(1)[ok] / (na[ok] * nb[ok])

def compare(cr, st):
    (cm, ck, cg), (sm, sk, sg) = cr, st
    ci = {k: i for i, k in enumerate(ck)}; si = {k: i for i, k in enumerate(sk)}
    shared = sorted(set(ci) & set(si)); union = len(set(ci) | set(si))
    cr_rows = np.array([ci[k] for k in shared]); st_rows = np.array([si[k] for k in shared])
    C = cm[cr_rows]; S = sm[st_rows]
    tot_r = ss.pearsonr(np.asarray(C.sum(1)).ravel(), np.asarray(S.sum(1)).ravel()).statistic
    cu = {g for g in cg if cg.count(g) == 1} if len(set(cg)) != len(cg) else set(cg)
    su = {g for g in sg if sg.count(g) == 1} if len(set(sg)) != len(sg) else set(sg)
    genes = sorted(cu & su); cgi = {g: i for i, g in enumerate(cg)}; sgi = {g: i for i, g in enumerate(sg)}
    C = C[:, [cgi[g] for g in genes]].tocsr(); S = S[:, [sgi[g] for g in genes]].tocsr()
    tc = np.asarray(C.sum(0)).ravel(); ts = np.asarray(S.sum(0)).ravel()
    g_s = ss.spearmanr(tc, ts).statistic; g_p = ss.pearsonr(tc, ts).statistic
    dc = np.asarray((C > 0).sum(0)).ravel(); ds = np.asarray((S > 0).sum(0)).ravel()
    keep = np.flatnonzero((tc >= 20) & (ts >= 20) & (dc >= 0.01 * len(shared)) & (ds >= 0.01 * len(shared)))
    Ck = C[:, keep].tocsr(); Sk = S[:, keep].tocsr(); rs, rh = [], []
    for i in range(0, len(shared), 4000):
        a = np.log1p(Ck[i:i + 4000].toarray()); b = np.log1p(Sk[i:i + 4000].toarray())
        rs.append(row_pearson(a, b))
        rh.append(row_pearson(ss.rankdata(a, axis=1), ss.rankdata(b, axis=1)))
    rs = np.concatenate(rs); rh = np.concatenate(rh)
    return dict(cr_cells=len(ck), star_cells=len(sk), shared=len(shared), jaccard=len(shared) / union,
                cell_total_r=tot_r, cell_profile_r=float(rs.mean()), cell_profile_r_median=float(np.median(rs)), cell_profile_rho=float(rh.mean()),
                profile_cells=len(rs), profile_genes=len(keep), genes=len(genes), gene_spearman=g_s, gene_pearson=g_p)

COLS = ["cr_cells", "star_cells", "shared", "jaccard", "cell_total_r", "cell_profile_r", "cell_profile_r_median", "cell_profile_rho",
        "profile_cells", "profile_genes", "genes", "gene_spearman", "gene_pearson"]
def emit(ds, level, name, r):
    print("\t".join([ds, level, name] + [f"{r[c]:.8f}" if isinstance(r[c], float) else str(r[c]) for c in COLS]), flush=True)

spec = json.load(open(sys.argv[1]))
print("\t".join(["dataset", "level", "sample"] + COLS), flush=True)
for ds in spec:
    tr = None
    if ds.get("translation"):
        tr = {}
        for l in lines(ds["translation"]):
            p = l.split("\t") if "\t" in l else l.split()
            if len(p) >= 2: tr[p[0]] = p[1]
    loaded, per = [], []
    for s in ds["samples"]:
        crm, crb, crg = read_mex(s["cr"], tr); crm, crb = collapse(crm, crb)
        stm, stb, stg = read_h5ad(s["star_h5ad"], crb, s["cr"]) if "star_h5ad" in s else read_mex(s["star"], tr)
        stm, stb = collapse(stm, stb)
        r = compare((crm, crb, crg), (stm, stb, stg))
        loaded.append(((crm, crb, crg), (stm, stb, stg))); per.append(r)
        emit(ds["name"], "per-sample", s["name"], r)
    if len(loaded) > 1:
        mean = {c: (float(np.mean([p[c] for p in per])) if isinstance(per[0][c], float)
                    else sum(p[c] for p in per) if c in ("cr_cells", "star_cells", "shared", "profile_cells") else "-")
                for c in COLS}
        emit(ds["name"], "per-sample", "MEAN", mean)
        def cat(side, keyf):
            ms, ks = [], []
            genes = loaded[0][side][2]
            for x in loaded:
                assert x[side][2] == genes, "feature lists differ between samples"
                ms.append(x[side][0]); ks += [keyf(k) for k in x[side][1]]
            return collapse(sp.vstack(ms).tocsr(), ks) + (genes,)
        emit(ds["name"], "pooled", "ALL", compare(cat(0, lambda k: k), cat(1, lambda k: k)))
        emit(ds["name"], "barcode", "ALL", compare(cat(0, lambda k: k[:16]), cat(1, lambda k: k[:16])))
