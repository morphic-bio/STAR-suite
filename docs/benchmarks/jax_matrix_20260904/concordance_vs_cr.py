#!/usr/bin/env python3
"""Concordance of STAR-Flex or cyto output against Cell Ranger 9.0.1.

The default mapping preserves the original four-sample JAX analysis. For the
public 320k benchmark, pass Cell Ranger's multi config: each public sample
groups two probe-barcode IDs, and those two output matrices are combined before
comparison with the corresponding Cell Ranger sample.

Examples:
  concordance_vs_cr.py STAR_RUN/per_sample star-bgzf \
    --cr-root CR_RUN/outs/per_sample_outs --cr-config cellranger.config.csv
  concordance_vs_cr.py CYTO_RUN/counts cyto-cbq --input-kind cyto \
    --cr-root CR_RUN/outs/per_sample_outs --cr-config cellranger.config.csv
"""
import argparse
import csv
import gzip
import os
from collections import Counter
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.stats

DEFAULT_CRROOT = Path("/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs")
DEFAULT_GROUPS = [
    ("WT-Day-7", ["BC004"]),
    ("PAX6-PTC-D9-Day7", ["BC006"]),
    ("WT-Day-8", ["BC007"]),
    ("PAX6-PTC-D9-Day8", ["BC008"]),
]


def read_groups(config):
    """Return [(Cell Ranger sample_id, [probe barcode IDs]), ...]."""
    if config is None:
        return DEFAULT_GROUPS
    rows = list(csv.reader(Path(config).open(newline="", encoding="utf-8")))
    try:
        start = next(i for i, row in enumerate(rows)
                     if row and row[0].strip() == "[samples]")
    except StopIteration as exc:
        raise SystemExit(f"no [samples] section in {config}") from exc
    if start + 1 >= len(rows) or rows[start + 1][:2] != [
            "sample_id", "probe_barcode_ids"]:
        raise SystemExit(f"unexpected [samples] header in {config}")
    groups = []
    for row in rows[start + 2:]:
        if not row or not row[0].strip() or row[0].lstrip().startswith("["):
            break
        groups.append((row[0].strip(),
                       [x.strip() for x in row[1].split("|") if x.strip()]))
    if not groups:
        raise SystemExit(f"no sample rows in {config}")
    return groups


def open_maybe_gz(directory, stem):
    plain, gz = directory / stem, directory / (stem + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    return plain.open("rt", encoding="utf-8")


DEFAULT_TAG_MAP = Path("/home/lhhung/jax_stage_20260903/ref/sample_whitelist_full_16.tsv")


def read_tag_map(path):
    """Probe-barcode whitelist: tag id (BC004) <-> listed 8-mer, both directions."""
    tag_map = {}
    if path is None or not Path(path).exists():
        return tag_map
    for line in Path(path).open():
        parts = line.split()
        if len(parts) >= 2:
            tag_map[parts[1]] = parts[0]
    return tag_map


def cell_key(raw_barcode, tag=None, tag_map=None):
    """A cell is the 16-base barcode AND its sample tag. Two tags in one GEM are two
    cells, so the key never truncates to CB16 when a tag is known. STAR and cyto
    per-tag outputs carry the tag in the directory/file name; Cell Ranger writes
    CB16 + tag8 (+ "-1") and the 8-mer is mapped back to the tag id."""
    raw = raw_barcode.strip().split("-")[0]
    cb16 = raw[:16]
    if tag is None:
        tag8 = raw[16:24]
        tag = (tag_map or {}).get(tag8, tag8) if tag8 else ""
    return f"{cb16}|{tag}" if tag else cb16


def read_mex(directory, tag=None, tag_map=None):
    """Return (cells x genes CSR, cell keys, gene ids). Keys are CB16|tag (see cell_key)."""
    directory = Path(directory)
    matrix_path = directory / "matrix.mtx"
    if not matrix_path.exists():
        matrix_path = directory / "matrix.mtx.gz"
    matrix = scipy.io.mmread(str(matrix_path))
    with open_maybe_gz(directory, "barcodes.tsv") as handle:
        barcodes = [cell_key(line, tag, tag_map) for line in handle]
    with open_maybe_gz(directory, "features.tsv") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    # Genes are keyed by Ensembl ID (column 0). STAR-Flex writes the ID in
    # both columns and Cell Ranger names all but 25 genes by symbol, so a
    # symbol match intersects only the 25 unnamed genes.
    symbols = [row[0] for row in rows]
    return sp.csr_matrix(matrix).T.tocsr(), barcodes, symbols


def collapse_duplicate_cells(matrix, barcodes):
    """Sum rows with the same cell key. With CB16|tag keys this only merges true
    duplicates; it must never merge cells of different tags (the pre-2026-09-06
    version truncated keys to CB16 and did exactly that for multi-tag samples)."""
    if len(set(barcodes)) == len(barcodes):
        return matrix, barcodes
    unique = sorted(set(barcodes))
    index = {barcode: i for i, barcode in enumerate(unique)}
    row_map = sp.csr_matrix(
        (np.ones(len(barcodes), dtype=np.int8),
         ([index[b] for b in barcodes], np.arange(len(barcodes)))),
        shape=(len(unique), len(barcodes)),
    )
    return (row_map @ matrix).tocsr(), unique


def combine_parts(parts):
    matrices, barcode_parts, symbols = [], [], None
    for matrix, barcodes, part_symbols in parts:
        if symbols is None:
            symbols = part_symbols
        elif symbols != part_symbols:
            raise SystemExit(
                "gene axes differ across grouped probe-barcode outputs")
        matrices.append(matrix)
        barcode_parts.extend(barcodes)
    if not matrices:
        raise SystemExit("no matrices to combine")
    combined = (matrices[0] if len(matrices) == 1
                else sp.vstack(matrices, format="csr"))
    combined, barcode_parts = collapse_duplicate_cells(
        combined, barcode_parts)
    return combined, barcode_parts, symbols


def read_star_group(root, barcode_ids):
    parts = []
    for bc in barcode_ids:
        directory = root / bc / "Gene" / "filtered"
        if not directory.exists():
            raise FileNotFoundError(directory)
        parts.append(read_mex(directory, tag=bc))
    return combine_parts(parts)


def read_cyto_group(root, barcode_ids):
    import anndata as ad

    parts = []
    for bc in barcode_ids:
        path = root / f"{bc}.filt.h5ad"
        if not path.exists():
            raise FileNotFoundError(path)
        data = ad.read_h5ad(path)
        parts.append((sp.csr_matrix(data.X),
                      [cell_key(str(x), tag=bc) for x in data.obs_names],
                      [str(x) for x in data.var_names]))
    return combine_parts(parts)


def unique_symbol_index(symbols):
    counts = Counter(symbols)
    return {symbol: i for i, symbol in enumerate(symbols)
            if counts[symbol] == 1}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "output_root", type=Path,
        help="STAR per_sample directory or cyto counts directory")
    parser.add_argument("label", nargs="?")
    parser.add_argument(
        "--input-kind", choices=("auto", "star", "cyto"), default="auto")
    parser.add_argument("--tag-map", type=Path, default=DEFAULT_TAG_MAP,
                        help="probe-barcode whitelist (tag id, 8-mer) used to key Cell Ranger cells by CB16|tag")
    parser.add_argument(
        "--cr-root", type=Path,
        default=Path(os.environ.get("CRROOT", DEFAULT_CRROOT)))
    parser.add_argument(
        "--cr-config", type=Path,
        help="Cell Ranger multi CSV; derives grouped probe-barcode mapping")
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.output_root
    label = args.label or root.name
    tag_map = read_tag_map(args.tag_map)
    kind = args.input_kind
    if kind == "auto":
        kind = "cyto" if any(root.glob("BC*.filt.h5ad")) else "star"
    groups = read_groups(args.cr_config)
    # cyto indexes genes by symbol; translate through Cell Ranger's feature
    # table (same reference for every sample) so all tools share the ID space.
    first_cr = (args.cr_root / groups[0][0] / "count" /
                "sample_filtered_feature_bc_matrix")
    with open_maybe_gz(first_cr, "features.tsv") as handle:
        cr_rows = [line.rstrip("\n").split("\t") for line in handle]
    sym_counts = Counter(r[1] for r in cr_rows if len(r) > 1)
    sym2id = {r[1]: r[0] for r in cr_rows
              if len(r) > 1 and sym_counts[r[1]] == 1}

    print(f"# concordance vs Cell Ranger 9.0.1 - {label}")
    print(f"# query_kind={kind} query_outputs={root}")
    print("# genes matched by Ensembl ID via Cell Ranger's feature table")
    print(f"# cr_root={args.cr_root} "
          f"cr_config={args.cr_config or 'JAX defaults'}")
    print(f"{'sample':28s} {'tags':13s} {'cells query/CR':>18s} "
          f"{'Jaccard':>8s} {'cellPear':>9s} {'r>=1-1e6':>9s} "
          f"{'p01':>8s} {'worst':>8s} {'cellSpear':>10s} "
          f"{'genePear':>9s}")

    jaccards, cell_pears, cell_spears, gene_pears = [], [], [], []
    fracs, p01s, worsts = [], [], []
    for sample, barcode_ids in groups:
        tags_label = "+".join(barcode_ids)
        try:
            if kind == "star":
                query_x, query_bcs, query_syms = read_star_group(
                    root, barcode_ids)
            else:
                query_x, query_bcs, query_syms = read_cyto_group(
                    root, barcode_ids)
                query_syms = [sym2id.get(sym, sym) for sym in query_syms]
        except FileNotFoundError as exc:
            print(f"{sample:28s} {tags_label:13s} MISSING {exc}")
            continue
        cr_dir = (args.cr_root / sample / "count" /
                  "sample_filtered_feature_bc_matrix")
        cr_x, cr_bcs, cr_syms = read_mex(cr_dir, tag_map=tag_map)
        cr_x, cr_bcs = collapse_duplicate_cells(cr_x, cr_bcs)

        query_sym_i = unique_symbol_index(query_syms)
        cr_sym_i = unique_symbol_index(cr_syms)
        common_g = sorted(set(query_sym_i) & set(cr_sym_i))
        gi_query = [query_sym_i[g] for g in common_g]
        gi_cr = [cr_sym_i[g] for g in common_g]

        set_query, set_cr = set(query_bcs), set(cr_bcs)
        common_c = sorted(set_query & set_cr)
        union_c = set_query | set_cr
        jaccard = (len(common_c) / len(union_c)
                   if union_c else float("nan"))
        ci_query = {barcode: i for i, barcode in enumerate(query_bcs)}
        ci_cr = {barcode: i for i, barcode in enumerate(cr_bcs)}
        query_common = query_x[
            [ci_query[b] for b in common_c]][:, gi_query]
        cr_common = cr_x[[ci_cr[b] for b in common_c]][:, gi_cr]

        pears, spears = [], []
        for i in range(len(common_c)):
            a = query_common.getrow(i).toarray().ravel()
            b = cr_common.getrow(i).toarray().ravel()
            if a.std() > 0 and b.std() > 0:
                pears.append(scipy.stats.pearsonr(a, b)[0])
                spears.append(scipy.stats.spearmanr(a, b)[0])
        pears_a = np.asarray(pears) if pears else np.array([np.nan])
        cell_pear = float(np.nanmedian(pears_a))
        cell_spear = (
            float(np.nanmedian(spears)) if spears else float("nan"))
        frac_exact = float(np.nanmean(pears_a >= 0.999999))
        p01 = float(np.nanpercentile(pears_a, 1))
        worst = float(np.nanmin(pears_a))

        gene_query = np.asarray(query_common.sum(axis=0)).ravel()
        gene_cr = np.asarray(cr_common.sum(axis=0)).ravel()
        gene_pear = (
            float(scipy.stats.pearsonr(gene_query, gene_cr)[0])
            if gene_query.std() > 0 and gene_cr.std() > 0
            else float("nan"))

        print(f"{sample:28s} {tags_label:13s} "
              f"{len(set_query):8d}/{len(set_cr):<9d} "
              f"{jaccard:8.4f} {cell_pear:9.6f} {frac_exact:9.4f} "
              f"{p01:8.4f} {worst:8.4f} {cell_spear:10.6f} "
              f"{gene_pear:9.6f}")
        jaccards.append(jaccard)
        cell_pears.append(cell_pear)
        cell_spears.append(cell_spear)
        gene_pears.append(gene_pear)
        fracs.append(frac_exact)
        p01s.append(p01)
        worsts.append(worst)

    if jaccards:
        print(f"{'median':28s} {'':13s} {'':18s} "
              f"{np.nanmedian(jaccards):8.4f} "
              f"{np.nanmedian(cell_pears):9.6f} "
              f"{np.nanmedian(fracs):9.4f} "
              f"{np.nanmedian(p01s):8.4f} {np.nanmin(worsts):8.4f} "
              f"{np.nanmedian(cell_spears):10.6f} "
              f"{np.nanmedian(gene_pears):9.6f}")
        print("# r>=1-1e6 is the fraction of common cells at Pearson "
              ">= 0.999999; p01 is the 1st percentile and worst is the minimum.")
        print("# Cells called by only one tool are excluded from correlation "
              "columns; the Jaccard column reports that disagreement.")


if __name__ == "__main__":
    main()
