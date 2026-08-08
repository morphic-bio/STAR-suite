#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path
from typing import Any


PEAK_RE = re.compile(r"^(?P<chrom>.+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a MuData object with RNA and ATAC modalities from a STAR "
            "downstream h5ad or GEX MEX plus an ATAC peak MEX. Optional "
            "feature-barcode modalities such as protein/ADT can be added from "
            "10x-style MEX inputs."
        )
    )
    rna = parser.add_mutually_exclusive_group(required=True)
    rna.add_argument("--rna-h5ad", help="RNA AnnData input, usually downstream GeneFull+Velocyto h5ad")
    rna.add_argument("--rna-mex-dir", help="RNA 10x-style MEX directory")
    parser.add_argument("--atac-mex-dir", required=True, help="ATAC peak 10x-style MEX directory")
    parser.add_argument("--protein-mex-dir", help="Optional protein/ADT 10x-style MEX directory")
    parser.add_argument(
        "--protein-feature-ref",
        help=(
            "Optional ADT feature reference CSV/TSV. Columns such as id, name, "
            "sequence, feature_type, target_gene, clone, vendor, and "
            "isotype_control are preserved when present."
        ),
    )
    parser.add_argument(
        "--protein-normalization",
        choices=["none", "clr"],
        default="none",
        help="Optional derived protein normalization layer to add (default: none)",
    )
    parser.add_argument("--guide-mex-dir", help="Optional CRISPR guide/GDO 10x-style MEX directory")
    parser.add_argument(
        "--guide-feature-ref",
        help="Optional guide feature reference CSV/TSV merged into guide.var when present",
    )
    parser.add_argument("--guide-source", default="", help="Provenance label for the guide MEX input")
    parser.add_argument("--hash-mex-dir", help="Optional HTO/hash 10x-style MEX directory")
    parser.add_argument(
        "--hash-feature-ref",
        help="Optional hash feature reference CSV/TSV merged into hash.var when present",
    )
    parser.add_argument("--hash-source", default="", help="Provenance label for the hash MEX input")
    parser.add_argument("--state-mex-dir", help="Optional viral/LARRY/custom state 10x-style MEX directory")
    parser.add_argument(
        "--state-feature-ref",
        help="Optional state feature reference CSV/TSV merged into state.var when present",
    )
    parser.add_argument("--state-source", default="", help="Provenance label for the state MEX input")
    parser.add_argument("--output-h5mu", required=True, help="Output .h5mu path")
    parser.add_argument("--per-barcode-metrics", help="ARC/STAR per-barcode metrics CSV/TSV")
    parser.add_argument("--metrics-barcode-column", default="barcode", help="Barcode column in metrics table")
    parser.add_argument(
        "--filtered-barcodes",
        help="Optional barcode list used to set is_cell=true for matching barcodes",
    )
    parser.add_argument(
        "--subset-to-filtered-barcodes",
        action="store_true",
        help=(
            "Restrict the joined MuData to barcodes present in --filtered-barcodes "
            "instead of only marking them as is_cell"
        ),
    )
    parser.add_argument(
        "--hash-demux-assignments",
        help="Optional hash demux assignments TSV to merge into mdata.obs",
    )
    parser.add_argument(
        "--all-barcodes-are-cells",
        action="store_true",
        help="Set is_cell=true for every barcode in the joined MuData object",
    )
    parser.add_argument(
        "--metrics-is-cell-column",
        default="is_cell",
        help="Metrics column used for is_cell when present",
    )
    parser.add_argument(
        "--strip-barcode-suffix",
        action="store_true",
        help="Strip a trailing -1 from RNA, ATAC, metrics, and filtered barcode names before joining",
    )
    parser.add_argument(
        "--require-rna-velocyto-layers",
        action="store_true",
        help="Require RNA layers spliced, unspliced, and ambiguous",
    )
    parser.add_argument(
        "--allow-empty-barcode-intersection",
        action="store_true",
        help=(
            "Permit writing a zero-observation MuData object when RNA and ATAC "
            "barcodes do not overlap. Intended for sparse smoke-test filtered outputs."
        ),
    )
    parser.add_argument(
        "--y-removal-enabled",
        choices=["true", "false", "unknown"],
        default="unknown",
        help="Y-removal mode recorded in uns metadata",
    )
    parser.add_argument("--cell-call-source", default="unknown")
    parser.add_argument("--rna-source", default="")
    parser.add_argument("--atac-source", default="")
    parser.add_argument("--protein-source", default="")
    parser.add_argument("--fragments-source", default="")
    parser.add_argument("--peaks-source", default="")
    parser.add_argument("--evidence-source", default="")
    parser.add_argument("--metadata-json", help="Optional JSON object merged into mdata.uns['multiome']")
    return parser.parse_args()


def import_deps():
    try:
        import anndata as ad
        import mudata as md
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp
        from scipy.io import mmread
    except ImportError as exc:
        raise SystemExit(
            "Missing Python dependency for MuData construction. Install anndata, "
            "mudata, pandas, numpy, and scipy in the execution environment. "
            f"Original error: {exc}"
        ) from exc
    return ad, md, np, pd, sp, mmread


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def resolve_required(base: Path, stem: str) -> Path:
    plain = base / stem
    gz = base / f"{stem}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing required file: {plain}(.gz)")


def normalize_barcode(barcode: str, strip_suffix: bool) -> str:
    barcode = str(barcode)
    if strip_suffix and barcode.endswith("-1"):
        return barcode[:-2]
    return barcode


def read_barcodes(path: Path, strip_suffix: bool) -> list[str]:
    barcodes: list[str] = []
    with open_text(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                barcodes.append(normalize_barcode(line.split("\t")[0], strip_suffix))
    return barcodes


def read_features(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with open_text(path, "r") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line:
                rows.append(line.split("\t"))
    return rows


def ensure_unique(values: list[str], label: str) -> None:
    seen: set[str] = set()
    dupes: list[str] = []
    for value in values:
        if value in seen:
            dupes.append(value)
            if len(dupes) >= 5:
                break
        seen.add(value)
    if dupes:
        preview = ", ".join(dupes)
        raise ValueError(f"Duplicate {label} after normalization: {preview}")


def build_var(features: list[list[str]], modality: str, pd: Any):
    ids = [row[0] for row in features]
    names = [row[1] if len(row) > 1 else row[0] for row in features]
    types = [row[2] if len(row) > 2 else modality for row in features]
    ensure_unique(ids, f"{modality} feature ids")

    var = pd.DataFrame(index=pd.Index(ids, name=None))
    var["feature_ids"] = ids
    if modality == "rna":
        var["gene_symbols"] = names
    elif modality == "atac":
        var["peak_ids"] = names
    else:
        var["feature_names"] = names
    var["feature_types"] = types

    chroms: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    have_any_coords = False
    for row in features:
        chrom = ""
        start = -1
        end = -1
        if len(row) >= 6:
            chrom = row[3]
            try:
                start = int(row[4])
                end = int(row[5])
                have_any_coords = True
            except ValueError:
                start = -1
                end = -1
        elif modality == "atac":
            match = PEAK_RE.match(row[0])
            if match:
                chrom = match.group("chrom")
                start = int(match.group("start"))
                end = int(match.group("end"))
                have_any_coords = True
        chroms.append(chrom)
        starts.append(start)
        ends.append(end)

    if modality == "atac" or have_any_coords:
        var["chrom"] = chroms
        var["chromStart"] = starts
        var["chromEnd"] = ends
    return var


def normalize_ref_column(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def load_feature_ref(path: Path, pd: Any):
    if path.name.lower().endswith((".csv", ".csv.gz")):
        frame = pd.read_csv(path)
    else:
        frame = pd.read_csv(path, sep="\t")
    frame = frame.rename(columns={column: normalize_ref_column(column) for column in frame.columns})
    if "id" not in frame.columns:
        for candidate in ("feature_id", "feature_ids", "gene_id"):
            if candidate in frame.columns:
                frame = frame.rename(columns={candidate: "id"})
                break
    if "name" not in frame.columns:
        for candidate in ("feature_name", "feature_names", "target_name"):
            if candidate in frame.columns:
                frame = frame.rename(columns={candidate: "name"})
                break
    if "id" in frame.columns:
        frame["id"] = frame["id"].astype(str)
    if "name" in frame.columns:
        frame["name"] = frame["name"].astype(str)
    return frame


def merge_feature_ref(var, feature_ref, pd: Any):
    if feature_ref is None:
        return var
    merged = var.copy()
    ref = feature_ref.copy()
    join_key: tuple[str, str] | None = None
    if "id" in ref.columns and "feature_ids" in merged.columns:
        join_key = ("feature_ids", "id")
    elif "name" in ref.columns and "feature_names" in merged.columns:
        join_key = ("feature_names", "name")
    if join_key is None:
        return merged

    ref_columns = [
        column
        for column in ref.columns
        if column
        not in {
            "read",
            "pattern",
        }
    ]
    local_column, ref_column = join_key
    ref = ref.drop_duplicates(subset=[ref_column], keep="first")
    ref = ref.set_index(ref_column, drop=False)
    mapped = ref.reindex(merged[local_column].astype(str))
    for column in ref_columns:
        if column in {"id", "name"}:
            continue
        if column in merged.columns:
            continue
        merged[column] = pd.Series(mapped[column].values, index=merged.index).fillna("")
    if "tag_sequence" not in merged.columns:
        for candidate in ("sequence", "barcode", "adt_sequence"):
            if candidate in merged.columns:
                merged["tag_sequence"] = merged[candidate]
                break
    return merged


def read_mex_as_anndata(mex_dir: Path, modality: str, strip_suffix: bool, ad: Any, pd: Any, sp: Any, mmread: Any):
    features_path = resolve_required(mex_dir, "features.tsv")
    barcodes_path = resolve_required(mex_dir, "barcodes.tsv")
    matrix_path = resolve_required(mex_dir, "matrix.mtx")

    features = read_features(features_path)
    barcodes = read_barcodes(barcodes_path, strip_suffix)
    ensure_unique(barcodes, f"{modality} barcodes")

    matrix = mmread(matrix_path)
    if not sp.issparse(matrix):
        matrix = sp.coo_matrix(matrix)
    matrix = matrix.tocsr()
    expected = (len(features), len(barcodes))
    if matrix.shape != expected:
        raise ValueError(
            f"{modality} matrix shape mismatch for {mex_dir}: "
            f"matrix={matrix.shape}, features/barcodes={expected}"
        )

    obs = pd.DataFrame(index=pd.Index(barcodes, name=None))
    obs["barcode_raw"] = barcodes
    obs["barcode_canonical"] = barcodes
    var = build_var(features, modality, pd)
    adata = ad.AnnData(X=matrix.T.tocsr(), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.uns[f"{modality}_mex_source"] = str(mex_dir)
    return adata


def row_sums(X, np: Any):
    return np.asarray(X.sum(axis=1)).ravel()


def row_nnz(X, np: Any):
    return np.asarray((X > 0).sum(axis=1)).ravel()


def reindex_anndata_obs(adata, target, modality: str, ad: Any, pd: Any, sp: Any, np: Any):
    target_index = pd.Index(target)
    indexer = adata.obs_names.get_indexer(target_index)
    present = indexer >= 0

    def reindex_matrix(matrix):
        matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
        selected = matrix[indexer[present], :].tocoo()
        target_rows = np.flatnonzero(present)
        rows = target_rows[selected.row]
        return sp.csr_matrix((selected.data, (rows, selected.col)), shape=(len(target_index), matrix.shape[1]))

    obs = adata.obs.reindex(target_index).copy()
    obs["barcode_raw"] = list(target_index)
    obs["barcode_canonical"] = list(target_index)
    obs[f"{modality}_barcode_present"] = present
    out = ad.AnnData(X=reindex_matrix(adata.X), obs=obs, var=adata.var.copy())
    for layer_name, layer in adata.layers.items():
        out.layers[layer_name] = reindex_matrix(layer)
    out.uns.update(adata.uns)
    return out


def add_feature_library_metrics(adata, modality: str, np: Any) -> None:
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    adata.obs[f"{modality}_umis"] = row_sums(counts, np)
    adata.obs[f"{modality}_features_detected"] = row_nnz(counts, np)
    if f"{modality}_barcode_present" not in adata.obs:
        adata.obs[f"{modality}_barcode_present"] = adata.obs[f"{modality}_umis"] > 0
    adata.obs[f"{modality}_module_call"] = adata.obs[f"{modality}_umis"] > 0
    top_fraction = np.zeros(adata.n_obs, dtype=float)
    if adata.n_vars > 0 and adata.n_obs > 0:
        max_counts = counts.max(axis=1)
        max_counts = np.asarray(max_counts.toarray() if hasattr(max_counts, "toarray") else max_counts).ravel()
        total = np.maximum(adata.obs[f"{modality}_umis"].to_numpy(dtype=float), 1.0)
        top_fraction = max_counts / total
    adata.obs[f"{modality}_top_feature_fraction"] = top_fraction


def add_protein_metrics(protein, np: Any) -> None:
    add_feature_library_metrics(protein, "protein", np)


def add_protein_clr_layer(protein, np: Any) -> None:
    if protein.n_vars == 0:
        protein.layers["clr"] = protein.X.copy()
        return
    counts = protein.layers["counts"] if "counts" in protein.layers else protein.X
    dense = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    dense = dense.astype("float32", copy=False)
    logged = np.log1p(dense)
    logged -= logged.mean(axis=1, keepdims=True)
    protein.layers["clr"] = logged


FEATURE_LIBRARY_MODALITIES = ("protein", "guide", "hash", "state")


def load_feature_library_modality(
    mex_dir: Path,
    modality: str,
    feature_ref_path: Path | None,
    source_label: str,
    strip_suffix: bool,
    ad: Any,
    pd: Any,
    sp: Any,
    mmread: Any,
):
    adata = read_mex_as_anndata(
        mex_dir,
        modality,
        strip_suffix,
        ad,
        pd,
        sp,
        mmread,
    )
    feature_ref = load_feature_ref(feature_ref_path, pd) if feature_ref_path else None
    adata.var = merge_feature_ref(adata.var, feature_ref, pd)
    adata.uns[f"{modality}_feature_ref"] = str(feature_ref_path) if feature_ref_path else ""
    adata.uns[f"{modality}_source"] = source_label
    adata.uns[f"{modality}_mex_dir"] = str(mex_dir)
    return adata


def build_feature_library_qc(
    adata,
    modality: str,
    mex_dir: Path,
    feature_ref_path: Path | None,
    source_label: str,
    n_source_barcodes: int,
    n_joined_cells: int,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    present = adata.obs[f"{modality}_barcode_present"].to_numpy(dtype=bool)
    module = adata.obs[f"{modality}_module_call"].to_numpy(dtype=bool)
    qc = {
        f"{modality}_mex_dir": str(mex_dir),
        f"{modality}_feature_ref": str(feature_ref_path) if feature_ref_path else "",
        f"{modality}_source": source_label,
        "n_source_barcodes": int(n_source_barcodes),
        "n_source_features": int(adata.n_vars),
        "n_joined_cells": int(n_joined_cells),
        "n_barcodes_present_on_joined_cells": int(present.sum()),
        "n_module_call_cells": int(module.sum()),
        "barcode_overlap_fraction": float(present.sum() / max(n_joined_cells, 1)),
        "module_call_rule": f"{modality}_umis > 0",
    }
    if extra:
        qc.update(extra)
    return qc


def attach_feature_library_modalities(
    modalities: dict[str, Any],
    common,
    args: argparse.Namespace,
    ad: Any,
    pd: Any,
    sp: Any,
    mmread: Any,
    np: Any,
) -> dict[str, dict[str, Any]]:
    specs: list[tuple[str, str, str | None, str, str | None]] = [
        ("protein", args.protein_mex_dir, args.protein_feature_ref, args.protein_source, "protein_normalization"),
        ("guide", args.guide_mex_dir, args.guide_feature_ref, args.guide_source, None),
        ("hash", args.hash_mex_dir, args.hash_feature_ref, args.hash_source, None),
        ("state", args.state_mex_dir, args.state_feature_ref, args.state_source, None),
    ]
    qc_summaries: dict[str, dict[str, Any]] = {}
    for modality, mex_dir_arg, feature_ref_arg, source_label, normalization_attr in specs:
        if not mex_dir_arg:
            continue
        mex_dir = Path(mex_dir_arg).resolve()
        feature_ref_path = Path(feature_ref_arg).resolve() if feature_ref_arg else None
        adata = load_feature_library_modality(
            mex_dir,
            modality,
            feature_ref_path,
            source_label,
            args.strip_barcode_suffix,
            ad,
            pd,
            sp,
            mmread,
        )
        source_barcodes = pd.Index(adata.obs_names)
        adata = reindex_anndata_obs(adata, common, modality, ad, pd, sp, np)
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        add_feature_library_metrics(adata, modality, np)
        if normalization_attr == "protein_normalization" and getattr(args, normalization_attr) == "clr":
            add_protein_clr_layer(adata, np)
        modalities[modality] = adata
        extra: dict[str, Any] = {}
        if normalization_attr == "protein_normalization":
            extra["normalization"] = args.protein_normalization
        qc_summaries[f"{modality}_qc"] = build_feature_library_qc(
            adata,
            modality,
            mex_dir,
            feature_ref_path,
            source_label,
            len(source_barcodes),
            len(common),
            extra=extra or None,
        )
    return qc_summaries


def summarize_modality_overlap(
    modalities: dict[str, Any],
    n_joined_cells: int,
    np: Any,
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "n_rna_cells": int(n_joined_cells),
        "n_atac_cells": int(n_joined_cells),
        "n_joined_cells": int(n_joined_cells),
    }
    present_masks: list[Any] = []
    for modality in FEATURE_LIBRARY_MODALITIES:
        if modality not in modalities:
            summary[f"n_{modality}_barcodes_present"] = 0
            continue
        present = modalities[modality].obs[f"{modality}_barcode_present"].to_numpy(dtype=bool)
        summary[f"n_{modality}_barcodes_present"] = int(present.sum())
        present_masks.append(present)

    guide_or_state_masks = []
    for modality in ("guide", "hash", "state"):
        if modality in modalities:
            guide_or_state_masks.append(
                modalities[modality].obs[f"{modality}_barcode_present"].to_numpy(dtype=bool)
            )
    if guide_or_state_masks:
        combined = guide_or_state_masks[0]
        for mask in guide_or_state_masks[1:]:
            combined = combined | mask
        summary["n_guide_or_state_cells"] = int(combined.sum())
    else:
        summary["n_guide_or_state_cells"] = 0

    if "protein" in modalities:
        protein_present = modalities["protein"].obs["protein_barcode_present"].to_numpy(dtype=bool)
        summary["n_rna_atac_protein_overlap"] = int(protein_present.sum())
    else:
        summary["n_rna_atac_protein_overlap"] = 0

    if present_masks:
        all_present = present_masks[0]
        for mask in present_masks[1:]:
            all_present = all_present & mask
        summary["n_feature_library_overlap"] = int(all_present.sum())
    else:
        summary["n_feature_library_overlap"] = 0

    if "protein" in modalities and guide_or_state_masks:
        protein_present = modalities["protein"].obs["protein_barcode_present"].to_numpy(dtype=bool)
        identity_present = guide_or_state_masks[0]
        for mask in guide_or_state_masks[1:]:
            identity_present = identity_present | mask
        summary["n_four_factor_overlap"] = int((protein_present & identity_present).sum())
    else:
        summary["n_four_factor_overlap"] = 0
    return summary


def load_rna_h5ad(path: Path, strip_suffix: bool, ad: Any, pd: Any, sp: Any):
    rna = ad.read_h5ad(path)
    new_names = [normalize_barcode(value, strip_suffix) for value in rna.obs_names.astype(str)]
    ensure_unique(new_names, "rna obs names")
    rna.obs_names = new_names
    if not sp.issparse(rna.X):
        rna.X = sp.csr_matrix(rna.X)
    else:
        rna.X = rna.X.tocsr()
    if "counts" not in rna.layers:
        rna.layers["counts"] = rna.X.copy()
    if "barcode_raw" not in rna.obs:
        rna.obs["barcode_raw"] = list(rna.obs_names)
    if "barcode_canonical" not in rna.obs:
        rna.obs["barcode_canonical"] = list(rna.obs_names)
    if "gene_symbols" not in rna.var:
        if "gene_ids" in rna.var:
            rna.var["gene_symbols"] = rna.var_names.astype(str)
        elif "_index" in rna.var:
            rna.var["gene_symbols"] = rna.var["_index"].astype(str)
        else:
            rna.var["gene_symbols"] = rna.var_names.astype(str)
    if "feature_types" not in rna.var:
        rna.var["feature_types"] = "Gene Expression"
    rna.uns["rna_h5ad_source"] = str(path)
    return rna


def read_table_auto(path: Path, pd: Any):
    name = path.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def read_barcode_set(path: Path, strip_suffix: bool) -> set[str]:
    values: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                values.add(normalize_barcode(line.split("\t")[0], strip_suffix))
    return values


def attach_hash_demux_assignments(obs, assignments_path: Path, strip_suffix: bool, pd: Any):
    frame = read_table_auto(assignments_path, pd)
    barcode_column = None
    for candidate in ("barcode", "barcode_raw", "barcode_canonical"):
        if candidate in frame.columns:
            barcode_column = candidate
            break
    if barcode_column is None:
        raise ValueError(
            f"Hash demux assignments table {assignments_path} is missing a barcode column"
        )
    frame[barcode_column] = frame[barcode_column].map(
        lambda value: normalize_barcode(value, strip_suffix)
    )
    frame = frame.drop_duplicates(subset=[barcode_column], keep="first")
    frame = frame.set_index(barcode_column).reindex(obs.index)
    for column in frame.columns:
        if column in obs.columns:
            continue
        obs[column] = frame[column].values
    return obs


def normalize_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def numeric_series(frame, column: str, pd: Any):
    if column not in frame.columns:
        return None
    return pd.to_numeric(frame[column], errors="coerce").fillna(0)


def attach_metrics(adata, metrics, columns: list[str]) -> None:
    for column in columns:
        if column not in metrics.columns:
            continue
        adata.obs[column] = metrics[column].values


def add_metrics_and_calls(
    rna,
    atac,
    args: argparse.Namespace,
    common,
    pd: Any,
    np: Any,
) -> Any:
    obs = pd.DataFrame(index=common)
    obs["barcode_raw"] = list(common)
    obs["barcode_canonical"] = list(common)

    metrics = None
    if args.per_barcode_metrics:
        metrics_path = Path(args.per_barcode_metrics).resolve()
        metrics = read_table_auto(metrics_path, pd)
        if args.metrics_barcode_column not in metrics.columns:
            raise ValueError(
                f"Missing metrics barcode column {args.metrics_barcode_column!r} in {metrics_path}"
            )
        metrics[args.metrics_barcode_column] = metrics[args.metrics_barcode_column].map(
            lambda value: normalize_barcode(value, args.strip_barcode_suffix)
        )
        metrics = metrics.drop_duplicates(subset=[args.metrics_barcode_column], keep="first")
        metrics = metrics.set_index(args.metrics_barcode_column).reindex(common)
        if (
            "atac_peak_fraction" not in metrics.columns
            and "atac_peak_region_fragments" in metrics.columns
            and "atac_fragments" in metrics.columns
        ):
            numerator = pd.to_numeric(metrics["atac_peak_region_fragments"], errors="coerce")
            denominator = pd.to_numeric(metrics["atac_fragments"], errors="coerce")
            metrics["atac_peak_fraction"] = np.where(denominator > 0, numerator / denominator, 0.0)
        for column in metrics.columns:
            obs[column] = metrics[column].values

    if args.all_barcodes_are_cells:
        obs["is_cell"] = True
    elif args.filtered_barcodes:
        filtered = read_barcode_set(Path(args.filtered_barcodes).resolve(), args.strip_barcode_suffix)
        obs["is_cell"] = [barcode in filtered for barcode in common]
    elif metrics is not None and args.metrics_is_cell_column in metrics.columns:
        obs["is_cell"] = [normalize_bool(value) for value in metrics[args.metrics_is_cell_column].values]
    elif "is_cell" in rna.obs:
        obs["is_cell"] = [normalize_bool(value) for value in rna.obs.loc[common, "is_cell"].values]
    else:
        obs["is_cell"] = False

    obs["cell_call_source"] = args.cell_call_source
    obs["call_provenance"] = args.cell_call_source
    if metrics is not None and args.metrics_is_cell_column in metrics.columns:
        obs["arc_is_cell"] = [
            normalize_bool(value) for value in metrics[args.metrics_is_cell_column].values
        ]
    elif "arc_is_cell" not in obs.columns:
        obs["arc_is_cell"] = False

    if "gex_module_call" not in obs.columns:
        gex_umis = numeric_series(obs, "gex_umis_count", pd)
        obs["gex_module_call"] = (gex_umis > 0).values if gex_umis is not None else obs["is_cell"].values
    if "atac_module_call" not in obs.columns:
        atac_cutsites = numeric_series(obs, "atac_peak_region_cutsites", pd)
        if atac_cutsites is not None:
            obs["atac_module_call"] = (atac_cutsites >= 1).values
        else:
            atac_fragments = numeric_series(obs, "atac_fragments", pd)
            obs["atac_module_call"] = (atac_fragments > 0).values if atac_fragments is not None else False
    if "atac_low_targeting" not in obs.columns:
        obs["atac_low_targeting"] = False
    if "effective_atac_module_call" not in obs.columns:
        obs["effective_atac_module_call"] = [
            normalize_bool(call) and not normalize_bool(low)
            for call, low in zip(obs["atac_module_call"].values, obs["atac_low_targeting"].values)
        ]
    if "gex_rescue_eligible" not in obs.columns:
        obs["gex_rescue_eligible"] = [
            normalize_bool(atac) and not normalize_bool(gex)
            for atac, gex in zip(obs["effective_atac_module_call"].values, obs["gex_module_call"].values)
        ]

    rna_columns = [
        "barcode_raw",
        "barcode_canonical",
        "is_cell",
        "cell_call_source",
        "call_provenance",
        "gex_module_call",
        "gex_rescue_eligible",
        "arc_is_cell",
        "gex_umis_count",
        "gex_genes_count",
        "gex_exonic_umis",
        "gex_intronic_umis",
    ]
    atac_columns = [
        "barcode_raw",
        "barcode_canonical",
        "is_cell",
        "cell_call_source",
        "call_provenance",
        "atac_fragments",
        "atac_TSS_fragments",
        "atac_peak_region_fragments",
        "atac_peak_region_cutsites",
        "atac_peak_fraction",
        "atac_module_call",
        "effective_atac_module_call",
        "gex_rescue_eligible",
        "atac_low_targeting",
        "arc_is_cell",
    ]
    attach_metrics(rna, obs, rna_columns)
    attach_metrics(atac, obs, atac_columns)
    return obs


def validate_velocyto_layers(rna) -> None:
    required = ["spliced", "unspliced", "ambiguous"]
    missing = [layer for layer in required if layer not in rna.layers]
    if missing:
        raise ValueError(f"RNA AnnData is missing required Velocyto layers: {', '.join(missing)}")
    for layer in required:
        if rna.layers[layer].shape != rna.X.shape:
            raise ValueError(
                f"RNA layer {layer!r} shape {rna.layers[layer].shape} does not match X {rna.X.shape}"
            )


def load_metadata_json(path: str | None) -> dict[str, Any]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("--metadata-json must contain a JSON object")
    return data


def main() -> None:
    args = parse_args()
    ad, md, np, pd, sp, mmread = import_deps()

    atac = read_mex_as_anndata(
        Path(args.atac_mex_dir).resolve(),
        "atac",
        args.strip_barcode_suffix,
        ad,
        pd,
        sp,
        mmread,
    )
    if args.rna_h5ad:
        rna = load_rna_h5ad(Path(args.rna_h5ad).resolve(), args.strip_barcode_suffix, ad, pd, sp)
    else:
        rna = read_mex_as_anndata(
            Path(args.rna_mex_dir).resolve(),
            "rna",
            args.strip_barcode_suffix,
            ad,
            pd,
            sp,
            mmread,
        )

    if args.require_rna_velocyto_layers:
        validate_velocyto_layers(rna)

    common = rna.obs_names.intersection(atac.obs_names)
    if args.subset_to_filtered_barcodes:
        if not args.filtered_barcodes:
            raise ValueError("--subset-to-filtered-barcodes requires --filtered-barcodes")
        filtered = read_barcode_set(Path(args.filtered_barcodes).resolve(), args.strip_barcode_suffix)
        common = common.intersection(filtered)
    if len(common) == 0 and not args.allow_empty_barcode_intersection:
        raise ValueError("RNA and ATAC modalities have no overlapping barcodes")
    if len(common) == 0:
        print(
            "WARNING: RNA and ATAC modalities have no overlapping barcodes; "
            "writing a zero-observation MuData object because "
            "--allow-empty-barcode-intersection was set.",
            file=sys.stderr,
        )
    rna = rna[common, :].copy()
    atac = atac[common, :].copy()

    obs = add_metrics_and_calls(rna, atac, args, common, pd, np)
    if args.hash_demux_assignments:
        obs = attach_hash_demux_assignments(
            obs,
            Path(args.hash_demux_assignments).resolve(),
            args.strip_barcode_suffix,
            pd,
        )
    modalities = {"rna": rna, "atac": atac}
    qc_summaries = attach_feature_library_modalities(
        modalities,
        common,
        args,
        ad,
        pd,
        sp,
        mmread,
        np,
    )

    mdata = md.MuData(modalities)
    mdata.update()
    for column in obs.columns:
        mdata.obs[column] = obs[column].values

    protein = modalities.get("protein")
    guide = modalities.get("guide")
    hash_mod = modalities.get("hash")
    state = modalities.get("state")
    overlap = summarize_modality_overlap(modalities, len(common), np)

    metadata = {
        "builder": "scripts/build_multiome_mudata.py",
        "rna_input": str(Path(args.rna_h5ad or args.rna_mex_dir).resolve()),
        "atac_mex_dir": str(Path(args.atac_mex_dir).resolve()),
        "per_barcode_metrics": str(Path(args.per_barcode_metrics).resolve()) if args.per_barcode_metrics else "",
        "rna_source": args.rna_source,
        "atac_source": args.atac_source,
        "protein_source": args.protein_source,
        "guide_source": args.guide_source,
        "hash_source": args.hash_source,
        "state_source": args.state_source,
        "hash_demux_assignments": (
            str(Path(args.hash_demux_assignments).resolve()) if args.hash_demux_assignments else ""
        ),
        "subset_to_filtered_barcodes": bool(args.subset_to_filtered_barcodes),
        "fragments_source": args.fragments_source,
        "peaks_source": args.peaks_source,
        "evidence_source": args.evidence_source,
        "cell_call_source": args.cell_call_source,
        "y_removal_enabled": args.y_removal_enabled,
        "n_obs": int(mdata.n_obs),
        "rna_n_vars": int(rna.n_vars),
        "atac_n_vars": int(atac.n_vars),
        "protein_n_vars": int(protein.n_vars) if protein is not None else 0,
        "guide_n_vars": int(guide.n_vars) if guide is not None else 0,
        "hash_n_vars": int(hash_mod.n_vars) if hash_mod is not None else 0,
        "state_n_vars": int(state.n_vars) if state is not None else 0,
    }
    metadata.update(overlap)
    metadata.update(load_metadata_json(args.metadata_json))
    mdata.uns["multiome"] = metadata
    for qc_key, qc_value in qc_summaries.items():
        mdata.uns[qc_key] = qc_value

    output = Path(args.output_h5mu).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(mdata, "write_h5mu"):
        mdata.write_h5mu(str(output))
    else:
        mdata.write(str(output))

    print(f"Wrote {output}")
    print(f"obs={mdata.n_obs}")
    print(f"mods={', '.join(mdata.mod.keys())}")
    print(f"rna_vars={rna.n_vars}")
    print(f"atac_vars={atac.n_vars}")
    for modality in FEATURE_LIBRARY_MODALITIES:
        adata = modalities.get(modality)
        if adata is None:
            continue
        qc_key = f"{modality}_qc"
        print(f"{modality}_vars={adata.n_vars}")
        if qc_key in mdata.uns:
            print(f"{modality}_barcodes_present={mdata.uns[qc_key]['n_barcodes_present_on_joined_cells']}")


if __name__ == "__main__":
    main()
