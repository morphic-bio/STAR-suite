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
            "downstream h5ad or GEX MEX plus an ATAC peak MEX."
        )
    )
    rna = parser.add_mutually_exclusive_group(required=True)
    rna.add_argument("--rna-h5ad", help="RNA AnnData input, usually downstream GeneFull+Velocyto h5ad")
    rna.add_argument("--rna-mex-dir", help="RNA 10x-style MEX directory")
    parser.add_argument("--atac-mex-dir", required=True, help="ATAC peak 10x-style MEX directory")
    parser.add_argument("--output-h5mu", required=True, help="Output .h5mu path")
    parser.add_argument("--per-barcode-metrics", help="ARC/STAR per-barcode metrics CSV/TSV")
    parser.add_argument("--metrics-barcode-column", default="barcode", help="Barcode column in metrics table")
    parser.add_argument(
        "--filtered-barcodes",
        help="Optional barcode list used to set is_cell=true for matching barcodes",
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
    else:
        var["peak_ids"] = names
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
    if path.suffix == ".csv":
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

    mdata = md.MuData({"rna": rna, "atac": atac})
    mdata.update()
    for column in obs.columns:
        mdata.obs[column] = obs[column].values

    metadata = {
        "builder": "scripts/build_multiome_mudata.py",
        "rna_input": str(Path(args.rna_h5ad or args.rna_mex_dir).resolve()),
        "atac_mex_dir": str(Path(args.atac_mex_dir).resolve()),
        "per_barcode_metrics": str(Path(args.per_barcode_metrics).resolve()) if args.per_barcode_metrics else "",
        "rna_source": args.rna_source,
        "atac_source": args.atac_source,
        "fragments_source": args.fragments_source,
        "peaks_source": args.peaks_source,
        "evidence_source": args.evidence_source,
        "cell_call_source": args.cell_call_source,
        "y_removal_enabled": args.y_removal_enabled,
        "n_obs": int(mdata.n_obs),
        "rna_n_vars": int(rna.n_vars),
        "atac_n_vars": int(atac.n_vars),
    }
    metadata.update(load_metadata_json(args.metadata_json))
    mdata.uns["multiome"] = metadata

    output = Path(args.output_h5mu).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(mdata, "write_h5mu"):
        mdata.write_h5mu(str(output))
    else:
        mdata.write(str(output))

    print(f"Wrote {output}")
    print(f"obs={mdata.n_obs}")
    print(f"rna_vars={rna.n_vars}")
    print(f"atac_vars={atac.n_vars}")


if __name__ == "__main__":
    main()
