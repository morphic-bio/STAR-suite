#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return slug or "feature_library"


def first_existing(directory: Path, candidates):
    for candidate in candidates:
        path = directory / candidate
        if path.exists():
            return path
    return None


def read_table(path: Path) -> pd.DataFrame:
    path_str = str(path)
    if path_str.endswith(".csv") or path_str.endswith(".csv.gz"):
        return pd.read_csv(path, sep=",", header=None, comment="#")
    if path_str.endswith(".tsv") or path_str.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", header=None, comment="#")

    opener = open
    if path_str.endswith(".gz"):
        import gzip

        opener = gzip.open

    with opener(path, "rt", encoding="utf-8") as handle:
        rows = [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]
    return pd.DataFrame(rows)


def read_provenance(library_dir: Path) -> dict:
    provenance_path = library_dir / "pf_library_provenance.tsv"
    if not provenance_path.exists():
        return {}
    provenance = pd.read_csv(provenance_path, sep="\t")
    return dict(zip(provenance["key"], provenance["value"]))


def load_feature_matrix(matrix_dir: Path) -> ad.AnnData:
    matrix_path = first_existing(matrix_dir, ["matrix.mtx", "matrix.mtx.gz", "features_matrix.mtx", "features_matrix.mtx.gz"])
    barcodes_path = first_existing(matrix_dir, ["barcodes.txt", "barcodes.tsv", "barcodes.tsv.gz", "barcodes.csv"])
    features_path = first_existing(matrix_dir, ["features.txt", "features.tsv", "features.tsv.gz", "features.csv"])
    if matrix_path is None or barcodes_path is None or features_path is None:
        raise FileNotFoundError(f"Missing matrix/barcodes/features in {matrix_dir}")

    matrix = scipy.io.mmread(matrix_path).tocsr()
    barcodes = read_table(barcodes_path).iloc[:, 0].astype(str).tolist()
    features_df = read_table(features_path)
    feature_names = features_df.iloc[:, 0].astype(str).tolist()

    if matrix.shape[0] != len(feature_names):
        raise ValueError(f"{matrix_dir}: matrix rows {matrix.shape[0]} do not match features {len(feature_names)}")
    if matrix.shape[1] != len(barcodes):
        raise ValueError(f"{matrix_dir}: matrix columns {matrix.shape[1]} do not match barcodes {len(barcodes)}")

    obs = pd.DataFrame(index=pd.Index(barcodes, dtype=str, name="barcode"))
    var = pd.DataFrame(index=pd.Index(feature_names, dtype=str, name="feature_name"))
    var["feature_name"] = var.index.astype(str)
    if features_df.shape[1] > 1:
        for col in range(1, features_df.shape[1]):
            var[f"feature_field_{col}"] = features_df.iloc[:, col].astype(str).to_numpy()

    adata = ad.AnnData(X=matrix.T, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def annotate_feature_barcodes(adata: ad.AnnData, matrix_dir: Path):
    feature_per_cell = matrix_dir / "feature_per_cell.csv"
    if not feature_per_cell.exists():
        return

    df = pd.read_csv(feature_per_cell)
    if "barcode" not in df.columns:
        return

    df["barcode"] = df["barcode"].astype(str)
    df = df.set_index("barcode")

    for source_col, target_col in [
        ("num_features", "num_features"),
        ("top_feature_index", "top_feature_index"),
        ("total_deduped_umi", "total_deduped_umi"),
    ]:
        if source_col in df.columns:
            series = df[source_col].reindex(adata.obs_names)
            if pd.api.types.is_numeric_dtype(series):
                fill_value = -1 if target_col == "top_feature_index" else 0
                adata.obs[target_col] = series.fillna(fill_value).astype(int)
            else:
                adata.obs[target_col] = series.astype(str).fillna("")

    if "num_features" in adata.obs.columns:
        adata.obs["is_featured"] = (adata.obs["num_features"] > 0).astype(bool)

    if "top_feature_index" in adata.obs.columns:
        top_names = []
        feature_names = list(adata.var_names)
        for value in adata.obs["top_feature_index"]:
            if pd.isna(value):
                top_names.append("")
                continue
            index = int(value)
            if 0 <= index < len(feature_names):
                top_names.append(feature_names[index])
            else:
                top_names.append("")
        adata.obs["top_feature_name"] = top_names


def build_feature_outputs(library_dir: Path, feature_output_dir: Path, provenance: dict) -> dict:
    outputs = {}

    for label, matrix_dir in [("raw", library_dir), ("filtered", library_dir / "filtered")]:
        if not matrix_dir.exists():
            continue
        try:
            adata = load_feature_matrix(matrix_dir)
        except FileNotFoundError:
            continue
        annotate_feature_barcodes(adata, matrix_dir)
        adata.uns["feature_library_provenance"] = provenance
        adata.uns["feature_library_label"] = label
        adata.uns["feature_library_source_dir"] = str(matrix_dir)
        output_path = feature_output_dir / f"{label}_feature_library.h5ad"
        adata.write(output_path)
        outputs[label] = str(output_path)

    return outputs


def integrate_calls(counts_path: Path, calls_path: Path, prefix: str, provenance: dict, generic_aliases: bool):
    counts_adata = ad.read_h5ad(counts_path)
    calls = pd.read_csv(calls_path)

    expected_columns = ["cell_barcode", "num_features", "feature_call", "num_umis"]
    missing = [column for column in expected_columns if column not in calls.columns]
    if missing:
        raise ValueError(f"{calls_path} is missing columns: {', '.join(missing)}")

    calls["cell_barcode"] = calls["cell_barcode"].astype(str)
    calls = calls.set_index("cell_barcode")

    num_features = calls["num_features"].reindex(counts_adata.obs_names).fillna(0).astype(int)
    num_umis = calls["num_umis"].reindex(counts_adata.obs_names).fillna(0).astype(int)
    feature_call = calls["feature_call"].reindex(counts_adata.obs_names).fillna("").astype(str)
    feature_call = feature_call.replace({"None": "", "nan": "", "<NA>": ""})
    feature_call = feature_call.where(num_features > 0, "")

    counts_adata.obs[f"{prefix}__num_features"] = num_features.astype(int)
    counts_adata.obs[f"{prefix}__num_umis"] = num_umis.astype(int)
    counts_adata.obs[f"{prefix}__feature_call"] = feature_call
    counts_adata.obs[f"{prefix}__is_featured"] = (num_features > 0).astype(bool)

    if generic_aliases:
        counts_adata.obs["is_featured"] = counts_adata.obs[f"{prefix}__is_featured"]
        counts_adata.obs["feature_call"] = counts_adata.obs[f"{prefix}__feature_call"]
        counts_adata.obs["feature_call_num_features"] = counts_adata.obs[f"{prefix}__num_features"]
        counts_adata.obs["feature_call_num_umis"] = counts_adata.obs[f"{prefix}__num_umis"]

    feature_libraries = dict(counts_adata.uns.get("feature_libraries", {}))
    feature_libraries[prefix] = {
        "library_id": provenance.get("library_id", ""),
        "sample": provenance.get("sample", ""),
        "feature_type": provenance.get("feature_type", ""),
        "call_source": str(calls_path),
        "obs_columns": [
            f"{prefix}__is_featured",
            f"{prefix}__feature_call",
            f"{prefix}__num_features",
            f"{prefix}__num_umis",
        ],
    }
    counts_adata.uns["feature_libraries"] = feature_libraries
    if generic_aliases:
        counts_adata.uns["feature_library_generic_alias"] = prefix

    counts_adata.write(counts_path)


def main():
    parser = argparse.ArgumentParser(description="Integrate one feature library into downstream h5ad outputs.")
    parser.add_argument("--library-dir", required=True, help="Feature library directory under run/cr_assign")
    parser.add_argument("--feature-output-root", required=True, help="Output root for per-library h5ad artifacts")
    parser.add_argument("--counts-h5ad", action="append", default=[], help="Counts h5ad to annotate; may be repeated")
    parser.add_argument("--calls-csv", help="Per-cell feature calls aligned to counts barcodes")
    parser.add_argument("--set-generic-aliases", action="store_true", help="Also set generic feature-call obs aliases")
    args = parser.parse_args()

    library_dir = Path(args.library_dir).resolve()
    provenance = read_provenance(library_dir)
    library_id = provenance.get("library_id", library_dir.parent.name)
    sample = provenance.get("sample", library_dir.name)
    feature_type = provenance.get("feature_type", library_dir.parent.parent.name if library_dir.parent.parent else "feature_library")

    library_slug = slugify(library_id)
    feature_output_dir = Path(args.feature_output_root).resolve() / library_slug
    feature_output_dir.mkdir(parents=True, exist_ok=True)

    feature_outputs = build_feature_outputs(library_dir, feature_output_dir, provenance)

    prefix = slugify(f"{feature_type}_{library_id}")
    if args.calls_csv:
        calls_path = Path(args.calls_csv).resolve()
        for counts_h5ad in args.counts_h5ad:
            integrate_calls(Path(counts_h5ad).resolve(), calls_path, prefix, provenance, args.set_generic_aliases)

    manifest = {
        "library_id": library_id,
        "sample": sample,
        "feature_type": feature_type,
        "source_dir": str(library_dir),
        "feature_outputs": feature_outputs,
        "calls_csv": str(Path(args.calls_csv).resolve()) if args.calls_csv else "",
        "counts_h5ads": [str(Path(path).resolve()) for path in args.counts_h5ad],
        "obs_prefix": prefix,
        "generic_aliases": bool(args.set_generic_aliases),
    }
    with open(feature_output_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
