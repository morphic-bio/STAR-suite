#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return slug or "feature_library"


def canonical_barcode(value: str) -> str:
    return re.sub(r"-[0-9]+$", "", str(value).strip())


def complement_base(base: str) -> str:
    return {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }.get(base.upper(), base)


def translate_nxt_middle_two_bases(barcode: str) -> str:
    barcode = canonical_barcode(barcode).upper()
    if len(barcode) >= 9:
        chars = list(barcode)
        chars[7] = complement_base(chars[7])
        chars[8] = complement_base(chars[8])
        return "".join(chars)
    return barcode


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
    # Prefer barcodes.tsv over barcodes.txt. In current pf outputs, barcodes.tsv
    # is the output-namespace surface that should align with GEX barcodes,
    # while barcodes.txt can remain in the assignment namespace.
    barcodes_path = first_existing(matrix_dir, ["barcodes.tsv", "barcodes.tsv.gz", "barcodes.txt", "barcodes.csv"])
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
    direct_index = pd.Index(df["barcode"].map(canonical_barcode), dtype=str, name="barcode")
    translated_index = pd.Index(df["barcode"].map(translate_nxt_middle_two_bases), dtype=str, name="barcode")
    target_index = pd.Index(adata.obs_names.map(canonical_barcode), dtype=str, name="barcode")

    direct_overlap = int(direct_index.isin(target_index).sum())
    translated_overlap = int(translated_index.isin(target_index).sum())
    if translated_overlap > direct_overlap:
        df.index = translated_index
        df["barcode_namespace_transform"] = "translated"
    else:
        df.index = direct_index
        df["barcode_namespace_transform"] = "direct"

    adata.obs["barcode_feature_namespace"] = df["barcode_namespace_transform"].iloc[0]

    for source_col, target_col in [
        ("num_features", "num_features"),
        ("top_feature_index", "top_feature_index"),
        ("total_deduped_umi", "total_deduped_umi"),
    ]:
        if source_col in df.columns:
            series = df[source_col].reindex(target_index)
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


def derive_feature_calls(adata: ad.AnnData) -> pd.DataFrame:
    if adata.n_obs == 0:
        return pd.DataFrame(
            columns=[
                "best_feature",
                "feature1_count",
                "feature2_count",
                "feature_call_category",
            ]
        )

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = sp.csr_matrix(X)

    best_feature = []
    feature1_count = np.zeros(adata.n_obs, dtype=np.int64)
    feature2_count = np.zeros(adata.n_obs, dtype=np.int64)

    for i in range(adata.n_obs):
        row = X.getrow(i)
        if row.nnz == 0:
            best_feature.append("")
            continue

        counts = row.data.astype(np.int64, copy=False)
        cols = row.indices
        order = np.argsort(counts)[::-1]
        top_idx = order[0]
        top_count = int(counts[top_idx])
        second_count = int(counts[order[1]]) if len(order) > 1 else 0

        feature1_count[i] = top_count
        feature2_count[i] = second_count
        if top_count > second_count:
            best_feature.append(str(adata.var_names[cols[top_idx]]))
        else:
            best_feature.append("")

    num_features = adata.obs.get("num_features", pd.Series(0, index=adata.obs_names)).astype(int)
    category = np.where(
        num_features.to_numpy() <= 0,
        "none",
        np.where(
            np.array(best_feature, dtype=object) == "",
            "ambiguous",
            np.where(num_features.to_numpy() > 1, "multi", "single"),
        ),
    )

    return pd.DataFrame(
        {
            "best_feature": pd.Series(best_feature, index=adata.obs_names, dtype="string"),
            "feature1_count": pd.Series(feature1_count, index=adata.obs_names, dtype="int64"),
            "feature2_count": pd.Series(feature2_count, index=adata.obs_names, dtype="int64"),
            "feature_call_category": pd.Series(category, index=adata.obs_names, dtype="string"),
        }
    )


def build_feature_obs_table(adata: ad.AnnData) -> pd.DataFrame:
    metrics = derive_feature_calls(adata)
    obs = adata.obs.copy()
    for col in metrics.columns:
        obs[col] = metrics[col]

    obs["barcode_raw"] = obs.index.astype(str)
    obs["barcode_canonical"] = obs.index.map(canonical_barcode)
    obs["barcode_translated"] = obs["barcode_raw"].map(translate_nxt_middle_two_bases)
    return obs


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


def integrate_calls(
    counts_path: Path,
    prefix: str,
    provenance: dict,
    generic_aliases: bool,
    feature_obs: pd.DataFrame,
    call_source: str,
):
    counts_adata = ad.read_h5ad(counts_path)
    counts_canonical = counts_adata.obs_names.map(canonical_barcode)
    counts_index = pd.Index(counts_canonical, dtype=str, name="barcode")
    if counts_index.duplicated().any():
        dupes = counts_canonical[counts_canonical.duplicated()].tolist()[:5]
        raise ValueError(f"{counts_path} produced duplicate canonical barcodes: {dupes}")

    direct_index = pd.Index(feature_obs["barcode_canonical"].astype(str), dtype=str, name="barcode")
    translated_index = pd.Index(feature_obs["barcode_translated"].astype(str), dtype=str, name="barcode")
    direct_overlap = int(direct_index.isin(counts_index).sum())
    translated_overlap = int(translated_index.isin(counts_index).sum())
    if translated_overlap > direct_overlap:
        barcode_key = "barcode_translated"
        barcode_transform = "translated"
    else:
        barcode_key = "barcode_canonical"
        barcode_transform = "direct"

    if feature_obs[barcode_key].duplicated().any():
        dupes = feature_obs.loc[feature_obs[barcode_key].duplicated(), barcode_key].astype(str).tolist()[:5]
        raise ValueError(f"{counts_path} produced duplicate feature barcodes after {barcode_transform} mapping: {dupes}")

    mapped = feature_obs.set_index(barcode_key, drop=False).reindex(counts_index)

    num_features = mapped["num_features"].fillna(0).astype(int)
    num_umis = mapped["total_deduped_umi"].fillna(0).astype(int)
    feature_call = mapped["best_feature"].fillna("").astype(str)
    feature1_count = mapped["feature1_count"].fillna(0).astype(int)
    feature2_count = mapped["feature2_count"].fillna(0).astype(int)
    feature_category = mapped["feature_call_category"].fillna("none").astype(str)
    is_featured = mapped["is_featured"].fillna(False).astype(bool)

    counts_adata.obs[f"{prefix}__num_features"] = num_features.to_numpy()
    counts_adata.obs[f"{prefix}__num_umis"] = num_umis.to_numpy()
    counts_adata.obs[f"{prefix}__feature_call"] = feature_call.to_numpy()
    counts_adata.obs[f"{prefix}__is_featured"] = is_featured.to_numpy()
    counts_adata.obs[f"{prefix}__feature1_count"] = feature1_count.to_numpy()
    counts_adata.obs[f"{prefix}__feature2_count"] = feature2_count.to_numpy()
    counts_adata.obs[f"{prefix}__feature_call_category"] = pd.Categorical(feature_category.to_numpy())

    if generic_aliases:
        counts_adata.obs["is_featured"] = counts_adata.obs[f"{prefix}__is_featured"]
        counts_adata.obs["feature_call"] = counts_adata.obs[f"{prefix}__feature_call"]
        counts_adata.obs["feature_call_num_features"] = counts_adata.obs[f"{prefix}__num_features"]
        counts_adata.obs["feature_call_num_umis"] = counts_adata.obs[f"{prefix}__num_umis"]
        counts_adata.obs["best_feature"] = counts_adata.obs[f"{prefix}__feature_call"]
        counts_adata.obs["feature1_count"] = counts_adata.obs[f"{prefix}__feature1_count"]
        counts_adata.obs["feature2_count"] = counts_adata.obs[f"{prefix}__feature2_count"]
        counts_adata.obs["feature_call_category"] = counts_adata.obs[f"{prefix}__feature_call_category"]

    feature_libraries = dict(counts_adata.uns.get("feature_libraries", {}))
    feature_libraries[prefix] = {
        "library_id": provenance.get("library_id", ""),
        "sample": provenance.get("sample", ""),
        "feature_type": provenance.get("feature_type", ""),
        "call_source": call_source,
        "barcode_transform": barcode_transform,
        "obs_columns": [
            f"{prefix}__is_featured",
            f"{prefix}__feature_call",
            f"{prefix}__num_features",
            f"{prefix}__num_umis",
            f"{prefix}__feature1_count",
            f"{prefix}__feature2_count",
            f"{prefix}__feature_call_category",
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
    raw_feature_dir = feature_output_dir / "raw_feature_library.h5ad"
    feature_obs = None
    feature_obs_source = ""
    if raw_feature_dir.exists():
        raw_feature_adata = ad.read_h5ad(raw_feature_dir)
        feature_obs = build_feature_obs_table(raw_feature_adata)
        feature_obs_source = str(raw_feature_dir)

    prefix = slugify(f"{feature_type}_{library_id}")
    if feature_obs is not None:
        for counts_h5ad in args.counts_h5ad:
            integrate_calls(
                Path(counts_h5ad).resolve(),
                prefix,
                provenance,
                args.set_generic_aliases,
                feature_obs,
                feature_obs_source,
            )

    manifest = {
        "library_id": library_id,
        "sample": sample,
        "feature_type": feature_type,
        "source_dir": str(library_dir),
        "feature_outputs": feature_outputs,
        "calls_csv": str(Path(args.calls_csv).resolve()) if args.calls_csv else "",
        "call_source": feature_obs_source,
        "counts_h5ads": [str(Path(path).resolve()) for path in args.counts_h5ad],
        "obs_prefix": prefix,
        "generic_aliases": bool(args.set_generic_aliases),
    }
    with open(feature_output_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
