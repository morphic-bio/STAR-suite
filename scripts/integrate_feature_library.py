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


def paired_barcode_index(matrix_dir: Path, barcodes: pd.Series) -> pd.Index | None:
    """Map assignment-namespace barcodes to the selected output namespace."""
    barcodes_txt = first_existing(matrix_dir, ["barcodes.txt", "barcodes.csv"])
    barcodes_tsv = first_existing(matrix_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    if barcodes_txt is None or barcodes_tsv is None:
        return None

    assignment = read_table(barcodes_txt).iloc[:, 0].map(canonical_barcode).astype(str)
    output = read_table(barcodes_tsv).iloc[:, 0].map(canonical_barcode).astype(str)
    if len(assignment) != len(output) or len(assignment) == 0:
        return None
    if assignment.duplicated().any() or output.duplicated().any():
        return None

    barcode_map = pd.Series(output.to_numpy(), index=assignment.to_numpy())
    mapped = barcodes.map(canonical_barcode).astype(str).map(barcode_map)
    if mapped.notna().sum() == 0:
        return None
    return pd.Index(mapped.fillna("").to_numpy(), dtype=str, name="barcode")


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
    paired_index = paired_barcode_index(matrix_dir, df["barcode"])
    target_index = pd.Index(adata.obs_names.map(canonical_barcode), dtype=str, name="barcode")

    candidates = [
        ("direct", direct_index, int(direct_index.isin(target_index).sum())),
        ("translated", translated_index, int(translated_index.isin(target_index).sum())),
    ]
    if paired_index is not None:
        candidates.append(
            (
                "paired_barcodes_txt_to_tsv",
                paired_index,
                int(paired_index.isin(target_index).sum()),
            )
        )
    transform, index, _overlap = max(candidates, key=lambda item: item[2])
    df.index = index
    df["barcode_namespace_transform"] = transform

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
            # process_features reserves 0 for ties/no unambiguous top feature;
            # positive ids are one-based MEX feature row ids.
            if 1 <= index <= len(feature_names):
                top_names.append(feature_names[index - 1])
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


def load_calls_csv(path: Path) -> pd.DataFrame:
    calls = pd.read_csv(path)
    barcode_col = "cell_barcode" if "cell_barcode" in calls.columns else "barcode"
    if barcode_col not in calls.columns:
        raise ValueError(f"{path} must contain cell_barcode or barcode")

    calls = calls.copy()
    calls["barcode_raw"] = calls[barcode_col].astype(str)
    calls["barcode_canonical"] = calls["barcode_raw"].map(canonical_barcode)
    calls["barcode_translated"] = calls["barcode_raw"].map(translate_nxt_middle_two_bases)
    return calls


def calls_csv_is_ambient_fdr(calls: pd.DataFrame) -> bool:
    ambient_columns = {"min_qvalue", "default_fdr", "call_status", "caller"}
    if not ambient_columns.intersection(calls.columns):
        return False
    if "caller" not in calls.columns:
        return True
    callers = calls["caller"].dropna().astype(str).str.lower().unique()
    return len(callers) == 0 or any("ambient" in caller and "fdr" in caller for caller in callers)


def map_calls_to_counts(
    counts_path: Path,
    counts_index: pd.Index,
    calls: pd.DataFrame,
) -> tuple[pd.DataFrame, str]:
    direct_index = pd.Index(calls["barcode_canonical"].astype(str), dtype=str, name="barcode")
    translated_index = pd.Index(calls["barcode_translated"].astype(str), dtype=str, name="barcode")
    direct_overlap = int(direct_index.isin(counts_index).sum())
    translated_overlap = int(translated_index.isin(counts_index).sum())
    if translated_overlap > direct_overlap:
        barcode_key = "barcode_translated"
        barcode_transform = "translated"
    else:
        barcode_key = "barcode_canonical"
        barcode_transform = "direct"

    if calls[barcode_key].duplicated().any():
        dupes = calls.loc[calls[barcode_key].duplicated(), barcode_key].astype(str).tolist()[:5]
        raise ValueError(
            f"{counts_path} produced duplicate call barcodes after {barcode_transform} mapping: {dupes}"
        )

    mapped = calls.set_index(barcode_key, drop=False).reindex(counts_index)
    return mapped, barcode_transform


def numeric_call_column(mapped: pd.DataFrame, column: str, default, dtype):
    if column not in mapped.columns:
        return pd.Series(default, index=mapped.index).astype(dtype)
    return pd.to_numeric(mapped[column], errors="coerce").fillna(default).astype(dtype)


def string_call_column(mapped: pd.DataFrame, column: str, default: str) -> pd.Series:
    if column not in mapped.columns:
        return pd.Series(default, index=mapped.index, dtype="string")
    series = mapped[column].astype("string").fillna(default)
    return series.replace({"None": "", "nan": "", "<NA>": ""})


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


def integrate_ambient_fdr_calls(
    counts_path: Path,
    prefix: str,
    provenance: dict,
    generic_aliases: bool,
    calls: pd.DataFrame,
    calls_csv: Path,
):
    counts_adata = ad.read_h5ad(counts_path)
    counts_canonical = counts_adata.obs_names.map(canonical_barcode)
    counts_index = pd.Index(counts_canonical, dtype=str, name="barcode")
    if counts_index.duplicated().any():
        dupes = counts_canonical[counts_canonical.duplicated()].tolist()[:5]
        raise ValueError(f"{counts_path} produced duplicate canonical barcodes: {dupes}")

    mapped, barcode_transform = map_calls_to_counts(counts_path, counts_index, calls)

    num_features = numeric_call_column(mapped, "num_features", 0, "int64")
    num_umis = numeric_call_column(mapped, "num_umis", 0, "int64")
    min_called_umi = numeric_call_column(mapped, "min_called_umi", 0, "int64")
    max_called_umi = numeric_call_column(mapped, "max_called_umi", 0, "int64")
    min_qvalue = numeric_call_column(mapped, "min_qvalue", 1.0, "float64")
    if "num_features_at_default_fdr" in mapped.columns:
        num_features_at_default_fdr = numeric_call_column(mapped, "num_features_at_default_fdr", 0, "int64")
    else:
        num_features_at_default_fdr = num_features
    default_fdr = numeric_call_column(mapped, "default_fdr", 0.01, "float64")
    feature_call = string_call_column(mapped, "feature_call", "")
    call_status = string_call_column(mapped, "call_status", "none")
    caller = string_call_column(mapped, "caller", "ambient-fdr").replace({"": "ambient-fdr"})
    is_called = (num_features > 0).astype(bool)

    prefixed_columns = {
        f"{prefix}__guide_fdr_num_features": num_features.to_numpy(),
        f"{prefix}__guide_fdr_feature_call": pd.Categorical(feature_call.to_numpy()),
        f"{prefix}__guide_fdr_num_umis": num_umis.to_numpy(),
        f"{prefix}__guide_fdr_min_called_umi": min_called_umi.to_numpy(),
        f"{prefix}__guide_fdr_max_called_umi": max_called_umi.to_numpy(),
        f"{prefix}__guide_fdr_min_qvalue": min_qvalue.to_numpy(),
        f"{prefix}__guide_fdr_num_features_at_default_fdr": num_features_at_default_fdr.to_numpy(),
        f"{prefix}__guide_fdr_call_status": pd.Categorical(call_status.to_numpy()),
        f"{prefix}__guide_fdr_default_fdr": default_fdr.to_numpy(),
        f"{prefix}__guide_fdr_caller": pd.Categorical(caller.to_numpy()),
        f"{prefix}__guide_fdr_is_called": is_called.to_numpy(),
    }
    for column, values in prefixed_columns.items():
        counts_adata.obs[column] = values

    generic_columns = {
        "guide_fdr_num_features": num_features.to_numpy(),
        "guide_fdr_feature_call": pd.Categorical(feature_call.to_numpy()),
        "guide_fdr_num_umis": num_umis.to_numpy(),
        "guide_fdr_min_called_umi": min_called_umi.to_numpy(),
        "guide_fdr_max_called_umi": max_called_umi.to_numpy(),
        "guide_fdr_min_qvalue": min_qvalue.to_numpy(),
        "guide_fdr_num_features_at_default_fdr": num_features_at_default_fdr.to_numpy(),
        "guide_fdr_call_status": pd.Categorical(call_status.to_numpy()),
        "guide_fdr_default_fdr": default_fdr.to_numpy(),
        "guide_fdr_caller": pd.Categorical(caller.to_numpy()),
        "guide_fdr_is_called": is_called.to_numpy(),
    }
    if generic_aliases:
        for column, values in generic_columns.items():
            counts_adata.obs[column] = values

    feature_libraries = dict(counts_adata.uns.get("feature_libraries", {}))
    entry = dict(feature_libraries.get(prefix, {}))
    entry.setdefault("library_id", provenance.get("library_id", ""))
    entry.setdefault("sample", provenance.get("sample", ""))
    entry.setdefault("feature_type", provenance.get("feature_type", ""))
    entry["ambient_fdr_call_source"] = str(calls_csv)
    entry["ambient_fdr_barcode_transform"] = barcode_transform
    obs_columns = list(entry.get("obs_columns", []))
    for column in prefixed_columns:
        if column not in obs_columns:
            obs_columns.append(column)
    entry["obs_columns"] = obs_columns
    entry["ambient_fdr_obs_columns"] = list(prefixed_columns.keys())
    feature_libraries[prefix] = entry
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
    parser.add_argument(
        "--ambient-fdr-calls-csv",
        help="Ambient-FDR guide_fdr_calls_per_cell.csv to annotate into h5ad obs",
    )
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

    ambient_fdr_calls_csv = Path(args.ambient_fdr_calls_csv).resolve() if args.ambient_fdr_calls_csv else None
    calls_csv_path = Path(args.calls_csv).resolve() if args.calls_csv else None
    calls_csv_integrated = ""
    if ambient_fdr_calls_csv is None and calls_csv_path is not None:
        candidate_calls = load_calls_csv(calls_csv_path)
        if calls_csv_is_ambient_fdr(candidate_calls):
            ambient_fdr_calls_csv = calls_csv_path
            calls_csv = candidate_calls
        else:
            sibling = calls_csv_path.parent / "ambient_fdr" / "guide_fdr_calls_per_cell.csv"
            if sibling.exists():
                ambient_fdr_calls_csv = sibling
                calls_csv = load_calls_csv(sibling)
            else:
                calls_csv = None
    elif ambient_fdr_calls_csv is not None:
        calls_csv = load_calls_csv(ambient_fdr_calls_csv)
    else:
        calls_csv = None

    if ambient_fdr_calls_csv is not None and calls_csv is not None:
        for counts_h5ad in args.counts_h5ad:
            integrate_ambient_fdr_calls(
                Path(counts_h5ad).resolve(),
                prefix,
                provenance,
                args.set_generic_aliases,
                calls_csv,
                ambient_fdr_calls_csv,
            )
        calls_csv_integrated = str(ambient_fdr_calls_csv)

    manifest = {
        "library_id": library_id,
        "sample": sample,
        "feature_type": feature_type,
        "source_dir": str(library_dir),
        "feature_outputs": feature_outputs,
        "calls_csv": str(Path(args.calls_csv).resolve()) if args.calls_csv else "",
        "ambient_fdr_calls_csv": str(ambient_fdr_calls_csv) if ambient_fdr_calls_csv else "",
        "call_source": feature_obs_source,
        "calls_csv_integrated": calls_csv_integrated,
        "counts_h5ads": [str(Path(path).resolve()) for path in args.counts_h5ad],
        "obs_prefix": prefix,
        "generic_aliases": bool(args.set_generic_aliases),
    }
    with open(feature_output_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
