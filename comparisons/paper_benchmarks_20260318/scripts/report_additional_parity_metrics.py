#!/usr/bin/env python3
"""
Additional parity metrics for STAR vs Cell Ranger runs:
1) GEX per-barcode Pearson/Spearman correlations
2) GEX per-gene Pearson/Spearman correlations
3) CRISPR feature-call parity (protospacer_calls_per_cell.csv)
"""

import argparse
import csv
import gzip
import importlib.util
import math
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_translation_domains(translation_path: Path) -> Tuple[set, set]:
    """
    Load barcode translation pairs as (left_domain, right_domain).
    Supports whitespace- or comma-separated two-column files.
    """
    left = set()
    right = set()
    with open_maybe_gz(translation_path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2 and "," in line:
                parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            left.add(parts[0])
            right.add(parts[1])
    return left, right


def detect_mixed_pf_libraries(star_run: Path) -> Tuple[bool, List[str]]:
    """
    Best-effort parse of pf_multi_config.csv [libraries] section.
    Returns (is_mixed, sorted_library_types).
    """
    config_path = star_run / "pf_multi_config.csv"
    if not config_path.exists():
        return False, []

    library_types = set()
    in_libraries = False
    header = []
    idx_type = -1

    with open(config_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                in_libraries = (line.lower() == "[libraries]")
                header = []
                idx_type = -1
                continue
            if not in_libraries:
                continue
            if line.startswith("#"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if not header:
                header = parts
                if "library_type" in header:
                    idx_type = header.index("library_type")
                continue
            if idx_type >= 0 and idx_type < len(parts):
                if parts[idx_type]:
                    library_types.add(parts[idx_type])

    types_sorted = sorted(library_types)
    return len(types_sorted) > 1, types_sorted


def import_compare_feature_mex(repo_root: Path):
    module_path = repo_root / "tests" / "compare_feature_mex.py"
    spec = importlib.util.spec_from_file_location("compare_feature_mex", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def safe_int(value: str) -> int:
    v = (value or "").strip()
    if not v:
        return 0
    try:
        return int(v)
    except ValueError:
        try:
            return int(float(v))
        except ValueError:
            return 0


def parse_num_umis(value: str) -> int:
    v = (value or "").strip()
    if not v:
        return 0
    total = 0
    for token in v.split("|"):
        total += safe_int(token)
    return total


def normalize_barcode(raw: str, strip_suffix: bool) -> str:
    bc = (raw or "").strip()
    if not bc:
        return ""
    if "," in bc:
        bc = bc.split(",")[-1].strip()
    if strip_suffix and "-" in bc:
        base, suffix = bc.rsplit("-", 1)
        if suffix.isdigit():
            bc = base
    return bc


def normalize_call(raw: str) -> str:
    value = (raw or "").strip()
    if not value or value.lower() == "none":
        return "None"
    parts = sorted([p.strip() for p in value.split("|") if p.strip()])
    return "|".join(parts) if parts else "None"


_RE_GUIDE_SUFFIX = re.compile(r"^(.*)_([AB])$")
_RE_PROMOTER_SUFFIX = re.compile(r"^(.*)_P\d+(?:P\d+)*$")


def call_to_set(call: str) -> set:
    if call == "None":
        return set()
    return set(x for x in call.split("|") if x)


def feature_to_target(feature: str) -> str:
    m = _RE_GUIDE_SUFFIX.match(feature)
    return m.group(1) if m else feature


def target_to_gene(target: str) -> str:
    m = _RE_PROMOTER_SUFFIX.match(target)
    return m.group(1) if m else target


def collapse_call_set(call_set: set, level: str) -> set:
    if level == "guide":
        return call_set
    if level == "target":
        return {feature_to_target(f) for f in call_set}
    if level == "gene":
        return {target_to_gene(feature_to_target(f)) for f in call_set}
    raise ValueError(f"Unknown collapse level: {level}")


def rankdata(values: List[float]) -> List[float]:
    pairs = sorted((v, i) for i, v in enumerate(values))
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        current = pairs[i][0]
        while j < len(values) and pairs[j][0] == current:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[pairs[k][1]] = avg_rank
        i = j
    return ranks


def pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n == 0:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0.0 or den_y == 0.0:
        return float("nan")
    return num / (den_x * den_y)


def spearman(xs: List[float], ys: List[float]) -> float:
    if len(xs) == 0:
        return float("nan")
    return pearson(rankdata(xs), rankdata(ys))


def fmt_float(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.6f}"


def maybe_translate_counts(
    counts: Dict[str, float],
    cfm,
    translation_path: Path,
    direction: str,
) -> Dict[str, float]:
    if not counts:
        return counts
    mapping = cfm.load_translation_map(translation_path, direction, set(counts.keys()))
    return cfm.apply_translation(counts, mapping)


def maybe_translate_set(values: set, cfm, translation_path: Path, direction: str) -> set:
    if not values:
        return values
    mapping = cfm.load_translation_map(translation_path, direction, set(values))
    out = set()
    for v in values:
        out.add(mapping.get(v, v))
    return out


def fmt_num(value: float) -> str:
    rounded = round(value)
    if abs(value - rounded) < 1e-9:
        return str(int(rounded))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def summarize_corr(label: str, left: Dict[str, float], right: Dict[str, float]) -> List[str]:
    common = sorted(set(left.keys()) & set(right.keys()))
    xs = [left[b] for b in common]
    ys = [right[b] for b in common]
    return [
        f"{label}:",
        f"  common_barcodes: {len(common)}",
        f"  pearson: {fmt_float(pearson(xs, ys))}",
        f"  spearman: {fmt_float(spearman(xs, ys))}",
        f"  left_sum_on_common: {fmt_num(sum(xs))}",
        f"  right_sum_on_common: {fmt_num(sum(ys))}",
    ]


def _resolve_features_file(cfm, mex_dir: Path) -> Path:
    for name in ("features.tsv", "features.txt"):
        try:
            return cfm.resolve_mex_file(mex_dir, name)
        except FileNotFoundError:
            continue
    raise FileNotFoundError(f"Missing features.tsv(.gz) or features.txt(.gz) in {mex_dir}")


def _resolve_barcodes_file(cfm, mex_dir: Path) -> Path:
    for name in ("barcodes.tsv", "barcodes.txt"):
        try:
            return cfm.resolve_mex_file(mex_dir, name)
        except FileNotFoundError:
            continue
    raise FileNotFoundError(f"Missing barcodes.tsv(.gz) or barcodes.txt(.gz) in {mex_dir}")


def build_gene_meta(
    mex_dir: Path,
    feature_types: set,
    strip_suffix: bool,
    cfm,
    matrix_basename: str,
    translation_path: Path | None,
    translation_direction: str,
    translate_enabled: bool,
) -> Dict[str, object]:
    features_path = _resolve_features_file(cfm, mex_dir)
    features, _ = cfm.read_features(features_path)

    selected_features = [
        feat for feat in features if (not feature_types or feat["type"] in feature_types)
    ]
    row_to_id = {feat["row"]: feat["id"] for feat in selected_features}
    feature_ids = {feat["id"] for feat in selected_features}

    barcodes_path = _resolve_barcodes_file(cfm, mex_dir)
    barcodes = cfm.read_barcodes(barcodes_path, strip_suffix)
    if translate_enabled and translation_path is not None:
        mapping = cfm.load_translation_map(
            translation_path, translation_direction, set(barcodes)
        )
        barcodes = [mapping.get(bc, bc) for bc in barcodes]

    matrix_path = cfm.resolve_mex_file(mex_dir, matrix_basename)
    barcode_collisions = len(set(barcodes)) < len(barcodes)

    return {
        "matrix_path": matrix_path,
        "row_to_id": row_to_id,
        "feature_ids": feature_ids,
        "barcodes": barcodes,
        "barcode_collisions": barcode_collisions,
    }


def aggregate_gene_totals_and_cells(
    meta: Dict[str, object],
    allowed_barcodes: set,
    cfm,
) -> Tuple[Dict[str, float], Dict[str, int]]:
    gene_totals: Dict[str, float] = defaultdict(float)
    gene_cells: Dict[str, int] = defaultdict(int)
    if not allowed_barcodes:
        return {}, {}

    matrix_path = meta["matrix_path"]
    row_to_id = meta["row_to_id"]
    barcodes = meta["barcodes"]
    barcode_collisions = bool(meta["barcode_collisions"])

    if barcode_collisions:
        # Rare path: if translation merges barcodes, avoid double-counting cells.
        seen: Dict[str, set] = defaultdict(set)
        for row, col, val in cfm.iter_matrix_entries(matrix_path):
            if val == 0 or row not in row_to_id:
                continue
            if col <= 0 or col - 1 >= len(barcodes):
                continue
            bc = barcodes[col - 1]
            if bc not in allowed_barcodes:
                continue
            gene_id = row_to_id[row]
            gene_totals[gene_id] += val
            if bc not in seen[gene_id]:
                seen[gene_id].add(bc)
                gene_cells[gene_id] += 1
    else:
        for row, col, val in cfm.iter_matrix_entries(matrix_path):
            if val == 0 or row not in row_to_id:
                continue
            if col <= 0 or col - 1 >= len(barcodes):
                continue
            bc = barcodes[col - 1]
            if bc not in allowed_barcodes:
                continue
            gene_id = row_to_id[row]
            gene_totals[gene_id] += val
            gene_cells[gene_id] += 1

    return dict(gene_totals), dict(gene_cells)


def summarize_gene_corr(
    label: str,
    left_meta: Dict[str, object],
    right_meta: Dict[str, object],
    left_counts: Dict[str, float],
    right_counts: Dict[str, float],
    min_counts: int,
    min_cells_pct: float,
    cfm,
) -> List[str]:
    common_barcodes = sorted(set(left_counts.keys()) & set(right_counts.keys()))
    common_barcode_set = set(common_barcodes)
    common_genes = sorted(
        set(left_meta["feature_ids"]) & set(right_meta["feature_ids"])
    )
    min_cells_abs = max(1, int(math.ceil(len(common_barcodes) * min_cells_pct)))

    left_totals, left_cells = aggregate_gene_totals_and_cells(
        left_meta, common_barcode_set, cfm
    )
    right_totals, right_cells = aggregate_gene_totals_and_cells(
        right_meta, common_barcode_set, cfm
    )

    xs_all: List[float] = []
    ys_all: List[float] = []
    xs: List[float] = []
    ys: List[float] = []
    for gene_id in common_genes:
        left_sum = left_totals.get(gene_id, 0)
        right_sum = right_totals.get(gene_id, 0)
        xs_all.append(left_sum)
        ys_all.append(right_sum)
        left_cell_n = left_cells.get(gene_id, 0)
        right_cell_n = right_cells.get(gene_id, 0)
        if (
            left_sum >= min_counts
            and right_sum >= min_counts
            and left_cell_n >= min_cells_abs
            and right_cell_n >= min_cells_abs
        ):
            xs.append(left_sum)
            ys.append(right_sum)

    lines = [
        f"{label}:",
        f"  common_barcodes: {len(common_barcodes)}",
        f"  common_genes: {len(common_genes)}",
        f"  min_counts_per_gene: {min_counts}",
        f"  min_cells_per_gene_pct: {min_cells_pct}",
        f"  min_cells_per_gene_abs: {min_cells_abs}",
        f"  pearson_all_genes: {fmt_float(pearson(xs_all, ys_all)) if len(xs_all) > 1 else 'NA'}",
        f"  spearman_all_genes: {fmt_float(spearman(xs_all, ys_all)) if len(xs_all) > 1 else 'NA'}",
        f"  filtered_genes: {len(xs)}",
    ]
    if len(xs) > 1:
        lines.extend(
            [
                f"  pearson_filtered_genes: {fmt_float(pearson(xs, ys))}",
                f"  spearman_filtered_genes: {fmt_float(spearman(xs, ys))}",
                # Backward-compatible aliases used by existing parsing scripts.
                f"  pearson: {fmt_float(pearson(xs, ys))}",
                f"  spearman: {fmt_float(spearman(xs, ys))}",
                f"  left_sum_on_filtered_genes: {fmt_num(sum(xs))}",
                f"  right_sum_on_filtered_genes: {fmt_num(sum(ys))}",
            ]
        )
    else:
        lines.extend(
            [
                "  pearson_filtered_genes: NA",
                "  spearman_filtered_genes: NA",
                "  pearson: NA",
                "  spearman: NA",
                f"  left_sum_on_filtered_genes: {fmt_num(sum(xs))}",
                f"  right_sum_on_filtered_genes: {fmt_num(sum(ys))}",
            ]
        )
    if len(xs) < 100:
        lines.append("  stability_warning: low_filtered_gene_count")
    return lines


def load_feature_calls(path: Path, strip_suffix: bool) -> Dict[str, Dict[str, int | str]]:
    out: Dict[str, Dict[str, int | str]] = {}
    with open_maybe_gz(path) as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return out
        key_col = "cell_barcode" if "cell_barcode" in reader.fieldnames else reader.fieldnames[0]
        call_col = "feature_call" if "feature_call" in reader.fieldnames else ""
        num_features_col = "num_features" if "num_features" in reader.fieldnames else ""
        num_umis_col = "num_umis" if "num_umis" in reader.fieldnames else ""

        for row in reader:
            bc = normalize_barcode(row.get(key_col, ""), strip_suffix)
            if not bc:
                continue
            record = {
                "call": normalize_call(row.get(call_col, "")) if call_col else "None",
                "num_features": safe_int(row.get(num_features_col, "")) if num_features_col else 0,
                "umis": parse_num_umis(row.get(num_umis_col, "")) if num_umis_col else 0,
            }
            existing = out.get(bc)
            if existing is None:
                out[bc] = record
            else:
                if int(record["umis"]) > int(existing["umis"]):
                    out[bc] = record
    return out


def maybe_translate_calls(
    calls: Dict[str, Dict[str, int | str]],
    cfm,
    translation_path: Path,
    direction: str,
) -> Dict[str, Dict[str, int | str]]:
    if not calls:
        return calls
    mapping = cfm.load_translation_map(translation_path, direction, set(calls.keys()))
    out: Dict[str, Dict[str, int | str]] = {}
    for bc, payload in calls.items():
        new_bc = mapping.get(bc, bc)
        existing = out.get(new_bc)
        if existing is None:
            out[new_bc] = payload
            continue
        if int(payload["umis"]) > int(existing["umis"]):
            out[new_bc] = payload
    return out


def main():
    parser = argparse.ArgumentParser(description="Report additional parity metrics.")
    parser.add_argument("--cr-run", required=True, help="Cell Ranger run directory (contains outs/)")
    parser.add_argument("--star-run", required=True, help="STAR run directory")
    parser.add_argument(
        "--cr-filtered-barcodes",
        default="",
        help="Optional CR filtered barcode file (defaults to outs/filtered_feature_bc_matrix/barcodes.tsv.gz)",
    )
    parser.add_argument(
        "--barcode-translation",
        default="",
        help="Optional two-column translation file (txt/.gz)",
    )
    parser.add_argument(
        "--translation-direction",
        default="left-to-right",
        choices=("left-to-right", "right-to-left"),
    )
    parser.add_argument(
        "--translate",
        default="both",
        choices=("star", "cr", "both", "none"),
    )
    parser.add_argument(
        "--keep-suffix",
        action="store_true",
        help="Keep trailing -<digits> in barcodes (default strips)",
    )
    parser.add_argument(
        "--include-star-none-in-call-parity",
        action="store_true",
        help=(
            "Include STAR rows with feature_call=None in call-parity statistics. "
            "Default excludes STAR None rows from concordance denominators."
        ),
    )
    parser.add_argument(
        "--gene-corr-min-counts",
        type=int,
        default=20,
        help="Minimum total counts per gene in both datasets (default: 20)",
    )
    parser.add_argument(
        "--gene-corr-min-cells-pct",
        type=float,
        default=0.01,
        help="Minimum fraction of common barcodes expressing a gene in both datasets (default: 0.01)",
    )
    parser.add_argument(
        "--gene-corr-feature-type",
        default="Gene Expression",
        help='Feature type for gene-level correlations (default: "Gene Expression"). Use "all" for no filter.',
    )
    parser.add_argument(
        "--cr-raw-matrix-basename",
        default="matrix.mtx",
        help="CR raw matrix basename (default: matrix.mtx)",
    )
    parser.add_argument(
        "--cr-filtered-matrix-basename",
        default="matrix.mtx",
        help="CR filtered matrix basename (default: matrix.mtx)",
    )
    parser.add_argument(
        "--star-raw-matrix-basename",
        default="matrix.mtx",
        help="STAR raw matrix basename (default: matrix.mtx)",
    )
    parser.add_argument(
        "--star-filtered-matrix-basename",
        default="matrix.mtx",
        help="STAR filtered matrix basename (default: matrix.mtx)",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    cfm = import_compare_feature_mex(repo_root)

    strip_suffix = not args.keep_suffix
    cr_run = Path(args.cr_run)
    star_run = Path(args.star_run)

    cr_raw_mex = cr_run / "outs" / "raw_feature_bc_matrix"
    cr_filt_mex = cr_run / "outs" / "filtered_feature_bc_matrix"
    star_raw_gex = star_run / "Solo.out" / "GeneFull" / "raw"
    star_filt_gex = star_run / "Solo.out" / "GeneFull" / "filtered"
    cr_calls_csv = cr_run / "outs" / "crispr_analysis" / "protospacer_calls_per_cell.csv"
    star_calls_csv = star_run / "outs" / "crispr_analysis" / "protospacer_calls_per_cell.csv"

    if not cr_raw_mex.exists() or not cr_filt_mex.exists():
        raise SystemExit(f"Missing CR MEX dirs under: {cr_run}")
    if not star_raw_gex.exists() or not star_filt_gex.exists():
        raise SystemExit(f"Missing STAR GeneFull MEX dirs under: {star_run}")
    if not cr_calls_csv.exists() or not star_calls_csv.exists():
        raise SystemExit("Missing protospacer_calls_per_cell.csv in CR or STAR run")

    cr_filtered_barcodes = (
        Path(args.cr_filtered_barcodes)
        if args.cr_filtered_barcodes
        else cr_filt_mex / "barcodes.tsv.gz"
    )
    if not cr_filtered_barcodes.exists():
        raise SystemExit(f"Missing CR filtered barcode file: {cr_filtered_barcodes}")

    feature_types = {"Gene Expression"}
    _, cr_raw_counts, _, _ = cfm.load_mex_counts(
        cr_raw_mex, feature_types, strip_suffix, args.cr_raw_matrix_basename
    )
    _, star_raw_counts, _, _ = cfm.load_mex_counts(
        star_raw_gex, feature_types, strip_suffix, args.star_raw_matrix_basename
    )
    _, cr_filt_counts, _, _ = cfm.load_mex_counts(
        cr_filt_mex, feature_types, strip_suffix, args.cr_filtered_matrix_basename
    )
    _, star_filt_counts, _, _ = cfm.load_mex_counts(
        star_filt_gex, feature_types, strip_suffix, args.star_filtered_matrix_basename
    )

    filter_set = cfm.load_barcode_filter(cr_filtered_barcodes, strip_suffix)

    translation_path = Path(args.barcode_translation) if args.barcode_translation else None
    if translation_path is not None and not translation_path.exists():
        raise SystemExit(f"Missing translation file: {translation_path}")

    translate_enabled = translation_path is not None and args.translate != "none"
    source_domain = set()
    target_domain = set()
    cr_raw_keys_pre = set(cr_raw_counts.keys())
    star_raw_keys_pre = set(star_raw_counts.keys())
    if translate_enabled:
        left_domain, right_domain = load_translation_domains(translation_path)
        if args.translation_direction == "left-to-right":
            source_domain, target_domain = left_domain, right_domain
        else:
            source_domain, target_domain = right_domain, left_domain

    if translation_path is not None and args.translate != "none":
        if args.translate in ("cr", "both"):
            cr_raw_counts = maybe_translate_counts(
                cr_raw_counts, cfm, translation_path, args.translation_direction
            )
            cr_filt_counts = maybe_translate_counts(
                cr_filt_counts, cfm, translation_path, args.translation_direction
            )
        if args.translate in ("star", "both"):
            star_raw_counts = maybe_translate_counts(
                star_raw_counts, cfm, translation_path, args.translation_direction
            )
            star_filt_counts = maybe_translate_counts(
                star_filt_counts, cfm, translation_path, args.translation_direction
            )
        filter_set = maybe_translate_set(filter_set, cfm, translation_path, args.translation_direction)

    if args.gene_corr_min_counts < 0:
        raise SystemExit("--gene-corr-min-counts must be >= 0")
    if args.gene_corr_min_cells_pct < 0.0 or args.gene_corr_min_cells_pct > 1.0:
        raise SystemExit("--gene-corr-min-cells-pct must be within [0, 1]")

    cr_raw_on_filter = {bc: c for bc, c in cr_raw_counts.items() if bc in filter_set}
    star_raw_on_filter = {bc: c for bc, c in star_raw_counts.items() if bc in filter_set}

    print("=" * 70)
    print("Parity Metrics Configuration")
    print("=" * 70)
    print(f"cr_run: {cr_run}")
    print(f"star_run: {star_run}")
    print(f"cr_filtered_barcodes: {cr_filtered_barcodes}")
    print(f"keep_suffix: {args.keep_suffix}")
    print(f"translate_enabled: {translate_enabled}")
    if translation_path is not None:
        print(f"translation_file: {translation_path}")
        print(f"translation_direction: {args.translation_direction}")
        print(f"translate_side: {args.translate}")
    else:
        print("translation_file: none")
        print("translation_direction: none")
        print("translate_side: none")
    if translate_enabled:
        cr_in_source = len(cr_raw_keys_pre & source_domain)
        cr_in_target = len(cr_raw_keys_pre & target_domain)
        star_in_source = len(star_raw_keys_pre & source_domain)
        star_in_target = len(star_raw_keys_pre & target_domain)
        cr_in_both = len(cr_raw_keys_pre & source_domain & target_domain)
        star_in_both = len(star_raw_keys_pre & source_domain & target_domain)
        print(f"translation_observed_raw_cr_in_source_domain: {cr_in_source}")
        print(f"translation_observed_raw_cr_in_target_domain: {cr_in_target}")
        print(f"translation_observed_raw_star_in_source_domain: {star_in_source}")
        print(f"translation_observed_raw_star_in_target_domain: {star_in_target}")
        print(f"translation_observed_raw_cr_in_both_domains: {cr_in_both}")
        print(f"translation_observed_raw_star_in_both_domains: {star_in_both}")
        if args.translate in ("cr", "star") and (cr_in_both > 0 or star_in_both > 0):
            print(
                "translation_warning: source/target domains overlap observed barcodes; "
                "single-side translation can shift overlap/correlation metrics. "
                "Prefer --translate both (or none) for stable comparisons."
            )
    print()

    print("=" * 70)
    print("GEX Correlations (Per-Barcode UMI Totals)")
    print("=" * 70)
    print(f"matrix_raw_cr: {args.cr_raw_matrix_basename}")
    print(f"matrix_raw_star: {args.star_raw_matrix_basename}")
    print(f"matrix_filtered_cr: {args.cr_filtered_matrix_basename}")
    print(f"matrix_filtered_star: {args.star_filtered_matrix_basename}")
    print()
    for line in summarize_corr("raw_all", cr_raw_counts, star_raw_counts):
        print(line)
    for line in summarize_corr("raw_restricted_to_cr_filtered", cr_raw_on_filter, star_raw_on_filter):
        print(line)
    for line in summarize_corr("filtered_vs_filtered", cr_filt_counts, star_filt_counts):
        print(line)
    print()

    translate_cr = translation_path is not None and args.translate in ("cr", "both")
    translate_star = translation_path is not None and args.translate in ("star", "both")
    gene_feature_types = set()
    if args.gene_corr_feature_type.lower() not in ("all", "*", "any", "none"):
        gene_feature_types = {args.gene_corr_feature_type}

    cr_raw_gene_meta = build_gene_meta(
        cr_raw_mex,
        gene_feature_types,
        strip_suffix,
        cfm,
        args.cr_raw_matrix_basename,
        translation_path,
        args.translation_direction,
        translate_cr,
    )
    star_raw_gene_meta = build_gene_meta(
        star_raw_gex,
        gene_feature_types,
        strip_suffix,
        cfm,
        args.star_raw_matrix_basename,
        translation_path,
        args.translation_direction,
        translate_star,
    )
    cr_filt_gene_meta = build_gene_meta(
        cr_filt_mex,
        gene_feature_types,
        strip_suffix,
        cfm,
        args.cr_filtered_matrix_basename,
        translation_path,
        args.translation_direction,
        translate_cr,
    )
    star_filt_gene_meta = build_gene_meta(
        star_filt_gex,
        gene_feature_types,
        strip_suffix,
        cfm,
        args.star_filtered_matrix_basename,
        translation_path,
        args.translation_direction,
        translate_star,
    )

    print("=" * 70)
    print("GEX Correlations (Per-Gene Totals Across Common Barcodes)")
    print("=" * 70)
    print(f"method_feature_type: {args.gene_corr_feature_type}")
    print("method_barcode_set: normalized/translated common barcodes per comparison")
    print(f"method_min_counts_per_gene: {args.gene_corr_min_counts}")
    print(f"method_min_cells_per_gene_pct: {args.gene_corr_min_cells_pct}")
    mixed_pf, pf_library_types = detect_mixed_pf_libraries(star_run)
    if mixed_pf:
        print(
            "method_warning: STAR run has mixed PF libraries "
            f"({', '.join(pf_library_types)}). "
            "GeneFull may include non-GEX-library reads in this run."
        )
        print(
            "method_recommendation: use a GEX-only STAR run for strict gene-level "
            "parity against Cell Ranger GEX."
        )
    print()
    for line in summarize_gene_corr(
        "raw_all",
        cr_raw_gene_meta,
        star_raw_gene_meta,
        cr_raw_counts,
        star_raw_counts,
        args.gene_corr_min_counts,
        args.gene_corr_min_cells_pct,
        cfm,
    ):
        print(line)
    for line in summarize_gene_corr(
        "raw_restricted_to_cr_filtered",
        cr_raw_gene_meta,
        star_raw_gene_meta,
        cr_raw_on_filter,
        star_raw_on_filter,
        args.gene_corr_min_counts,
        args.gene_corr_min_cells_pct,
        cfm,
    ):
        print(line)
    for line in summarize_gene_corr(
        "filtered_vs_filtered",
        cr_filt_gene_meta,
        star_filt_gene_meta,
        cr_filt_counts,
        star_filt_counts,
        args.gene_corr_min_counts,
        args.gene_corr_min_cells_pct,
        cfm,
    ):
        print(line)
    print()

    cr_calls = load_feature_calls(cr_calls_csv, strip_suffix)
    star_calls = load_feature_calls(star_calls_csv, strip_suffix)

    if translation_path is not None and args.translate != "none":
        if args.translate in ("cr", "both"):
            cr_calls = maybe_translate_calls(
                cr_calls, cfm, translation_path, args.translation_direction
            )
        if args.translate in ("star", "both"):
            star_calls = maybe_translate_calls(
                star_calls, cfm, translation_path, args.translation_direction
            )

    common_calls_all = sorted(set(cr_calls.keys()) & set(star_calls.keys()))
    star_none_rows = sum(1 for payload in star_calls.values() if str(payload["call"]) == "None")
    star_non_none_rows = len(star_calls) - star_none_rows

    if args.include_star_none_in_call_parity:
        common_calls_eval = common_calls_all
    else:
        common_calls_eval = [
            bc for bc in common_calls_all if str(star_calls[bc]["call"]) != "None"
        ]

    excluded_common_star_none = len(common_calls_all) - len(common_calls_eval)

    set_equiv = sum(
        1 for bc in common_calls_eval if str(cr_calls[bc]["call"]) == str(star_calls[bc]["call"])
    )
    set_equiv_target = 0
    set_equiv_gene = 0
    for bc in common_calls_eval:
        cr_call_set = call_to_set(str(cr_calls[bc]["call"]))
        star_call_set = call_to_set(str(star_calls[bc]["call"]))
        if collapse_call_set(cr_call_set, "target") == collapse_call_set(star_call_set, "target"):
            set_equiv_target += 1
        if collapse_call_set(cr_call_set, "gene") == collapse_call_set(star_call_set, "gene"):
            set_equiv_gene += 1
    num_features_match = sum(
        1
        for bc in common_calls_eval
        if int(cr_calls[bc]["num_features"]) == int(star_calls[bc]["num_features"])
    )
    cr_umis_common = [int(cr_calls[bc]["umis"]) for bc in common_calls_eval]
    star_umis_common = [int(star_calls[bc]["umis"]) for bc in common_calls_eval]

    print("=" * 70)
    print("Feature-Call Parity (protospacer_calls_per_cell.csv)")
    print("=" * 70)
    if len(star_calls) == 0:
        print(
            "warning: STAR protospacer_calls_per_cell.csv has no data rows; "
            "call-parity statistics are not informative."
        )
    if len(cr_calls) == 0:
        print(
            "warning: CR protospacer_calls_per_cell.csv has no data rows; "
            "call-parity statistics are not informative."
        )
    print(f"rows_cr: {len(cr_calls)}")
    print(f"rows_star: {len(star_calls)}")
    print(f"rows_star_non_none: {star_non_none_rows}")
    print(f"rows_star_none: {star_none_rows}")
    print(f"common_rows_all: {len(common_calls_all)}")
    print(f"common_rows_eval: {len(common_calls_eval)}")
    print(f"excluded_common_rows_star_none: {excluded_common_star_none}")
    print(
        "call_parity_mode: "
        + (
            "include_star_none"
            if args.include_star_none_in_call_parity
            else "exclude_star_none_default"
        )
    )
    print(f"only_cr_rows: {len(set(cr_calls.keys()) - set(star_calls.keys()))}")
    print(f"only_star_rows: {len(set(star_calls.keys()) - set(cr_calls.keys()))}")
    print(f"num_features_exact_match: {num_features_match}")
    print(
        "num_features_exact_match_pct: "
        + (f"{(100.0 * num_features_match / len(common_calls_eval)):.4f}" if common_calls_eval else "NA")
    )
    print(f"set_equivalent_calls: {set_equiv}")
    print(
        "set_equivalent_calls_pct: "
        + (f"{(100.0 * set_equiv / len(common_calls_eval)):.4f}" if common_calls_eval else "NA")
    )
    print(f"set_equivalent_calls_target: {set_equiv_target}")
    print(
        "set_equivalent_calls_target_pct: "
        + (f"{(100.0 * set_equiv_target / len(common_calls_eval)):.4f}" if common_calls_eval else "NA")
    )
    print(f"set_equivalent_calls_gene: {set_equiv_gene}")
    print(
        "set_equivalent_calls_gene_pct: "
        + (f"{(100.0 * set_equiv_gene / len(common_calls_eval)):.4f}" if common_calls_eval else "NA")
    )
    print(f"umi_sum_pearson_on_common: {fmt_float(pearson(cr_umis_common, star_umis_common))}")
    print(f"umi_sum_spearman_on_common: {fmt_float(spearman(cr_umis_common, star_umis_common))}")
    print(f"total_assigned_umis_cr: {sum(int(v['umis']) for v in cr_calls.values())}")
    print(f"total_assigned_umis_star: {sum(int(v['umis']) for v in star_calls.values())}")
    print(
        "total_assigned_umis_delta_star_minus_cr: "
        + str(sum(int(v["umis"]) for v in star_calls.values()) - sum(int(v["umis"]) for v in cr_calls.values()))
    )


if __name__ == "__main__":
    main()
