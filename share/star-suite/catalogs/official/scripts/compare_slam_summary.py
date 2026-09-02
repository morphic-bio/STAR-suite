#!/usr/bin/env python3
"""
Compare two STAR-SLAM gene summaries.

Primary mode (runbook / Phase 2): structured diff of two tabular SlamQuant.out files:
  - Fails if any gene ID appears in only one file (strict identity of gene sets)
  - |ΔNTR| vs --max-abs-delta on **all** shared genes (not read-count filtered)
  - ASCII histogram of |ΔNTR| over all shared genes
  - Pearson correlation of NTR only for shared genes with readcount >= --min-read-count on both sides

  python3 scripts/compare_slam_summary.py ref.SlamQuant.out cand.SlamQuant.out --max-abs-delta 0

Backward-compatible QC JSON mode (keyed comparison, dotted paths):
  python3 scripts/compare_slam_summary.py --json --a run1.slam_qc.json --b run2.slam_qc.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

FLOAT_TOL = 1e-12


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def find_col(header: List[str], names: List[str]) -> Optional[int]:
    lower = {name.lower(): i for i, name in enumerate(header)}
    for name in names:
        idx = lower.get(name.lower())
        if idx is not None:
            return idx
    return None


def find_suffix(header: List[str], suffixes: List[str]) -> Optional[int]:
    for suffix in suffixes:
        for i, name in enumerate(header):
            if name.endswith(" " + suffix) or name.endswith(suffix):
                return i
    return None


def parse_slamquant_out(path: str) -> Dict[str, Dict[str, float]]:
    """Parse STAR SlamQuant.out-like TSV (Gene, ReadCount, Conversions, Coverage, NTR)."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Empty file: {path}")

        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_col(header, ["ReadCount", "Readcount"])
        conv_idx = find_col(header, ["Conversions", "ConversionCount", "MismatchCount"])
        cov_idx = find_col(header, ["Coverage", "TCount"])
        ntr_idx = find_col(header, ["NTR", "NTR_MAP", "MAP"])

        if read_idx is None:
            read_idx = find_suffix(header, ["Readcount", "ReadCount"])
        if conv_idx is None:
            conv_idx = find_suffix(header, ["Conversions", "ConversionCount", "MismatchCount"])
        if cov_idx is None:
            cov_idx = find_suffix(header, ["Coverage", "TCount"])
        if ntr_idx is None:
            ntr_idx = find_suffix(header, ["MAP", "NTR", "NTR_MAP"])

        missing = [
            name
            for name, idx in [
                ("Gene", gene_idx),
                ("ReadCount", read_idx),
                ("Conversions", conv_idx),
                ("Coverage", cov_idx),
                ("NTR", ntr_idx),
            ]
            if idx is None
        ]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}; header={header!r}")

        data: Dict[str, Dict[str, float]] = {}
        for row in reader:
            if not row or len(row) <= max(gene_idx, read_idx, conv_idx, cov_idx, ntr_idx):
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "conversions": float(row[conv_idx]),
                "coverage": float(row[cov_idx]),
                "ntr": float(row[ntr_idx]),
            }
        return data


def pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denx = sum((x - mx) ** 2 for x in xs)
    deny = sum((y - my) ** 2 for y in ys)
    den = math.sqrt(denx * deny)
    if den == 0.0:
        return float("nan")
    return num / den


def ascii_histogram(values: List[float], edges: List[float], title: str) -> None:
    if not values:
        print(f"{title}: (no samples)")
        return
    bins = [0] * (len(edges) - 1)
    last_i = len(bins) - 1
    for v in values:
        placed = False
        for i in range(len(bins)):
            lo, hi = edges[i], edges[i + 1]
            if i == last_i:
                if v >= lo:
                    bins[i] += 1
                    placed = True
                    break
            elif lo <= v < hi:
                bins[i] += 1
                placed = True
                break
        if not placed:
            bins[last_i] += 1
    print(title)
    max_cnt = max(bins) if bins else 1
    scale = 40.0 / max_cnt if max_cnt else 1.0
    for i in range(len(bins)):
        lo, hi = edges[i], edges[i + 1]
        bar = "#" * max(1, int(bins[i] * scale)) if bins[i] else ""
        if i == last_i:
            label = f"  [{lo:g}, +inf)"
        else:
            label = f"  [{lo:g}, {hi:g})"
        print(f"{label}: {bins[i]:5d}  {bar}")


def compare_slamquant(
    ref_path: str,
    cand_path: str,
    max_abs_delta: float,
    min_read_count: float,
    min_ntr_pearson: Optional[float],
    max_report: int,
    allow_gene_set_mismatch: bool,
) -> int:
    ref = parse_slamquant_out(ref_path)
    cand = parse_slamquant_out(cand_path)

    ref_ids = set(ref.keys())
    cand_ids = set(cand.keys())
    only_ref = sorted(ref_ids - cand_ids)
    only_cand = sorted(cand_ids - ref_ids)
    shared = sorted(ref_ids & cand_ids)

    print("=== compare_slam_summary (SlamQuant.out) ===")
    print(f"Reference:  {ref_path}  ({len(ref_ids)} genes)")
    print(f"Candidate:  {cand_path}  ({len(cand_ids)} genes)")
    print(f"Shared:     {len(shared)}")
    print(f"Only ref:   {len(only_ref)}")
    print(f"Only cand: {len(only_cand)}")
    if only_ref[:max_report]:
        print("  sample only-in-ref:", ", ".join(only_ref[:max_report]))
    if only_cand[:max_report]:
        print("  sample only-in-cand:", ", ".join(only_cand[:max_report]))

    ok = True
    if only_ref or only_cand:
        if allow_gene_set_mismatch:
            print(
                "WARNING: gene ID set mismatch (only-in-ref or only-in-cand); "
                "continuing due to --allow-gene-set-mismatch",
                file=sys.stderr,
            )
        else:
            print(
                "FAIL: gene ID set mismatch - comparisons require identical gene IDs in both outputs",
                file=sys.stderr,
            )
            ok = False

    if not shared:
        print("FAIL: no overlapping genes", file=sys.stderr)
        return 1

    ntr_deltas_strict: List[float] = []
    ntr_bad: List[Tuple[str, float, float, float]] = []
    lim_nd = max_abs_delta + FLOAT_TOL if max_abs_delta == 0 else max_abs_delta

    # Strict Phase 2 SE gate: every shared gene counts for |ΔNTR|, regardless of read depth.
    for gid in shared:
        r0 = ref[gid]["ntr"]
        c0 = cand[gid]["ntr"]
        d = abs(r0 - c0)
        ntr_deltas_strict.append(d)
        if d > lim_nd:
            ntr_bad.append((gid, r0, c0, d))

    # Pearson/higher-stability view: optionally restrict to adequately covered genes only.
    ntr_r: List[float] = []
    ntr_c: List[float] = []
    for gid in shared:
        rr = ref[gid]["readcount"]
        cc = cand[gid]["readcount"]
        if rr < min_read_count or cc < min_read_count:
            continue
        ntr_r.append(ref[gid]["ntr"])
        ntr_c.append(cand[gid]["ntr"])

    if ntr_deltas_strict:
        edges = [0.0, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1.0, 10.0, 1e30]
        ascii_histogram(ntr_deltas_strict, edges, "|ΔNTR| histogram (all shared genes, strict gate)")

    corr = pearson(ntr_r, ntr_c)
    print(f"NTR Pearson (shared, readcount >= {min_read_count} both sides): {corr:.6f}")

    if min_ntr_pearson is not None:
        if math.isnan(corr) or corr < min_ntr_pearson:
            print(
                f"FAIL: NTR Pearson {corr:.6f} < --min-ntr-pearson {min_ntr_pearson}",
                file=sys.stderr,
            )
            ok = False

    if ntr_bad:
        ntr_bad.sort(key=lambda t: -t[3])
        print(f"FAIL: |ΔNTR| > max_abs_delta ({max_abs_delta:g}): {len(ntr_bad)} genes", file=sys.stderr)
        for gid, r0, c0, d in ntr_bad[:max_report]:
            print(f"  {gid}\tref_ntr={r0:.8f}\tcand_ntr={c0:.8f}\t|Δ|={d:.8g}", file=sys.stderr)
        ok = False

    if ok:
        print("compare_slam_summary: OK")
    else:
        print("compare_slam_summary: MISMATCH", file=sys.stderr)
    return 0 if ok else 1


# --- QC JSON legacy (optional) ---

def load_json(path: str) -> Dict[str, Any]:
    import json

    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def get_path(obj: Any, dotted: str) -> Any:
    cur: Any = obj
    for part in dotted.split("."):
        if not isinstance(cur, dict):
            return None
        cur = cur.get(part)
    return cur


def compare_vals(a: Any, b: Any, tol: float, key: str) -> Tuple[bool, str]:
    if type(a) != type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        return False, f"{key}: type mismatch {type(a).__name__} vs {type(b).__name__}"
    if isinstance(a, float) or isinstance(b, float):
        af, bf = float(a), float(b)
        if tol <= 0:
            return af == bf, f"{key}: {af} vs {bf}"
        return abs(af - bf) <= tol, f"{key}: {af} vs {bf}"
    if a != b:
        return False, f"{key}: {a!r} vs {b!r}"
    return True, ""


def compare_json_mode(a_path: str, b_path: str, keys_csv: str, tolerance: float) -> int:
    try:
        ja = load_json(a_path)
        jb = load_json(b_path)
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    keys = [k.strip() for k in keys_csv.split(",") if k.strip()]
    failed: List[str] = []
    for key in keys:
        va = get_path(ja, key) if "." in key else ja.get(key)
        vb = get_path(jb, key) if "." in key else jb.get(key)
        if va is None and vb is None:
            continue
        if va is None or vb is None:
            failed.append(f"{key}: missing in one file ({va!r} vs {vb!r})")
            continue
        ok, msg = compare_vals(va, vb, tolerance, key)
        if not ok:
            failed.append(msg)
    if failed:
        print("compare_slam_summary (json): MISMATCH", file=sys.stderr)
        for line in failed:
            print(f"  {line}", file=sys.stderr)
        return 1
    print("compare_slam_summary (json): OK")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    # QC JSON: --a/--b (optional --json for clarity)
    if "--a" in argv and "--b" in argv:
        p = argparse.ArgumentParser(description="Compare SLAM QC JSON keys")
        p.add_argument("--json", action="store_true", help="Compare QC JSON (same as --a/--b without this flag)")
        p.add_argument("--a", required=True)
        p.add_argument("--b", required=True)
        p.add_argument(
            "--keys",
            default=(
                "version,variance_histogram_mode,trim5p,trim3p,"
                "trim5p_mate1,trim3p_mate1,trim5p_mate2,trim3p_mate2,reads_analyzed"
            ),
        )
        p.add_argument("--tolerance", type=float, default=1e-9)
        args = p.parse_args(argv)
        return compare_json_mode(args.a, args.b, args.keys, args.tolerance)

    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("reference", help="Baseline SlamQuant.out (or .gz)")
    p.add_argument("candidate", help="Candidate SlamQuant.out (or .gz)")
    p.add_argument(
        "--max-abs-delta",
        type=float,
        default=1e-9,
        help="Maximum allowed |ΔNTR| per shared gene (use 0 for ~exact match; adds tiny float tol)",
    )
    p.add_argument(
        "--min-read-count",
        type=float,
        default=20.0,
        help="Minimum read count on BOTH sides for NTR Pearson and histogram",
    )
    p.add_argument(
        "--min-ntr-pearson",
        type=float,
        default=None,
        help="If set, fail when Pearson(NTR) on filtered genes is below this threshold",
    )
    p.add_argument("--max-report", type=int, default=20, help="Max per-gene mismatches to print")
    p.add_argument(
        "--allow-gene-set-mismatch",
        action="store_true",
        help="Warn instead of failing when gene ID sets differ between reference and candidate",
    )
    args = p.parse_args(argv)

    if not os.path.isfile(args.reference):
        print(f"ERROR: not a file: {args.reference}", file=sys.stderr)
        return 1
    if not os.path.isfile(args.candidate):
        print(f"ERROR: not a file: {args.candidate}", file=sys.stderr)
        return 1

    return compare_slamquant(
        args.reference,
        args.candidate,
        args.max_abs_delta,
        args.min_read_count,
        args.min_ntr_pearson,
        args.max_report,
        args.allow_gene_set_mismatch,
    )


if __name__ == "__main__":
    raise SystemExit(main())
