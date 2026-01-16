#!/usr/bin/env python3
"""
Compare STAR-Slam output against GRAND-SLAM reference fixture.

Expected STAR-Slam columns (tab-separated, header required):
  Gene (or GeneID), ReadCount, Conversions, Coverage, NTR (or MAP)
"""

import argparse
import csv
import gzip
import math
import sys
try:
    from scipy.stats import spearmanr
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def find_col(header, names):
    lower = {name.lower(): i for i, name in enumerate(header)}
    for name in names:
        idx = lower.get(name.lower())
        if idx is not None:
            return idx
    return None


def find_suffix(header, suffixes):
    for suffix in suffixes:
        for i, name in enumerate(header):
            if name.endswith(" " + suffix) or name.endswith(suffix):
                return i
    return None


def parse_reference(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Reference file is empty")

        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_suffix(header, ["Readcount"])
        conv_idx = find_suffix(header, ["Conversions"])
        cov_idx = find_suffix(header, ["Coverage"])
        map_idx = find_suffix(header, ["MAP"])

        missing = [name for name, idx in [
            ("Gene", gene_idx),
            ("Readcount", read_idx),
            ("MAP", map_idx),
        ] if idx is None]
        if missing:
            raise ValueError("Reference missing columns: " + ", ".join(missing))

        data = {}
        has_conv_cov = conv_idx is not None and cov_idx is not None
        for row in reader:
            if not row or len(row) <= map_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            conversions = None
            coverage = None
            if has_conv_cov:
                conversions = float(row[conv_idx])
                coverage = float(row[cov_idx])
            data[gene] = {
                "readcount": float(row[read_idx]),
                "conversions": conversions,
                "coverage": coverage,
                "ntr": float(row[map_idx]),
            }
        return data


def parse_test(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Test file is empty")

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

        missing = [name for name, idx in [
            ("Gene", gene_idx),
            ("ReadCount", read_idx),
            ("Conversions", conv_idx),
            ("Coverage", cov_idx),
            ("NTR", ntr_idx),
        ] if idx is None]
        if missing:
            raise ValueError("Test output missing columns: " + ", ".join(missing))

        data = {}
        for row in reader:
            if not row or len(row) <= ntr_idx:
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


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = sum((x - mean_x) ** 2 for x in xs)
    den_y = sum((y - mean_y) ** 2 for y in ys)
    den = math.sqrt(den_x * den_y)
    if den == 0.0:
        return float("nan")
    return num / den


def spearman(xs, ys):
    """Compute Spearman rank correlation."""
    if len(xs) < 2 or len(ys) < 2 or len(xs) != len(ys):
        return float("nan")
    if not HAS_SCIPY:
        # Manual implementation if scipy not available
        n = len(xs)
        rank_x = sorted(range(n), key=lambda i: xs[i])
        rank_y = sorted(range(n), key=lambda i: ys[i])
        rank_x_dict = {rank_x[i]: i + 1 for i in range(n)}
        rank_y_dict = {rank_y[i]: i + 1 for i in range(n)}
        d_sq = sum((rank_x_dict[i] - rank_y_dict[i]) ** 2 for i in range(n))
        return 1.0 - (6.0 * d_sq) / (n * (n * n - 1))
    try:
        corr, _ = spearmanr(xs, ys)
        return corr if not math.isnan(corr) else float("nan")
    except:
        return float("nan")


def main():
    parser = argparse.ArgumentParser(description="Compare STAR-Slam output to GRAND-SLAM fixture.")
    parser.add_argument("--reference", required=True, help="GRAND-SLAM fixture TSV(.gz)")
    parser.add_argument("--test", required=True, help="STAR-Slam output TSV")
    parser.add_argument("--min-read-count", type=float, default=None,
                        help="Minimum read count threshold (single threshold mode, overrides --thresholds)")
    parser.add_argument("--min-coverage", type=float, default=None,
                        help="Minimum coverage (nT) threshold (single threshold mode, overrides --thresholds)")
    parser.add_argument("--thresholds", type=str, default="20,50,100",
                        help="Comma-separated thresholds for correlation table (default: 20,50,100; ignored if --min-read-count or --min-coverage set)")
    parser.add_argument("--count-tol", type=float, default=0.0,
                        help="Tolerance for read/conversion/coverage deltas")
    parser.add_argument("--ntr-abs-max", type=float, default=1e-3,
                        help="Max absolute NTR difference (MAP vs NTR)")
    parser.add_argument("--corr-min", type=float, default=0.999,
                        help="Minimum Pearson correlation for NTR")
    parser.add_argument("--max-report", type=int, default=10,
                        help="Max mismatches to print")
    args = parser.parse_args()

    ref = parse_reference(args.reference)
    test = parse_test(args.test)

    ref_missing_conv_cov = any(
        ref_data.get("conversions") is None or ref_data.get("coverage") is None
        for ref_data in ref.values()
    )
    if ref_missing_conv_cov:
        print("WARNING: Reference missing Conversions/Coverage; conversion metrics will be skipped.")

    shared = sorted(set(ref.keys()) & set(test.keys()))
    if not shared:
        print("FAIL: no overlapping genes between reference and test")
        return 1

    def collect_deltas(field):
        deltas = []
        for gene in shared:
            ref_val = ref[gene].get(field)
            test_val = test[gene].get(field)
            if ref_val is None or test_val is None:
                continue
            delta = abs(ref_val - test_val)
            if delta > args.count_tol:
                deltas.append((delta, gene, ref_val, test_val))
        deltas.sort(reverse=True)
        return deltas

    read_deltas = collect_deltas("readcount")
    conv_deltas = collect_deltas("conversions")
    cov_deltas = collect_deltas("coverage")

    # Determine thresholds to use
    if args.min_read_count is not None or args.min_coverage is not None:
        # Single threshold mode (backward compatibility)
        thresholds = []
        if args.min_read_count is not None:
            thresholds = [(args.min_read_count, "readcount")]
        elif args.min_coverage is not None:
            thresholds = [(args.min_coverage, "coverage")]
    else:
        # Multi-threshold mode
        threshold_values = [float(t.strip()) for t in args.thresholds.split(",")]
        thresholds = [(t, "readcount") for t in threshold_values]

    # Collect pairs for correlation analysis
    def get_filtered_pairs(threshold, filter_field):
        pairs = []
        for gene in shared:
            if ref[gene][filter_field] < threshold:
                continue
            pairs.append((gene, ref[gene], test[gene]))
        return pairs

    # Compute correlations for each threshold
    correlation_results = []
    for threshold, filter_field in thresholds:
        pairs = get_filtered_pairs(threshold, filter_field)
        if len(pairs) < 2:
            continue
        
        # NTR pairs
        ntr_ref = [ref_data["ntr"] for _, ref_data, _ in pairs]
        ntr_test = [test_data["ntr"] for _, _, test_data in pairs]
        ntr_pearson = pearson(ntr_ref, ntr_test)
        ntr_spearman = spearman(ntr_ref, ntr_test)
        
        # Conversion fraction (k/nT) pairs
        conv_frac_pearson = float("nan")
        conv_frac_spearman = float("nan")
        conv_frac_ref = []
        conv_frac_test = []
        for _, ref_data, test_data in pairs:
            if ref_data["coverage"] is None or ref_data["conversions"] is None:
                continue
            if test_data["coverage"] is None or test_data["conversions"] is None:
                continue
            if ref_data["coverage"] > 0:
                conv_frac_ref.append(ref_data["conversions"] / ref_data["coverage"])
            else:
                conv_frac_ref.append(0.0)
            if test_data["coverage"] > 0:
                conv_frac_test.append(test_data["conversions"] / test_data["coverage"])
            else:
                conv_frac_test.append(0.0)
        
        if conv_frac_ref and conv_frac_test:
            conv_frac_pearson = pearson(conv_frac_ref, conv_frac_test)
            conv_frac_spearman = spearman(conv_frac_ref, conv_frac_test)
        
        correlation_results.append({
            "threshold": threshold,
            "filter_field": filter_field,
            "n_genes": len(pairs),
            "ntr_pearson": ntr_pearson,
            "ntr_spearman": ntr_spearman,
            "conv_frac_pearson": conv_frac_pearson,
            "conv_frac_spearman": conv_frac_spearman,
        })

    # For backward compatibility, use first threshold for NTR bad check
    ntr_pairs = []
    ntr_bad = []
    if thresholds:
        threshold, filter_field = thresholds[0]
        for gene in shared:
            if ref[gene][filter_field] < threshold:
                continue
            r = ref[gene]["ntr"]
            t = test[gene]["ntr"]
            ntr_pairs.append((r, t, gene))
            if abs(r - t) > args.ntr_abs_max:
                ntr_bad.append((abs(r - t), gene, r, t))
        ntr_bad.sort(reverse=True)
    
    # Use first threshold result for backward compatibility
    corr = correlation_results[0]["ntr_pearson"] if correlation_results else float("nan")

    ok = True
    if read_deltas:
        if ref_missing_conv_cov:
            print(f"WARNING: ReadCount mismatches: {len(read_deltas)} (reference format lacks Conversions/Coverage)")
        else:
            ok = False
            print(f"FAIL: ReadCount mismatches: {len(read_deltas)}")
        for delta, gene, r, t in read_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")
    if conv_deltas:
        if ref_missing_conv_cov:
            print(f"WARNING: Conversions mismatches: {len(conv_deltas)} (reference format lacks Conversions/Coverage)")
        else:
            ok = False
            print(f"FAIL: Conversions mismatches: {len(conv_deltas)}")
        for delta, gene, r, t in conv_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")
    if cov_deltas:
        if ref_missing_conv_cov:
            print(f"WARNING: Coverage mismatches: {len(cov_deltas)} (reference format lacks Conversions/Coverage)")
        else:
            ok = False
            print(f"FAIL: Coverage mismatches: {len(cov_deltas)}")
        for delta, gene, r, t in cov_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")

    if math.isnan(corr) or corr < args.corr_min:
        if ref_missing_conv_cov:
            print(f"WARNING: NTR correlation {corr:.6f} (min {args.corr_min})")
        else:
            ok = False
            print(f"FAIL: NTR correlation {corr:.6f} (min {args.corr_min})")
    if ntr_bad:
        if ref_missing_conv_cov:
            print(f"WARNING: NTR abs diff > {args.ntr_abs_max}: {len(ntr_bad)}")
        else:
            ok = False
            print(f"FAIL: NTR abs diff > {args.ntr_abs_max}: {len(ntr_bad)}")
        for delta, gene, r, t in ntr_bad[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")

    # Print correlation table
    if correlation_results:
        print("\n=== Correlation Metrics (Pearson + Spearman) ===")
        print(f"{'Threshold':<12} {'Filter':<12} {'N Genes':<10} {'NTR Pearson':<15} {'NTR Spearman':<15} {'k/nT Pearson':<15} {'k/nT Spearman':<15}")
        print("-" * 100)
        for result in correlation_results:
            threshold_label = f">={result['threshold']:.0f}" if result['threshold'] >= 1 else f">={result['threshold']}"
            print(f"{threshold_label:<12} {result['filter_field']:<12} {result['n_genes']:<10} "
                  f"{result['ntr_pearson']:>14.6f} {result['ntr_spearman']:>14.6f} "
                  f"{result['conv_frac_pearson']:>14.6f} {result['conv_frac_spearman']:>14.6f}")
        print()

    if ok:
        print("PASS: STAR-Slam parity checks")
        print(f"  Genes compared: {len(shared)}")
        if correlation_results:
            first = correlation_results[0]
            print(f"  NTR correlation (Pearson): {first['ntr_pearson']:.6f} (min {args.corr_min})")
            print(f"  NTR correlation (Spearman): {first['ntr_spearman']:.6f}")
            if ref_missing_conv_cov or math.isnan(first["conv_frac_pearson"]):
                print("  Conversion fraction correlation: N/A (reference missing Conversions/Coverage)")
            else:
                print(f"  Conversion fraction correlation (Pearson): {first['conv_frac_pearson']:.6f}")
                print(f"  Conversion fraction correlation (Spearman): {first['conv_frac_spearman']:.6f}")
            if ntr_pairs:
                threshold, filter_field = thresholds[0]
                print(f"  NTR genes ({filter_field} >= {threshold}): {len(ntr_pairs)}")
        else:
            print(f"  NTR correlation: {corr:.6f} (min {args.corr_min})")
            if ntr_pairs:
                threshold, filter_field = thresholds[0] if thresholds else (args.min_read_count or 50.0, "readcount")
                print(f"  NTR genes ({filter_field} >= {threshold}): {len(ntr_pairs)}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
