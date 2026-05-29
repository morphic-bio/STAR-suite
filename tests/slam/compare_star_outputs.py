#!/usr/bin/env python3
import argparse
import csv
import math

def read_star(path):
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["Gene"]
            data[gene] = {
                "readcount": float(row["ReadCount"]),
                "conversions": float(row["Conversions"]),
                "coverage": float(row["Coverage"]),
                "ntr": float(row["NTR"]),
            }
    return data

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denx = sum((x - mx) ** 2 for x in xs)
    deny = sum((y - my) ** 2 for y in ys)
    den = math.sqrt(denx * deny)
    return num / den if den else float("nan")

def spearman(xs, ys):
    def rank(vals):
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(vals):
            j = i
            while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + 1 + j + 1) / 2.0
            for k in range(i, j + 1):
                ranks[order[k]] = avg
            i = j + 1
        return ranks
    return pearson(rank(xs), rank(ys))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--thresholds", default="20,50,100")
    ap.add_argument("--min-pearson", type=float, default=None,
                    help="Fail if any reported Pearson correlation is below this value")
    ap.add_argument("--max-abs-delta", type=float, default=None,
                    help="Fail if any per-gene ReadCount/Conversions/Coverage/NTR delta exceeds this value")
    ap.add_argument("--require-same-genes", action="store_true",
                    help="Fail if the two outputs do not contain the same gene set")
    args = ap.parse_args()

    ref = read_star(args.reference)
    test = read_star(args.test)
    common = sorted(set(ref) & set(test))
    thresholds = [int(x) for x in args.thresholds.split(",") if x]

    print("=== STAR vs re-quant correlations ===")
    failed = False
    if args.require_same_genes:
        only_ref = sorted(set(ref) - set(test))
        only_test = sorted(set(test) - set(ref))
        if only_ref or only_test:
            failed = True
            print(f"FAIL: gene set mismatch only_ref={len(only_ref)} only_test={len(only_test)}")
            if only_ref:
                print("  sample only_ref:", ", ".join(only_ref[:10]))
            if only_test:
                print("  sample only_test:", ", ".join(only_test[:10]))

    if args.max_abs_delta is not None:
        for key, label in [
            ("readcount", "ReadCount"),
            ("conversions", "Conversions"),
            ("coverage", "Coverage"),
            ("ntr", "NTR"),
        ]:
            deltas = [
                (abs(ref[g][key] - test[g][key]), g, ref[g][key], test[g][key])
                for g in common
            ]
            deltas.sort(reverse=True)
            max_delta = deltas[0][0] if deltas else float("nan")
            bad = [d for d in deltas if d[0] > args.max_abs_delta]
            print(f"Exact {label}: max_abs_delta={max_delta:.12g} bad>{args.max_abs_delta}={len(bad)}")
            if bad:
                failed = True
                for delta, gene, r, t in bad[:10]:
                    print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")

    for thr in thresholds:
        genes = [g for g in common if ref[g]["readcount"] >= thr and test[g]["readcount"] >= thr]
        if len(genes) < 3:
            print(f"thr>={thr}: too few genes ({len(genes)})")
            continue
        for key, label in [("ntr", "NTR"), ("conversions", "Conversions"), ("coverage", "Coverage")]:
            xs = [ref[g][key] for g in genes]
            ys = [test[g][key] for g in genes]
            p = pearson(xs, ys)
            s = spearman(xs, ys)
            print(f"thr>={thr} {label}: Pearson={p:.6f} Spearman={s:.6f} (n={len(genes)})")
            if args.min_pearson is not None and (math.isnan(p) or p < args.min_pearson):
                failed = True
                print(f"FAIL: thr>={thr} {label} Pearson {p:.6f} < {args.min_pearson:.6f}")

    if failed:
        return 1
    if args.min_pearson is not None:
        print(f"PASS: all reported Pearson correlations >= {args.min_pearson:.6f}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
