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
    args = ap.parse_args()

    ref = read_star(args.reference)
    test = read_star(args.test)
    common = sorted(set(ref) & set(test))
    thresholds = [int(x) for x in args.thresholds.split(",") if x]

    print("=== STAR vs re-quant correlations ===")
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

if __name__ == "__main__":
    main()
