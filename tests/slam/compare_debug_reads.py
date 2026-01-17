#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict


def normalize_loc(loc):
    if loc.startswith("chr"):
        return loc[3:]
    return loc


def load_star_reads(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ReadLoc" not in reader.fieldnames:
            raise RuntimeError("STAR debug reads file missing ReadLoc column.")
        data = {}
        for row in reader:
            gene = row.get("Gene", "")
            loc = normalize_loc(row.get("ReadLoc", ""))
            if not gene or not loc:
                continue
            key = (gene, loc)
            status = row.get("Status", "")
            try:
                weight = float(row.get("Weight", "0") or 0.0)
            except ValueError:
                weight = 0.0
            entry = data.setdefault(
                key,
                {
                    "pass_weight": 0.0,
                    "total_weight": 0.0,
                    "status_counts": Counter(),
                },
            )
            entry["status_counts"][status] += 1
            entry["total_weight"] += weight
            if status == "PASS":
                entry["pass_weight"] += weight
        return data


def load_gedi_debug(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Gene", "Read", "GeneConsistent", "Weight"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"GEDI debug file missing columns: {sorted(missing)}")
        data = {}
        for row in reader:
            gene = row.get("Gene", "")
            loc = normalize_loc(row.get("Read", ""))
            if not gene or not loc:
                continue
            key = (gene, loc)
            try:
                weight = float(row.get("Weight", "0") or 0.0)
            except ValueError:
                weight = 0.0
            try:
                consistent = int(row.get("GeneConsistent", "0") or 0)
            except ValueError:
                consistent = 0
            entry = data.setdefault(
                key,
                {
                    "consistent_weight": 0.0,
                    "inconsistent_weight": 0.0,
                    "consistent_count": 0,
                    "inconsistent_count": 0,
                },
            )
            if consistent:
                entry["consistent_weight"] += weight
                entry["consistent_count"] += 1
            else:
                entry["inconsistent_weight"] += weight
                entry["inconsistent_count"] += 1
        return data


def format_status_counts(counter):
    if not counter:
        return ""
    return ",".join(f"{k}:{v}" for k, v in counter.most_common())


def main():
    parser = argparse.ArgumentParser(
        description="Compare STAR vs GEDI debug read outputs by gene + read location."
    )
    parser.add_argument("--star", required=True, help="STAR <prefix>.reads.tsv")
    parser.add_argument("--gedi", required=True, help="GEDI slam_debug.tsv")
    parser.add_argument(
        "--delta-min",
        type=float,
        default=0.1,
        help="Minimum absolute weight delta to report (default: 0.1)",
    )
    parser.add_argument(
        "--max-report",
        type=int,
        default=20,
        help="Max rows to show in top mismatches (default: 20)",
    )
    parser.add_argument("--out", help="Optional TSV to write all mismatches")
    args = parser.parse_args()

    star = load_star_reads(args.star)
    gedi = load_gedi_debug(args.gedi)

    all_keys = set(star.keys()) | set(gedi.keys())
    rows = []
    for key in all_keys:
        star_entry = star.get(
            key, {"pass_weight": 0.0, "total_weight": 0.0, "status_counts": Counter()}
        )
        gedi_entry = gedi.get(
            key,
            {
                "consistent_weight": 0.0,
                "inconsistent_weight": 0.0,
                "consistent_count": 0,
                "inconsistent_count": 0,
            },
        )
        star_weight = star_entry["pass_weight"]
        gedi_weight = gedi_entry["consistent_weight"]
        delta = star_weight - gedi_weight
        if abs(delta) < args.delta_min:
            continue
        gene, loc = key
        rows.append(
            {
                "Gene": gene,
                "ReadLoc": loc,
                "STAR_PassWeight": star_weight,
                "GEDI_ConsistentWeight": gedi_weight,
                "Delta": delta,
                "STAR_StatusCounts": format_status_counts(star_entry["status_counts"]),
                "GEDI_ConsistentCount": gedi_entry["consistent_count"],
                "GEDI_InconsistentCount": gedi_entry["inconsistent_count"],
                "GEDI_InconsistentWeight": gedi_entry["inconsistent_weight"],
            }
        )

    rows.sort(key=lambda r: abs(r["Delta"]), reverse=True)

    star_only = sum(
        1
        for key in all_keys
        if star.get(key, {}).get("pass_weight", 0.0) > 0.0
        and gedi.get(key, {}).get("consistent_weight", 0.0) == 0.0
    )
    gedi_only = sum(
        1
        for key in all_keys
        if gedi.get(key, {}).get("consistent_weight", 0.0) > 0.0
        and star.get(key, {}).get("pass_weight", 0.0) == 0.0
    )

    print("=== Debug Read Comparison Summary ===")
    print(f"STAR keys: {len(star)}")
    print(f"GEDI keys: {len(gedi)}")
    print(f"STAR-only (pass weight > 0): {star_only}")
    print(f"GEDI-only (consistent weight > 0): {gedi_only}")
    print(f"Mismatches (|delta| >= {args.delta_min}): {len(rows)}")
    print("")
    print("=== Top Mismatches ===")
    for row in rows[: args.max_report]:
        print(
            f"{row['Gene']}\t{row['ReadLoc']}\t"
            f"STAR={row['STAR_PassWeight']:.6f}\tGEDI={row['GEDI_ConsistentWeight']:.6f}\t"
            f"Delta={row['Delta']:.6f}\tSTAR_Status={row['STAR_StatusCounts']}\t"
            f"GEDI_consistent={row['GEDI_ConsistentCount']}\tGEDI_inconsistent={row['GEDI_InconsistentCount']}"
        )

    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="") as handle:
            fieldnames = [
                "Gene",
                "ReadLoc",
                "STAR_PassWeight",
                "GEDI_ConsistentWeight",
                "Delta",
                "STAR_StatusCounts",
                "GEDI_ConsistentCount",
                "GEDI_InconsistentCount",
                "GEDI_InconsistentWeight",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)


if __name__ == "__main__":
    main()
