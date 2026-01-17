#!/usr/bin/env python3
import argparse
import gzip
import os


def normalize_chrom(chrom):
    if chrom.startswith("chr"):
        if chrom in ("chrMT", "chrM"):
            return "chrM"
        return chrom
    if chrom in ("M", "MT"):
        return "chrM"
    if chrom in {"X", "Y"}:
        return "chr" + chrom
    if chrom.isdigit():
        val = int(chrom)
        if 1 <= val <= 22:
            return "chr" + chrom
    return chrom


def open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_bed_positions(path, max_span):
    positions = set()
    total_lines = 0
    used_lines = 0
    with open_any(path) as handle:
        for line in handle:
            total_lines += 1
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = normalize_chrom(parts[0])
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            span = end - start
            if span > max_span:
                continue
            used_lines += 1
            if span == 1:
                positions.add((chrom, start))
            else:
                for pos in range(start, end):
                    positions.add((chrom, pos))
    return positions, total_lines, used_lines


def summarize_overlap(name_a, set_a, name_b, set_b):
    inter = set_a & set_b
    union = set_a | set_b
    a_count = len(set_a)
    b_count = len(set_b)
    inter_count = len(inter)
    union_count = len(union)
    a_pct = (inter_count / a_count) if a_count else 0.0
    b_pct = (inter_count / b_count) if b_count else 0.0
    jacc = (inter_count / union_count) if union_count else 0.0
    subset = ""
    if inter_count == a_count and a_count <= b_count:
        subset = f"{name_a} subset of {name_b}"
    elif inter_count == b_count and b_count <= a_count:
        subset = f"{name_b} subset of {name_a}"
    return {
        "a_count": a_count,
        "b_count": b_count,
        "inter": inter_count,
        "a_pct": a_pct,
        "b_pct": b_pct,
        "jacc": jacc,
        "subset": subset,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare SNP mask BED overlap.")
    parser.add_argument("--vcf-mask", required=True, help="VCF-derived mask BED (.bed or .bed.gz)")
    parser.add_argument("--mask-a", required=True, help="Reference mask A BED")
    parser.add_argument("--mask-b", required=True, help="Reference mask B BED")
    parser.add_argument("--max-span", type=int, default=1000, help="Max span per BED line to expand")
    parser.add_argument("--out", required=True, help="Output report path")
    args = parser.parse_args()

    vcf_set, vcf_lines, vcf_used = load_bed_positions(args.vcf_mask, args.max_span)
    a_set, a_lines, a_used = load_bed_positions(args.mask_a, args.max_span)
    b_set, b_lines, b_used = load_bed_positions(args.mask_b, args.max_span)

    stats_a = summarize_overlap("vcf", vcf_set, "mask_a", a_set)
    stats_b = summarize_overlap("vcf", vcf_set, "mask_b", b_set)

    with open(args.out, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"vcf_mask\t{args.vcf_mask}\n")
        out.write(f"mask_a\t{args.mask_a}\n")
        out.write(f"mask_b\t{args.mask_b}\n")
        out.write(f"vcf_lines_total\t{vcf_lines}\n")
        out.write(f"vcf_lines_used\t{vcf_used}\n")
        out.write(f"mask_a_lines_total\t{a_lines}\n")
        out.write(f"mask_a_lines_used\t{a_used}\n")
        out.write(f"mask_b_lines_total\t{b_lines}\n")
        out.write(f"mask_b_lines_used\t{b_used}\n")
        out.write(f"vcf_sites\t{stats_a['a_count']}\n")
        out.write(f"mask_a_sites\t{stats_a['b_count']}\n")
        out.write(f"mask_b_sites\t{stats_b['b_count']}\n")
        out.write("\n")
        out.write("comparison\tvalue\n")
        out.write(f"vcf_vs_mask_a_overlap\t{stats_a['inter']}\n")
        out.write(f"vcf_vs_mask_a_pct\t{stats_a['a_pct']:.6f}\n")
        out.write(f"mask_a_vs_vcf_pct\t{stats_a['b_pct']:.6f}\n")
        out.write(f"vcf_vs_mask_a_jaccard\t{stats_a['jacc']:.6f}\n")
        if stats_a["subset"]:
            out.write(f"vcf_vs_mask_a_subset\t{stats_a['subset']}\n")
        out.write(f"vcf_vs_mask_b_overlap\t{stats_b['inter']}\n")
        out.write(f"vcf_vs_mask_b_pct\t{stats_b['a_pct']:.6f}\n")
        out.write(f"mask_b_vs_vcf_pct\t{stats_b['b_pct']:.6f}\n")
        out.write(f"vcf_vs_mask_b_jaccard\t{stats_b['jacc']:.6f}\n")
        if stats_b["subset"]:
            out.write(f"vcf_vs_mask_b_subset\t{stats_b['subset']}\n")

    print(f"Wrote overlap report to {args.out}")


if __name__ == "__main__":
    main()
