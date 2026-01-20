#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import sys


def normalize_header(value):
    return value.strip().lower()


def load_feature_csv(path):
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("feature CSV is empty")

        header_map = {normalize_header(name): idx for idx, name in enumerate(header)}
        name_idx = header_map.get("name")
        id_idx = header_map.get("id")
        ftype_idx = header_map.get("feature_type")
        if ftype_idx is None:
            ftype_idx = header_map.get("type")

        if name_idx is None and id_idx is None:
            raise ValueError("feature CSV header must include 'name' or 'id'")

        rows = []
        for row in reader:
            if not row:
                continue
            name = row[name_idx].strip() if name_idx is not None and name_idx < len(row) else ""
            fid = row[id_idx].strip() if id_idx is not None and id_idx < len(row) else ""
            ftype = row[ftype_idx].strip() if ftype_idx is not None and ftype_idx < len(row) else ""
            rows.append((fid, name, ftype))

    if not rows:
        raise ValueError("feature CSV has no data rows")
    return rows


def read_features_txt(path):
    with open(path, "r") as handle:
        return [line.strip() for line in handle if line.strip()]


def ensure_dir(path):
    if not os.path.isdir(path):
        raise ValueError(f"missing output directory: {path}")


def write_features_tsv(out_path, feature_rows, default_type, force):
    if os.path.exists(out_path) and not force:
        return False
    with open(out_path, "w", newline="") as handle:
        for fid, name, ftype in feature_rows:
            if not fid:
                fid = name
            if not name:
                name = fid
            if not ftype:
                ftype = default_type
            handle.write(f"{fid}\t{name}\t{ftype}\n")
    return True


def copy_barcodes(barcodes_txt, barcodes_tsv, force):
    if not os.path.exists(barcodes_txt):
        return False
    if os.path.exists(barcodes_tsv) and not force:
        return False
    shutil.copyfile(barcodes_txt, barcodes_tsv)
    return True


def compare_feature_names(feature_rows, features_txt):
    names = [name for _, name, _ in feature_rows]
    if len(names) != len(features_txt):
        return f"feature count mismatch: csv={len(names)} features.txt={len(features_txt)}"
    for idx, (a, b) in enumerate(zip(names, features_txt), start=1):
        if a != b:
            return f"name mismatch at row {idx}: csv='{a}' features.txt='{b}'"
    return ""


def run(args):
    feature_rows = load_feature_csv(args.feature_csv)
    outdirs = [args.assign_out]
    filtered_dir = os.path.join(args.assign_out, "filtered")
    if os.path.isdir(filtered_dir):
        outdirs.append(filtered_dir)

    warnings = []
    wrote_any = False

    for outdir in outdirs:
        ensure_dir(outdir)
        features_txt = os.path.join(outdir, "features.txt")
        barcodes_txt = os.path.join(outdir, "barcodes.txt")
        barcodes_tsv = os.path.join(outdir, "barcodes.tsv")
        features_tsv = os.path.join(outdir, "features.tsv")

        if os.path.exists(features_txt):
            features_txt_rows = read_features_txt(features_txt)
            warn = compare_feature_names(feature_rows, features_txt_rows)
            if warn:
                warnings.append(f"{outdir}: {warn}")
        else:
            warnings.append(f"{outdir}: features.txt not found")

        if write_features_tsv(features_tsv, feature_rows, args.default_feature_type, args.force):
            wrote_any = True
        if copy_barcodes(barcodes_txt, barcodes_tsv, args.force):
            wrote_any = True

    for warn in warnings:
        print(f"WARNING: {warn}", file=sys.stderr)

    if not wrote_any:
        print("No outputs written (files may already exist).", file=sys.stderr)
        return 1

    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Create 10x-style features.tsv/barcodes.tsv from assignBarcodes outputs"
    )
    parser.add_argument("--assign-out", required=True, help="assignBarcodes output directory")
    parser.add_argument("--feature-csv", required=True, help="feature reference CSV used by assignBarcodes")
    parser.add_argument(
        "--default-feature-type",
        default="Custom",
        help="feature_type fallback when missing (default: Custom)",
    )
    parser.add_argument("--force", action="store_true", help="overwrite existing TSVs")
    args = parser.parse_args()

    sys.exit(run(args))


if __name__ == "__main__":
    main()
