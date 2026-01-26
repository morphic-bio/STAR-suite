#!/usr/bin/env python3
"""Validate perturb-seq compatibility output layout and MEX format.

Checks:
- CR9-style directory layout (per-sample MEX and crispr_analysis/)
- MEX files are gzip-compressed (optional override)
- Matrix Market header and dimension line
- features.tsv has >=3 columns
- barcodes.tsv non-empty

Layout is defined in a JSON file for easy extension.
"""

import argparse
import gzip
import json
import sys
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def first_nonempty_line(handle):
    for line in handle:
        line = line.strip()
        if line:
            return line
    return ""


def resolve_file(base_dir: Path, name: str, require_gzip: bool, allow_legacy: bool) -> Path:
    legacy_map = {
        "barcodes.tsv": "barcodes.txt",
        "features.tsv": "features.txt",
    }

    if name.endswith(".gz"):
        candidate = base_dir / name
        if candidate.exists():
            return candidate
        return Path()

    if require_gzip:
        candidate = base_dir / (name + ".gz")
        if candidate.exists():
            return candidate
        return Path()

    # Allow either gz or plain
    candidate = base_dir / name
    if candidate.exists():
        return candidate
    candidate_gz = base_dir / (name + ".gz")
    if candidate_gz.exists():
        return candidate_gz

    if allow_legacy and name in legacy_map:
        legacy_candidate = base_dir / legacy_map[name]
        if legacy_candidate.exists():
            return legacy_candidate

    return Path()


def validate_matrix(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"missing: {path}"
    if path.stat().st_size == 0:
        return False, f"empty: {path}"

    with open_text(path) as handle:
        header = handle.readline().strip()
        if header != "%%MatrixMarket matrix coordinate integer general":
            return False, f"bad header: {header}"

        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        parts = line.strip().split()
        if len(parts) != 3:
            return False, f"bad dimensions line: {line.strip()}"
        try:
            rows, cols, entries = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            return False, f"non-integer dimensions: {line.strip()}"

    return True, f"OK {rows}x{cols}, entries={entries}"


def validate_features(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"missing: {path}"
    if path.stat().st_size == 0:
        return False, f"empty: {path}"

    with open_text(path) as handle:
        line = first_nonempty_line(handle)
        if not line:
            return False, "no data lines"
        cols = line.split("\t")
        if len(cols) < 3:
            return False, f"expected >=3 columns, got {len(cols)}"
    return True, "OK"


def validate_barcodes(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"missing: {path}"
    if path.stat().st_size == 0:
        return False, f"empty: {path}"

    with open_text(path) as handle:
        line = first_nonempty_line(handle)
        if not line:
            return False, "no barcodes"
    return True, "OK"


def validate_generic(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"missing: {path}"
    if path.stat().st_size == 0:
        return False, f"empty: {path}"
    return True, "OK"


def load_layout(path: Path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def discover_samples(root: Path, mex_dir_names, crispr_dir_name):
    samples = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        if (entry / crispr_dir_name).exists():
            samples.append(entry.name)
            continue
        for mex_dir in mex_dir_names:
            if (entry / mex_dir).exists():
                samples.append(entry.name)
                break
    return sorted(set(samples))


def validate_sample(root: Path, sample: str, layout: dict, allow_legacy: bool, require_gzip_override):
    mex_conf = layout["mex"]
    mex_dirs = list(mex_conf.get("dir_names", []))
    if allow_legacy:
        mex_dirs.extend(mex_conf.get("legacy_dir_names", []))

    require_gzip = mex_conf.get("require_gzip", True)
    if require_gzip_override is not None:
        require_gzip = require_gzip_override

    sample_root = root / sample

    mex_dir = None
    for name in mex_dirs:
        candidate = sample_root / name
        if candidate.exists():
            mex_dir = candidate
            break

    if mex_dir is None:
        return False, [f"missing MEX dir under {sample_root} (expected one of: {', '.join(mex_dirs)})"]

    errors = []
    for entry in mex_conf.get("files", []):
        file_path = resolve_file(mex_dir, entry["name"], require_gzip, allow_legacy)
        if not file_path.exists():
            errors.append(f"missing {entry['name']} (gz required={require_gzip})")
            continue

        if entry.get("type") == "matrix_market":
            ok, msg = validate_matrix(file_path)
        elif entry.get("type") == "features":
            ok, msg = validate_features(file_path)
        elif entry.get("type") == "barcodes":
            ok, msg = validate_barcodes(file_path)
        else:
            ok, msg = validate_generic(file_path)

        if not ok:
            errors.append(f"{entry['name']}: {msg}")

    crispr_conf = layout.get("crispr_analysis")
    if crispr_conf:
        crispr_dir = sample_root / crispr_conf["dir_name"]
        if not crispr_dir.exists():
            errors.append(f"missing {crispr_conf['dir_name']}/")
        else:
            for entry in crispr_conf.get("files", []):
                file_path = crispr_dir / entry["name"]
                ok, msg = validate_generic(file_path)
                if not ok:
                    errors.append(f"{entry['name']}: {msg}")

    return len(errors) == 0, errors


def main():
    parser = argparse.ArgumentParser(description="Validate perturb-seq output layout and MEX format")
    parser.add_argument("--root", required=True, type=Path, help="Root output directory")
    parser.add_argument("--layout", default=Path(__file__).parent / "perturb_seq_layout.json", type=Path,
                        help="Layout JSON file")
    parser.add_argument("--samples", nargs="+", help="Sample directory names (defaults to auto-discovery)")
    parser.add_argument("--allow-legacy-filtered", action="store_true",
                        help="Allow legacy 'filtered' MEX dir name")
    parser.add_argument("--require-gzip", action="store_true", help="Require .gz MEX files (default from layout)")
    parser.add_argument("--allow-plain", action="store_true", help="Allow non-gz MEX files")

    args = parser.parse_args()

    layout = load_layout(args.layout)

    mex_dir_names = layout["mex"].get("dir_names", [])
    if args.allow_legacy_filtered:
        mex_dir_names.extend(layout["mex"].get("legacy_dir_names", []))

    samples = args.samples or discover_samples(args.root, mex_dir_names, layout["crispr_analysis"]["dir_name"])
    if not samples:
        print("ERROR: No samples found under root", file=sys.stderr)
        sys.exit(2)

    require_gzip_override = None
    if args.require_gzip and args.allow_plain:
        print("ERROR: --require-gzip and --allow-plain are mutually exclusive", file=sys.stderr)
        sys.exit(2)
    if args.require_gzip:
        require_gzip_override = True
    elif args.allow_plain:
        require_gzip_override = False

    all_ok = True
    for sample in samples:
        ok, errors = validate_sample(args.root, sample, layout, args.allow_legacy_filtered, require_gzip_override)
        if ok:
            print(f"{sample}: OK")
        else:
            all_ok = False
            print(f"{sample}: FAIL")
            for err in errors:
                print(f"  - {err}")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
