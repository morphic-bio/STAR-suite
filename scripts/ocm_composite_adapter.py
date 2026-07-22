#!/usr/bin/env python3
"""
Utilities for the JAX OCM composite-barcode smoke runs.

The adapter is intentionally explicit and file-oriented. It supports both the
old helper-only smoke path and the native STAR OCM-Flex barcode mode:

  legacy helper: raw CB16 + observed OCM code -> effective CB17
  native STAR:   raw CB16 -> effective CB16+OCM_TAG8 inside STAR
  materialize:   effective composite CB -> per-OCM sample MEX with CB16 labels
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import json
import math
import os
import resource
import shutil
import subprocess
import sys
import time
from array import array
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


OVERHANG_TO_OCM = {
    "GT": ("OB1", "A"),
    "CA": ("OB2", "C"),
    "TC": ("OB3", "G"),
    "AG": ("OB4", "T"),
}
CODE_TO_OCM = {code: ob for ob, code in (v for v in OVERHANG_TO_OCM.values())}
OCM_TO_CODE = {ob: code for _overhang, (ob, code) in OVERHANG_TO_OCM.items()}
OCM_TO_TAG8 = {
    "OB1": "GTGTGTGT",
    "OB2": "CACACACA",
    "OB3": "TCTCTCTC",
    "OB4": "AGAGAGAG",
}
TAG8_TO_OCM = {tag8: ob for ob, tag8 in OCM_TO_TAG8.items()}
OCM_ORDER = ["OB1", "OB2", "OB3", "OB4"]


def log(message: str) -> None:
    print(f"[ocm_composite] {message}", file=sys.stderr, flush=True)


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def open_text(path: Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz" or str(path).endswith(".gz"):
        kwargs = {"encoding": "utf-8"}
        if "w" in mode:
            kwargs["compresslevel"] = 1
        return gzip.open(path, mode, **kwargs)
    return open(path, mode, encoding="utf-8")


def resolve_file(directory: Path, basename: str) -> Path:
    direct = directory / basename
    gz = directory / f"{basename}.gz"
    if direct.exists():
        return direct
    if gz.exists():
        return gz
    die(f"Missing {basename}(.gz) in {directory}")


def read_lines(path: Path) -> List[str]:
    with open_text(path, "rt") as handle:
        return [line.rstrip("\n\r") for line in handle]


def write_lines_gz(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=1) as handle:
        for line in lines:
            handle.write(line)
            handle.write("\n")


def normalize_barcode(barcode: str) -> str:
    barcode = barcode.strip()
    if "\t" in barcode:
        barcode = barcode.split("\t", 1)[0]
    if "," in barcode:
        barcode = barcode.split(",")[-1].strip()
    if "-" in barcode:
        base, suffix = barcode.rsplit("-", 1)
        if suffix.isdigit():
            barcode = base
    return barcode


def output_cb16_from_composite(barcode: str) -> str:
    base = normalize_barcode(barcode)
    if len(base) >= 24 and base[-8:] in TAG8_TO_OCM:
        return f"{base[:-8][:16]}-1"
    if len(base) >= 17 and base[16] in CODE_TO_OCM:
        return f"{base[:16]}-1"
    return f"{base[:16]}-1"


def classify_cb16(cb16: str) -> Tuple[Optional[str], str]:
    if len(cb16) < 9:
        return None, "N"
    entry = OVERHANG_TO_OCM.get(cb16[7:9])
    if entry is None:
        return None, "N"
    return entry


def classify_composite_barcode(barcode: str) -> Optional[str]:
    base = normalize_barcode(barcode)
    if len(base) >= 24 and base[-8:] in TAG8_TO_OCM:
        return TAG8_TO_OCM[base[-8:]]
    if len(base) >= 17 and base[16] in CODE_TO_OCM:
        return CODE_TO_OCM[base[16]]
    ob, _code = classify_cb16(base[:16])
    return ob


def parse_ocm_config(config_path: Path) -> List[Dict[str, object]]:
    samples: List[Dict[str, object]] = []
    section = ""
    header: List[str] = []
    with open(config_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line.strip("[]").lower()
                header = []
                continue
            if section != "samples" or line.startswith("#"):
                continue
            row = next(csv.reader([line]))
            if not header:
                header = [x.strip() for x in row]
                continue
            data = {header[i]: row[i].strip() if i < len(row) else "" for i in range(len(header))}
            sample_id = data.get("sample_id", "")
            tags = [x.strip() for x in data.get("ocm_barcode_ids", "").split("|") if x.strip()]
            if not sample_id or not tags:
                die(f"Invalid [samples] row in {config_path}: {line}")
            bad = [tag for tag in tags if tag not in OCM_TO_CODE]
            if bad:
                die(f"Unsupported OCM tag(s) for {sample_id}: {','.join(bad)}")
            samples.append(
                {
                    "sample_id": sample_id,
                    "ocm_ids": tags,
                    "description": data.get("description", ""),
                }
            )
    if not samples:
        die(f"No [samples] entries found in {config_path}")
    return samples


def source_fastq(raw_dir: Path, sample_stem: str, lane: str, read: str) -> Path:
    return raw_dir / f"{sample_stem}_{lane}_{read}_001.fastq.gz"


def hardlink_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def build_composite_whitelist(source: Path, output: Path) -> Dict[str, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = {ob: 0 for ob in OCM_ORDER}
    counts["unknown"] = 0
    with open_text(source, "rt") as src, open(output, "w", encoding="utf-8") as out:
        for raw in src:
            cb = raw.strip().split()[0]
            ob, code = classify_cb16(cb)
            if ob is None:
                counts["unknown"] += 1
                continue
            out.write(f"{cb}{code}\n")
            counts[ob] += 1
    return counts


def count_whitelist_ocm(source: Path) -> Dict[str, int]:
    counts = {ob: 0 for ob in OCM_ORDER}
    counts["unknown"] = 0
    with open_text(source, "rt") as src:
        for raw in src:
            cb = raw.strip().split()[0]
            ob, _code = classify_cb16(cb)
            counts[ob or "unknown"] += 1
    return counts


def convert_r1_record(seq: str, qual: str) -> Tuple[str, str, Optional[str]]:
    seq = seq.rstrip("\n\r")
    qual = qual.rstrip("\n\r")
    cb16 = seq[:16]
    ob, code = classify_cb16(cb16)
    return f"{cb16}{code}{seq[16:]}", f"{qual[:16]}I{qual[16:]}", ob


def prepare(args: argparse.Namespace) -> None:
    raw_dir = Path(args.raw_dir)
    cr_fastq_dir = Path(args.cr_fastq_dir)
    star_fastq_dir = Path(args.star_fastq_dir)
    lanes = [x.strip() for x in args.lanes.split(",") if x.strip()]
    if not lanes:
        die("--lanes must include at least one lane")

    read_pairs = int(args.read_pairs)
    base = read_pairs // len(lanes)
    remainder = read_pairs % len(lanes)

    if args.force:
        for path in (cr_fastq_dir, star_fastq_dir):
            if path.exists():
                shutil.rmtree(path)
    cr_fastq_dir.mkdir(parents=True, exist_ok=True)
    star_fastq_dir.mkdir(parents=True, exist_ok=True)

    if args.star_mode == "legacy17":
        whitelist_stats = build_composite_whitelist(Path(args.whitelist), Path(args.composite_whitelist))
    else:
        whitelist_stats = count_whitelist_ocm(Path(args.whitelist))
    lane_rows = []
    total_tag_counts = {ob: 0 for ob in OCM_ORDER}
    total_tag_counts["unknown"] = 0

    for lane_index, lane in enumerate(lanes):
        lane_pairs = base + (1 if lane_index < remainder else 0)
        src_r1 = source_fastq(raw_dir, args.sample_stem, lane, "R1")
        src_r2 = source_fastq(raw_dir, args.sample_stem, lane, "R2")
        if not src_r1.exists() or not src_r2.exists():
            die(f"Missing source FASTQs for {lane}: {src_r1}, {src_r2}")

        cr_r1 = cr_fastq_dir / src_r1.name
        cr_r2 = cr_fastq_dir / src_r2.name
        star_r1 = star_fastq_dir / src_r1.name
        star_r2 = star_fastq_dir / src_r2.name

        log(f"staging {lane}: {lane_pairs} read pairs")
        tag_counts = {ob: 0 for ob in OCM_ORDER}
        tag_counts["unknown"] = 0
        written = 0
        with ExitStack() as stack:
            r1_in = stack.enter_context(gzip.open(src_r1, "rt", encoding="utf-8"))
            r2_in = stack.enter_context(gzip.open(src_r2, "rt", encoding="utf-8"))
            cr_r1_out = stack.enter_context(gzip.open(cr_r1, "wt", encoding="utf-8", compresslevel=1))
            cr_r2_out = stack.enter_context(gzip.open(cr_r2, "wt", encoding="utf-8", compresslevel=1))
            star_r1_out = None
            if args.star_mode == "legacy17":
                star_r1_out = stack.enter_context(gzip.open(star_r1, "wt", encoding="utf-8", compresslevel=1))
            for _ in range(lane_pairs):
                rec1 = [r1_in.readline() for _ in range(4)]
                rec2 = [r2_in.readline() for _ in range(4)]
                if not all(rec1) or not all(rec2):
                    die(f"{lane} ended before {lane_pairs} read pairs; wrote {written}")
                cr_r1_out.writelines(rec1)
                cr_r2_out.writelines(rec2)
                if args.star_mode == "legacy17":
                    assert star_r1_out is not None
                    new_seq, new_qual, ob = convert_r1_record(rec1[1], rec1[3])
                    star_r1_out.write(rec1[0])
                    star_r1_out.write(new_seq + "\n")
                    star_r1_out.write(rec1[2])
                    star_r1_out.write(new_qual + "\n")
                else:
                    ob, _code = classify_cb16(rec1[1].strip()[:16])
                tag_counts[ob or "unknown"] += 1
                written += 1

        if args.star_mode == "legacy17":
            hardlink_or_copy(cr_r2, star_r2)
        else:
            hardlink_or_copy(cr_r1, star_r1)
            hardlink_or_copy(cr_r2, star_r2)
        for key, value in tag_counts.items():
            total_tag_counts[key] += value
        lane_rows.append(
            {
                "lane": lane,
                "read_pairs": written,
                "cr_r1": str(cr_r1),
                "cr_r2": str(cr_r2),
                "star_r1": str(star_r1),
                "star_r2": str(star_r2),
                **{f"reads_{key}": value for key, value in tag_counts.items()},
            }
        )

    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "lane",
            "read_pairs",
            "cr_r1",
            "cr_r2",
            "star_r1",
            "star_r2",
            *(f"reads_{key}" for key in [*OCM_ORDER, "unknown"]),
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in lane_rows:
            writer.writerow(row)

    stats = {
        "sample_id": args.sample_id,
        "sample_stem": args.sample_stem,
        "read_pairs": read_pairs,
        "lanes": lane_rows,
        "read_ocm_counts": total_tag_counts,
        "whitelist_ocm_counts": whitelist_stats,
        "star_mode": args.star_mode,
        "composite_whitelist": str(Path(args.composite_whitelist)) if args.star_mode == "legacy17" else None,
    }
    stats_path = Path(args.stats_json)
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.star_mode == "legacy17":
        log(f"prepared composite FASTQs and whitelist: {stats_path}")
    else:
        log(f"prepared native STAR FASTQ hardlinks and original-whitelist stats: {stats_path}")


def read_mtx_shape(matrix_path: Path) -> Tuple[int, int, int]:
    with open_text(matrix_path, "rt") as handle:
        for raw in handle:
            if raw.startswith("%"):
                continue
            fields = raw.strip().split()
            if len(fields) >= 3:
                return int(fields[0]), int(fields[1]), int(fields[2])
    die(f"Could not read MatrixMarket shape from {matrix_path}")


def iter_mtx_entries(matrix_path: Path):
    with open_text(matrix_path, "rt") as handle:
        saw_shape = False
        for raw in handle:
            if raw.startswith("%"):
                continue
            if not saw_shape:
                saw_shape = True
                continue
            fields = raw.strip().split()
            if len(fields) < 3:
                continue
            yield int(fields[0]), int(fields[1]), fields[2]


def write_features(out_dir: Path, features_lines: Sequence[str]) -> None:
    write_lines_gz(out_dir / "features.tsv.gz", features_lines)


def build_column_maps(
    barcodes_path: Path, samples: Sequence[Dict[str, object]]
) -> Tuple[array, List[array], List[List[str]], Dict[str, int]]:
    tag_to_samples: Dict[str, List[int]] = {tag: [] for tag in OCM_ORDER}
    for sample_idx, sample in enumerate(samples):
        for tag in sample["ocm_ids"]:  # type: ignore[index]
            tag_to_samples[str(tag)].append(sample_idx)

    raw_barcodes = read_lines(barcodes_path)
    ncols = len(raw_barcodes)
    col_tag = array("b", [-1]) * (ncols + 1)
    col_new_by_sample = [array("i", [0]) * (ncols + 1) for _sample in samples]
    sample_barcodes: List[List[str]] = [[] for _sample in samples]
    tag_counts = {tag: 0 for tag in OCM_ORDER}
    tag_counts["unknown"] = 0

    tag_index = {tag: i for i, tag in enumerate(OCM_ORDER)}
    for col, barcode in enumerate(raw_barcodes, start=1):
        tag = classify_composite_barcode(barcode)
        if tag is None:
            tag_counts["unknown"] += 1
            continue
        col_tag[col] = tag_index[tag]
        tag_counts[tag] += 1
        output_bc = output_cb16_from_composite(barcode)
        for sample_idx in tag_to_samples[tag]:
            sample_barcodes[sample_idx].append(output_bc)
            col_new_by_sample[sample_idx][col] = len(sample_barcodes[sample_idx])

    return col_tag, col_new_by_sample, sample_barcodes, tag_counts


def split_matrix_layers_by_maps(
    input_dir: Path,
    output_dirs: Sequence[Path],
    samples: Sequence[Dict[str, object]],
    matrix_names: Sequence[str],
) -> Dict[str, object]:
    barcodes_path = resolve_file(input_dir, "barcodes.tsv")
    features_path = resolve_file(input_dir, "features.tsv")
    features_lines = read_lines(features_path)
    col_tag, col_new_by_sample, sample_barcodes, tag_counts = build_column_maps(barcodes_path, samples)

    tag_to_samples: Dict[int, List[int]] = {i: [] for i in range(len(OCM_ORDER))}
    for sample_idx, sample in enumerate(samples):
        for tag in sample["ocm_ids"]:  # type: ignore[index]
            tag_to_samples[OCM_ORDER.index(str(tag))].append(sample_idx)

    for out_dir, barcodes in zip(output_dirs, sample_barcodes):
        out_dir.mkdir(parents=True, exist_ok=True)
        write_features(out_dir, features_lines)
        write_lines_gz(out_dir / "barcodes.tsv.gz", barcodes)

    layer_stats = {}
    for matrix_name in matrix_names:
        matrix_path = resolve_file(input_dir, matrix_name)
        nrows, ncols, _nnz = read_mtx_shape(matrix_path)
        if ncols + 1 != len(col_tag):
            die(f"{matrix_path} columns ({ncols}) do not match barcode count ({len(col_tag) - 1})")
        nnz_counts = [0 for _sample in samples]
        for _row, col, _val in iter_mtx_entries(matrix_path):
            if col <= 0 or col >= len(col_tag):
                continue
            tag_idx = col_tag[col]
            if tag_idx < 0:
                continue
            for sample_idx in tag_to_samples[int(tag_idx)]:
                if col_new_by_sample[sample_idx][col] > 0:
                    nnz_counts[sample_idx] += 1

        writers = []
        try:
            for sample_idx, out_dir in enumerate(output_dirs):
                layer_out = out_dir / f"{matrix_name}.gz"
                handle = gzip.open(layer_out, "wt", encoding="utf-8", compresslevel=1)
                handle.write("%%MatrixMarket matrix coordinate integer general\n")
                handle.write("%\n")
                handle.write(f"{nrows} {len(sample_barcodes[sample_idx])} {nnz_counts[sample_idx]}\n")
                writers.append(handle)
            for row, col, val in iter_mtx_entries(matrix_path):
                if col <= 0 or col >= len(col_tag):
                    continue
                tag_idx = col_tag[col]
                if tag_idx < 0:
                    continue
                for sample_idx in tag_to_samples[int(tag_idx)]:
                    new_col = col_new_by_sample[sample_idx][col]
                    if new_col > 0:
                        writers[sample_idx].write(f"{row} {new_col} {val}\n")
        finally:
            for writer in writers:
                writer.close()
        layer_stats[matrix_name] = {
            samples[sample_idx]["sample_id"]: nnz_counts[sample_idx] for sample_idx in range(len(samples))
        }

    velocity_layers = [layer for layer in ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx") if layer in matrix_names]
    if "matrix.mtx" not in matrix_names and velocity_layers:
        total_maps: List[Dict[Tuple[int, int], int]] = [dict() for _sample in samples]
        nrows_total: Optional[int] = None
        for layer in velocity_layers:
            matrix_path = resolve_file(input_dir, layer)
            nrows, ncols, _nnz = read_mtx_shape(matrix_path)
            if nrows_total is None:
                nrows_total = nrows
            elif nrows_total != nrows:
                die(f"Velocity layer row mismatch while synthesizing matrix.mtx: {layer}")
            if ncols + 1 != len(col_tag):
                die(f"{matrix_path} columns ({ncols}) do not match barcode count ({len(col_tag) - 1})")
            for row, col, val_s in iter_mtx_entries(matrix_path):
                if col <= 0 or col >= len(col_tag):
                    continue
                tag_idx = col_tag[col]
                if tag_idx < 0:
                    continue
                val = int(float(val_s))
                if val == 0:
                    continue
                for sample_idx in tag_to_samples[int(tag_idx)]:
                    new_col = col_new_by_sample[sample_idx][col]
                    if new_col > 0:
                        key = (row, int(new_col))
                        total_maps[sample_idx][key] = total_maps[sample_idx].get(key, 0) + val
        nnz_counts = []
        for sample_idx, out_dir in enumerate(output_dirs):
            entries = sorted(total_maps[sample_idx].items())
            nnz_counts.append(len(entries))
            with gzip.open(out_dir / "matrix.mtx.gz", "wt", encoding="utf-8", compresslevel=1) as out:
                out.write("%%MatrixMarket matrix coordinate integer general\n")
                out.write("%\n")
                out.write(f"{nrows_total or 0} {len(sample_barcodes[sample_idx])} {len(entries)}\n")
                for (row, col), val in entries:
                    out.write(f"{row} {col} {val}\n")
        layer_stats["matrix.mtx"] = {
            samples[sample_idx]["sample_id"]: nnz_counts[sample_idx] for sample_idx in range(len(samples))
        }

    return {
        "input_dir": str(input_dir),
        "matrix_names": list(matrix_names),
        "tag_counts": tag_counts,
        "sample_columns": {
            str(samples[i]["sample_id"]): len(sample_barcodes[i]) for i in range(len(samples))
        },
        "layer_nnz": layer_stats,
    }


def peak_rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports kilobytes; macOS reports bytes. This host is Linux, but keep
    # the conversion defensive for local developer runs.
    if usage > 10_000_000_000:
        return usage / (1024.0 * 1024.0)
    return usage / 1024.0


def now_s() -> float:
    return time.perf_counter()


def log_timing(timings: Dict[str, Dict[str, float]], name: str, started: float) -> None:
    elapsed = time.perf_counter() - started
    timings[name] = {"seconds": elapsed, "peak_rss_mb": peak_rss_mb()}
    log(f"{name} completed in {elapsed:.2f}s; peak_rss_mb={peak_rss_mb():.1f}")


def make_tag_to_sample_indices(samples: Sequence[Dict[str, object]]) -> Dict[str, List[int]]:
    tag_to_samples: Dict[str, List[int]] = {tag: [] for tag in OCM_ORDER}
    for sample_idx, sample in enumerate(samples):
        for tag in sample["ocm_ids"]:  # type: ignore[index]
            tag_to_samples[str(tag)].append(sample_idx)
    return tag_to_samples


def build_filtered_column_maps(
    barcodes_path: Path,
    samples: Sequence[Dict[str, object]],
    desired_barcodes_by_sample: Sequence[Sequence[str]],
) -> Tuple[List[array], Dict[str, Dict[str, int]]]:
    raw_barcodes = read_lines(barcodes_path)
    ncols = len(raw_barcodes)
    col_new_by_sample = [array("i", [0]) * (ncols + 1) for _sample in samples]
    desired_by_sample = [
        {normalize_barcode(barcode): idx + 1 for idx, barcode in enumerate(desired)}
        for desired in desired_barcodes_by_sample
    ]
    tag_to_samples = make_tag_to_sample_indices(samples)
    matched = {str(sample["sample_id"]): 0 for sample in samples}
    missing = {str(sample["sample_id"]): len(desired_barcodes_by_sample[i]) for i, sample in enumerate(samples)}

    for col, barcode in enumerate(raw_barcodes, start=1):
        tag = classify_composite_barcode(barcode)
        if tag is None:
            continue
        output_bc = normalize_barcode(output_cb16_from_composite(barcode))
        for sample_idx in tag_to_samples.get(tag, []):
            new_col = desired_by_sample[sample_idx].get(output_bc, 0)
            if new_col:
                col_new_by_sample[sample_idx][col] = new_col
                sample_id = str(samples[sample_idx]["sample_id"])
                matched[sample_id] += 1

    for sample in samples:
        sample_id = str(sample["sample_id"])
        missing[sample_id] -= matched[sample_id]
    return col_new_by_sample, {"matched_barcodes": matched, "missing_barcodes": missing}


def prepare_mex_output_dirs(
    output_dirs: Sequence[Path],
    features_lines: Sequence[str],
    barcodes_by_sample: Sequence[Sequence[str]],
) -> None:
    for out_dir, barcodes in zip(output_dirs, barcodes_by_sample):
        out_dir.mkdir(parents=True, exist_ok=True)
        write_features(out_dir, features_lines)
        write_lines_gz(out_dir / "barcodes.tsv.gz", barcodes)


def finalize_matrix_body(
    body_path: Path,
    output_path: Path,
    nrows: int,
    ncols: int,
    nnz: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wb", compresslevel=1) as out:
        out.write(b"%%MatrixMarket matrix coordinate integer general\n")
        out.write(b"%\n")
        out.write(f"{nrows} {ncols} {nnz}\n".encode("ascii"))
        with open(body_path, "rb") as body:
            shutil.copyfileobj(body, out, length=1024 * 1024)
    body_path.unlink(missing_ok=True)


def stream_matrix_layers_to_groups(
    input_dir: Path,
    groups: Sequence[Dict[str, object]],
    matrix_names: Sequence[str],
    temp_root: Path,
    threads: int,
    label: str,
    synthesize_matrix_from_layers: bool = False,
) -> Dict[str, object]:
    features_lines = read_lines(resolve_file(input_dir, "features.tsv"))
    for group in groups:
        prepare_mex_output_dirs(
            group["output_dirs"],  # type: ignore[arg-type]
            features_lines,
            group["barcodes_by_sample"],  # type: ignore[arg-type]
        )

    temp_root.mkdir(parents=True, exist_ok=True)
    worker_count = max(1, int(threads))
    stats: Dict[str, object] = {
        "input_dir": str(input_dir),
        "matrix_names": list(matrix_names),
        "synthetic_matrix": bool(synthesize_matrix_from_layers),
        "groups": {},
    }

    def finalize_many(tasks: Sequence[Tuple[Path, Path, int, int, int]]) -> None:
        if worker_count <= 1 or len(tasks) <= 1:
            for task in tasks:
                finalize_matrix_body(*task)
            return
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [executor.submit(finalize_matrix_body, *task) for task in tasks]
            for future in as_completed(futures):
                future.result()

    synthetic_handles: List[Dict[str, object]] = []
    if synthesize_matrix_from_layers and "matrix.mtx" not in matrix_names:
        for group_idx, group in enumerate(groups):
            group_name = str(group["name"])
            output_dirs = group["output_dirs"]  # type: ignore[assignment]
            barcodes_by_sample = group["barcodes_by_sample"]  # type: ignore[assignment]
            handles = []
            counts = [0 for _ in output_dirs]
            body_paths = []
            for sample_idx, _out_dir in enumerate(output_dirs):
                body_path = temp_root / f"{label}.{group_name}.matrix.sample{sample_idx}.body"
                body_paths.append(body_path)
                handles.append(open(body_path, "wb"))
            synthetic_handles.append(
                {
                    "group_idx": group_idx,
                    "group_name": group_name,
                    "handles": handles,
                    "counts": counts,
                    "body_paths": body_paths,
                    "nrows": None,
                    "ncols": [len(barcodes_by_sample[i]) for i in range(len(output_dirs))],
                }
            )

    try:
        for matrix_name in matrix_names:
            matrix_path = resolve_file(input_dir, matrix_name)
            nrows, ncols, _nnz = read_mtx_shape(matrix_path)
            layer_handles: List[Dict[str, object]] = []
            for group_idx, group in enumerate(groups):
                group_name = str(group["name"])
                output_dirs = group["output_dirs"]  # type: ignore[assignment]
                barcodes_by_sample = group["barcodes_by_sample"]  # type: ignore[assignment]
                col_new_by_sample = group["col_new_by_sample"]  # type: ignore[assignment]
                for col_map in col_new_by_sample:
                    if ncols + 1 != len(col_map):
                        die(
                            f"{matrix_path} columns ({ncols}) do not match barcode map "
                            f"for group {group_name} ({len(col_map) - 1})"
                        )
                handles = []
                counts = [0 for _ in output_dirs]
                body_paths = []
                for sample_idx, _out_dir in enumerate(output_dirs):
                    safe_matrix = matrix_name.replace("/", "_").replace(".", "_")
                    body_path = temp_root / f"{label}.{group_name}.{safe_matrix}.sample{sample_idx}.body"
                    body_paths.append(body_path)
                    handles.append(open(body_path, "wb"))
                layer_handles.append(
                    {
                        "group_idx": group_idx,
                        "group_name": group_name,
                        "handles": handles,
                        "counts": counts,
                        "body_paths": body_paths,
                        "nrows": nrows,
                        "ncols": [len(barcodes_by_sample[i]) for i in range(len(output_dirs))],
                        "col_new_by_sample": col_new_by_sample,
                        "output_dirs": output_dirs,
                    }
                )

            for synthetic in synthetic_handles:
                if synthetic["nrows"] is None:
                    synthetic["nrows"] = nrows
                elif int(synthetic["nrows"]) != nrows:
                    die(f"Velocity layer row mismatch while synthesizing matrix.mtx: {matrix_name}")

            for row, col, val in iter_mtx_entries(matrix_path):
                if col <= 0 or col > ncols:
                    continue
                for group_state in layer_handles:
                    col_new_by_sample = group_state["col_new_by_sample"]  # type: ignore[assignment]
                    handles = group_state["handles"]  # type: ignore[assignment]
                    counts = group_state["counts"]  # type: ignore[assignment]
                    synthetic = synthetic_handles[group_state["group_idx"]] if synthetic_handles else None  # type: ignore[index]
                    for sample_idx, col_map in enumerate(col_new_by_sample):
                        new_col = col_map[col]
                        if new_col <= 0:
                            continue
                        line = f"{row} {new_col} {val}\n".encode("ascii")
                        handles[sample_idx].write(line)
                        counts[sample_idx] += 1
                        if synthetic is not None:
                            synthetic["handles"][sample_idx].write(line)  # type: ignore[index]
                            synthetic["counts"][sample_idx] += 1  # type: ignore[index]

            finalize_tasks = []
            for group_state in layer_handles:
                handles = group_state["handles"]  # type: ignore[assignment]
                for handle in handles:
                    handle.close()
                output_dirs = group_state["output_dirs"]  # type: ignore[assignment]
                counts = group_state["counts"]  # type: ignore[assignment]
                body_paths = group_state["body_paths"]  # type: ignore[assignment]
                ncols_by_sample = group_state["ncols"]  # type: ignore[assignment]
                group_name = str(group_state["group_name"])
                group_stats = stats["groups"].setdefault(group_name, {"layer_nnz": {}})  # type: ignore[index]
                layer_stats = {}
                for sample_idx, out_dir in enumerate(output_dirs):
                    finalize_tasks.append(
                        (
                            body_paths[sample_idx],
                            out_dir / f"{matrix_name}.gz",
                            int(group_state["nrows"]),
                            int(ncols_by_sample[sample_idx]),
                            int(counts[sample_idx]),
                        )
                    )
                    layer_stats[sample_idx] = int(counts[sample_idx])
                group_stats["layer_nnz"][matrix_name] = layer_stats  # type: ignore[index]
            finalize_many(finalize_tasks)

        synthetic_tasks = []
        for synthetic in synthetic_handles:
            handles = synthetic["handles"]  # type: ignore[assignment]
            for handle in handles:
                handle.close()
            group_idx = int(synthetic["group_idx"])
            group_name = str(synthetic["group_name"])
            output_dirs = groups[group_idx]["output_dirs"]  # type: ignore[assignment]
            body_paths = synthetic["body_paths"]  # type: ignore[assignment]
            counts = synthetic["counts"]  # type: ignore[assignment]
            ncols_by_sample = synthetic["ncols"]  # type: ignore[assignment]
            group_stats = stats["groups"].setdefault(group_name, {"layer_nnz": {}})  # type: ignore[index]
            matrix_stats = {}
            for sample_idx, out_dir in enumerate(output_dirs):
                synthetic_tasks.append(
                    (
                        body_paths[sample_idx],
                        out_dir / "matrix.mtx.gz",
                        int(synthetic["nrows"] or 0),
                        int(ncols_by_sample[sample_idx]),
                        int(counts[sample_idx]),
                    )
                )
                matrix_stats[sample_idx] = int(counts[sample_idx])
            group_stats["layer_nnz"]["matrix.mtx"] = matrix_stats  # type: ignore[index]
        finalize_many(synthetic_tasks)
    finally:
        for synthetic in synthetic_handles:
            for handle in synthetic.get("handles", []):
                try:
                    if not handle.closed:
                        handle.close()
                except Exception:
                    pass

    return stats


def split_matrix_layers_by_maps_streaming(
    input_dir: Path,
    output_dirs: Sequence[Path],
    samples: Sequence[Dict[str, object]],
    matrix_names: Sequence[str],
    temp_root: Path,
    threads: int,
    label: str,
) -> Tuple[Dict[str, object], List[array], List[List[str]], Dict[str, int]]:
    barcodes_path = resolve_file(input_dir, "barcodes.tsv")
    _col_tag, col_new_by_sample, sample_barcodes, tag_counts = build_column_maps(barcodes_path, samples)
    stats = stream_matrix_layers_to_groups(
        input_dir,
        [
            {
                "name": "raw",
                "output_dirs": list(output_dirs),
                "barcodes_by_sample": sample_barcodes,
                "col_new_by_sample": col_new_by_sample,
            }
        ],
        matrix_names,
        temp_root,
        threads,
        label,
    )
    stats["tag_counts"] = tag_counts
    stats["sample_columns"] = {
        str(samples[i]["sample_id"]): len(sample_barcodes[i]) for i in range(len(samples))
    }
    return stats, col_new_by_sample, sample_barcodes, tag_counts


def link_tree(src_dir: Path, dst_dir: Path, force: bool = False) -> None:
    if force and dst_dir.exists():
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    for src in src_dir.iterdir():
        if not src.is_file():
            continue
        dst = dst_dir / src.name
        hardlink_or_copy(src, dst)


def write_stripped_pool_mex(input_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    features = read_lines(resolve_file(input_dir, "features.tsv"))
    barcodes = [output_cb16_from_composite(line) for line in read_lines(resolve_file(input_dir, "barcodes.tsv"))]
    write_lines_gz(out_dir / "features.tsv.gz", features)
    write_lines_gz(out_dir / "barcodes.tsv.gz", barcodes)
    matrix_path = resolve_file(input_dir, "matrix.mtx")
    with open_text(matrix_path, "rt") as src, gzip.open(
        out_dir / "matrix.mtx.gz", "wt", encoding="utf-8", compresslevel=1
    ) as dst:
        shutil.copyfileobj(src, dst)


def run_simpleed(
    repo_root: Path,
    raw_mex: Path,
    filtered_barcodes: Path,
    ed_out_dir: Path,
    mode: str,
    sim_n: int,
) -> None:
    cmd = [
        "bash",
        str(repo_root / "scripts" / "run_simpleed_fallback.sh"),
        "--raw-mex",
        str(raw_mex),
        "--filtered-barcodes",
        str(filtered_barcodes),
        "--force",
        "--mode",
        mode,
        "--use-bootstrap",
        "--use-fdr-gate",
        "--apply-bh-correction",
        "--out-dir",
        str(ed_out_dir),
    ]
    if sim_n > 0:
        cmd.extend(["--sim-n", str(sim_n)])
    log("running " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_simpleed_binary(repo_root: Path) -> None:
    simpleed_bin = repo_root / "core" / "features" / "libscrna" / "bin" / "scrna_simpleed"
    if simpleed_bin.exists() and os.access(simpleed_bin, os.X_OK):
        return
    log(f"building scrna_simpleed once before parallel ED workers: {simpleed_bin}")
    subprocess.run(["make", "-C", str(repo_root / "core" / "features" / "libscrna"), "tools"], check=True)


def subset_mex_by_barcodes(input_dir: Path, desired_barcodes: Sequence[str], output_dir: Path) -> Dict[str, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    features = read_lines(resolve_file(input_dir, "features.tsv"))
    raw_barcodes = read_lines(resolve_file(input_dir, "barcodes.tsv"))
    desired_by_norm = {normalize_barcode(bc): i + 1 for i, bc in enumerate(desired_barcodes)}
    col_new = array("i", [0]) * (len(raw_barcodes) + 1)
    matched = 0
    for col, barcode in enumerate(raw_barcodes, start=1):
        new_col = desired_by_norm.get(normalize_barcode(barcode), 0)
        if new_col:
            col_new[col] = new_col
            matched += 1

    matrix_names = ["matrix.mtx"]
    for layer in ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx"):
        if (input_dir / layer).exists() or (input_dir / f"{layer}.gz").exists():
            matrix_names.append(layer)

    write_lines_gz(output_dir / "features.tsv.gz", features)
    write_lines_gz(output_dir / "barcodes.tsv.gz", desired_barcodes)
    stats = {"desired_barcodes": len(desired_barcodes), "matched_barcodes": matched}

    for matrix_name in matrix_names:
        matrix_path = resolve_file(input_dir, matrix_name)
        nrows, ncols, _nnz = read_mtx_shape(matrix_path)
        if ncols + 1 != len(col_new):
            die(f"{matrix_path} columns ({ncols}) do not match barcode count ({len(col_new) - 1})")
        nnz = 0
        for _row, col, _val in iter_mtx_entries(matrix_path):
            if 0 < col < len(col_new) and col_new[col] > 0:
                nnz += 1
        with gzip.open(output_dir / f"{matrix_name}.gz", "wt", encoding="utf-8", compresslevel=1) as out:
            out.write("%%MatrixMarket matrix coordinate integer general\n")
            out.write("%\n")
            out.write(f"{nrows} {len(desired_barcodes)} {nnz}\n")
            for row, col, val in iter_mtx_entries(matrix_path):
                if 0 < col < len(col_new) and col_new[col] > 0:
                    out.write(f"{row} {col_new[col]} {val}\n")
        stats[f"{matrix_name}_nnz"] = nnz
    return stats


def materialize(args: argparse.Namespace) -> None:
    repo_root = Path(args.repo_root).resolve()
    star_run_dir = Path(args.star_run_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    samples = parse_ocm_config(Path(args.config))
    threads = max(1, int(getattr(args, "threads", 1)))
    if args.force and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_root = out_dir / ".tmp_ocm_materialize"
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_root.mkdir(parents=True, exist_ok=True)
    timings: Dict[str, Dict[str, float]] = {}

    raw_in = star_run_dir / "Solo.out" / args.feature / "raw"
    if not raw_in.exists():
        die(f"Missing STAR raw MEX for {args.feature}: {raw_in}")

    log(f"materializing {args.feature} raw MEX by OCM with streaming router")
    primary_raw_dirs = [
        out_dir / "outs" / "per_sample_outs" / str(sample["sample_id"]) / "count" / "sample_raw_feature_bc_matrix"
        for sample in samples
    ]
    started = now_s()
    raw_stats, _raw_col_maps, _raw_sample_barcodes, _tag_counts = split_matrix_layers_by_maps_streaming(
        raw_in,
        primary_raw_dirs,
        samples,
        ["matrix.mtx"],
        temp_root,
        threads,
        f"{args.feature}.raw",
    )
    log_timing(timings, f"{args.feature}.raw_split", started)

    started = now_s()
    write_stripped_pool_mex(raw_in, out_dir / "outs" / "multi" / "count" / "raw_feature_bc_matrix")
    for sample, raw_dir in zip(samples, primary_raw_dirs):
        downstream_raw = out_dir / "samples" / str(sample["sample_id"]) / "run" / "outs" / "raw_feature_bc_matrix"
        link_tree(raw_dir, downstream_raw, force=True)
    log_timing(timings, f"{args.feature}.raw_mirrors", started)

    cells_per_tag: Dict[str, List[str]] = {tag: [] for tag in OCM_ORDER}
    filtered_stats: Dict[str, object] = {}

    def run_one_ed(sample_and_dir: Tuple[Dict[str, object], Path]) -> Tuple[str, List[str]]:
        sample, raw_dir = sample_and_dir
        sample_id = str(sample["sample_id"])
        ed_dir = out_dir / "ed" / sample_id
        filtered_barcodes = ed_dir / "filtered_barcodes.tsv"
        run_simpleed(repo_root, raw_dir, filtered_barcodes, ed_dir, args.simpleed_mode, int(args.sim_n))
        return sample_id, read_lines(filtered_barcodes)

    started = now_s()
    ensure_simpleed_binary(repo_root)
    desired_by_sample_id: Dict[str, List[str]] = {}
    ed_workers = min(threads, max(1, len(samples)))
    if ed_workers <= 1:
        for pair in zip(samples, primary_raw_dirs):
            sample_id, desired = run_one_ed(pair)
            desired_by_sample_id[sample_id] = desired
    else:
        with ThreadPoolExecutor(max_workers=ed_workers) as executor:
            futures = [executor.submit(run_one_ed, pair) for pair in zip(samples, primary_raw_dirs)]
            for future in as_completed(futures):
                sample_id, desired = future.result()
                desired_by_sample_id[sample_id] = desired
    log_timing(timings, f"{args.feature}.emptydrops", started)

    desired_barcodes_by_sample = [desired_by_sample_id[str(sample["sample_id"])] for sample in samples]
    for sample, desired in zip(samples, desired_barcodes_by_sample):
        for tag in sample["ocm_ids"]:  # type: ignore[index]
            cells_per_tag[str(tag)].extend(desired)

    started = now_s()
    filtered_col_maps, filtered_match_stats = build_filtered_column_maps(
        resolve_file(raw_in, "barcodes.tsv"),
        samples,
        desired_barcodes_by_sample,
    )
    primary_filtered_dirs = [
        out_dir / "outs" / "per_sample_outs" / str(sample["sample_id"]) / "count" / "sample_filtered_feature_bc_matrix"
        for sample in samples
    ]
    filtered_stream_stats = stream_matrix_layers_to_groups(
        raw_in,
        [
            {
                "name": "filtered",
                "output_dirs": primary_filtered_dirs,
                "barcodes_by_sample": desired_barcodes_by_sample,
                "col_new_by_sample": filtered_col_maps,
            }
        ],
        ["matrix.mtx"],
        temp_root,
        threads,
        f"{args.feature}.filtered",
    )
    log_timing(timings, f"{args.feature}.filtered_split", started)

    for sample, desired, primary_filtered in zip(samples, desired_barcodes_by_sample, primary_filtered_dirs):
        sample_id = str(sample["sample_id"])
        filtered_stats[sample_id] = {
            "desired_barcodes": len(desired),
            "matched_barcodes": filtered_match_stats["matched_barcodes"][sample_id],
            "missing_barcodes": filtered_match_stats["missing_barcodes"][sample_id],
        }
        primary_filtered = (
            out_dir
            / "outs"
            / "per_sample_outs"
            / sample_id
            / "count"
            / "sample_filtered_feature_bc_matrix"
        )
        downstream_filtered = out_dir / "samples" / sample_id / "run" / "outs" / "filtered_feature_bc_matrix"
        link_tree(primary_filtered, downstream_filtered, force=True)
        csv_path = out_dir / "outs" / "per_sample_outs" / sample_id / "count" / "sample_filtered_barcodes.csv"
        with open(csv_path, "w", encoding="utf-8") as handle:
            for barcode in desired:
                handle.write(f"GRCh38,{barcode}\n")

    multi_mux = out_dir / "outs" / "multi" / "multiplexing_analysis"
    multi_mux.mkdir(parents=True, exist_ok=True)
    (multi_mux / "cells_per_tag.json").write_text(
        json.dumps(cells_per_tag, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for sample in samples:
        sample_id = str(sample["sample_id"])
        mux_dir = out_dir / "samples" / sample_id / "run" / "outs" / "multiplexing_analysis"
        mux_dir.mkdir(parents=True, exist_ok=True)
        (mux_dir / "cells_per_tag.json").write_text(
            json.dumps(cells_per_tag, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    velo_stats = {}
    if args.include_velocyto:
        velo_in = star_run_dir / "Solo.out" / "Velocyto" / "raw"
        if velo_in.exists():
            matrix_names = [
                layer
                for layer in ("matrix.mtx", "spliced.mtx", "unspliced.mtx", "ambiguous.mtx")
                if (velo_in / layer).exists() or (velo_in / f"{layer}.gz").exists()
            ]
            if matrix_names:
                started = now_s()
                log("materializing Velocyto raw and filtered MEX by OCM with one stream per layer")
                velo_raw_dirs = [
                    out_dir / "samples" / str(sample["sample_id"]) / "run" / "outs" / "raw_velocyto_feature_bc_matrix"
                    for sample in samples
                ]
                velo_filtered_dirs = [
                    out_dir
                    / "samples"
                    / str(sample["sample_id"])
                    / "run"
                    / "outs"
                    / "filtered_velocyto_feature_bc_matrix"
                    for sample in samples
                ]
                _velo_col_tag, velo_raw_col_maps, velo_sample_barcodes, velo_tag_counts = build_column_maps(
                    resolve_file(velo_in, "barcodes.tsv"),
                    samples,
                )
                velo_filtered_col_maps, velo_filtered_match_stats = build_filtered_column_maps(
                    resolve_file(velo_in, "barcodes.tsv"),
                    samples,
                    desired_barcodes_by_sample,
                )
                velo_stats = stream_matrix_layers_to_groups(
                    velo_in,
                    [
                        {
                            "name": "raw",
                            "output_dirs": velo_raw_dirs,
                            "barcodes_by_sample": velo_sample_barcodes,
                            "col_new_by_sample": velo_raw_col_maps,
                        },
                        {
                            "name": "filtered",
                            "output_dirs": velo_filtered_dirs,
                            "barcodes_by_sample": desired_barcodes_by_sample,
                            "col_new_by_sample": velo_filtered_col_maps,
                        },
                    ],
                    matrix_names,
                    temp_root,
                    threads,
                    "Velocyto",
                    synthesize_matrix_from_layers=("matrix.mtx" not in matrix_names),
                )
                velo_stats["tag_counts"] = velo_tag_counts
                velo_stats["filtered_match_stats"] = velo_filtered_match_stats
                log_timing(timings, "Velocyto.raw_filtered_split", started)
            else:
                log(f"Velocyto raw dir has no matrix layers: {velo_in}")
        else:
            log(f"Velocyto raw dir not present; skipping: {velo_in}")

    summary = {
        "feature": args.feature,
        "star_run_dir": str(star_run_dir),
        "out_dir": str(out_dir),
        "samples": samples,
        "raw": raw_stats,
        "filtered": filtered_stats,
        "filtered_stream": filtered_stream_stats,
        "velocyto": velo_stats,
        "cells_per_tag_counts": {tag: len(values) for tag, values in cells_per_tag.items()},
        "timings": timings,
        "threads": threads,
    }
    (out_dir / "materialization_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    shutil.rmtree(temp_root, ignore_errors=True)
    log(f"materialization complete: {out_dir}")


def load_mex_summary(mex_dir: Path) -> Dict[str, object]:
    barcodes = read_lines(resolve_file(mex_dir, "barcodes.tsv"))
    features = [line.split("\t", 1)[0] for line in read_lines(resolve_file(mex_dir, "features.tsv"))]
    matrix = resolve_file(mex_dir, "matrix.mtx")
    bc_totals: Dict[str, int] = {normalize_barcode(bc): 0 for bc in barcodes}
    feature_totals: Dict[str, int] = {feature_id: 0 for feature_id in features}
    total = 0
    for row, col, val_s in iter_mtx_entries(matrix):
        val = int(float(val_s))
        total += val
        if 1 <= col <= len(barcodes):
            bc_totals[normalize_barcode(barcodes[col - 1])] = bc_totals.get(normalize_barcode(barcodes[col - 1]), 0) + val
        if 1 <= row <= len(features):
            feature_totals[features[row - 1]] = feature_totals.get(features[row - 1], 0) + val
    return {
        "barcodes": set(bc_totals.keys()),
        "bc_totals": bc_totals,
        "feature_totals": feature_totals,
        "total_umis": total,
        "n_features": len(features),
    }


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    if not xs or len(xs) != len(ys):
        return float("nan")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0.0 or den_y == 0.0:
        return float("nan")
    return num / (den_x * den_y)


def fmt_float(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.6f}"


def find_cr_sample_mex(cr_run_dir: Path, sample_id: str) -> Path:
    candidates = [
        cr_run_dir / "outs" / "per_sample_outs" / sample_id / "count" / "sample_filtered_feature_bc_matrix",
        cr_run_dir / "outs" / "per_sample_outs" / sample_id / "count" / "filtered_feature_bc_matrix",
        cr_run_dir / "outs" / "per_sample_outs" / sample_id / "count" / "sample_feature_bc_matrix",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    die(f"Could not locate CR per-sample filtered MEX for {sample_id} under {cr_run_dir}")


def compare(args: argparse.Namespace) -> None:
    samples = parse_ocm_config(Path(args.config))
    star_root = Path(args.star_materialized_dir)
    cr_run_dir = Path(args.cr_run_dir)
    rows = []
    for sample in samples:
        sample_id = str(sample["sample_id"])
        star_mex = star_root / "samples" / sample_id / "run" / "outs" / "filtered_feature_bc_matrix"
        cr_mex = find_cr_sample_mex(cr_run_dir, sample_id)
        if not star_mex.exists():
            die(f"Missing STAR materialized MEX for {sample_id}: {star_mex}")
        log(f"comparing {sample_id}")
        star = load_mex_summary(star_mex)
        cr = load_mex_summary(cr_mex)
        star_bcs = star["barcodes"]  # type: ignore[assignment]
        cr_bcs = cr["barcodes"]  # type: ignore[assignment]
        common_bcs = sorted(star_bcs & cr_bcs)  # type: ignore[operator]
        union_bcs = star_bcs | cr_bcs  # type: ignore[operator]
        star_bc_totals = star["bc_totals"]  # type: ignore[assignment]
        cr_bc_totals = cr["bc_totals"]  # type: ignore[assignment]
        common_features = sorted(
            set(star["feature_totals"].keys()) & set(cr["feature_totals"].keys())  # type: ignore[union-attr]
        )
        star_feature_totals = star["feature_totals"]  # type: ignore[assignment]
        cr_feature_totals = cr["feature_totals"]  # type: ignore[assignment]
        row = {
            "sample_id": sample_id,
            "star_cells": len(star_bcs),  # type: ignore[arg-type]
            "cr_cells": len(cr_bcs),  # type: ignore[arg-type]
            "common_cells": len(common_bcs),
            "barcode_jaccard": (len(common_bcs) / len(union_bcs)) if union_bcs else 1.0,
            "star_total_umis": star["total_umis"],
            "cr_total_umis": cr["total_umis"],
            "common_feature_count": len(common_features),
            "barcode_umi_pearson": pearson(
                [float(star_bc_totals[bc]) for bc in common_bcs],
                [float(cr_bc_totals[bc]) for bc in common_bcs],
            ),
            "feature_umi_pearson": pearson(
                [float(star_feature_totals[f]) for f in common_features],
                [float(cr_feature_totals[f]) for f in common_features],
            ),
            "star_mex": str(star_mex),
            "cr_mex": str(cr_mex),
        }
        rows.append(row)

    out_json = Path(args.out_json)
    out_tsv = Path(args.out_tsv)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with open(out_tsv, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "sample_id",
            "star_cells",
            "cr_cells",
            "common_cells",
            "barcode_jaccard",
            "star_total_umis",
            "cr_total_umis",
            "common_feature_count",
            "barcode_umi_pearson",
            "feature_umi_pearson",
            "star_mex",
            "cr_mex",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            formatted["barcode_jaccard"] = fmt_float(float(row["barcode_jaccard"]))
            formatted["barcode_umi_pearson"] = fmt_float(float(row["barcode_umi_pearson"]))
            formatted["feature_umi_pearson"] = fmt_float(float(row["feature_umi_pearson"]))
            writer.writerow(formatted)
    log(f"comparison complete: {out_tsv}")


def extract_sam_tag(fields: Sequence[str], tag: str) -> Optional[str]:
    prefix = f"{tag}:Z:"
    for field in fields[11:]:
        if field.startswith(prefix):
            return field[len(prefix) :]
    return None


def normalize_cb_tag_value(cb: str) -> str:
    return output_cb16_from_composite(cb)


def normalized_record_fields(fields: List[str], ocm_id: str) -> List[str]:
    out = fields[:]
    saw_cb = False
    saw_zt = False
    for idx in range(11, len(out)):
        if out[idx].startswith("CB:Z:"):
            out[idx] = "CB:Z:" + normalize_cb_tag_value(out[idx][5:])
            saw_cb = True
        elif out[idx].startswith("ZT:Z:"):
            out[idx] = "ZT:Z:" + ocm_id
            saw_zt = True
    if not saw_cb:
        raw = extract_sam_tag(out, "CR")
        if raw:
            out.append("CB:Z:" + normalize_cb_tag_value(raw))
    if not saw_zt:
        out.append("ZT:Z:" + ocm_id)
    return out


def split_bam(args: argparse.Namespace) -> None:
    bam_path = Path(args.bam).resolve()
    out_dir = Path(args.out_dir).resolve()
    samples = parse_ocm_config(Path(args.config))
    if not bam_path.exists():
        die(f"Missing input BAM: {bam_path}")
    if args.force and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pooled_dir = out_dir / "pooled"
    per_sample_dir = out_dir / "per_sample"
    pooled_dir.mkdir(parents=True, exist_ok=True)
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    hardlink_or_copy(bam_path, pooled_dir / "pooled_composite_cb.bam")

    tag_to_samples: Dict[str, List[str]] = {tag: [] for tag in OCM_ORDER}
    for sample in samples:
        sample_id = str(sample["sample_id"])
        for tag in sample["ocm_ids"]:  # type: ignore[index]
            tag_to_samples[str(tag)].append(sample_id)

    output_paths = {
        sample_id: per_sample_dir / f"{sample_id}.bam"
        for sample_id in (str(sample["sample_id"]) for sample in samples)
    }
    output_paths["__pooled_cb16__"] = pooled_dir / "pooled_cb16_with_ocm_tag.bam"
    output_paths["__unassigned__"] = out_dir / "unassigned.bam"

    writers = {
        key: subprocess.Popen(
            ["samtools", "view", "-b", "-o", str(path), "-"],
            stdin=subprocess.PIPE,
            text=True,
        )
        for key, path in output_paths.items()
    }

    counts = {
        "input_records": 0,
        "assigned_records": 0,
        "unassigned_records": 0,
        "records_by_ocm": {tag: 0 for tag in OCM_ORDER},
        "records_by_sample": {str(sample["sample_id"]): 0 for sample in samples},
    }

    reader = subprocess.Popen(
        ["samtools", "view", "-h", str(bam_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert reader.stdout is not None
    try:
        for line in reader.stdout:
            if line.startswith("@"):
                for proc in writers.values():
                    assert proc.stdin is not None
                    proc.stdin.write(line)
                continue
            counts["input_records"] += 1
            fields = line.rstrip("\n").split("\t")
            cb = extract_sam_tag(fields, "CB")
            ocm_id = classify_composite_barcode(cb) if cb else None
            if ocm_id is None:
                assert writers["__unassigned__"].stdin is not None
                writers["__unassigned__"].stdin.write(line)
                counts["unassigned_records"] += 1
                continue
            norm_fields = normalized_record_fields(fields, ocm_id)
            norm_line = "\t".join(norm_fields) + "\n"
            assert writers["__pooled_cb16__"].stdin is not None
            writers["__pooled_cb16__"].stdin.write(norm_line)
            counts["assigned_records"] += 1
            counts["records_by_ocm"][ocm_id] += 1  # type: ignore[index]
            for sample_id in tag_to_samples.get(ocm_id, []):
                assert writers[sample_id].stdin is not None
                writers[sample_id].stdin.write(norm_line)
                counts["records_by_sample"][sample_id] += 1  # type: ignore[index]
    finally:
        if reader.stdout:
            reader.stdout.close()
        stderr = reader.stderr.read() if reader.stderr else ""
        rc = reader.wait()
        for proc in writers.values():
            if proc.stdin:
                proc.stdin.close()
        writer_failures = {}
        for key, proc in writers.items():
            writer_rc = proc.wait()
            if writer_rc != 0:
                writer_failures[key] = writer_rc
        if rc != 0:
            die(f"samtools view failed for {bam_path}: rc={rc} {stderr.strip()}")
        if writer_failures:
            die(f"samtools BAM writers failed: {writer_failures}")

    summary = {
        "input_bam": str(bam_path),
        "out_dir": str(out_dir),
        "pooled_composite_bam": str(pooled_dir / "pooled_composite_cb.bam"),
        "pooled_cb16_bam": str(pooled_dir / "pooled_cb16_with_ocm_tag.bam"),
        "per_sample_bams": {key: str(path) for key, path in output_paths.items() if not key.startswith("__")},
        "unassigned_bam": str(output_paths["__unassigned__"]),
        **counts,
    }
    (out_dir / "split_bam_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    log(f"BAM split complete: {out_dir}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    prep = sub.add_parser("prepare", help="downsample FASTQs and build composite R1/whitelist")
    prep.add_argument("--raw-dir", required=True)
    prep.add_argument("--sample-id", required=True)
    prep.add_argument("--sample-stem", required=True)
    prep.add_argument("--lanes", default="L007,L008")
    prep.add_argument("--read-pairs", required=True, type=int)
    prep.add_argument("--cr-fastq-dir", required=True)
    prep.add_argument("--star-fastq-dir", required=True)
    prep.add_argument(
        "--star-mode",
        choices=("legacy17", "native"),
        default="legacy17",
        help="legacy17 rewrites STAR R1 as CB16+1bp; native hardlinks original FASTQs for STAR --ocmMultiBarcodeMode flex",
    )
    prep.add_argument("--whitelist", required=True)
    prep.add_argument("--composite-whitelist", required=True)
    prep.add_argument("--manifest", required=True)
    prep.add_argument("--stats-json", required=True)
    prep.add_argument("--force", action="store_true")
    prep.set_defaults(func=prepare)

    mat = sub.add_parser("materialize", help="split composite STAR MEX into per-OCM sample outputs")
    mat.add_argument("--repo-root", required=True)
    mat.add_argument("--star-run-dir", required=True)
    mat.add_argument("--config", required=True)
    mat.add_argument("--feature", default="GeneFull")
    mat.add_argument("--out-dir", required=True)
    mat.add_argument("--simpleed-mode", default="full", choices=["simple", "full"])
    mat.add_argument("--sim-n", default=100000, type=int)
    mat.add_argument("--include-velocyto", action="store_true")
    mat.add_argument("--threads", default=1, type=int, help="worker threads for ED and MEX finalization")
    mat.add_argument("--force", action="store_true")
    mat.set_defaults(func=materialize)

    prod = sub.add_parser(
        "materialize-production",
        help="optimized production materialization: GeneFull plus Velocyto, no Gene comparator",
    )
    prod.add_argument("--repo-root", required=True)
    prod.add_argument("--star-run-dir", required=True)
    prod.add_argument("--config", required=True)
    prod.add_argument("--out-dir", required=True)
    prod.add_argument("--simpleed-mode", default="full", choices=["simple", "full"])
    prod.add_argument("--sim-n", default=100000, type=int)
    prod.add_argument("--threads", default=1, type=int, help="worker threads for ED and MEX finalization")
    prod.add_argument("--force", action="store_true")
    prod.set_defaults(func=materialize, feature="GeneFull", include_velocyto=True)

    cmp_p = sub.add_parser("compare", help="compare materialized STAR per-sample MEX with CR multi")
    cmp_p.add_argument("--star-materialized-dir", required=True)
    cmp_p.add_argument("--cr-run-dir", required=True)
    cmp_p.add_argument("--config", required=True)
    cmp_p.add_argument("--out-json", required=True)
    cmp_p.add_argument("--out-tsv", required=True)
    cmp_p.set_defaults(func=compare)

    bam_p = sub.add_parser("split-bam", help="split pooled composite STAR BAM into per-OCM BAMs")
    bam_p.add_argument("--bam", required=True)
    bam_p.add_argument("--config", required=True)
    bam_p.add_argument("--out-dir", required=True)
    bam_p.add_argument("--force", action="store_true")
    bam_p.set_defaults(func=split_bam)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
