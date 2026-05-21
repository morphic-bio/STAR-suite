#!/usr/bin/env python3
"""FASTQ preflight: detect the 10x whitelist family/namespace from R1 barcodes.

This is a cheap guard before STAR/pf-multi production runs. It samples cell
barcodes from FASTQ R1, scores exact membership against one or more whitelist
families, and reports whether the observed barcode space matches the expected
manifest whitelist.

Typical use:

  scripts/preflight_whitelist_family.py \
    --manifest docs/MSK_30KO_FASTQ_MANIFEST.tsv \
    --whitelist feb2018:TRU:/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt \
    --whitelist feb2018:NXT:/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt \
    --whitelist may2023_gemx:TRU:/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt \
    --whitelist may2023_gemx:NXT:/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_NXT.txt \
    --outdir /tmp/msk_whitelist_preflight
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

log = logging.getLogger("preflight_whitelist_family")

DEFAULT_SAMPLE_READS = 200_000
DEFAULT_BARCODE_LEN = 16
DEFAULT_MIN_HIT_RATE = 0.05
DEFAULT_DOMINANCE_RATIO = 3.0
DEFAULT_MAX_FASTQS_PER_MANIFEST_ROW = 4


@dataclass(frozen=True)
class WhitelistSpec:
    family: str
    chemistry: str
    path: str
    column: int = 1

    @property
    def key(self) -> str:
        return f"{self.family}:{self.chemistry}"


@dataclass
class FastqCheck:
    label: str
    fastq_paths: List[str]
    expected_chemistry: str = ""
    expected_whitelist: str = ""
    expected_family: str = ""
    expected_key: str = ""
    manifest_row: Dict[str, str] = field(default_factory=dict)


@dataclass
class WhitelistLoadReport:
    family: str
    chemistry: str
    path: str
    column: int
    barcode_count: int


@dataclass
class FastqReport:
    label: str
    fastq_paths: List[str]
    reads_sampled: int
    barcode_len: int
    expected_chemistry: str
    expected_whitelist: str
    expected_family: str
    expected_key: str
    best_key: str
    best_family: str
    best_chemistry: str
    best_hits: int
    best_rate: float
    second_key: str
    second_rate: float
    sn_ratio: float
    status: str
    rates: Dict[str, float]
    hits: Dict[str, int]


def open_text(path: Path):
    """Open plain or gzip-compressed text by magic bytes."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_whitelist_spec(raw: str) -> WhitelistSpec:
    """Parse FAMILY:CHEMISTRY:PATH[:COLUMN]."""
    parts = raw.split(":")
    if len(parts) < 3:
        raise ValueError(f"Expected FAMILY:CHEMISTRY:PATH[:COLUMN], got {raw!r}")

    family = parts[0].strip()
    chemistry = parts[1].strip().upper()
    if chemistry not in {"TRU", "NXT"}:
        raise ValueError(f"Chemistry must be TRU or NXT in {raw!r}")

    column = 1
    path_parts = parts[2:]
    if len(path_parts) > 1 and path_parts[-1].isdigit():
        column = int(path_parts[-1])
        path_parts = path_parts[:-1]
    path = ":".join(path_parts)

    if not family:
        raise ValueError(f"Family is empty in {raw!r}")
    if not path:
        raise ValueError(f"Path is empty in {raw!r}")
    if column < 1:
        raise ValueError(f"Column must be 1-based in {raw!r}")
    return WhitelistSpec(family=family, chemistry=chemistry, path=path, column=column)


def load_whitelist(spec: WhitelistSpec, barcode_len: int) -> Tuple[Set[str], WhitelistLoadReport]:
    barcodes: Set[str] = set()
    path = Path(spec.path)
    with open_text(path) as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < spec.column:
                continue
            bc = parts[spec.column - 1].upper()
            if len(bc) == barcode_len and _is_acgt(bc):
                barcodes.add(bc)

    if not barcodes:
        raise ValueError(
            f"Whitelist {spec.key} loaded zero {barcode_len}bp ACGT barcodes from {spec.path}"
        )

    return barcodes, WhitelistLoadReport(
        family=spec.family,
        chemistry=spec.chemistry,
        path=spec.path,
        column=spec.column,
        barcode_count=len(barcodes),
    )


def _is_acgt(seq: str) -> bool:
    return all(base in "ACGT" for base in seq)


def sample_barcodes(paths: Sequence[Path], sample_reads: int, barcode_len: int) -> List[str]:
    """Sample up to sample_reads total R1 barcodes across paths."""
    if not paths:
        return []
    per_file = max(1, math.ceil(sample_reads / len(paths)))
    out: List[str] = []
    for path in paths:
        with open_text(path) as handle:
            sampled_here = 0
            while sampled_here < per_file and len(out) < sample_reads:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().strip().upper()
                handle.readline()
                handle.readline()
                if len(seq) >= barcode_len:
                    out.append(seq[:barcode_len])
                    sampled_here += 1
    return out


def score_barcodes(
    barcodes: Sequence[str],
    whitelist_sets: Dict[str, Set[str]],
) -> Tuple[Dict[str, int], Dict[str, float]]:
    hits = {key: 0 for key in whitelist_sets}
    scored = 0
    for bc in barcodes:
        if len(bc) == 0:
            continue
        scored += 1
        for key, wl in whitelist_sets.items():
            if bc in wl:
                hits[key] += 1
    rates = {key: (count / scored if scored else 0.0) for key, count in hits.items()}
    return hits, rates


def call_status(
    expected_key: str,
    best_key: str,
    best_rate: float,
    second_rate: float,
    min_hit_rate: float,
    dominance_ratio: float,
) -> str:
    if not best_key or best_rate < min_hit_rate:
        return "FAIL_LOW_HIT_RATE"

    sn_ratio = best_rate / second_rate if second_rate > 0 else float("inf")
    low_conf = sn_ratio < dominance_ratio

    if expected_key:
        if best_key != expected_key:
            return "FAIL_EXPECTED_MISMATCH"
        return "WARN_LOW_CONFIDENCE" if low_conf else "PASS"

    return "WARN_LOW_CONFIDENCE" if low_conf else "INFO_BEST_CALL"


def build_fastq_check_from_arg(raw: str) -> FastqCheck:
    """Parse LABEL:PATH[,PATH...] for ad-hoc checks."""
    parts = raw.split(":", 1)
    if len(parts) != 2:
        raise ValueError(f"Expected LABEL:PATH[,PATH...], got {raw!r}")
    label, paths = parts
    fastq_paths = [p for p in paths.split(",") if p]
    if not label or not fastq_paths:
        raise ValueError(f"Expected LABEL:PATH[,PATH...], got {raw!r}")
    return FastqCheck(label=label, fastq_paths=fastq_paths)


def load_manifest_checks(
    manifest: Path,
    max_fastqs_per_row: int,
    whitelist_path_to_key: Dict[str, str],
    whitelist_path_to_family: Dict[str, str],
) -> List[FastqCheck]:
    checks: List[FastqCheck] = []
    with open(manifest, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "provider_group",
            "library",
            "chemistry",
            "whitelist",
            "fastq_root",
            "fastq_sample_ids",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest {manifest} is missing required columns: {sorted(missing)}")

        for row in reader:
            root = Path(row["fastq_root"])
            sample_ids = [sid for sid in row["fastq_sample_ids"].split(";") if sid]
            paths: List[Path] = []
            for sample_id in sample_ids:
                matches = sorted(root.glob(f"{sample_id}_*_R1_001.fastq*"))
                paths.extend(matches)
            paths = sorted(paths)[:max_fastqs_per_row]

            label = f"{row['provider_group']}:{row['library']}"
            expected_whitelist = row.get("whitelist", "")
            expected_key = whitelist_path_to_key.get(expected_whitelist, "")
            expected_family = whitelist_path_to_family.get(expected_whitelist, "")
            checks.append(FastqCheck(
                label=label,
                fastq_paths=[str(path) for path in paths],
                expected_chemistry=row.get("chemistry", ""),
                expected_whitelist=expected_whitelist,
                expected_family=expected_family,
                expected_key=expected_key,
                manifest_row=dict(row),
            ))
    return checks


def write_reports(
    reports: List[FastqReport],
    whitelist_reports: List[WhitelistLoadReport],
    outdir: Path,
    parameters: Dict[str, object],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "whitelist_family_summary.tsv"
    rates_path = outdir / "whitelist_family_rates.tsv"
    json_path = outdir / "whitelist_family_report.json"

    summary_fields = [
        "label",
        "reads_sampled",
        "expected_key",
        "expected_chemistry",
        "expected_whitelist",
        "best_key",
        "best_rate",
        "second_key",
        "second_rate",
        "sn_ratio",
        "status",
        "fastq_paths",
    ]
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields)
        writer.writeheader()
        for report in reports:
            writer.writerow({
                "label": report.label,
                "reads_sampled": report.reads_sampled,
                "expected_key": report.expected_key,
                "expected_chemistry": report.expected_chemistry,
                "expected_whitelist": report.expected_whitelist,
                "best_key": report.best_key,
                "best_rate": f"{report.best_rate:.6f}",
                "second_key": report.second_key,
                "second_rate": f"{report.second_rate:.6f}",
                "sn_ratio": "inf" if math.isinf(report.sn_ratio) else f"{report.sn_ratio:.3f}",
                "status": report.status,
                "fastq_paths": ";".join(report.fastq_paths),
            })

    with open(rates_path, "w", newline="") as handle:
        fields = ["label", "whitelist_key", "hits", "rate"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for report in reports:
            for key in sorted(report.rates):
                writer.writerow({
                    "label": report.label,
                    "whitelist_key": key,
                    "hits": report.hits[key],
                    "rate": f"{report.rates[key]:.6f}",
                })

    payload = {
        "parameters": parameters,
        "whitelists": [asdict(w) for w in whitelist_reports],
        "reports": [asdict(r) for r in reports],
        "version": "0.1.0",
    }
    with open(json_path, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    log.info("Wrote %s", summary_path)
    log.info("Wrote %s", rates_path)
    log.info("Wrote %s", json_path)


def run(args: argparse.Namespace) -> int:
    specs = [parse_whitelist_spec(raw) for raw in args.whitelist]
    duplicate_keys = [spec.key for spec in specs if [s.key for s in specs].count(spec.key) > 1]
    if duplicate_keys:
        raise ValueError(f"Duplicate whitelist keys: {sorted(set(duplicate_keys))}")

    whitelist_sets: Dict[str, Set[str]] = {}
    whitelist_reports: List[WhitelistLoadReport] = []
    path_to_key: Dict[str, str] = {}
    path_to_family: Dict[str, str] = {}
    for spec in specs:
        path = Path(spec.path)
        if not path.exists():
            raise FileNotFoundError(f"Whitelist not found: {path}")
        log.info("Loading whitelist %s from %s", spec.key, path)
        wl, load_report = load_whitelist(spec, args.barcode_len)
        whitelist_sets[spec.key] = wl
        whitelist_reports.append(load_report)
        path_to_key[str(path)] = spec.key
        path_to_key[str(path.resolve())] = spec.key
        path_to_family[str(path)] = spec.family
        path_to_family[str(path.resolve())] = spec.family
        log.info("  %d barcodes", len(wl))

    checks: List[FastqCheck] = []
    if args.manifest:
        checks.extend(load_manifest_checks(
            args.manifest,
            args.max_fastqs_per_manifest_row,
            path_to_key,
            path_to_family,
        ))
    for raw_fastq in args.fastq or []:
        checks.append(build_fastq_check_from_arg(raw_fastq))

    if not checks:
        raise ValueError("Provide --manifest and/or one or more --fastq LABEL:PATH[,PATH...]")

    reports: List[FastqReport] = []
    for check in checks:
        paths = [Path(path) for path in check.fastq_paths]
        missing = [str(path) for path in paths if not path.exists()]
        if missing:
            raise FileNotFoundError(f"{check.label}: FASTQ path(s) not found: {missing}")
        if not paths:
            reports.append(FastqReport(
                label=check.label,
                fastq_paths=[],
                reads_sampled=0,
                barcode_len=args.barcode_len,
                expected_chemistry=check.expected_chemistry,
                expected_whitelist=check.expected_whitelist,
                expected_family=check.expected_family,
                expected_key=check.expected_key,
                best_key="",
                best_family="",
                best_chemistry="",
                best_hits=0,
                best_rate=0.0,
                second_key="",
                second_rate=0.0,
                sn_ratio=0.0,
                status="FAIL_NO_FASTQ",
                rates={key: 0.0 for key in whitelist_sets},
                hits={key: 0 for key in whitelist_sets},
            ))
            continue

        log.info("Sampling %s from %d FASTQ(s)", check.label, len(paths))
        barcodes = sample_barcodes(paths, args.sample_reads, args.barcode_len)
        hits, rates = score_barcodes(barcodes, whitelist_sets)
        ranked = sorted(rates.items(), key=lambda item: item[1], reverse=True)
        best_key, best_rate = ranked[0] if ranked else ("", 0.0)
        second_key, second_rate = ranked[1] if len(ranked) > 1 else ("", 0.0)
        sn_ratio = best_rate / second_rate if second_rate > 0 else float("inf")
        best_family, best_chemistry = ("", "")
        if best_key:
            best_family, best_chemistry = best_key.split(":", 1)
        status = call_status(
            check.expected_key,
            best_key,
            best_rate,
            second_rate,
            args.min_hit_rate,
            args.dominance_ratio,
        )
        if check.expected_key and check.expected_chemistry and best_chemistry != check.expected_chemistry:
            status = "FAIL_EXPECTED_MISMATCH"

        log.info(
            "  %s best=%s %.4f second=%s %.4f status=%s",
            check.label,
            best_key,
            best_rate,
            second_key,
            second_rate,
            status,
        )
        reports.append(FastqReport(
            label=check.label,
            fastq_paths=[str(path) for path in paths],
            reads_sampled=len(barcodes),
            barcode_len=args.barcode_len,
            expected_chemistry=check.expected_chemistry,
            expected_whitelist=check.expected_whitelist,
            expected_family=check.expected_family,
            expected_key=check.expected_key,
            best_key=best_key,
            best_family=best_family,
            best_chemistry=best_chemistry,
            best_hits=hits.get(best_key, 0),
            best_rate=best_rate,
            second_key=second_key,
            second_rate=second_rate,
            sn_ratio=sn_ratio,
            status=status,
            rates=rates,
            hits=hits,
        ))

    params = {
        "sample_reads": args.sample_reads,
        "barcode_len": args.barcode_len,
        "min_hit_rate": args.min_hit_rate,
        "dominance_ratio": args.dominance_ratio,
        "manifest": str(args.manifest) if args.manifest else "",
        "max_fastqs_per_manifest_row": args.max_fastqs_per_manifest_row,
    }
    write_reports(reports, whitelist_reports, args.outdir, params)

    print()
    print("=" * 72)
    print("WHITELIST FAMILY PREFLIGHT")
    print("=" * 72)
    for report in reports:
        sn = "inf" if math.isinf(report.sn_ratio) else f"{report.sn_ratio:.1f}x"
        expected = report.expected_key or report.expected_chemistry or "not-set"
        print(
            f"{report.status:24s} {report.label:28s} "
            f"expected={expected:18s} best={report.best_key:18s} "
            f"rate={report.best_rate:.3f} S/N={sn}"
        )
    print(f"Full report: {args.outdir / 'whitelist_family_report.json'}")
    print("=" * 72)

    has_failures = any(report.status.startswith("FAIL") for report in reports)
    if has_failures and not args.no_fail:
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Detect FASTQ barcode whitelist family/namespace before STAR runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Optional TSV manifest with provider_group/library/chemistry/whitelist/fastq_root/fastq_sample_ids columns.",
    )
    parser.add_argument(
        "--fastq",
        action="append",
        metavar="LABEL:PATH[,PATH...]",
        help="Ad-hoc R1 FASTQ set to inspect. May be repeated.",
    )
    parser.add_argument(
        "--whitelist",
        action="append",
        required=True,
        metavar="FAMILY:CHEMISTRY:PATH[:COLUMN]",
        help="Whitelist to score against. CHEMISTRY is TRU or NXT. COLUMN is 1-based and defaults to 1.",
    )
    parser.add_argument(
        "--sample-reads",
        type=int,
        default=DEFAULT_SAMPLE_READS,
        help=f"Total reads to sample per manifest row or --fastq set (default: {DEFAULT_SAMPLE_READS:,}).",
    )
    parser.add_argument(
        "--barcode-len",
        type=int,
        default=DEFAULT_BARCODE_LEN,
        help=f"Cell barcode length at R1 start (default: {DEFAULT_BARCODE_LEN}).",
    )
    parser.add_argument(
        "--min-hit-rate",
        type=float,
        default=DEFAULT_MIN_HIT_RATE,
        help=f"Minimum best whitelist exact-hit rate to avoid low-hit failure (default: {DEFAULT_MIN_HIT_RATE}).",
    )
    parser.add_argument(
        "--dominance-ratio",
        type=float,
        default=DEFAULT_DOMINANCE_RATIO,
        help=f"Minimum best/second hit-rate ratio for high confidence (default: {DEFAULT_DOMINANCE_RATIO}).",
    )
    parser.add_argument(
        "--max-fastqs-per-manifest-row",
        type=int,
        default=DEFAULT_MAX_FASTQS_PER_MANIFEST_ROW,
        help="Maximum R1 FASTQ files to sample per manifest row (default: 4).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for whitelist-family preflight reports.",
    )
    parser.add_argument(
        "--no-fail",
        action="store_true",
        help="Always exit 0 even when rows fail. Reports still record FAIL_* statuses.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    try:
        return run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
