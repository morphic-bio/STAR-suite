#!/usr/bin/env python3
"""Verify the pinned SLAM fixture and the effective STAR command."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shlex
import sys
from pathlib import Path
from typing import Iterable


def sha256_bytes(chunks: Iterable[bytes]) -> str:
    digest = hashlib.sha256()
    for chunk in chunks:
        digest.update(chunk)
    return digest.hexdigest()


def file_chunks(path: Path, opener=open) -> Iterable[bytes]:
    with opener(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                return
            yield chunk


def artifact_digest(path: Path, mode: str) -> str:
    if mode == "raw":
        return sha256_bytes(file_chunks(path))
    if mode == "gzip_content":
        return sha256_bytes(file_chunks(path, gzip.open))
    if mode == "sorted_lines":
        lines = path.read_bytes().splitlines(keepends=True)
        return sha256_bytes(sorted(lines))
    raise ValueError(f"unsupported hash mode: {mode}")


def count_gzip_lines(path: Path) -> int:
    with gzip.open(path, "rb") as handle:
        return sum(1 for _ in handle)


def read_index_parameters(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("###"):
                continue
            fields = line.split()
            if len(fields) >= 2:
                values[fields[0]] = " ".join(fields[1:])
    return values


def command_from_file(path: Path) -> list[str]:
    commands = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line and not line.startswith("#") and line != "set -euo pipefail":
            commands.append(line)
    if len(commands) != 1:
        raise ValueError(f"expected one command in {path}, found {len(commands)}")
    return shlex.split(commands[0])


def option_value(command: list[str], flag: str) -> str | None:
    positions = [i for i, token in enumerate(command) if token == flag]
    if len(positions) != 1:
        return None
    position = positions[0]
    if position + 1 >= len(command):
        return None
    return command[position + 1]


def verify(
    manifest: dict,
    fastq: Path,
    snp_bed: Path,
    reference: Path,
    genome_dir: Path,
    command: list[str],
    thresholds: dict[str, float],
) -> list[str]:
    errors: list[str] = []
    paths = {"fastq": fastq, "snp_bed": snp_bed, "reference": reference}

    for label, spec in manifest["artifacts"].items():
        path = paths[label]
        if not path.is_file():
            errors.append(f"{label}: missing file: {path}")
            continue
        actual = artifact_digest(path, spec["hash_mode"])
        if actual != spec["sha256"]:
            errors.append(
                f"{label}: sha256 mismatch for {path}: expected {spec['sha256']}, got {actual}"
            )
        if "line_count" in spec:
            with path.open("rb") as handle:
                actual_lines = sum(1 for _ in handle)
            if actual_lines != spec["line_count"]:
                errors.append(
                    f"{label}: line count mismatch: expected {spec['line_count']}, got {actual_lines}"
                )
        if "fastq_records" in spec:
            line_count = count_gzip_lines(path)
            if line_count % 4 != 0 or line_count // 4 != spec["fastq_records"]:
                errors.append(
                    f"{label}: FASTQ record mismatch: expected {spec['fastq_records']}, "
                    f"got {line_count / 4:g}"
                )

    if not genome_dir.is_dir():
        errors.append(f"index: missing directory: {genome_dir}")
    else:
        for relative, spec in manifest["index"]["files"].items():
            path = genome_dir / relative
            if not path.is_file():
                errors.append(f"index: missing file: {path}")
                continue
            if "size" in spec and path.stat().st_size != spec["size"]:
                errors.append(
                    f"index: {relative} size mismatch: expected {spec['size']}, "
                    f"got {path.stat().st_size}"
                )
            if "sha256" in spec:
                actual = artifact_digest(path, "raw")
                if actual != spec["sha256"]:
                    errors.append(
                        f"index: {relative} sha256 mismatch: expected {spec['sha256']}, got {actual}"
                    )

        parameters_path = genome_dir / "genomeParameters.txt"
        if not parameters_path.is_file():
            errors.append(f"index: missing file: {parameters_path}")
        else:
            actual_parameters = read_index_parameters(parameters_path)
            for key, expected in manifest["index"]["parameters"].items():
                actual = actual_parameters.get(key)
                if actual != expected:
                    errors.append(
                        f"index: parameter {key} mismatch: expected {expected!r}, got {actual!r}"
                    )

    contract = manifest["alignment_contract"]
    for flag, expected in contract["required_options"].items():
        actual = option_value(command, flag)
        if actual != expected:
            errors.append(
                f"command: {flag} must occur once with value {expected!r}; got {actual!r}"
            )
    for flag in contract["required_flags"]:
        if command.count(flag) != 1:
            errors.append(f"command: {flag} must occur exactly once")
    for flag in contract["forbidden_options"]:
        if flag in command:
            errors.append(f"command: forbidden fixture option present: {flag}")
    for flag, forbidden_values in contract.get("forbidden_option_values", {}).items():
        actual = option_value(command, flag)
        if actual in forbidden_values:
            errors.append(f"command: {flag} has forbidden value {actual!r}")

    dynamic_options = {
        "--genomeDir": str(genome_dir),
        "--readFilesIn": str(fastq),
        "--slamSnpMaskIn": str(snp_bed),
    }
    for flag, expected in dynamic_options.items():
        actual = option_value(command, flag)
        if actual != expected:
            errors.append(f"command: {flag} expected {expected!r}, got {actual!r}")

    parity_contract = manifest["parity_contract"]
    minimum_fields = ("ref_corr_min", "exact_min_pearson", "align_min_pearson")
    for field in minimum_fields:
        if thresholds[field] < parity_contract[field]:
            errors.append(
                f"threshold: {field} cannot be relaxed below {parity_contract[field]}; "
                f"got {thresholds[field]}"
            )
    maximum_field = "exact_max_abs_delta"
    if thresholds[maximum_field] > parity_contract[maximum_field]:
        errors.append(
            f"threshold: {maximum_field} cannot exceed {parity_contract[maximum_field]}; "
            f"got {thresholds[maximum_field]}"
        )

    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify pinned SLAM fixture artifacts, index, and STAR arguments."
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--snp-bed", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--genome-dir", required=True)
    parser.add_argument("--ref-corr-min", required=True, type=float)
    parser.add_argument("--exact-min-pearson", required=True, type=float)
    parser.add_argument("--exact-max-abs-delta", required=True, type=float)
    parser.add_argument("--align-min-pearson", required=True, type=float)
    parser.add_argument("--command-file")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command_file and args.command:
        parser.error("use either --command-file or a command after --, not both")
    if not args.command_file and not args.command:
        parser.error("provide --command-file or the effective STAR command after --")
    return args


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.manifest)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        command = (
            command_from_file(Path(args.command_file))
            if args.command_file
            else list(args.command)
        )
        errors = verify(
            manifest,
            Path(args.fastq),
            Path(args.snp_bed),
            Path(args.reference),
            Path(args.genome_dir),
            command,
            {
                "ref_corr_min": args.ref_corr_min,
                "exact_min_pearson": args.exact_min_pearson,
                "exact_max_abs_delta": args.exact_max_abs_delta,
                "align_min_pearson": args.align_min_pearson,
            },
        )
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"SLAM fixture manifest: ERROR: {exc}", file=sys.stderr)
        return 2

    if errors:
        print("SLAM fixture manifest: FAIL", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(f"SLAM fixture manifest: PASS ({manifest['fixture']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
