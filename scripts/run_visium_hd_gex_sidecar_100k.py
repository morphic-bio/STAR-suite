#!/usr/bin/env python3
"""Run the frozen 100K Visium HD 3-prime GEX sidecar clean-room workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO


SUITE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = Path(
    "/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/"
    "downsample_100k_v1"
)
DEFAULT_GENOME = Path(
    "/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star"
)
DEFAULT_COMPANION = Path(
    "/mnt/pikachu/visium-hd-processing-ovarian-replication-20260722"
)
DEFAULT_CONTRACT = Path(
    "/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract"
)
DEFAULT_BC1 = Path(
    "/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt"
)
DEFAULT_BC2 = Path(
    "/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt"
)
DEFAULT_OUTPUT = Path(
    "/mnt/pikachu/star-spatial/gex_sidecar_tests/20260722_ovarian_100k_v1"
)
PROHIBITED_RUN_TREE = Path("/mnt/pikachu/star-spatial/runs")


def producer_thread_budgets(
    total: int,
    mode: str,
    r1_threads: int | None,
    star_threads: int | None,
) -> tuple[int, int]:
    if total < 1:
        raise ValueError("total thread budget must be positive")
    if mode not in {"serial", "concurrent"}:
        raise ValueError(f"invalid producer mode: {mode}")
    for label, value in (("R1", r1_threads), ("STAR", star_threads)):
        if value is not None and value < 1:
            raise ValueError(f"{label} thread budget must be positive")

    if mode == "serial":
        resolved_r1 = total if r1_threads is None else r1_threads
        resolved_star = total if star_threads is None else star_threads
        if resolved_r1 > total or resolved_star > total:
            raise ValueError("a serial producer thread budget exceeds --threads")
        return resolved_r1, resolved_star

    if total < 2:
        raise ValueError("concurrent producer mode requires --threads >= 2")
    if r1_threads is None and star_threads is None:
        resolved_r1 = total // 2
        resolved_star = total - resolved_r1
    elif r1_threads is None:
        assert star_threads is not None
        resolved_star = star_threads
        resolved_r1 = total - resolved_star
    elif star_threads is None:
        resolved_r1 = r1_threads
        resolved_star = total - resolved_r1
    else:
        resolved_r1 = r1_threads
        resolved_star = star_threads
    if resolved_r1 < 1 or resolved_star < 1:
        raise ValueError("concurrent R1 and STAR producers each require at least one thread")
    if resolved_r1 + resolved_star > total:
        raise ValueError("concurrent R1 + STAR thread budgets exceed --threads")
    return resolved_r1, resolved_star


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--genome-dir", type=Path, default=DEFAULT_GENOME)
    parser.add_argument("--companion-root", type=Path, default=DEFAULT_COMPANION)
    parser.add_argument("--barcode-contract-dir", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--bc1-oligos", type=Path, default=DEFAULT_BC1)
    parser.add_argument("--bc2-oligos", type=Path, default=DEFAULT_BC2)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--threads", type=int, default=16,
        help="total build/downstream thread budget and concurrent producer cap",
    )
    parser.add_argument(
        "--producer-mode", choices=("concurrent", "serial"), default="serial",
        help="run R1 decode/H0 and STAR sidecar as concurrent branches or serial controls",
    )
    parser.add_argument(
        "--r1-threads", type=int,
        help="R1 decoder threads (default: half of --threads concurrently, all serially)",
    )
    parser.add_argument(
        "--star-threads", type=int,
        help="STAR threads (default: remaining --threads concurrently, all serially)",
    )
    parser.add_argument("--sort-memory-mb", type=int, default=2048)
    result = parser.parse_args()
    if result.threads < 1 or result.sort_memory_mb < 1:
        parser.error("thread and sort-memory values must be positive")
    try:
        result.r1_threads, result.star_threads = producer_thread_budgets(
            result.threads, result.producer_mode, result.r1_threads, result.star_threads
        )
    except ValueError as error:
        parser.error(str(error))
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def git_identity(root: Path) -> dict[str, object]:
    def git(*values: str) -> str:
        return subprocess.check_output(["git", *values], cwd=root, text=True).strip()

    status = git("status", "--porcelain")
    return {
        "root": str(root.resolve()),
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "dirty": bool(status),
        "status": status.splitlines(),
    }


@dataclass
class RunningCommand:
    label: str
    process: subprocess.Popen[bytes]
    stdout: BinaryIO
    stderr: BinaryIO
    record: dict[str, object]
    started: float
    finished: bool = False


class Driver:
    def __init__(self, out_dir: Path) -> None:
        self.out_dir = out_dir
        self.logs = out_dir / "logs"
        self.commands: list[dict[str, object]] = []

    def save_commands(self) -> None:
        temporary = self.out_dir / "commands.json.tmp"
        temporary.write_text(
            json.dumps(self.commands, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(self.out_dir / "commands.json")

    def start(
        self, label: str, command: list[str], cwd: Path | None = None
    ) -> RunningCommand:
        record = {
            "label": label,
            "argv": command,
            "cwd": None if cwd is None else str(cwd.resolve()),
            "state": "running",
            "started_unix_seconds": round(time.time(), 6),
        }
        self.commands.append(record)
        self.save_commands()
        stdout = (self.logs / f"{label}.stdout.log").open("wb")
        stderr = (self.logs / f"{label}.stderr.log").open("wb")
        try:
            process = subprocess.Popen(
                command,
                cwd=cwd,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
        except BaseException as error:
            stdout.close()
            stderr.close()
            record["state"] = "launch_failed"
            record["launch_error"] = str(error)
            self.save_commands()
            raise
        return RunningCommand(label, process, stdout, stderr, record, time.monotonic())

    def _finish(self, running: RunningCommand, returncode: int) -> int:
        if running.finished:
            return int(running.record["exit_code"])
        running.stdout.close()
        running.stderr.close()
        running.record["exit_code"] = returncode
        running.record["finished_unix_seconds"] = round(time.time(), 6)
        running.record["wall_seconds"] = round(time.monotonic() - running.started, 6)
        if running.record.get("terminated_by_driver"):
            running.record["state"] = "terminated"
        else:
            running.record["state"] = "complete" if returncode == 0 else "failed"
        running.finished = True
        self.save_commands()
        return returncode

    def poll(self, running: RunningCommand) -> int | None:
        if running.finished:
            return int(running.record["exit_code"])
        returncode = running.process.poll()
        if returncode is None:
            return None
        return self._finish(running, returncode)

    def stop(self, running: RunningCommand) -> None:
        if self.poll(running) is not None:
            return
        running.record["terminated_by_driver"] = True
        try:
            os.killpg(running.process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            returncode = running.process.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(running.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            returncode = running.process.wait()
        self._finish(running, returncode)

    def ensure_success(self, running: RunningCommand, returncode: int) -> None:
        if returncode:
            raise RuntimeError(
                f"{running.label} failed with exit code {returncode}; see {self.logs}"
            )

    def run(self, label: str, command: list[str], cwd: Path | None = None) -> None:
        running = self.start(label, command, cwd)
        try:
            while True:
                returncode = self.poll(running)
                if returncode is not None:
                    self.ensure_success(running, returncode)
                    return
                time.sleep(0.05)
        except BaseException:
            self.stop(running)
            raise


def run_producer_stages(
    driver: Driver,
    mode: str,
    decoder_command: list[str],
    h0_command: list[str],
    star_command: list[str],
) -> float:
    started = time.monotonic()
    if mode == "serial":
        driver.run("r1_decode", decoder_command)
        driver.run("h0_prior", h0_command)
        driver.run("star_sidecar", star_command)
        return time.monotonic() - started
    if mode != "concurrent":
        raise ValueError(f"invalid producer mode: {mode}")

    running: dict[str, RunningCommand] = {}
    completed: set[str] = set()
    try:
        running["r1_decode"] = driver.start("r1_decode", decoder_command)
        running["star_sidecar"] = driver.start("star_sidecar", star_command)
        while completed != {"r1_decode", "h0_prior", "star_sidecar"}:
            for label, command in list(running.items()):
                if label in completed:
                    continue
                returncode = driver.poll(command)
                if returncode is None:
                    continue
                driver.ensure_success(command, returncode)
                completed.add(label)
                if label == "r1_decode":
                    running["h0_prior"] = driver.start("h0_prior", h0_command)
            if completed != {"r1_decode", "h0_prior", "star_sidecar"}:
                time.sleep(0.05)
    except BaseException:
        for command in running.values():
            driver.stop(command)
        raise
    return time.monotonic() - started


def require_sources(args: argparse.Namespace) -> tuple[list[Path], list[Path]]:
    fixture_summary = args.fixture / "summary.json"
    checksums = args.fixture / "checksums.sha256"
    required_files = [
        fixture_summary,
        checksums,
        args.bc1_oligos,
        args.bc2_oligos,
        args.barcode_contract_dir / "barcode_coords.tsv",
        args.barcode_contract_dir / "raw28_barcode_coords.tsv",
        args.barcode_contract_dir / "bc1_whitelist.txt",
        args.barcode_contract_dir / "bc2_whitelist.txt",
        args.barcode_contract_dir / "summary.json",
        args.companion_root / "native/hd_r1_anchored_decode.cpp",
        args.companion_root / "scripts/build_hd_h0_read_prior_from_oligo_stats.py",
    ]
    for path in required_files:
        if not path.is_file():
            raise FileNotFoundError(path)
    if not args.genome_dir.is_dir():
        raise FileNotFoundError(args.genome_dir)
    summary = json.loads(fixture_summary.read_text(encoding="utf-8"))
    if summary.get("counts") != {"output_fastqs": 8, "read_pairs": 100000}:
        raise ValueError("fixture summary does not declare exactly 100,000 pairs/eight FASTQs")
    lanes = summary.get("lanes")
    if not isinstance(lanes, list) or len(lanes) != 4:
        raise ValueError("fixture summary does not contain four lanes")
    r1: list[Path] = []
    r2: list[Path] = []
    for lane_index, lane in enumerate(lanes):
        if lane.get("lane_order") != lane_index or lane.get("read_pairs") != 25000:
            raise ValueError("fixture lane order/count contract failed")
        if lane.get("read_lengths") != {"R1": {"43": 25000}, "R2": {"75": 25000}}:
            raise ValueError("fixture read-length contract failed")
        r1_path = args.fixture / "fastqs" / Path(lane["output_r1"]["path"]).name
        r2_path = args.fixture / "fastqs" / Path(lane["output_r2"]["path"]).name
        if not r1_path.is_file() or not r2_path.is_file():
            raise FileNotFoundError("frozen lane FASTQ is missing")
        r1.append(r1_path)
        r2.append(r2_path)
    return r1, r2


def validate_output_root(path: Path) -> None:
    resolved = path.resolve()
    prohibited = PROHIBITED_RUN_TREE.resolve()
    if resolved == prohibited or prohibited in resolved.parents:
        raise ValueError(f"output may not be under prohibited run tree: {prohibited}")
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            raise FileExistsError(f"output directory must be absent or empty: {path}")
    else:
        path.mkdir(parents=True)


def parse_tsv(
    path: Path, expected_header: list[str] | None = None
) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty TSV: {path}")
    if expected_header is not None:
        matching = [
            index
            for index, line in enumerate(lines)
            if line.split("\t") == expected_header
        ]
        if len(matching) != 1:
            raise ValueError(f"expected exactly one table header in TSV: {path}")
        lines = lines[matching[0] :]
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"invalid TSV row: {path}")
        rows.append(dict(zip(header, fields, strict=True)))
    return rows


def main() -> int:
    args = arguments()
    validate_output_root(args.out_dir)
    for directory in ("logs", "bin", "decoder", "star", "join", "resolver_a", "materialized"):
        (args.out_dir / directory).mkdir()
    driver = Driver(args.out_dir)
    r1, r2 = require_sources(args)

    source_provenance = {
        "star_suite": git_identity(SUITE_ROOT),
        "visium_hd_processing": git_identity(args.companion_root),
    }
    (args.out_dir / "source_provenance.json").write_text(
        json.dumps(source_provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    dirty_sources = [
        name for name, provenance in source_provenance.items() if provenance["dirty"]
    ]
    if dirty_sources:
        raise RuntimeError(
            "source-only run requires committed, clean worktrees: "
            + ", ".join(dirty_sources)
        )
    driver.run("fixture_checksums", ["sha256sum", "-c", "checksums.sha256"], args.fixture)

    source = SUITE_ROOT / "core/legacy/source"
    driver.run("star_clean", ["make", "-C", str(source), "clean"])
    driver.run(
        "star_build", ["make", "-C", str(source), f"-j{args.threads}", "STAR"]
    )
    shutil.copy2(source / "STAR", args.out_dir / "bin/STAR")

    hts_cflags = subprocess.check_output(
        ["pkg-config", "--cflags", "htslib"], text=True
    ).split()
    hts_libs = subprocess.check_output(
        ["pkg-config", "--libs", "htslib"], text=True
    ).split()
    driver.run(
        "decoder_build",
        [
            "g++", "-std=c++17", "-O3", "-pthread", *hts_cflags,
            str(args.companion_root / "native/hd_r1_anchored_decode.cpp"),
            "-o", str(args.out_dir / "bin/hd_r1_anchored_decode"),
            *hts_libs, "-lz",
        ],
    )
    resolver_tool = SUITE_ROOT / "flex/tools/molecule_first_resolver"
    flex_source = SUITE_ROOT / "flex/source/solo/MoleculeFirstResolver.cpp"
    multigene_source = source / "MultiGeneUmiCr.cpp"
    driver.run(
        "resolver_build",
        [
            "g++", "-std=c++11", "-O3", "-Wall", "-Wextra",
            f"-I{SUITE_ROOT / 'flex/source'}", f"-I{source}",
            str(resolver_tool / "molecule_first_resolver.cpp"), str(flex_source),
            str(multigene_source), "-o", str(args.out_dir / "bin/molecule_first_resolver"),
        ],
    )
    driver.run(
        "join_build",
        [
            "g++", "-std=c++11", "-O3", "-Wall", "-Wextra", f"-I{source}",
            str(resolver_tool / "spatial_feature_sidecar_join.cpp"),
            str(source / "SpatialFeatureSidecar.cpp"), "-o",
            str(args.out_dir / "bin/spatial_feature_sidecar_join"), "-lz", "-lcrypto",
        ],
    )
    driver.run(
        "materializer_build",
        [
            "g++", "-std=c++11", "-O3", "-Wall", "-Wextra",
            str(resolver_tool / "molecule_first_materialize.cpp"), "-o",
            str(args.out_dir / "bin/molecule_first_materialize"),
        ],
    )

    candidate_path = args.out_dir / "decoder/raw_r1_candidates.tsv"
    oligo_stats = args.out_dir / "decoder/oligo_mutation_stats.tsv"
    decode_reads = args.out_dir / "decoder/decode_reads.tsv"
    decoder_command = [str(args.out_dir / "bin/hd_r1_anchored_decode")]
    for path in r1:
        decoder_command.extend(("--r1-fastq", str(path)))
    decoder_command.extend(
        (
            "--bc1-oligos", str(args.bc1_oligos), "--bc2-oligos", str(args.bc2_oligos),
            "--out", str(decode_reads), "--candidate-preserving-out", str(candidate_path),
            "--oligo-mutation-stats-out", str(oligo_stats), "--direct-tiered-h2-decode",
            "--full-start-min", "8", "--full-start-max", "12", "--grid-rows", "3350",
            "--grid-cols", "3350", "--threads", str(args.r1_threads),
            "--progress-log-every-batches", "1",
        )
    )

    h0_prior = args.out_dir / "decoder/h0_read_prior.tsv"
    h0_summary = args.out_dir / "decoder/h0_read_prior_summary.json"
    h0_command = [
        sys.executable,
        str(args.companion_root / "scripts/build_hd_h0_read_prior_from_oligo_stats.py"),
        "--oligo-stats", str(oligo_stats), "--out-tsv", str(h0_prior),
        "--summary-json", str(h0_summary),
    ]

    sidecar_prefix = args.out_dir / "star/gex_features"
    star_command = [
        str(args.out_dir / "bin/STAR"), "--runThreadN", str(args.star_threads),
        "--genomeDir", str(args.genome_dir), "--readFilesIn",
        ",".join(map(str, r2)), ",".join(map(str, r1)), "--readFilesCommand", "zcat",
        "--outFileNamePrefix", str(args.out_dir / "star/") ,
        "--clipAdapterType", "CellRanger4", "--outFilterScoreMin", "30",
        "--soloType", "None", "--soloFeatures", "GeneFull",
        "--soloCrGexFeature", "GeneFull", "--soloCrMultimapRescue", "yes",
        "--soloCrMultimapRescueIntronic", "auto", "--soloUMIdedup", "1MM_CR",
        "--soloUMIfiltering", "MultiGeneUMI_CR", "--soloMultiMappers", "Unique",
        "--soloStrand", "Forward", "--soloCellFilter", "None", "--outSAMtype", "None",
        "--soloSpatialFeatureSidecar", str(sidecar_prefix),
    ]
    prohibited_tokens = {"GX", "GN", "UR", "UB", "CB", "CR", "SAM", "BAM"}
    if any(token in prohibited_tokens for token in star_command):
        raise AssertionError("rendered STAR command contains a prohibited tag or output token")
    producer_wall_seconds = run_producer_stages(
        driver,
        args.producer_mode,
        decoder_command,
        h0_command,
        star_command,
    )

    normalized = args.out_dir / "join/normalized_evidence.tsv"
    join_summary = args.out_dir / "join/summary.json"
    join_command = [
        str(args.out_dir / "bin/spatial_feature_sidecar_join"),
        "--sidecar", str(sidecar_prefix) + ".bin",
        "--features", str(sidecar_prefix) + ".features.tsv",
        "--read-name-digests", str(sidecar_prefix) + ".read_name_digests.tsv",
        "--candidates", str(candidate_path), "--h0-prior", str(h0_prior),
        "--barcode-contract", str(args.barcode_contract_dir / "barcode_coords.tsv"),
        "--expected-reads", "100000", "--r1-length", "43", "--r2-length", "75",
        "--output", str(normalized), "--summary", str(join_summary),
    ]
    for path in r1:
        join_command.extend(("--r1-fastq", str(path)))
    for path in r2:
        join_command.extend(("--r2-fastq", str(path)))
    driver.run("sidecar_join", join_command)

    resolver_command = [
        str(args.out_dir / "bin/molecule_first_resolver"), "--input", str(normalized),
        "--gex-multigene-umi-cr", "--out-dir", str(args.out_dir / "resolver_a"),
    ]
    driver.run("resolver_a", resolver_command)
    resolver_b = args.out_dir / "resolver_b"
    driver.run(
        "resolver_b",
        [
            str(args.out_dir / "bin/molecule_first_resolver"), "--input", str(normalized),
            "--gex-multigene-umi-cr", "--out-dir", str(resolver_b),
        ],
    )
    deterministic_files = (
        "strict_molecules.tsv", "hard_molecules.tsv", "gated_hard_molecules.tsv"
    )
    for name in deterministic_files:
        if sha256(args.out_dir / "resolver_a" / name) != sha256(resolver_b / name):
            raise RuntimeError(f"integer resolver product is nondeterministic: {name}")

    driver.run(
        "materialize",
        [
            str(args.out_dir / "bin/molecule_first_materialize"),
            "--resolved-dir", str(args.out_dir / "resolver_a"),
            "--out-dir", str(args.out_dir / "materialized"), "--assay", "visium-hd",
            "--umi-mode", "1mm_cr", "--sort-memory-mb", str(args.sort_memory_mb),
            "--tmp-dir", str(args.out_dir / "materialized"),
        ],
    )

    sidecar_summary = json.loads(
        (Path(str(sidecar_prefix) + ".summary.json")).read_text(encoding="utf-8")
    )
    joined = json.loads(join_summary.read_text(encoding="utf-8"))
    if sidecar_summary.get("status") != "complete" or sidecar_summary.get("total_reads") != 100000:
        raise RuntimeError("sidecar completion/read-count gate failed")
    if joined.get("status") != "complete" or joined.get("total_reads") != 100000:
        raise RuntimeError("ordinal join completion/read-count gate failed")
    resolver_config_rows = parse_tsv(args.out_dir / "resolver_a/resolved_config.tsv")
    if any(set(row) != {"key", "value"} for row in resolver_config_rows):
        raise RuntimeError("unexpected resolver config schema")
    resolver_config = {row["key"]: row["value"] for row in resolver_config_rows}
    if resolver_config.get("gex_multigene_umi_cr") != "1":
        raise RuntimeError("GEX MultiGeneUMI_CR resolver gate failed")
    matrices = parse_tsv(
        args.out_dir / "materialized/summary.tsv",
        [
            "product", "scale", "features", "barcodes", "nnz", "mass",
            "matrix_field",
        ],
    )
    if len(matrices) != 12:
        raise RuntimeError("materializer did not produce four policies at three scales")
    expected_products = {"strict", "soft_expected", "hard", "gated_hard"}
    expected_scales = {"square_002um", "square_008um", "square_016um"}
    if {row["product"] for row in matrices} != expected_products or {
        row["scale"] for row in matrices
    } != expected_scales:
        raise RuntimeError("materialized product/scale axes are incomplete")
    for row in matrices:
        mass = float(row["mass"])
        if not (mass >= 0.0 and mass < float("inf")):
            raise RuntimeError("matrix mass is non-finite or negative")
        if int(row["features"]) < 1 or int(row["barcodes"]) < 1 or int(row["nnz"]) < 0:
            raise RuntimeError("invalid MEX dimensions")

    completion = {
        "schema": "star_suite.visium_hd_gex_sidecar_100k.v1",
        "status": "complete",
        "fixture": identity(args.fixture / "summary.json"),
        "source_provenance": source_provenance,
        "producer_execution": {
            "mode": args.producer_mode,
            "total_thread_budget": args.threads,
            "r1_threads": args.r1_threads,
            "star_threads": args.star_threads,
            "wall_seconds": round(producer_wall_seconds, 6),
        },
        "counts": {
            "read_pairs": 100000,
            "sidecar_records": sidecar_summary["total_reads"],
            "candidate_reads": joined["candidate_reads"],
            "candidate_rows": joined["candidate_rows"],
            "emitted_reads": joined["emitted_reads"],
            "emitted_rows": joined["emitted_rows"],
            "matrices": len(matrices),
        },
        "invariants": {
            "fixture_checksums_passed": True,
            "producer_join_barrier_complete": True,
            "concurrent_thread_budget_respected": (
                args.producer_mode != "concurrent"
                or args.r1_threads + args.star_threads <= args.threads
            ),
            "lane_order_and_name_digests_match": True,
            "raw_umi_from_r1": True,
            "coordinate_contract_validated": True,
            "gex_multigene_umi_cr_enabled": True,
            "integer_products_repeat_deterministically": True,
            "all_scale_masses_conserved_by_materializer": True,
            "no_bam_or_alignment_tag_carrier": True,
        },
        "outputs": {
            "sidecar": identity(Path(str(sidecar_prefix) + ".bin")),
            "normalized_evidence": identity(normalized),
            "resolver_summary": identity(args.out_dir / "resolver_a/summary.tsv"),
            "materializer_summary": identity(args.out_dir / "materialized/summary.tsv"),
            "commands": identity(args.out_dir / "commands.json"),
        },
    }
    temporary = args.out_dir / "RUN_COMPLETE.json.tmp"
    temporary.write_text(
        json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(args.out_dir / "RUN_COMPLETE.json")
    print(json.dumps(completion["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
