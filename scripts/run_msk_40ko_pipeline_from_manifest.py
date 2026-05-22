#!/usr/bin/env python3
"""Run the MSK 40KO three-library STAR/CR-compatible pipeline from a manifest.

The input FASTQs are flat in /mnt/pikachu/scRNAseq_40KO. This wrapper stages
per-sample symlink directories for GEX, PolyIII/gRNA, and LARRY, writes the
pf-multi config, runs STAR with the same perturb-seq surface used for the MSK
30KO production runs, and can optionally run remote downstream/CellBender and
Globus BAM transfer.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "docs" / "MSK_40KO_FASTQ_MANIFEST.tsv"
DEFAULT_STAR_BIN = REPO_ROOT / "core" / "legacy" / "source" / "STAR"
DEFAULT_GENOME_DIR = Path("/storage/autoindex_110_44/bulk_index")
DEFAULT_OUT_ROOT = Path(
    "/storage/MSK-perturb-comparison/"
    f"msk40ko_prod_{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%d_%H%M%S')}"
)
DEFAULT_REMOTE_ROOT = Path("/tmp/msk40ko_cellbender")
DEFAULT_DOWNSTREAM = REPO_ROOT / "scripts" / "run_remote_scrna_downstream_rsync.sh"
DEFAULT_DOWNSTREAM_NAME = "downstream_genefull_velocyto_cellbender"

LIBRARY_ORDER = ("GEX", "CRISPR_PolyIII", "LARRY")
STAGE_SUBDIR = {
    "GEX": "mRNA",
    "CRISPR_PolyIII": "PolyIII",
    "LARRY": "LARRY",
}
PF_LIBRARY_TYPE = {
    "GEX": "Gene Expression",
    "CRISPR_PolyIII": "CRISPR Guide Capture",
    "LARRY": "Custom",
}
PF_FEATURE_TYPES = {
    "GEX": "Gene Expression",
    "CRISPR_PolyIII": "CRISPR Guide Capture",
    "LARRY": "Custom",
}


@dataclass
class DownstreamJob:
    group: str
    sample_root: Path
    process: subprocess.Popen
    submit_log: Path
    pid_file: Path


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(message: str) -> None:
    print(f"[{utc_now()}] {message}", flush=True)


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sanitize_id(raw: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", raw.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "sample"


def split_semicolon(value: str) -> List[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def require_file(path: Path, label: str = "file") -> None:
    if not path.is_file():
        die(f"missing {label}: {path}")


def require_executable(path: Path, label: str = "executable") -> None:
    require_file(path, label)
    if not os.access(path, os.X_OK):
        die(f"{label} is not executable: {path}")


def require_dir(path: Path, label: str = "directory") -> None:
    if not path.is_dir():
        die(f"missing {label}: {path}")


def read_manifest(path: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
    required = {
        "provider_group",
        "library",
        "chemistry",
        "whitelist",
        "fastq_root",
        "fastq_sample_ids",
    }
    groups: Dict[str, Dict[str, Dict[str, str]]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            die(f"manifest {path} is missing required columns: {sorted(missing)}")
        for row in reader:
            group = row["provider_group"].strip()
            library = row["library"].strip()
            if library not in LIBRARY_ORDER:
                die(f"{group}: unsupported library {library!r}")
            groups.setdefault(group, {})[library] = row

    for group, rows in groups.items():
        missing_libs = [lib for lib in LIBRARY_ORDER if lib not in rows]
        if missing_libs:
            die(f"{group}: manifest missing libraries: {', '.join(missing_libs)}")
    return groups


def select_groups(all_groups: Sequence[str], samples: str, exclude_samples: str) -> List[str]:
    selected = list(all_groups)
    if samples:
        requested = split_semicolon(samples.replace(",", ";"))
        unknown = sorted(set(requested) - set(all_groups))
        if unknown:
            die(f"unknown samples in --samples: {', '.join(unknown)}")
        selected = requested
    if exclude_samples:
        excluded = set(split_semicolon(exclude_samples.replace(",", ";")))
        selected = [sample for sample in selected if sample not in excluded]
    if not selected:
        die("no samples selected")
    return selected


def row_fastqs(row: Dict[str, str]) -> List[Path]:
    root = Path(row["fastq_root"])
    require_dir(root, "FASTQ root")
    paths: List[Path] = []
    for sample_id in split_semicolon(row["fastq_sample_ids"]):
        matches = sorted(root.glob(f"{sample_id}_L*_*.fastq.gz"))
        if not matches:
            die(f"{row['provider_group']} {row['library']}: no FASTQs for sample id {sample_id} under {root}")
        paths.extend(matches)
    unique = sorted(dict.fromkeys(paths))
    expected = row.get("file_count", "").strip()
    if expected:
        try:
            expected_n = int(expected)
        except ValueError:
            expected_n = -1
        if expected_n >= 0 and len(unique) != expected_n:
            die(
                f"{row['provider_group']} {row['library']}: expected {expected_n} FASTQs "
                f"from manifest, found {len(unique)}"
            )
    return unique


def stage_fastqs(paths: Sequence[Path], stage_dir: Path, force: bool) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    for source in paths:
        dest = stage_dir / source.name
        if dest.exists() or dest.is_symlink():
            if dest.is_symlink() and Path(os.readlink(dest)) == source:
                continue
            if force:
                dest.unlink()
            else:
                die(f"staged FASTQ already exists and differs: {dest}")
        os.symlink(source, dest)


def list_role_fastqs(stage_dir: Path, role: str) -> List[Path]:
    return sorted(stage_dir.glob(f"*_L*_{role}_001.fastq.gz"))


def validate_pairs(stage_dir: Path) -> None:
    r1_files = list_role_fastqs(stage_dir, "R1")
    if not r1_files:
        die(f"no R1 FASTQs in {stage_dir}")
    for r1 in r1_files:
        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.exists():
            die(f"missing R2 pair for {r1}")


def comma_join(paths: Iterable[Path]) -> str:
    return ",".join(str(path) for path in paths)


def star_whitelist(row: Dict[str, str]) -> str:
    return (row.get("star_whitelist") or row.get("whitelist") or "").strip()


def feature_ref(row: Dict[str, str]) -> str:
    value = (row.get("feature_ref") or "").strip()
    return "" if value in {"", "NA", "na"} else value


def write_pf_config(group: str, rows: Dict[str, Dict[str, str]], stage_dirs: Dict[str, Path], out_path: Path) -> None:
    safe = sanitize_id(group).lower()
    with out_path.open("w", newline="") as handle:
        handle.write("[libraries]\n")
        writer = csv.writer(handle)
        writer.writerow([
            "fastqs",
            "sample",
            "library_type",
            "feature_types",
            "star_chemistry",
            "star_whitelist",
            "star_feature_ref",
            "star_library_id",
            "star_max_hamming",
        ])
        writer.writerow([
            stage_dirs["GEX"],
            group,
            PF_LIBRARY_TYPE["GEX"],
            PF_FEATURE_TYPES["GEX"],
            rows["GEX"]["chemistry"],
            "",
            "",
            f"gex_{safe}",
            "",
        ])
        writer.writerow([
            stage_dirs["CRISPR_PolyIII"],
            group,
            PF_LIBRARY_TYPE["CRISPR_PolyIII"],
            PF_FEATURE_TYPES["CRISPR_PolyIII"],
            rows["CRISPR_PolyIII"]["chemistry"],
            star_whitelist(rows["CRISPR_PolyIII"]),
            feature_ref(rows["CRISPR_PolyIII"]),
            f"grna_{safe}",
            "1",
        ])
        writer.writerow([
            stage_dirs["LARRY"],
            group,
            PF_LIBRARY_TYPE["LARRY"],
            PF_FEATURE_TYPES["LARRY"],
            rows["LARRY"]["chemistry"],
            star_whitelist(rows["LARRY"]),
            feature_ref(rows["LARRY"]),
            f"larry_{safe}",
            "1",
        ])


def build_star_cmd(args: argparse.Namespace, group: str, rows: Dict[str, Dict[str, str]], sample_root: Path) -> List[str]:
    run_dir = sample_root / "run"
    tmp_dir = sample_root / "tmp_STAR"
    pf_config = sample_root / "pf_multi_config.csv"
    gex_stage = sample_root / "staged_fastqs" / STAGE_SUBDIR["GEX"]
    r1_files = list_role_fastqs(gex_stage, "R1")
    r2_files = [Path(str(path).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")) for path in r1_files]
    solo_cb_whitelist = star_whitelist(rows["GEX"])
    if not solo_cb_whitelist:
        die(f"{group}: GEX row has no whitelist/star_whitelist")

    cmd = [
        str(args.star_bin),
        "--runThreadN", str(args.threads),
        "--genomeDir", str(args.genome_dir),
        "--readFilesIn", comma_join(r2_files), comma_join(r1_files),
        "--readFilesCommand", "zcat",
        "--outFileNamePrefix", str(run_dir) + "/",
        "--outTmpDir", str(tmp_dir),
        "--clipAdapterType", "CellRanger4",
        "--alignEndsType", "Local",
        "--chimSegmentMin", "1000000",
        "--clip3pPolyG", "yes",
        "--soloType", "CB_UMI_Simple",
        "--soloCBstart", "1",
        "--soloCBlen", "16",
        "--soloUMIstart", "17",
        "--soloUMIlen", "12",
        "--soloBarcodeReadLength", "0",
        "--soloInlineHashMode", "no",
        "--soloCBwhitelist", solo_cb_whitelist,
        "--soloStrand", "Forward",
        "--soloFeatures", "GeneFull", "Velocyto",
        "--soloUMIdedup", "1MM_CR",
        "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
        "--soloCellFilter", "EmptyDrops_CR",
        "--soloUMIfiltering", "MultiGeneUMI_CR",
        "--soloMultiMappers", "Unique",
        "--soloCbUbRequireTogether", "no",
        "--soloCrGexFeature", "genefull",
        "--soloCrMultimapRescue", "yes",
        "--pfMultiConfig", str(pf_config),
        "--crChemistry", "auto",
        "--crOutputChemistry", "TRU",
        "--crMinUmi", str(args.cr_min_umi),
        "--crAssignMaxHamming", "1",
        "--crAssignFeatureOffset", "0",
        "--crAssignLimitSearch", "-1",
        "--crAssignMinCounts", "0",
        "--crAssignMaxBarcodeMismatches", "5",
        "--crAssignFeatureN", "0",
        "--crAssignBarcodeN", "1",
        "--crAssignConsumerThreads", "-1",
        "--crAssignSearchThreads", "1",
        "--crAssignSkipQcOutputs", "1",
        "--defaultCrCompat", "yes",
        "--dynamicThreadInterface", "1",
        "--dynamicThreadConstMapPermits", str(args.threads),
        "--dynamicThreadTelemetry", "1",
    ]

    out_sam = args.out_samtype.strip()
    if out_sam.lower() == "none":
        cmd.extend(["--outSAMtype", "None"])
    else:
        cmd.append("--outSAMtype")
        cmd.extend(shlex.split(out_sam))
    return cmd


def command_script_text(cmd: Sequence[str], env: Dict[str, str]) -> str:
    lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
    for key, value in env.items():
        lines.append(f"export {shlex.quote(key)}={shlex.quote(value)}")
    lines.append(" ".join(shlex.quote(part) for part in cmd) + ' "$@"')
    lines.append("")
    return "\n".join(lines)


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if executable:
        path.chmod(path.stat().st_mode | 0o111)


def run_subprocess(cmd: Sequence[str], env: Optional[Dict[str, str]] = None) -> None:
    subprocess.run(cmd, check=True, env={**os.environ, **(env or {})})


def build_downstream_cmd(args: argparse.Namespace, sample_root: Path) -> List[str]:
    if not args.remote_host:
        die("--run-downstream requires --remote-host")
    require_file(args.downstream_wrapper, "downstream wrapper")
    cmd = [
        str(args.downstream_wrapper),
        "--sample-dir", str(sample_root),
        "--remote-host", args.remote_host,
        "--remote-root", str(args.remote_root),
        "--output-name", args.downstream_output_name,
        "--adaptive-filter",
        "--cellbender-cpu-cores", str(args.cellbender_cpu_cores),
        "--cellbender-layer", args.cellbender_layer,
    ]
    if args.run_cellbender:
        cmd.append("--run-cellbender")
    if args.cellbender_gpu:
        cmd.append("--cellbender-gpu")
    return cmd


def run_downstream(args: argparse.Namespace, group: str, sample_root: Path) -> Optional[DownstreamJob]:
    if not args.run_downstream:
        return None
    cmd = build_downstream_cmd(args, sample_root)
    log(f"{sample_root.name}: remote downstream start")
    if args.downstream_async:
        submit_log = sample_root / f"{args.downstream_output_name}.submit.log"
        pid_file = sample_root / f"{args.downstream_output_name}.pid"
        cmd_file = sample_root / f"{args.downstream_output_name}.command.sh"
        write_text(
            cmd_file,
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            + " ".join(shlex.quote(part) for part in cmd)
            + "\n",
            executable=True,
        )
        with submit_log.open("ab") as handle:
            handle.write(f"[{utc_now()}] async downstream launch\n".encode())
            handle.write((" ".join(shlex.quote(part) for part in cmd) + "\n").encode())
            process = subprocess.Popen(
                cmd,
                stdout=handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                env=os.environ.copy(),
            )
        pid_file.write_text(f"{process.pid}\n")
        log(f"{sample_root.name}: remote downstream submitted pid={process.pid}")
        return DownstreamJob(group=group, sample_root=sample_root, process=process, submit_log=submit_log, pid_file=pid_file)
    run_subprocess(cmd)
    log(f"{sample_root.name}: remote downstream complete")
    return None


def wait_downstream_jobs(jobs: Sequence[DownstreamJob], run_tsv: Path) -> None:
    for job in jobs:
        log(f"{job.group}: waiting for remote downstream pid={job.process.pid}")
        return_code = job.process.wait()
        finish_downstream_job(job, run_tsv, return_code)


def finish_downstream_job(job: DownstreamJob, run_tsv: Path, return_code: int) -> None:
    if return_code != 0:
        append_status(run_tsv, job.group, f"downstream_failed_rc_{return_code}", job.sample_root)
        die(f"{job.group}: remote downstream failed with rc={return_code}; see {job.submit_log}")
    append_status(run_tsv, job.group, "downstream_done", job.sample_root)
    write_text(job.sample_root / "RUN_COMPLETE.ok", f"completed_utc={utc_now()}\n")
    append_status(run_tsv, job.group, "done", job.sample_root)
    log(f"{job.group}: remote downstream complete")


def reap_finished_downstream_jobs(jobs: Sequence[DownstreamJob], run_tsv: Path) -> List[DownstreamJob]:
    active: List[DownstreamJob] = []
    for job in jobs:
        return_code = job.process.poll()
        if return_code is None:
            active.append(job)
        else:
            finish_downstream_job(job, run_tsv, return_code)
    return active


def wait_for_downstream_slot(jobs: Sequence[DownstreamJob], run_tsv: Path, max_active: int) -> List[DownstreamJob]:
    if max_active <= 0:
        return list(jobs)
    active = reap_finished_downstream_jobs(jobs, run_tsv)
    while len(active) >= max_active:
        job = active.pop(0)
        log(
            f"{job.group}: waiting for remote downstream pid={job.process.pid} "
            f"to keep <= {max_active} active remote jobs"
        )
        finish_downstream_job(job, run_tsv, job.process.wait())
        active = reap_finished_downstream_jobs(active, run_tsv)
    return active


def globus_enabled(args: argparse.Namespace) -> bool:
    return bool(args.globus_src_endpoint and args.globus_dst_endpoint and args.globus_dst_root)


def archive_bam_if_needed(args: argparse.Namespace, group: str, bam: Path) -> Path:
    if not args.bam_archive_root:
        return bam
    archive_dir = args.bam_archive_root / group
    archive_dir.mkdir(parents=True, exist_ok=True)
    archive_bam = archive_dir / bam.name
    if not archive_bam.exists() or archive_bam.stat().st_size != bam.stat().st_size:
        log(f"{group}: copying BAM to archive path for Globus source: {archive_bam}")
        shutil.copy2(bam, archive_bam)
    return archive_bam


def transfer_bam(args: argparse.Namespace, group: str, sample_root: Path) -> None:
    if not globus_enabled(args):
        return
    if shutil.which("globus") is None:
        die("globus CLI not found on PATH")
    bam = sample_root / "run" / "Aligned.out.bam"
    if not bam.exists():
        log(f"{group}: no BAM found for Globus transfer: {bam}")
        return
    source_bam = archive_bam_if_needed(args, group, bam)
    batch = sample_root / "globus_batch.tsv"
    task_json = sample_root / "globus_transfer_task.json"
    task_final = sample_root / "globus_transfer_task_final.json"
    dst_root = args.globus_dst_root.rstrip("/")
    dst_path = f"{dst_root}/{group}/{source_bam.name}"
    write_text(batch, f"{source_bam}\t{dst_path}\n")
    label = f"MSK 40KO BAM {group} {sample_root.parent.parent.name}"
    log(f"{group}: Globus transfer start -> {dst_path}")
    cmd = [
        "globus", "transfer",
        args.globus_src_endpoint,
        args.globus_dst_endpoint,
        "--batch", str(batch),
        "--sync-level", "checksum",
        "--label", label,
        "--notify", "off",
        "--format", "json",
    ]
    result = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    task_json.write_text(result.stdout)
    try:
        task_id = json.loads(result.stdout)["task_id"]
    except (KeyError, json.JSONDecodeError) as exc:
        die(f"{group}: could not parse Globus task id from {task_json}: {exc}")
    (sample_root / "globus_transfer_task_id.txt").write_text(task_id + "\n")
    if args.globus_wait:
        run_subprocess(["globus", "task", "wait", "--polling-interval", str(args.globus_poll_seconds), task_id])
        final = subprocess.run(
            ["globus", "task", "show", task_id, "--format", "json"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
        task_final.write_text(final.stdout)
        status = json.loads(final.stdout).get("status")
        if status != "SUCCEEDED":
            die(f"{group}: Globus task {task_id} ended with status {status}")
        if args.delete_local_bam_after_transfer:
            deleted = sample_root / "deleted_generated_large_files.tsv"
            with deleted.open("w") as handle:
                handle.write("deleted_utc\tsource_path\tsize_bytes\n")
                for path in (bam, source_bam):
                    if path.exists():
                        size = path.stat().st_size
                        path.unlink()
                        handle.write(f"{utc_now()}\t{path}\t{size}\n")
            log(f"{group}: deleted local BAM copies after successful Globus transfer")
    log(f"{group}: Globus transfer submitted task_id={task_id}")


def append_status(run_tsv: Path, group: str, status: str, sample_root: Path) -> None:
    new_file = not run_tsv.exists()
    with run_tsv.open("a") as handle:
        if new_file:
            handle.write("timestamp_utc\tsample\tstatus\tsample_root\n")
        handle.write(f"{utc_now()}\t{group}\t{status}\t{sample_root}\n")


def record_launch_command(out_root: Path, argv: Sequence[str]) -> None:
    text = " ".join(shlex.quote(part) for part in argv) + "\n"
    primary = out_root / "LAUNCH_COMMAND.txt"
    if not primary.exists():
        write_text(primary, text)
    else:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        write_text(out_root / f"LAUNCH_COMMAND_{stamp}.txt", text)
    log_path = out_root / "LAUNCH_COMMANDS.tsv"
    new_file = not log_path.exists()
    with log_path.open("a") as handle:
        if new_file:
            handle.write("timestamp_utc\tcommand\n")
        handle.write(f"{utc_now()}\t{text.strip()}\n")


def process_group(
    args: argparse.Namespace,
    group: str,
    rows: Dict[str, Dict[str, str]],
    run_tsv: Path,
    downstream_jobs: List[DownstreamJob],
) -> None:
    sample_root = args.out_root / "samples" / group
    run_dir = sample_root / "run"
    done_file = sample_root / "RUN_COMPLETE.ok"
    if done_file.exists() and args.skip_existing:
        log(f"{group}: already complete, skipping")
        append_status(run_tsv, group, "skipped_existing", sample_root)
        return

    sample_root.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    stage_dirs: Dict[str, Path] = {}
    for library in LIBRARY_ORDER:
        row = rows[library]
        paths = row_fastqs(row)
        stage_dir = sample_root / "staged_fastqs" / STAGE_SUBDIR[library]
        stage_fastqs(paths, stage_dir, args.force_stage)
        validate_pairs(stage_dir)
        stage_dirs[library] = stage_dir

    write_pf_config(group, rows, stage_dirs, sample_root / "pf_multi_config.csv")
    cmd = build_star_cmd(args, group, rows, sample_root)
    env = {}
    if args.low_mem:
        env["STAR_VELOCYTO_LOW_MEM"] = "1"
        env["STAR_SOLO_BINARY_SPOOL"] = "1"
    write_text(sample_root / "RUN_COMMAND.sh", command_script_text(cmd, env), executable=True)

    manifest_lines = [
        f"created_utc={utc_now()}",
        f"sample={group}",
        f"repo={REPO_ROOT}",
        f"manifest={args.manifest}",
        f"star_bin={args.star_bin}",
        f"genome_dir={args.genome_dir}",
        f"threads={args.threads}",
        f"cr_min_umi={args.cr_min_umi}",
        f"out_samtype={args.out_samtype}",
        f"run_downstream={int(args.run_downstream)}",
        f"downstream_async={int(args.downstream_async)}",
        f"downstream_max_active={args.downstream_max_active}",
        f"globus_enabled={int(globus_enabled(args))}",
    ]
    write_text(sample_root / "RUN_MANIFEST.txt", "\n".join(manifest_lines) + "\n")

    if args.dry_run:
        log(f"{group}: dry-run command written to {sample_root / 'RUN_COMMAND.sh'}")
        append_status(run_tsv, group, "dry_run", sample_root)
        return

    append_status(run_tsv, group, "star_start", sample_root)
    log(f"{group}: STAR start")
    run_subprocess(cmd, env)
    append_status(run_tsv, group, "star_done", sample_root)
    log(f"{group}: STAR complete")

    if args.run_downstream and args.downstream_async:
        downstream_jobs[:] = wait_for_downstream_slot(downstream_jobs, run_tsv, args.downstream_max_active)
    downstream_job = run_downstream(args, group, sample_root)
    if args.run_downstream and downstream_job is None:
        append_status(run_tsv, group, "downstream_done", sample_root)
    elif downstream_job is not None:
        downstream_jobs.append(downstream_job)
        append_status(run_tsv, group, "downstream_submitted", sample_root)

    transfer_bam(args, group, sample_root)
    if globus_enabled(args):
        append_status(run_tsv, group, "globus_done", sample_root)

    if downstream_job is None:
        write_text(done_file, f"completed_utc={utc_now()}\n")
        append_status(run_tsv, group, "done", sample_root)
    else:
        append_status(run_tsv, group, "stream_done_pending_downstream", sample_root)
    return


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run MSK 40KO STAR/CR-compatible production samples from manifest.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--samples", default="", help="Comma-separated provider_group values to run")
    parser.add_argument("--exclude-samples", default="", help="Comma-separated provider_group values to skip")
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--star-bin", type=Path, default=DEFAULT_STAR_BIN)
    parser.add_argument("--genome-dir", type=Path, default=DEFAULT_GENOME_DIR)
    parser.add_argument("--out-samtype", default="BAM Unsorted", help='STAR --outSAMtype value, e.g. "BAM Unsorted" or "None"')
    parser.add_argument("--cr-min-umi", type=int, default=2)
    parser.add_argument("--low-mem", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--force-stage", action="store_true", help="Replace existing staged FASTQ symlinks if they differ")
    parser.add_argument("--dry-run", action="store_true")

    parser.add_argument("--run-downstream", action="store_true")
    parser.add_argument("--downstream-wrapper", type=Path, default=DEFAULT_DOWNSTREAM)
    parser.add_argument("--downstream-output-name", default=DEFAULT_DOWNSTREAM_NAME)
    parser.add_argument("--remote-host", default="")
    parser.add_argument("--remote-root", type=Path, default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-cellbender", action="store_true")
    parser.add_argument("--cellbender-gpu", action="store_true")
    parser.add_argument("--cellbender-cpu-cores", type=int, default=8)
    parser.add_argument("--cellbender-layer", default="denoised")
    parser.add_argument(
        "--downstream-async",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Launch remote downstream in the background so local sample processing can continue",
    )
    parser.add_argument(
        "--downstream-wait-end",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="With --downstream-async, wait for submitted remote downstream jobs after local sample processing finishes",
    )
    parser.add_argument(
        "--downstream-max-active",
        type=int,
        default=4,
        help="Maximum concurrent async remote downstream jobs; 0 means no cap",
    )

    parser.add_argument("--globus-src-endpoint", default="")
    parser.add_argument("--globus-dst-endpoint", default="")
    parser.add_argument("--globus-dst-root", default="")
    parser.add_argument("--globus-poll-seconds", type=int, default=60)
    parser.add_argument("--globus-wait", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--bam-archive-root", type=Path, default=None)
    parser.add_argument("--delete-local-bam-after-transfer", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    require_file(args.manifest, "manifest")
    require_executable(args.star_bin, "STAR binary")
    require_dir(args.genome_dir, "STAR genomeDir")
    groups = read_manifest(args.manifest)
    selected = select_groups(list(groups), args.samples, args.exclude_samples)
    args.out_root.mkdir(parents=True, exist_ok=True)
    run_tsv = args.out_root / "RUNS.tsv"
    record_launch_command(args.out_root, sys.argv)
    log(f"Output root: {args.out_root}")
    log(f"Selected samples: {', '.join(selected)}")
    downstream_jobs: List[DownstreamJob] = []
    for group in selected:
        process_group(args, group, groups[group], run_tsv, downstream_jobs)
    if downstream_jobs and args.downstream_wait_end:
        wait_downstream_jobs(downstream_jobs, run_tsv)
    elif downstream_jobs:
        log(f"Submitted {len(downstream_jobs)} remote downstream jobs without waiting at end")
    log("All selected samples complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
