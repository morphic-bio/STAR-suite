"""Local, serialized recipe jobs. Only the native Multiome entry point can run."""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import uuid

from .multiome import WORKFLOW_IDS


class JobConflict(ValueError):
    pass


def now():
    return datetime.now(timezone.utc).isoformat()


def tail(path, size=16000):
    try:
        with Path(path).open("rb") as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - size))
            return stream.read(size).decode("utf-8", errors="replace")
    except OSError:
        return ""


class JobManager:
    def __init__(self):
        self.lock = threading.RLock()
        self.jobs = {}
        self.processes = {}
        self.watchers = {}
        self.active = None

    def _save(self, job):
        path = Path(job["record_path"])
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(job, indent=2) + "\n")
        tmp.replace(path)

    def _save_status(self, job):
        try:
            self._save(job)
        except OSError as exc:
            # Disk errors must not prevent cancellation or process cleanup.
            job["message"] = f"Could not update run record: {exc}"

    def start(self, spec, rendered, identity, artifact_root):
        root = Path(spec["source_repo"]).resolve()
        entry = Path(spec["entry_script"]).resolve()
        if spec["id"] not in WORKFLOW_IDS or Path(rendered["argv"][0]).resolve() != entry:
            raise ValueError("Only the bundled Multiome recipe supports Run.")
        params = rendered["normalized_params"]
        output = Path(params["out_dir"])
        with self.lock:
            if self.active:
                raise JobConflict("A Multiome job is already running. Wait for it or cancel it first.")
            # Reserve a new directory atomically; never reuse another run's files.
            try:
                output.mkdir(parents=True, exist_ok=False)
            except FileExistsError as exc:
                raise JobConflict("Choose a new output directory; that path already exists.") from exc
            job_id = uuid.uuid4().hex
            folder = Path(artifact_root) / "launchpad" / job_id
            job = None
            try:
                folder.mkdir(parents=True)
                env = {k: v for k, v in os.environ.items() if not k.startswith("STAR_MULTIOME_")}
                env.update(rendered.get("env_overrides", {}))
                # The legacy recipe variable must select the binary just verified.
                env["STAR_BIN"] = identity["binary"]
                argv = ["nice", "-n", "10", *rendered["argv"]]
                job = {"id": job_id, "recipe_id": spec["id"], "status": "running",
                       "started_at": now(), "finished_at": None, "exit_code": None,
                       "dry_run": bool(params.get("dry_run")), "output_dir": str(output),
                       "argv": argv, "params": params, "runtime": identity,
                       "source_metadata": spec.get("source_metadata", {}),
                       "log_path": str(folder / "run.log"), "record_path": str(folder / "job.json")}
                self._save(job)
                with Path(job["log_path"]).open("wb") as log:
                    process = subprocess.Popen(argv, cwd=root, env=env, stdin=subprocess.DEVNULL,
                                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            except Exception as exc:
                if job is not None:
                    job.update(status="failed", finished_at=now(), message=str(exc))
                    self._save_status(job)
                # No process was created; release only our empty reservation.
                output.rmdir()
                raise
            self.jobs[job_id] = job
            self.processes[job_id] = process
            self.active = job_id
            watcher = threading.Thread(target=self._watch, args=(job_id,), daemon=True)
            self.watchers[job_id] = watcher
            watcher.start()
            return dict(job)

    def _watch(self, job_id):
        process = self.processes[job_id]
        code = process.wait()
        with self.lock:
            job = self.jobs[job_id]
            cancelled = job["status"] == "cancelling"
            if cancelled:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            marker = "DRY_RUN_PREVIEW.txt" if job["dry_run"] else "LOCAL_MEX_READY.txt"
            complete = (Path(job["output_dir"]) / marker).is_file()
            job.update(exit_code=code, finished_at=now(),
                       status="cancelled" if cancelled else "succeeded" if code == 0 and complete else "failed")
            if code == 0 and not complete and not cancelled:
                job["message"] = f"Recipe exited without {marker}; outputs are incomplete."
            self._save_status(job)
            self.active = None

    def get(self, job_id):
        with self.lock:
            if job_id not in self.jobs:
                raise KeyError(job_id)
            job = dict(self.jobs[job_id])
        job["log_tail"] = tail(job["log_path"])
        logs = Path(job["output_dir"]) / "logs"
        job["stage_logs"] = {name: tail(logs / name, 4000) for name in (
            "star_multiome.stdout.log", "star_multiome.stderr.log",
            "package_star_genefull_mex.log", "prepare_velocyto_mex.log", "build_atac_peak_matrix.log")
            if (logs / name).is_file()}
        return job

    def cancel(self, job_id):
        with self.lock:
            job = self.jobs[job_id]
            if job["status"] != "running":
                return self.get(job_id)
            job["status"] = "cancelling"
            self._save_status(job)
            process = self.processes[job_id]
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            threading.Thread(target=self._kill_after_grace, args=(job_id, process), daemon=True).start()
        return self.get(job_id)

    def _kill_after_grace(self, job_id, process):
        # Descendants can outlive the shell. Keep the process group bounded even
        # when the shell exits on TERM before an aligner does.
        threading.Event().wait(3)
        with self.lock:
            if self.active == job_id and self.jobs[job_id]["status"] == "cancelling":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    def shutdown(self):
        with self.lock:
            active = self.active
        if active:
            self.cancel(active)
            self.processes[active].wait(timeout=5)
            self.watchers[active].join(timeout=5)
