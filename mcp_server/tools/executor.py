"""Script execution for MCP server."""

import json
import os
import signal
import subprocess
import threading
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Optional

from ..config import get_config
from ..schemas.responses import (
    CollectOutputsResponse,
    ErrorResponse,
    OutputFile,
    RunCompletionResponse,
    RunScriptResponse,
)
from ..schemas.run_config import RunConfig
from .preflight import run_preflight
from .utils import ensure_run_log_dir, generate_run_id, get_run_log_dir


class RunStatus(str, Enum):
    """Status of a script run."""

    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    TIMEOUT = "timeout"
    CANCELLED = "cancelled"


@dataclass
class RunState:
    """State of a script run."""

    run_id: str
    run_config: RunConfig
    status: RunStatus
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    exit_code: Optional[int] = None
    process: Optional[subprocess.Popen] = None
    log_dir: Optional[Path] = None
    outputs: list[str] = field(default_factory=list)
    error_message: Optional[str] = None
    # Track actual output paths for collect_outputs
    output_roots: list[str] = field(default_factory=list)


class ExecutorQueue:
    """Manages script execution with queuing and true concurrency.
    
    Uses a ThreadPoolExecutor to run multiple jobs in parallel up to
    max_concurrent. The dispatcher thread monitors the queue and submits
    jobs to the pool as slots become available.
    """

    def __init__(
        self,
        max_concurrent: int = 1,
        max_queue_size: int = 10,
    ):
        # ThreadPoolExecutor requires at least 1 worker
        self.max_concurrent = max(1, max_concurrent)
        self.max_queue_size = max_queue_size

        self._lock = threading.Lock()
        self._queue: deque[RunState] = deque()
        self._running: dict[str, RunState] = {}
        self._completed: dict[str, RunState] = {}

        # Thread pool for parallel execution
        self._executor: Optional[ThreadPoolExecutor] = None
        
        # Background thread for dispatching from queue to pool
        self._stop_event = threading.Event()
        self._pause_event = threading.Event()  # For testing: pause dispatch
        self._dispatcher_thread: Optional[threading.Thread] = None

    def start(self):
        """Start the queue dispatcher and thread pool."""
        if self._executor is None:
            self._executor = ThreadPoolExecutor(
                max_workers=self.max_concurrent,
                thread_name_prefix="mcp_executor_"
            )
        
        if self._dispatcher_thread is None or not self._dispatcher_thread.is_alive():
            self._stop_event.clear()
            self._dispatcher_thread = threading.Thread(
                target=self._dispatch_queue, daemon=True, name="mcp_dispatcher"
            )
            self._dispatcher_thread.start()

    def stop(self):
        """Stop the queue dispatcher and thread pool."""
        self._stop_event.set()
        if self._dispatcher_thread:
            self._dispatcher_thread.join(timeout=5)
        if self._executor:
            self._executor.shutdown(wait=True, cancel_futures=False)
            self._executor = None

    def update_settings(self, max_concurrent: int, max_queue_size: int):
        """Update executor settings (called after config reload).
        
        Note: Changes to max_concurrent only affect new jobs. Running jobs
        continue until completion. The thread pool is resized on next start().
        """
        with self._lock:
            self.max_queue_size = max_queue_size
            if max_concurrent != self.max_concurrent:
                self.max_concurrent = max_concurrent
                # Recreate executor with new worker count
                if self._executor:
                    old_executor = self._executor
                    self._executor = ThreadPoolExecutor(
                        max_workers=max_concurrent,
                        thread_name_prefix="mcp_executor_"
                    )
                    # Shutdown old executor gracefully (running tasks continue)
                    old_executor.shutdown(wait=False)

    def submit(self, run_config: RunConfig) -> tuple[str, RunStatus]:
        """Submit a run to the queue.

        Args:
            run_config: Run configuration.

        Returns:
            Tuple of (run_id, initial_status).

        Raises:
            ValueError: If queue is full.
        """
        with self._lock:
            if len(self._queue) >= self.max_queue_size:
                raise ValueError(
                    f"Queue is full ({self.max_queue_size} jobs). Try again later."
                )

            run_id = generate_run_id()
            state = RunState(
                run_id=run_id,
                run_config=run_config,
                status=RunStatus.QUEUED,
            )
            self._queue.append(state)

            # Start processor if not running
            self.start()

            return run_id, RunStatus.QUEUED

    def get_state(self, run_id: str) -> Optional[RunState]:
        """Get the state of a run."""
        with self._lock:
            # Check running
            if run_id in self._running:
                return self._running[run_id]
            # Check completed
            if run_id in self._completed:
                return self._completed[run_id]
            # Check queue
            for state in self._queue:
                if state.run_id == run_id:
                    return state
            return None

    def cancel(self, run_id: str) -> bool:
        """Cancel a queued or running job.

        Returns:
            True if cancelled, False if not found or already completed.
        """
        with self._lock:
            # Check queue first
            for i, state in enumerate(self._queue):
                if state.run_id == run_id:
                    state.status = RunStatus.CANCELLED
                    state.completed_at = datetime.now(timezone.utc)
                    self._completed[run_id] = state
                    del self._queue[i]
                    return True

            # Check running - terminate process
            if run_id in self._running:
                state = self._running[run_id]
                if state.process:
                    state.process.terminate()
                    try:
                        state.process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        state.process.kill()
                state.status = RunStatus.CANCELLED
                state.completed_at = datetime.now(timezone.utc)
                self._completed[run_id] = state
                del self._running[run_id]
                return True

            return False

    def _dispatch_queue(self):
        """Background thread that dispatches jobs from queue to thread pool.
        
        This enables true parallel execution up to max_concurrent workers.
        """
        while not self._stop_event.is_set():
            # Check if paused (for testing)
            if self._pause_event.is_set():
                time.sleep(0.1)
                continue
                
            state_to_run: Optional[RunState] = None

            with self._lock:
                # Check if we can dispatch more jobs (pool slots available)
                if (
                    len(self._running) < self.max_concurrent
                    and len(self._queue) > 0
                    and self._executor is not None
                ):
                    state_to_run = self._queue.popleft()
                    state_to_run.status = RunStatus.RUNNING
                    state_to_run.started_at = datetime.now(timezone.utc)
                    self._running[state_to_run.run_id] = state_to_run

            if state_to_run:
                # Submit to thread pool for parallel execution
                self._executor.submit(self._execute_run, state_to_run)
            else:
                # Wait a bit before checking again
                time.sleep(0.1)
    
    def pause_dispatch(self):
        """Pause the dispatcher (for testing)."""
        self._pause_event.set()
    
    def resume_dispatch(self):
        """Resume the dispatcher (for testing)."""
        self._pause_event.clear()

    def _execute_run(self, state: RunState):
        """Execute a single run."""
        config = get_config()
        script_config = config.get_script(state.run_config.script)

        if not script_config:
            state.status = RunStatus.FAILED
            state.error_message = f"Script not found: {state.run_config.script}"
            state.completed_at = datetime.now(timezone.utc)
            with self._lock:
                if state.run_id in self._running:
                    del self._running[state.run_id]
                self._completed[state.run_id] = state
            return

        # Create log directory
        log_dir = ensure_run_log_dir(state.run_id)
        state.log_dir = log_dir

        # Determine output roots for collect_outputs
        # Always include the log directory, and out_dir if specified
        output_roots = [str(log_dir)]
        if state.run_config.out_dir:
            output_roots.append(state.run_config.out_dir)
        else:
            # Default output location when out_dir not specified
            default_out = config.paths.temp_root / "mcp_runs" / state.run_id
            output_roots.append(str(default_out))
        state.output_roots = output_roots

        # Write run.json (includes output_roots for later retrieval)
        run_json = {
            "run_id": state.run_id,
            "script": state.run_config.script,
            "args": state.run_config.args or [],
            "dataset_id": state.run_config.dataset_id,
            "out_dir": state.run_config.out_dir,
            "output_roots": output_roots,  # For collect_outputs when out_dir is None
            "env_overrides": state.run_config.env_overrides or {},
            "started_at": state.started_at.isoformat() if state.started_at else None,
        }
        with open(log_dir / "run.json", "w") as f:
            json.dump(run_json, f, indent=2)

        # Build command
        script_path = Path(script_config.path)
        if not script_path.is_absolute():
            script_path = config.paths.repo_root / script_path

        # For build commands, use the command directly
        if script_config.path in ("make", "cmake", "ninja"):
            cmd = [script_config.path]
            if script_config.args:
                cmd.extend(script_config.args)
        else:
            cmd = [str(script_path)]

        # Add user args
        if state.run_config.args:
            cmd.extend(state.run_config.args)

        # Determine working directory
        if script_config.working_dir:
            working_dir = Path(script_config.working_dir)
            if not working_dir.is_absolute():
                working_dir = config.paths.repo_root / working_dir
        else:
            working_dir = config.paths.repo_root

        # Build environment
        env = os.environ.copy()
        if state.run_config.env_overrides:
            env.update(state.run_config.env_overrides)

        # Determine timeout
        timeout = script_config.timeout_seconds or config.execution.default_timeout_seconds

        # Open log files
        stdout_path = log_dir / "stdout.log"
        stderr_path = log_dir / "stderr.log"
        combined_path = log_dir / "combined.log"

        try:
            with (
                open(stdout_path, "w") as stdout_file,
                open(stderr_path, "w") as stderr_file,
                open(combined_path, "w") as combined_file,
            ):
                # Start process in a new process group so we can kill all children
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    cwd=working_dir,
                    env=env,
                    text=True,
                    start_new_session=True,  # Creates new process group
                )
                state.process = process

                # Read output with timeout
                try:
                    stdout, stderr = process.communicate(timeout=timeout)

                    # Write to log files
                    stdout_file.write(stdout)
                    stderr_file.write(stderr)
                    combined_file.write(stdout)
                    if stderr:
                        combined_file.write("\n--- STDERR ---\n")
                        combined_file.write(stderr)

                    state.exit_code = process.returncode
                    if process.returncode == 0:
                        state.status = RunStatus.COMPLETED
                    else:
                        state.status = RunStatus.FAILED
                        state.error_message = f"Exit code: {process.returncode}"

                except subprocess.TimeoutExpired:
                    # Kill the entire process group (script + all children)
                    try:
                        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                        # Give processes a moment to terminate gracefully
                        try:
                            process.wait(timeout=2)
                        except subprocess.TimeoutExpired:
                            # Force kill if SIGTERM didn't work
                            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                    except (ProcessLookupError, OSError):
                        # Process already dead
                        pass

                    stdout, stderr = process.communicate()
                    stdout_file.write(stdout)
                    stderr_file.write(stderr)
                    combined_file.write(stdout)
                    if stderr:
                        combined_file.write("\n--- STDERR ---\n")
                        combined_file.write(stderr)

                    state.status = RunStatus.TIMEOUT
                    state.error_message = f"Timeout after {timeout} seconds"
                    state.exit_code = -1

        except Exception as e:
            state.status = RunStatus.FAILED
            state.error_message = str(e)
            state.exit_code = -1

        state.completed_at = datetime.now(timezone.utc)
        state.process = None

        # Write summary.json
        summary = {
            "run_id": state.run_id,
            "status": state.status.value,
            "exit_code": state.exit_code,
            "started_at": state.started_at.isoformat() if state.started_at else None,
            "completed_at": state.completed_at.isoformat() if state.completed_at else None,
            "duration_seconds": (
                (state.completed_at - state.started_at).total_seconds()
                if state.started_at and state.completed_at
                else None
            ),
            "error_message": state.error_message,
            "outputs": state.outputs,
        }
        with open(log_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # Move from running to completed
        with self._lock:
            if state.run_id in self._running:
                del self._running[state.run_id]
            self._completed[state.run_id] = state


# Global executor instance
_executor: Optional[ExecutorQueue] = None


def get_executor() -> ExecutorQueue:
    """Get or create the global executor instance."""
    global _executor
    if _executor is None:
        config = get_config()
        _executor = ExecutorQueue(
            max_concurrent=config.execution.max_concurrent_jobs,
            max_queue_size=config.execution.queue_max_size,
        )
        _executor.start()
    return _executor


def run_script(run_config: RunConfig) -> RunScriptResponse | ErrorResponse:
    """Submit a script for execution.

    Args:
        run_config: Run configuration.

    Returns:
        RunScriptResponse with run_id and status, or ErrorResponse on failure.
    """
    # Run preflight first
    preflight_result = run_preflight(run_config)
    if not preflight_result.valid:
        return ErrorResponse(
            code="PREFLIGHT_FAILED",
            message="Preflight validation failed",
            details={
                "errors": preflight_result.errors,
                "checks": [c.model_dump() for c in preflight_result.checks if not c.passed],
            },
        )

    # Submit to executor
    executor = get_executor()
    try:
        run_id, status = executor.submit(run_config)
    except ValueError as e:
        return ErrorResponse(
            code="QUEUE_FULL",
            message=str(e),
        )

    log_dir = get_run_log_dir(run_id)

    return RunScriptResponse(
        run_id=run_id,
        status=status.value,
        log_path=str(log_dir),
    )


def collect_outputs(run_id: str) -> CollectOutputsResponse | ErrorResponse:
    """Collect outputs and logs for a completed run.

    Args:
        run_id: Run identifier.

    Returns:
        CollectOutputsResponse with run details, or ErrorResponse on failure.
    """
    executor = get_executor()
    state = executor.get_state(run_id)

    if state is None:
        # Check if logs exist on disk (for runs from previous server instances)
        log_dir = get_run_log_dir(run_id)
        if log_dir.exists() and (log_dir / "summary.json").exists():
            return _collect_from_disk(run_id, log_dir)

        return ErrorResponse(
            code="RUN_NOT_FOUND",
            message=f"Run not found: {run_id}",
        )

    log_dir = state.log_dir or get_run_log_dir(run_id)

    # Build log files dict
    log_files = {}
    for name in ["stdout.log", "stderr.log", "combined.log"]:
        path = log_dir / name
        if path.exists():
            log_files[name.replace(".log", "")] = str(path)

    # Build outputs list from all output roots
    outputs: list[OutputFile] = []
    seen_paths: set[str] = set()  # Avoid duplicates
    
    # Collect from output_roots (includes log_dir and out_dir/default)
    for root in state.output_roots:
        root_path = Path(root)
        if root_path.exists():
            for f in root_path.rglob("*"):
                if f.is_file():
                    # Skip log files (already in log_files dict)
                    if f.suffix == ".log" and f.parent == log_dir:
                        continue
                    # Skip run metadata files
                    if f.name in ("run.json", "summary.json"):
                        continue
                    str_path = str(f)
                    if str_path in seen_paths:
                        continue
                    seen_paths.add(str_path)
                    try:
                        stat = f.stat()
                        outputs.append(
                            OutputFile(
                                path=str_path,
                                size_bytes=stat.st_size,
                                mtime=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc),
                            )
                        )
                    except OSError:
                        pass

    return CollectOutputsResponse(
        run_id=state.run_id,
        status=state.status.value,
        exit_code=state.exit_code,
        duration_seconds=(
            (state.completed_at - state.started_at).total_seconds()
            if state.started_at and state.completed_at
            else None
        ),
        started_at=state.started_at,
        completed_at=state.completed_at,
        script=state.run_config.script,
        args=state.run_config.args or [],
        log_files=log_files,
        outputs=outputs,
    )


def _collect_from_disk(run_id: str, log_dir: Path) -> CollectOutputsResponse | ErrorResponse:
    """Collect outputs from disk for a run not in memory."""
    try:
        with open(log_dir / "summary.json") as f:
            summary = json.load(f)

        with open(log_dir / "run.json") as f:
            run_json = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        return ErrorResponse(
            code="RUN_NOT_FOUND",
            message=f"Could not read run data: {e}",
        )

    # Build log files dict
    log_files = {}
    for name in ["stdout.log", "stderr.log", "combined.log"]:
        path = log_dir / name
        if path.exists():
            log_files[name.replace(".log", "")] = str(path)

    # Parse timestamps
    started_at = None
    completed_at = None
    if summary.get("started_at"):
        started_at = datetime.fromisoformat(summary["started_at"])
    if summary.get("completed_at"):
        completed_at = datetime.fromisoformat(summary["completed_at"])

    # Collect outputs from output_roots stored in run.json
    outputs: list[OutputFile] = []
    output_roots = run_json.get("output_roots", [])
    seen_paths: set[str] = set()
    
    for root in output_roots:
        root_path = Path(root)
        if root_path.exists():
            for f in root_path.rglob("*"):
                if f.is_file():
                    # Skip log files and metadata
                    if f.suffix == ".log" and f.parent == log_dir:
                        continue
                    if f.name in ("run.json", "summary.json"):
                        continue
                    str_path = str(f)
                    if str_path in seen_paths:
                        continue
                    seen_paths.add(str_path)
                    try:
                        stat = f.stat()
                        outputs.append(
                            OutputFile(
                                path=str_path,
                                size_bytes=stat.st_size,
                                mtime=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc),
                            )
                        )
                    except OSError:
                        pass

    return CollectOutputsResponse(
        run_id=run_id,
        status=summary.get("status", "unknown"),
        exit_code=summary.get("exit_code"),
        duration_seconds=summary.get("duration_seconds"),
        started_at=started_at,
        completed_at=completed_at,
        script=run_json.get("script", ""),
        args=run_json.get("args", []),
        log_files=log_files,
        outputs=outputs,
    )


def get_log_tail(run_id: str, lines: int = 50) -> Optional[str]:
    """Get the last N lines of the combined log.

    Args:
        run_id: Run identifier.
        lines: Number of lines to return.

    Returns:
        Last N lines of log, or None if not found.
    """
    log_dir = get_run_log_dir(run_id)
    combined_log = log_dir / "combined.log"

    if not combined_log.exists():
        return None

    try:
        with open(combined_log) as f:
            all_lines = f.readlines()
            return "".join(all_lines[-lines:])
    except OSError:
        return None
