"""Tests for script execution."""

import time
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from mcp_server.config import load_config
from mcp_server.schemas.run_config import RunConfig
from mcp_server.tools.executor import (
    ExecutorQueue,
    RunStatus,
    collect_outputs,
    get_log_tail,
    run_script,
)
import mcp_server.config as config_module
import mcp_server.tools.executor as executor_module


@pytest.fixture
def executor_config(temp_dir: Path) -> dict:
    """Create a config for executor tests."""
    # Create a simple test script
    scripts_dir = temp_dir / "tests"
    scripts_dir.mkdir(parents=True)

    # Script that succeeds
    success_script = scripts_dir / "success.sh"
    success_script.write_text("#!/bin/bash\necho 'Hello World'\nexit 0\n")
    success_script.chmod(0o755)

    # Script that fails
    fail_script = scripts_dir / "fail.sh"
    fail_script.write_text("#!/bin/bash\necho 'Error' >&2\nexit 1\n")
    fail_script.chmod(0o755)

    # Script that times out
    slow_script = scripts_dir / "slow.sh"
    slow_script.write_text("#!/bin/bash\nsleep 60\n")
    slow_script.chmod(0o755)

    # Script that outputs to a file
    output_script = scripts_dir / "output.sh"
    output_script.write_text(
        '#!/bin/bash\necho "result" > "${OUT_DIR}/result.txt"\n'
    )
    output_script.chmod(0o755)

    return {
        "server": {"auth_token": ""},
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
            "temp_root": str(temp_dir / "tmp"),
        },
        "trusted_roots": [str(temp_dir), "/tmp", "/usr/bin", "/bin"],
        "datasets": [],
        "execution": {
            "default_timeout_seconds": 5,  # Short timeout for tests
            "max_concurrent_jobs": 1,
            "queue_max_size": 3,
        },
        "disk_space": {"default_min_gb": 1, "full_depth_min_gb": 50},
        "scripts": [
            {
                "name": "success",
                "path": "tests/success.sh",
                "module": "test",
                "description": "Script that succeeds",
            },
            {
                "name": "fail",
                "path": "tests/fail.sh",
                "module": "test",
                "description": "Script that fails",
            },
            {
                "name": "slow",
                "path": "tests/slow.sh",
                "module": "test",
                "description": "Script that times out",
                "timeout_seconds": 1,  # Very short timeout
            },
            {
                "name": "output",
                "path": "tests/output.sh",
                "module": "test",
                "description": "Script that outputs files",
            },
        ],
        "test_suites": [
            {
                "module": "test",
                "description": "Test suite",
                "scripts": ["success", "fail", "slow", "output"],
                "disk_space_min_gb": 1,
            }
        ],
        "required_binaries": [],
    }


@pytest.fixture
def setup_executor_config(temp_dir: Path, executor_config: dict):
    """Load executor config."""
    # Cleanup any existing executor first
    if executor_module._executor:
        executor_module._executor.stop()
        executor_module._executor = None

    # Create artifact directory
    (temp_dir / "artifacts").mkdir(parents=True)
    (temp_dir / "tmp").mkdir(parents=True)

    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(executor_config, f)

    config = load_config(config_path)

    # Reset the global executor
    executor_module._executor = None

    yield config

    # Cleanup
    if executor_module._executor:
        executor_module._executor.stop()
        executor_module._executor = None

    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


class TestExecutorQueue:
    """Tests for the ExecutorQueue class."""

    def test_submit_returns_run_id(self, setup_executor_config):
        """Test that submit returns a run_id."""
        queue = ExecutorQueue(max_concurrent=1, max_queue_size=5)
        queue.start()

        try:
            run_config = RunConfig(script="success")
            run_id, status = queue.submit(run_config)

            assert run_id is not None
            assert len(run_id) > 0
            assert status == RunStatus.QUEUED
        finally:
            queue.stop()

    def test_queue_full_raises(self, setup_executor_config):
        """Test that submitting to a full queue raises."""
        queue = ExecutorQueue(max_concurrent=1, max_queue_size=2)
        queue.pause_dispatch()  # Prevent jobs from running
        queue.start()

        try:
            # Fill the queue (jobs stay queued since dispatch is paused)
            queue.submit(RunConfig(script="slow"))
            queue.submit(RunConfig(script="slow"))

            # Should raise
            with pytest.raises(ValueError, match="Queue is full"):
                queue.submit(RunConfig(script="slow"))
        finally:
            queue.resume_dispatch()
            queue.stop()

    def test_get_state_returns_state(self, setup_executor_config):
        """Test that get_state returns the run state."""
        queue = ExecutorQueue(max_concurrent=1, max_queue_size=5)
        queue.start()

        try:
            run_config = RunConfig(script="success")
            run_id, _ = queue.submit(run_config)

            state = queue.get_state(run_id)
            assert state is not None
            assert state.run_id == run_id
        finally:
            queue.stop()

    def test_cancel_queued_job(self, setup_executor_config):
        """Test cancelling a queued job."""
        queue = ExecutorQueue(max_concurrent=1, max_queue_size=5)
        queue.pause_dispatch()  # Prevent jobs from running
        queue.start()

        try:
            run_config = RunConfig(script="slow")
            run_id, _ = queue.submit(run_config)

            # Cancel it (stays queued since dispatch is paused)
            result = queue.cancel(run_id)
            assert result is True

            # Check state
            state = queue.get_state(run_id)
            assert state.status == RunStatus.CANCELLED
        finally:
            queue.resume_dispatch()
            queue.stop()


class TestRunScript:
    """Tests for run_script function."""

    def test_run_script_success(self, setup_executor_config, temp_dir: Path):
        """Test running a successful script."""
        run_config = RunConfig(script="success")
        result = run_script(run_config)

        # Should return a response, not an error
        assert hasattr(result, "run_id") or "run_id" in result.model_dump()

        run_id = result.run_id if hasattr(result, "run_id") else result.model_dump()["run_id"]

        # Wait for completion
        time.sleep(1)

        # Check status
        output = collect_outputs(run_id)
        assert output.status in ("completed", "running", "queued")

    def test_run_script_preflight_failure(self, setup_executor_config):
        """Test that preflight failures are reported."""
        run_config = RunConfig(script="nonexistent_script")
        result = run_script(run_config)

        assert hasattr(result, "code") or result.model_dump().get("error")
        if hasattr(result, "code"):
            assert result.code == "PREFLIGHT_FAILED"

    def test_run_script_creates_logs(self, setup_executor_config, temp_dir: Path):
        """Test that running a script creates log files."""
        run_config = RunConfig(script="success")
        result = run_script(run_config)

        run_id = result.run_id

        # Wait for completion
        time.sleep(1)

        # Check log directory exists
        log_dir = temp_dir / "artifacts" / f"mcp_runs_{run_id.split('_')[0]}" / run_id
        assert log_dir.exists()
        assert (log_dir / "run.json").exists()


class TestCollectOutputs:
    """Tests for collect_outputs function."""

    def test_collect_outputs_not_found(self, setup_executor_config):
        """Test collecting outputs for non-existent run."""
        result = collect_outputs("nonexistent_run_id")

        assert hasattr(result, "code")
        assert result.code == "RUN_NOT_FOUND"

    def test_collect_outputs_after_completion(self, setup_executor_config, temp_dir: Path):
        """Test collecting outputs after script completes."""
        run_config = RunConfig(script="success")
        result = run_script(run_config)
        run_id = result.run_id

        # Wait for completion
        time.sleep(2)

        # Collect outputs
        output = collect_outputs(run_id)

        assert output.run_id == run_id
        assert output.status in ("completed", "running")
        assert output.script == "success"


class TestTimeout:
    """Tests for script timeout handling."""

    def test_script_timeout(self, setup_executor_config, temp_dir: Path):
        """Test that slow scripts are timed out."""
        run_config = RunConfig(script="slow")
        result = run_script(run_config)
        run_id = result.run_id

        # Poll for completion (script has 1s timeout)
        # Use short sleep intervals for faster test
        final_status = None
        for _ in range(30):  # Max 15 seconds
            time.sleep(0.5)
            output = collect_outputs(run_id)
            final_status = output.status
            if final_status not in ("running", "queued"):
                break

        # Check status
        assert final_status == "timeout", f"Expected timeout, got {final_status}"


class TestGetLogTail:
    """Tests for get_log_tail function."""

    def test_get_log_tail_not_found(self, setup_executor_config):
        """Test getting log tail for non-existent run."""
        result = get_log_tail("nonexistent_run_id")
        assert result is None

    def test_get_log_tail_after_run(self, setup_executor_config, temp_dir: Path):
        """Test getting log tail after script runs."""
        run_config = RunConfig(script="success")
        result = run_script(run_config)
        run_id = result.run_id

        # Wait for completion
        time.sleep(2)

        # Get log tail
        tail = get_log_tail(run_id)
        assert tail is not None
        assert "Hello World" in tail
