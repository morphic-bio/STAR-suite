"""Pytest fixtures for MCP server tests."""

import os
import tempfile
from pathlib import Path
from typing import Generator

import pytest
import yaml

from mcp_server.config import load_config, _config, _config_path, _config_loaded_at
import mcp_server.config as config_module


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_config_dict(temp_dir: Path) -> dict:
    """Create a sample configuration dictionary."""
    return {
        "server": {
            "host": "127.0.0.1",
            "port": 9999,
            "auth_token": "test-token",
            "transport": "http",
        },
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
            "temp_root": str(temp_dir / "tmp"),
        },
        "trusted_roots": [
            str(temp_dir),
            "/tmp",
            "/usr/bin",
            "/bin",
        ],
        "datasets": [
            {
                "id": "test-dataset",
                "path": str(temp_dir / "datasets" / "test"),
                "description": "Test dataset",
            }
        ],
        "execution": {
            "default_timeout_seconds": 3600,
            "max_concurrent_jobs": 1,
            "queue_max_size": 5,
        },
        "disk_space": {
            "default_min_gb": 10,
            "full_depth_min_gb": 50,
        },
        "scripts": [
            {
                "name": "test_script",
                "path": "tests/test_script.sh",
                "module": "test",
                "description": "A test script",
            },
            {
                "name": "script_with_fixtures",
                "path": "tests/fixture_script.sh",
                "module": "test",
                "description": "Script with fixtures",
                "fixtures": [str(temp_dir / "fixtures" / "required")],
            },
        ],
        "test_suites": [
            {
                "module": "test",
                "description": "Test suite",
                "scripts": ["test_script", "script_with_fixtures"],
                "disk_space_min_gb": 10,
            }
        ],
        "required_binaries": [
            {
                "name": "echo",
                "paths": ["/bin/echo", "/usr/bin/echo"],
            }
        ],
    }


@pytest.fixture
def sample_config_file(temp_dir: Path, sample_config_dict: dict) -> Path:
    """Create a sample config file."""
    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(sample_config_dict, f)
    return config_path


@pytest.fixture
def loaded_config(sample_config_file: Path):
    """Load config and yield, then clean up."""
    config = load_config(sample_config_file)
    yield config
    # Reset global state
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


@pytest.fixture
def mock_repo_structure(temp_dir: Path) -> Path:
    """Create a mock repository structure."""
    # Create directories
    (temp_dir / "docs").mkdir(parents=True)
    (temp_dir / "tests").mkdir(parents=True)
    (temp_dir / "datasets" / "test").mkdir(parents=True)

    # Create some doc files
    (temp_dir / "docs" / "README.md").write_text("# Test Documentation\n\nThis is a test.")
    (temp_dir / "docs" / "CBUB_GUIDE.md").write_text("# CB/UB Guide\n\nHow to use CB/UB tags.")

    # Create a test script
    (temp_dir / "tests" / "test_script.sh").write_text("#!/bin/bash\necho 'test'\n")

    return temp_dir
