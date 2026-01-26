"""Tests for preflight validation."""

from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from mcp_server.config import load_config
from mcp_server.schemas.run_config import RunConfig
from mcp_server.tools.preflight import (
    check_binaries,
    check_dataset,
    check_disk_space,
    check_fixtures,
    check_output_directory,
    check_paths_valid,
    check_script_allowed,
    check_script_path,
    run_preflight,
)
import mcp_server.config as config_module


@pytest.fixture
def preflight_config(temp_dir: Path) -> dict:
    """Create a config for preflight tests."""
    # Create test script
    scripts_dir = temp_dir / "tests"
    scripts_dir.mkdir(parents=True)
    test_script = scripts_dir / "test_script.sh"
    test_script.write_text("#!/bin/bash\necho 'test'\n")
    test_script.chmod(0o755)

    # Create dataset dir
    dataset_dir = temp_dir / "datasets" / "test"
    dataset_dir.mkdir(parents=True)

    # Create fixtures dir (one exists, one doesn't)
    fixtures_dir = temp_dir / "fixtures"
    fixtures_dir.mkdir(parents=True)
    (fixtures_dir / "existing").mkdir()

    return {
        "server": {"auth_token": ""},
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
            "temp_root": str(temp_dir / "tmp"),
        },
        "trusted_roots": [str(temp_dir), "/tmp", "/usr/bin"],
        "datasets": [
            {
                "id": "test-dataset",
                "path": str(dataset_dir),
                "description": "Test dataset",
            }
        ],
        "execution": {
            "default_timeout_seconds": 3600,
            "max_concurrent_jobs": 1,
            "queue_max_size": 5,
        },
        "disk_space": {
            "default_min_gb": 1,  # Low for testing
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
                "path": "tests/test_script.sh",
                "module": "test",
                "description": "Script with fixtures",
                "fixtures": [
                    str(fixtures_dir / "existing"),
                    str(fixtures_dir / "missing"),
                ],
            },
            {
                "name": "not_runnable",
                "path": "tests/test_script.sh",
                "module": "test",
                "description": "Not runnable",
                "runnable": False,
            },
            {
                "name": "build_core",
                "path": "make",
                "args": ["core"],
                "module": "build",
                "description": "Build command",
            },
        ],
        "test_suites": [
            {
                "module": "test",
                "description": "Test suite",
                "scripts": ["test_script", "script_with_fixtures"],
                "disk_space_min_gb": 1,
            }
        ],
        "required_binaries": [
            {"name": "echo", "paths": ["/bin/echo", "/usr/bin/echo"]},
            {"name": "nonexistent_binary_xyz", "paths": ["/nonexistent/path"]},
        ],
    }


@pytest.fixture
def setup_preflight_config(temp_dir: Path, preflight_config: dict):
    """Load preflight config."""
    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(preflight_config, f)

    config = load_config(config_path)
    yield config

    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


class TestCheckScriptAllowed:
    """Tests for script allowlist check."""

    def test_allowed_script(self, setup_preflight_config):
        """Test that allowed scripts pass."""
        run_config = RunConfig(script="test_script")
        result = check_script_allowed(run_config)
        assert result.passed
        assert result.details["script"] == "test_script"

    def test_unknown_script(self, setup_preflight_config):
        """Test that unknown scripts fail."""
        run_config = RunConfig(script="unknown_script")
        result = check_script_allowed(run_config)
        assert not result.passed
        assert "not in the allowlist" in result.message
        assert "unknown_script" in result.details["requested"]

    def test_not_runnable_script(self, setup_preflight_config):
        """Test that non-runnable scripts fail."""
        run_config = RunConfig(script="not_runnable")
        result = check_script_allowed(run_config)
        assert not result.passed
        assert "not currently runnable" in result.message


class TestCheckScriptPath:
    """Tests for script path check."""

    def test_existing_script(self, setup_preflight_config):
        """Test that existing script paths pass."""
        run_config = RunConfig(script="test_script")
        result = check_script_path(run_config)
        assert result.passed
        assert "test_script.sh" in result.details["path"]

    def test_build_command(self, setup_preflight_config):
        """Test that build commands check PATH."""
        run_config = RunConfig(script="build_core")
        result = check_script_path(run_config)
        # 'make' should be found on most systems
        if result.passed:
            assert result.details["type"] == "build_command"

    def test_directory_as_script_fails(self, setup_preflight_config, temp_dir: Path):
        """Test that a directory path fails the script check."""
        # Create a directory where the script should be
        script_dir = temp_dir / "tests" / "fake_dir_script"
        script_dir.mkdir(parents=True)

        # Add a script config pointing to the directory
        config = config_module._config
        from mcp_server.schemas.config import ScriptConfig
        config.scripts.append(
            ScriptConfig(
                name="dir_script",
                path="tests/fake_dir_script",
                module="test",
                description="Directory not a file",
            )
        )

        run_config = RunConfig(script="dir_script")
        result = check_script_path(run_config)
        assert not result.passed
        assert "not a file" in result.message


class TestCheckDataset:
    """Tests for dataset check."""

    def test_no_dataset_specified(self, setup_preflight_config):
        """Test that no dataset is OK."""
        run_config = RunConfig(script="test_script")
        result = check_dataset(run_config)
        assert result.passed

    def test_valid_dataset(self, setup_preflight_config):
        """Test that valid dataset passes."""
        run_config = RunConfig(script="test_script", dataset_id="test-dataset")
        result = check_dataset(run_config)
        assert result.passed
        assert result.details["dataset_id"] == "test-dataset"

    def test_unknown_dataset(self, setup_preflight_config):
        """Test that unknown dataset fails."""
        run_config = RunConfig(script="test_script", dataset_id="nonexistent")
        result = check_dataset(run_config)
        assert not result.passed
        assert "not found" in result.message


class TestCheckOutputDirectory:
    """Tests for output directory check."""

    def test_default_output_dir(self, setup_preflight_config, temp_dir: Path):
        """Test that default output dir passes."""
        run_config = RunConfig(script="test_script")
        result = check_output_directory(run_config)
        assert result.passed

    def test_valid_output_dir(self, setup_preflight_config, temp_dir: Path):
        """Test that valid output dir passes."""
        out_dir = temp_dir / "output"
        run_config = RunConfig(script="test_script", out_dir=str(out_dir))
        result = check_output_directory(run_config)
        assert result.passed

    def test_output_dir_outside_roots(self, setup_preflight_config):
        """Test that output dir outside trusted roots fails."""
        run_config = RunConfig(script="test_script", out_dir="/etc/output")
        result = check_output_directory(run_config)
        assert not result.passed
        assert "outside trusted roots" in result.message


class TestCheckBinaries:
    """Tests for binary check."""

    def test_some_binaries_missing(self, setup_preflight_config):
        """Test that missing binaries are reported."""
        run_config = RunConfig(script="test_script")
        result = check_binaries(run_config)
        # echo should be found, nonexistent_binary_xyz should not
        assert not result.passed
        assert "nonexistent_binary_xyz" in result.details["missing"]
        assert "echo" in result.details["found"]

    def test_build_scripts_skip_binary_check(self, setup_preflight_config):
        """Test that build module scripts skip the binary check."""
        run_config = RunConfig(script="build_core")
        result = check_binaries(run_config)
        # Build scripts should pass even if binaries are missing
        assert result.passed
        assert "Skipped for build scripts" in result.message


class TestCheckDiskSpace:
    """Tests for disk space check."""

    def test_sufficient_disk_space(self, setup_preflight_config, temp_dir: Path):
        """Test that sufficient disk space passes."""
        run_config = RunConfig(script="test_script", out_dir=str(temp_dir))
        result = check_disk_space(run_config)
        # Should pass since we set requirement to 1GB
        assert result.passed
        assert result.details["required_gb"] == 1

    def test_insufficient_disk_space(self, setup_preflight_config, temp_dir: Path):
        """Test that insufficient disk space fails."""
        run_config = RunConfig(script="test_script", out_dir=str(temp_dir))

        # Mock to return very low disk space
        with patch("mcp_server.tools.preflight.get_disk_space_gb", return_value=0.1):
            result = check_disk_space(run_config)
            assert not result.passed
            assert "Insufficient disk space" in result.message


class TestCheckFixtures:
    """Tests for fixture check."""

    def test_no_fixtures_required(self, setup_preflight_config):
        """Test that scripts without fixtures pass."""
        run_config = RunConfig(script="test_script")
        result = check_fixtures(run_config)
        assert result.passed

    def test_missing_fixtures(self, setup_preflight_config):
        """Test that missing fixtures fail."""
        run_config = RunConfig(script="script_with_fixtures")
        result = check_fixtures(run_config)
        assert not result.passed
        assert "missing" in result.details
        assert len(result.details["missing"]) == 1
        assert "missing" in result.details["missing"][0]


class TestCheckPathsValid:
    """Tests for path validation."""

    def test_valid_paths(self, setup_preflight_config, temp_dir: Path):
        """Test that valid paths pass."""
        run_config = RunConfig(script="test_script", out_dir=str(temp_dir / "output"))
        result = check_paths_valid(run_config)
        assert result.passed

    def test_invalid_paths(self, setup_preflight_config):
        """Test that invalid paths fail."""
        run_config = RunConfig(script="test_script", out_dir="/etc/bad_path")
        result = check_paths_valid(run_config)
        assert not result.passed
        assert "/etc/bad_path" in result.details["invalid_paths"]


class TestRunPreflight:
    """Tests for full preflight run."""

    def test_preflight_valid_config(self, setup_preflight_config, temp_dir: Path):
        """Test preflight with valid config."""
        # Remove the nonexistent binary requirement for this test
        config = config_module._config
        config.required_binaries = [b for b in config.required_binaries if b.name != "nonexistent_binary_xyz"]

        run_config = RunConfig(
            script="test_script",
            out_dir=str(temp_dir / "output"),
        )
        result = run_preflight(run_config)

        # Check that most checks passed
        passed_checks = [c for c in result.checks if c.passed]
        assert len(passed_checks) >= 6

    def test_preflight_invalid_script(self, setup_preflight_config):
        """Test preflight with invalid script."""
        run_config = RunConfig(script="nonexistent_script")
        result = run_preflight(run_config)
        assert not result.valid
        assert any("script_allowed" in c.name and not c.passed for c in result.checks)

    def test_preflight_returns_all_checks(self, setup_preflight_config):
        """Test that preflight returns all check types."""
        run_config = RunConfig(script="test_script")
        result = run_preflight(run_config)

        check_names = {c.name for c in result.checks}
        expected_checks = {
            "script_allowed",
            "script_path",
            "working_dir",
            "dataset",
            "output_directory",
            "binaries_present",
            "disk_space",
            "fixtures_present",
            "paths_valid",
        }
        assert check_names == expected_checks


# --- Trusted Root Enforcement Tests ---


class TestTrustedRootEnforcement:
    """Tests for trusted root enforcement across all path types."""

    @pytest.fixture
    def config_with_untrusted_paths(self, temp_dir: Path) -> dict:
        """Create a config with paths outside trusted roots."""
        # Create test script within trusted root
        scripts_dir = temp_dir / "tests"
        scripts_dir.mkdir(parents=True)
        test_script = scripts_dir / "test_script.sh"
        test_script.write_text("#!/bin/bash\necho 'test'\n")
        test_script.chmod(0o755)

        return {
            "server": {"auth_token": ""},
            "paths": {
                "repo_root": str(temp_dir),
                "artifact_log_root": str(temp_dir / "artifacts"),
                "temp_root": str(temp_dir / "tmp"),
            },
            # Intentionally limited trusted roots - only temp_dir
            "trusted_roots": [str(temp_dir)],
            "datasets": [
                {
                    "id": "trusted-dataset",
                    "path": str(temp_dir / "datasets" / "trusted"),
                    "description": "Dataset in trusted root",
                },
                {
                    "id": "untrusted-dataset",
                    "path": "/etc/untrusted_dataset",  # Outside trusted roots
                    "description": "Dataset outside trusted root",
                },
            ],
            "execution": {
                "default_timeout_seconds": 3600,
                "max_concurrent_jobs": 1,
                "queue_max_size": 5,
            },
            "disk_space": {"default_min_gb": 1, "full_depth_min_gb": 50},
            "scripts": [
                {
                    "name": "trusted_script",
                    "path": "tests/test_script.sh",
                    "module": "test",
                    "description": "Script in trusted root",
                },
                {
                    "name": "untrusted_script",
                    "path": "/etc/passwd",  # Outside trusted roots
                    "module": "test",
                    "description": "Script outside trusted root",
                },
                {
                    "name": "script_untrusted_fixtures",
                    "path": "tests/test_script.sh",
                    "module": "test",
                    "description": "Script with untrusted fixtures",
                    "fixtures": [
                        str(temp_dir / "fixtures" / "trusted"),
                        "/etc/untrusted_fixture",  # Outside trusted roots
                    ],
                },
            ],
            "test_suites": [],
            "required_binaries": [],
        }

    @pytest.fixture
    def setup_untrusted_config(self, temp_dir: Path, config_with_untrusted_paths: dict):
        """Load config with untrusted paths."""
        # Create the trusted dataset and fixture dirs
        (temp_dir / "datasets" / "trusted").mkdir(parents=True)
        (temp_dir / "fixtures" / "trusted").mkdir(parents=True)

        config_path = temp_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_with_untrusted_paths, f)

        config = load_config(config_path)
        yield config

        config_module._config = None
        config_module._config_path = None
        config_module._config_loaded_at = None

    def test_script_path_outside_trusted_roots_fails(self, setup_untrusted_config):
        """Test that script paths outside trusted roots fail."""
        run_config = RunConfig(script="untrusted_script")
        result = check_script_path(run_config)
        assert not result.passed
        assert "outside trusted roots" in result.message

    def test_script_path_inside_trusted_roots_passes(self, setup_untrusted_config):
        """Test that script paths inside trusted roots pass."""
        run_config = RunConfig(script="trusted_script")
        result = check_script_path(run_config)
        assert result.passed

    def test_dataset_path_outside_trusted_roots_fails(self, setup_untrusted_config):
        """Test that dataset paths outside trusted roots fail."""
        run_config = RunConfig(script="trusted_script", dataset_id="untrusted-dataset")
        result = check_dataset(run_config)
        assert not result.passed
        assert "outside trusted roots" in result.message

    def test_dataset_path_inside_trusted_roots_passes(self, setup_untrusted_config):
        """Test that dataset paths inside trusted roots pass."""
        run_config = RunConfig(script="trusted_script", dataset_id="trusted-dataset")
        result = check_dataset(run_config)
        assert result.passed

    def test_fixture_path_outside_trusted_roots_fails(self, setup_untrusted_config):
        """Test that fixture paths outside trusted roots fail."""
        run_config = RunConfig(script="script_untrusted_fixtures")
        result = check_fixtures(run_config)
        assert not result.passed
        assert "outside_roots" in result.details
        assert "/etc/untrusted_fixture" in result.details["outside_roots"]


class TestBinaryTrustedRoots:
    """Tests for binary trusted root enforcement."""

    @pytest.fixture
    def config_binary_roots(self, temp_dir: Path) -> dict:
        """Create a config for binary trusted root testing."""
        scripts_dir = temp_dir / "tests"
        scripts_dir.mkdir(parents=True)
        test_script = scripts_dir / "test_script.sh"
        test_script.write_text("#!/bin/bash\necho 'test'\n")
        test_script.chmod(0o755)

        return {
            "server": {"auth_token": ""},
            "paths": {
                "repo_root": str(temp_dir),
                "artifact_log_root": str(temp_dir / "artifacts"),
                "temp_root": str(temp_dir / "tmp"),
            },
            # Only temp_dir is trusted - /usr/bin is NOT trusted
            "trusted_roots": [str(temp_dir)],
            "datasets": [],
            "execution": {
                "default_timeout_seconds": 3600,
                "max_concurrent_jobs": 1,
                "queue_max_size": 5,
            },
            "disk_space": {"default_min_gb": 1, "full_depth_min_gb": 50},
            "scripts": [
                {
                    "name": "test_script",
                    "path": "tests/test_script.sh",
                    "module": "test",
                    "description": "Test script",
                },
                {
                    "name": "build_make",
                    "path": "make",
                    "args": ["core"],
                    "module": "build",
                    "description": "Make command",
                },
            ],
            "test_suites": [],
            # Binary in untrusted path
            "required_binaries": [
                {"name": "echo", "paths": ["/bin/echo", "/usr/bin/echo"]},
            ],
        }

    @pytest.fixture
    def setup_binary_config(self, temp_dir: Path, config_binary_roots: dict):
        """Load config for binary tests."""
        config_path = temp_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_binary_roots, f)

        config = load_config(config_path)
        yield config

        config_module._config = None
        config_module._config_path = None
        config_module._config_loaded_at = None

    def test_binary_outside_trusted_roots_not_found(self, setup_binary_config):
        """Test that binaries outside trusted roots are not found."""
        from mcp_server.tools.utils import find_binary

        # echo exists in /usr/bin but that's not a trusted root
        result = find_binary("echo", enforce_trusted_roots=True)
        # Should be None because /usr/bin is not in trusted roots
        assert result is None

    def test_binary_found_when_trusted_roots_not_enforced(self, setup_binary_config):
        """Test that binaries are found when trusted roots not enforced."""
        from mcp_server.tools.utils import find_binary

        # echo should be found when we don't enforce trusted roots
        result = find_binary("echo", enforce_trusted_roots=False)
        assert result is not None
        assert "echo" in str(result)

    def test_build_command_fails_when_binary_outside_trusted_roots(
        self, setup_binary_config
    ):
        """Test that build commands fail when binary is outside trusted roots."""
        run_config = RunConfig(script="build_make")
        result = check_script_path(run_config)
        # make is outside trusted roots (only temp_dir is trusted)
        assert not result.passed
        assert "not found in trusted roots" in result.message
