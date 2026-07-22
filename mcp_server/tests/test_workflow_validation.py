"""Tests for validate_workflow_parameters."""

import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import validate_workflow_parameters


WORKFLOW_SCHEMA = {
    "id": "val_wf",
    "title": "Validation Test Workflow",
    "summary": "Test validation.",
    "entry_script": "scripts/test.sh",
    "parameters": [
        {
            "name": "input_dir",
            "cli_flag": "--input-dir",
            "type": "directory",
            "required": True,
            "path_must_exist": True,
            "category": "inputs",
        },
        {
            "name": "ref_file",
            "cli_flag": "--ref-file",
            "type": "file",
            "path_must_exist": True,
            "category": "inputs",
        },
        {
            "name": "threads",
            "cli_flag": "--threads",
            "type": "int",
            "default": 4,
            "category": "performance",
        },
        {
            "name": "ratio",
            "cli_flag": "--ratio",
            "type": "float",
            "default": 0.5,
            "min_value": 0,
            "max_value": 1,
            "category": "performance",
        },
        {
            "name": "mode",
            "cli_flag": "--mode",
            "type": "enum",
            "choices": ["fast", "normal", "slow"],
            "default": "normal",
            "category": "mode",
        },
        {
            "name": "samples",
            "cli_flag": "--samples",
            "type": "string_list",
            "category": "selection",
        },
        {
            "name": "all_samples",
            "cli_flag": "--all-samples",
            "type": "bool",
            "default": False,
            "category": "selection",
        },
        {
            "name": "dry_run",
            "cli_flag": "--dry-run",
            "type": "bool",
            "default": False,
            "category": "mode",
        },
        {
            "name": "poll_seconds",
            "cli_flag": "--poll-seconds",
            "type": "int",
            "default": 30,
            "category": "transfer",
        },
        {
            "name": "downsample",
            "cli_flag": "--downsample",
            "type": "int",
            "default": 0,
            "category": "performance",
        },
        {
            "name": "glob_src",
            "cli_flag": "--glob-src",
            "type": "string",
            "category": "transfer",
        },
        {
            "name": "glob_dst",
            "cli_flag": "--glob-dst",
            "type": "string",
            "category": "transfer",
        },
        {
            "name": "glob_root",
            "cli_flag": "--glob-root",
            "type": "string",
            "category": "transfer",
        },
    ],
    "constraints": [
        {
            "kind": "mutual_exclusion",
            "params": ["samples", "all_samples"],
            "message": "Mutually exclusive.",
        },
        {
            "kind": "group_required",
            "params": ["glob_src", "glob_dst", "glob_root"],
            "message": "All three Globus params needed.",
        },
        {
            "kind": "positive",
            "params": ["threads", "poll_seconds"],
            "message": "Must be positive.",
        },
        {
            "kind": "non_negative",
            "params": ["downsample"],
            "message": "Must be non-negative.",
        },
        {
            "kind": "dependency",
            "params": ["dry_run", "all_samples"],
            "message": "dry_run is independent.",
        },
    ],
    "rendering": {"flag_order": []},
}


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def workflow_env():
    """Set up a config with the validation test workflow."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        # Create real dirs/files for path checks
        input_dir = tmp / "data"
        input_dir.mkdir()
        ref_file = tmp / "ref.csv"
        ref_file.write_text("header\n")

        schema_rel = "workflows/val_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(WORKFLOW_SCHEMA, f)

        cfg = {
            "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
            "paths": {
                "repo_root": str(tmp),
                "artifact_log_root": str(tmp / "artifacts"),
                "temp_root": str(tmp / "tmp"),
            },
            "trusted_roots": [str(tmp), "/tmp"],
            "workflows": [
                {
                    "id": "val_wf",
                    "title": "Validation Test",
                    "entry_script": "scripts/test.sh",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield {
            "tmp": tmp,
            "input_dir": str(input_dir),
            "ref_file": str(ref_file),
        }


class TestRequiredParams:
    def test_missing_required_param_fails(self, workflow_env):
        result = validate_workflow_parameters("val_wf", {}, check_paths=False)
        assert not result.valid
        assert any("input_dir" in e for e in result.errors)

    def test_required_param_present_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"]},
        )
        assert result.valid


class TestUnknownParams:
    def test_unknown_param_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "bogus": "value"},
            check_paths=False,
        )
        assert not result.valid
        assert any("Unknown parameter" in e for e in result.errors)


class TestTypeChecks:
    def test_wrong_int_type(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "threads": "abc"},
            check_paths=False,
        )
        assert not result.valid
        assert any("integer" in e for e in result.errors)

    def test_wrong_float_type(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "ratio": "xyz"},
            check_paths=False,
        )
        assert not result.valid
        assert any("float" in e for e in result.errors)

    def test_valid_int_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "threads": 8},
            check_paths=False,
        )
        assert result.valid

    def test_string_int_coercion(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "threads": "8"},
            check_paths=False,
        )
        assert result.valid


class TestNumericRanges:
    def test_min_value_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "ratio": -0.1},
            check_paths=False,
        )
        assert not result.valid
        assert any("ratio" in e and ">= 0" in e for e in result.errors)

    def test_max_value_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "ratio": 1.1},
            check_paths=False,
        )
        assert not result.valid
        assert any("ratio" in e and "<= 1" in e for e in result.errors)

    def test_non_finite_value_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "ratio": "nan"},
            check_paths=False,
        )
        assert not result.valid
        assert any("ratio" in e and "finite" in e for e in result.errors)


class TestEnumValidation:
    def test_invalid_enum_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "mode": "turbo"},
            check_paths=False,
        )
        assert not result.valid
        assert any("mode" in e and "turbo" in e for e in result.errors)

    def test_valid_enum_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "mode": "fast"},
            check_paths=False,
        )
        assert result.valid


class TestPathValidation:
    def test_existing_dir_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf", {"input_dir": workflow_env["input_dir"]}
        )
        assert result.valid

    def test_missing_dir_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf", {"input_dir": str(Path(workflow_env["tmp"]) / "nope")}
        )
        assert not result.valid
        assert any("does not exist" in e for e in result.errors)

    def test_file_param_not_file_fails(self, workflow_env):
        # Pass a directory where a file is expected
        result = validate_workflow_parameters(
            "val_wf",
            {
                "input_dir": workflow_env["input_dir"],
                "ref_file": workflow_env["input_dir"],  # dir, not file
            },
        )
        assert not result.valid
        assert any("not a file" in e for e in result.errors)

    def test_path_outside_trusted_roots_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": "/etc/passwd"},
        )
        assert not result.valid
        assert any("trusted roots" in e for e in result.errors)

    def test_skip_path_check(self, workflow_env):
        """check_paths=False skips path existence validation."""
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": "/some/fake/path"},
            check_paths=False,
        )
        assert result.valid


class TestMutualExclusion:
    def test_both_present_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {
                "input_dir": workflow_env["input_dir"],
                "samples": "s1,s2",
                "all_samples": True,
            },
            check_paths=False,
        )
        assert not result.valid
        assert any("Mutually exclusive" in e for e in result.errors)

    def test_one_present_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "samples": "s1,s2"},
            check_paths=False,
        )
        assert result.valid


class TestAtLeastOne:
    def test_missing_group_fails(self, tmp_path):
        schema = {
            "id": "oneof_wf",
            "title": "One of test",
            "summary": "Validation test.",
            "entry_script": "scripts/test.sh",
            "parameters": [
                {"name": "read_files_cbq", "cli_flag": "--readFilesIn", "type": "string"},
                {"name": "read_files_manifest", "cli_flag": "--readFilesManifest", "type": "file"},
            ],
            "constraints": [
                {
                    "kind": "at_least_one",
                    "params": ["read_files_cbq", "read_files_manifest"],
                    "message": "Provide CBQ list or manifest.",
                }
            ],
            "rendering": {"flag_order": []},
        }
        (tmp_path / "workflows").mkdir()
        (tmp_path / "scripts").mkdir()
        (tmp_path / "scripts" / "test.sh").write_text("#!/bin/bash\n")
        schema_rel = "workflows/oneof_wf.yaml"
        (tmp_path / schema_rel).write_text(yaml.dump(schema))
        cfg = {
            "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
            "paths": {
                "repo_root": str(tmp_path),
                "artifact_log_root": str(tmp_path / "artifacts"),
                "temp_root": str(tmp_path / "tmp"),
            },
            "trusted_roots": [str(tmp_path), "/tmp"],
            "workflows": [
                {
                    "id": "oneof_wf",
                    "title": "One of test",
                    "entry_script": "scripts/test.sh",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        load_config(cfg_path)

        missing = validate_workflow_parameters("oneof_wf", {}, check_paths=False)
        assert not missing.valid
        assert any("At least one" in e for e in missing.errors)

        present = validate_workflow_parameters(
            "oneof_wf",
            {"read_files_cbq": "sampleA.cbq"},
            check_paths=False,
        )
        assert present.valid


class TestGroupRequired:
    def test_partial_group_warns(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "glob_src": "uuid1"},
            check_paths=False,
        )
        assert result.valid  # warning, not error
        assert len(result.warnings) > 0
        assert any("Incomplete" in w for w in result.warnings)

    def test_full_group_no_warning(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {
                "input_dir": workflow_env["input_dir"],
                "glob_src": "a",
                "glob_dst": "b",
                "glob_root": "c",
            },
            check_paths=False,
        )
        assert result.valid
        assert not any("Incomplete" in w for w in result.warnings)


class TestPositiveNonNegative:
    def test_zero_threads_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "threads": 0},
            check_paths=False,
        )
        assert not result.valid
        assert any("positive" in e for e in result.errors)

    def test_negative_downsample_fails(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "downsample": -1},
            check_paths=False,
        )
        assert not result.valid
        assert any("non-negative" in e for e in result.errors)

    def test_zero_downsample_passes(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"], "downsample": 0},
            check_paths=False,
        )
        assert result.valid


class TestNormalization:
    def test_defaults_applied(self, workflow_env):
        result = validate_workflow_parameters(
            "val_wf",
            {"input_dir": workflow_env["input_dir"]},
            check_paths=False,
        )
        assert result.valid
        assert result.normalized_params["threads"] == 4
        assert result.normalized_params["mode"] == "normal"
        assert result.normalized_params["dry_run"] is False


class TestUnknownWorkflow:
    def test_raises(self, workflow_env):
        with pytest.raises(ValueError, match="Unknown workflow"):
            validate_workflow_parameters("nonexistent", {})
