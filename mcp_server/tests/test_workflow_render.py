"""Tests for render_workflow_command."""

import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import render_workflow_command


RENDER_SCHEMA = {
    "id": "render_wf",
    "title": "Render Test Workflow",
    "summary": "Test rendering.",
    "entry_script": "scripts/entry.sh",
    "parameters": [
        {
            "name": "out_root",
            "cli_flag": "--out-root",
            "type": "directory",
            "category": "output",
            "is_output_root": True,
        },
        {
            "name": "dataset_root",
            "cli_flag": "--dataset-root",
            "type": "directory",
            "default": "/data/default",
            "category": "inputs",
        },
        {
            "name": "threads",
            "cli_flag": "--threads",
            "type": "int",
            "default": 24,
            "category": "performance",
        },
        {
            "name": "downsample_reads",
            "cli_flag": "--downsample-reads",
            "type": "int",
            "default": 0,
            "category": "performance",
            "skip_when_default": True,
        },
        {
            "name": "downsample_seed",
            "cli_flag": "--downsample-seed",
            "type": "int",
            "default": 1,
            "category": "performance",
            "env_var": "DOWNSAMPLE_SEED",
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
            "name": "star_only",
            "cli_flag": "--star-only",
            "type": "bool",
            "default": False,
            "category": "mode",
        },
        {
            "name": "dry_run",
            "cli_flag": "--dry-run",
            "type": "bool",
            "default": False,
            "category": "mode",
        },
        {
            "name": "globus_src",
            "cli_flag": "--globus-src-endpoint",
            "type": "string",
            "default": "",
            "category": "transfer",
        },
    ],
    "constraints": [],
    "rendering": {
        "bool_style": "flag_only",
        "omit_absent_optionals": True,
        "flag_order": [
            "out_root",
            "dataset_root",
            "threads",
            "downsample_reads",
            "downsample_seed",
            "samples",
            "all_samples",
            "globus_src",
            "star_only",
            "dry_run",
        ],
    },
}


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def render_env():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        schema_rel = "workflows/render_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(RENDER_SCHEMA, f)

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
                    "id": "render_wf",
                    "title": "Render Test",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield tmp


class TestBasicRendering:
    def test_entry_script_is_first_argv(self, render_env):
        result = render_workflow_command("render_wf", {})
        assert result.argv[0] == str(render_env / "scripts" / "entry.sh")

    def test_defaults_rendered(self, render_env):
        result = render_workflow_command("render_wf", {})
        assert "--dataset-root" in result.argv
        assert "/data/default" in result.argv
        assert "--threads" in result.argv
        assert "24" in result.argv

    def test_workflow_id_in_response(self, render_env):
        result = render_workflow_command("render_wf", {})
        assert result.workflow_id == "render_wf"


class TestBooleanFlags:
    def test_true_bool_emits_flag(self, render_env):
        result = render_workflow_command("render_wf", {"dry_run": True})
        assert "--dry-run" in result.argv

    def test_false_bool_omitted(self, render_env):
        result = render_workflow_command("render_wf", {"dry_run": False})
        assert "--dry-run" not in result.argv

    def test_default_false_bool_omitted(self, render_env):
        result = render_workflow_command("render_wf", {})
        assert "--dry-run" not in result.argv
        assert "--star-only" not in result.argv
        assert "--all-samples" not in result.argv


class TestStringListParams:
    def test_csv_string(self, render_env):
        result = render_workflow_command("render_wf", {"samples": "s1,s2,s3"})
        idx = result.argv.index("--samples")
        assert result.argv[idx + 1] == "s1,s2,s3"

    def test_list_joined(self, render_env):
        result = render_workflow_command("render_wf", {"samples": ["s1", "s2"]})
        idx = result.argv.index("--samples")
        assert result.argv[idx + 1] == "s1,s2"


class TestIntParams:
    def test_int_rendered(self, render_env):
        result = render_workflow_command("render_wf", {"threads": 8})
        idx = result.argv.index("--threads")
        assert result.argv[idx + 1] == "8"

    def test_zero_downsample_omitted(self, render_env):
        """downsample_reads=0 should be omitted (per script behavior)."""
        result = render_workflow_command("render_wf", {"downsample_reads": 0})
        assert "--downsample-reads" not in result.argv

    def test_positive_downsample_emitted(self, render_env):
        result = render_workflow_command("render_wf", {"downsample_reads": 1000})
        assert "--downsample-reads" in result.argv
        idx = result.argv.index("--downsample-reads")
        assert result.argv[idx + 1] == "1000"


class TestEmptyStringParams:
    def test_empty_string_omitted(self, render_env):
        result = render_workflow_command("render_wf", {"globus_src": ""})
        assert "--globus-src-endpoint" not in result.argv


class TestOutputRoot:
    def test_output_root_tracked(self, render_env):
        result = render_workflow_command("render_wf", {"out_root": "/tmp/myout"})
        assert result.output_root == "/tmp/myout"

    def test_no_output_root(self, render_env):
        result = render_workflow_command("render_wf", {})
        assert result.output_root is None


class TestEnvOverrides:
    def test_downsample_seed_in_env(self, render_env):
        result = render_workflow_command("render_wf", {"downsample_seed": 42})
        assert result.env_overrides.get("DOWNSAMPLE_SEED") == "42"


class TestShellPreview:
    def test_shell_preview_is_string(self, render_env):
        result = render_workflow_command("render_wf", {"dry_run": True})
        assert isinstance(result.shell_preview, str)
        assert "entry.sh" in result.shell_preview
        assert "--dry-run" in result.shell_preview


class TestDeterministicOrder:
    def test_flag_order_matches_schema(self, render_env):
        result = render_workflow_command(
            "render_wf",
            {
                "out_root": "/tmp/out",
                "threads": 16,
                "dry_run": True,
                "all_samples": True,
            },
        )
        argv = result.argv[1:]  # skip entry script
        # Extract flag indices
        flags = [a for a in argv if a.startswith("--")]
        expected_order = ["--out-root", "--threads", "--all-samples", "--dry-run"]
        # All flags should appear in the expected relative order
        indices = [flags.index(f) for f in expected_order]
        assert indices == sorted(indices)


class TestDeclarativeEnvVar:
    """env_var field causes the value to be emitted as an env override."""

    def test_env_var_populated(self, render_env):
        result = render_workflow_command("render_wf", {"downsample_seed": 99})
        assert result.env_overrides["DOWNSAMPLE_SEED"] == "99"

    def test_env_var_not_set_when_absent(self, render_env):
        result = render_workflow_command("render_wf", {})
        # downsample_seed has a default of 1, so it IS set
        assert result.env_overrides["DOWNSAMPLE_SEED"] == "1"


class TestDeclarativeSkipWhenDefault:
    """skip_when_default field omits the flag when value equals default."""

    def test_default_value_omitted(self, render_env):
        result = render_workflow_command("render_wf", {"downsample_reads": 0})
        assert "--downsample-reads" not in result.argv

    def test_non_default_value_emitted(self, render_env):
        result = render_workflow_command("render_wf", {"downsample_reads": 500})
        assert "--downsample-reads" in result.argv
        idx = result.argv.index("--downsample-reads")
        assert result.argv[idx + 1] == "500"


class TestDeclarativeIsOutputRoot:
    """is_output_root field designates which param is the output root."""

    def test_output_root_from_schema_flag(self, render_env):
        result = render_workflow_command("render_wf", {"out_root": "/tmp/test_out"})
        assert result.output_root == "/tmp/test_out"

    def test_no_output_root_when_not_set(self, render_env):
        result = render_workflow_command("render_wf", {"threads": 8})
        assert result.output_root is None


OPERAND_SCHEMA = {
    "id": "operand_wf",
    "title": "Operand group test",
    "summary": "Test.",
    "entry_script": "bin/tool",
    "parameters": [
        {
            "name": "out_root",
            "cli_flag": "--out",
            "type": "string",
            "required": True,
        },
        {
            "name": "mate1",
            "cli_flag": "--readFilesIn",
            "type": "string",
            "required": True,
            "operand_group": "pe",
            "category": "reads",
        },
        {
            "name": "mate2",
            "cli_flag": "--readFilesIn",
            "type": "string",
            "required": True,
            "operand_group": "pe",
            "category": "reads",
        },
        {
            "name": "sam_kind",
            "cli_flag": "--outSAMtype",
            "type": "string",
            "required": True,
            "default": "BAM",
            "operand_group": "sam",
            "category": "sam",
        },
        {
            "name": "sam_sort",
            "cli_flag": "--outSAMtype",
            "type": "string",
            "required": True,
            "default": "SortedByCoordinate",
            "operand_group": "sam",
            "category": "sam",
        },
    ],
    "constraints": [],
    "rendering": {
        "omit_absent_optionals": True,
        "flag_order": ["out_root", "mate1", "mate2", "sam_kind", "sam_sort"],
    },
}


@pytest.fixture
def operand_render_env():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        schema_rel = "workflows/operand_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(OPERAND_SCHEMA, f)
        (tmp / "bin").mkdir()
        (tmp / "bin" / "tool").write_text("#!/bin/sh\necho\n")
        (tmp / "bin" / "tool").chmod(0o755)

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
                    "id": "operand_wf",
                    "title": "Operand",
                    "entry_script": "bin/tool",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield tmp


class TestOperandGroupRendering:
    def test_readfiles_emits_single_flag_two_operands(self, operand_render_env):
        result = render_workflow_command(
            "operand_wf",
            {
                "out_root": "/tmp/o",
                "mate1": "A_1.fq,B_1.fq",
                "mate2": "A_2.fq,B_2.fq",
            },
        )
        tool = str(operand_render_env / "bin" / "tool")
        assert result.argv[:4] == [tool, "--out", "/tmp/o", "--readFilesIn"]
        assert result.argv[4:6] == ["A_1.fq,B_1.fq", "A_2.fq,B_2.fq"]
        assert result.argv[6:9] == ["--outSAMtype", "BAM", "SortedByCoordinate"]

    def test_outsamtype_operand_group(self, operand_render_env):
        result = render_workflow_command(
            "operand_wf",
            {
                "out_root": "/x",
                "mate1": "m1",
                "mate2": "m2",
            },
        )
        assert "--outSAMtype" in result.argv
        i = result.argv.index("--outSAMtype")
        assert result.argv[i + 1 : i + 3] == ["BAM", "SortedByCoordinate"]


class TestUnknownWorkflow:
    def test_raises(self, render_env):
        with pytest.raises(ValueError, match="Unknown workflow"):
            render_workflow_command("nonexistent", {})
