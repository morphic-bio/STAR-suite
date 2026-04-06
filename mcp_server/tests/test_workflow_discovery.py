"""Tests for workflow discovery tools (list_workflows, describe_workflow, get_workflow_parameter_schema)."""

import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import (
    describe_workflow,
    get_workflow_parameter_schema,
    list_workflows,
)


SAMPLE_WORKFLOW_SCHEMA = {
    "id": "test_wf",
    "title": "Test Workflow",
    "summary": "A test.",
    "kind": "shell_workflow",
    "entry_script": "scripts/test.sh",
    "supported_modes": ["local", "dry-run-capable"],
    "caveats": ["Test caveat"],
    "default_output_layout": "/tmp/out",
    "stages": [
        {
            "name": "stage1",
            "title": "Stage One",
            "script": "scripts/stage1.sh",
            "description": "First stage",
        },
        {
            "name": "stage2",
            "title": "Stage Two",
            "gated_by": "skip_stage2",
            "description": "Gated stage",
        },
    ],
    "parameter_groups": [
        {"name": "inputs", "title": "Inputs", "parameters": ["input_dir", "threads"]},
        {"name": "mode", "title": "Mode", "parameters": ["dry_run"]},
    ],
    "parameters": [
        {
            "name": "input_dir",
            "cli_flag": "--input-dir",
            "type": "directory",
            "required": True,
            "description": "Input directory",
            "category": "inputs",
            "path_must_exist": True,
        },
        {
            "name": "threads",
            "cli_flag": "--threads",
            "type": "int",
            "default": 4,
            "description": "Thread count",
            "category": "inputs",
        },
        {
            "name": "dry_run",
            "cli_flag": "--dry-run",
            "type": "bool",
            "default": False,
            "description": "Dry run mode",
            "category": "mode",
        },
    ],
    "constraints": [
        {"kind": "positive", "params": ["threads"], "message": "Must be positive."},
    ],
    "rendering": {"bool_style": "flag_only", "flag_order": ["input_dir", "threads", "dry_run"]},
}


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def loaded_workflow_config():
    """Load a config with a test workflow schema."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        schema_rel = "workflows/test_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(SAMPLE_WORKFLOW_SCHEMA, f)

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
                    "id": "test_wf",
                    "title": "Test Workflow",
                    "summary": "A test.",
                    "entry_script": "scripts/test.sh",
                    "kind": "shell_workflow",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield tmp


class TestListWorkflows:
    def test_returns_configured_workflows(self, loaded_workflow_config):
        result = list_workflows()
        assert len(result.workflows) == 1
        wf = result.workflows[0]
        assert wf.id == "test_wf"
        assert wf.title == "Test Workflow"
        assert wf.kind == "shell_workflow"
        assert "local" in wf.supported_modes

    def test_empty_when_no_workflows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cfg = {
                "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
                "paths": {
                    "repo_root": str(tmp),
                    "artifact_log_root": str(tmp / "artifacts"),
                    "temp_root": str(tmp / "tmp"),
                },
                "trusted_roots": [str(tmp)],
            }
            cfg_path = tmp / "config.yaml"
            with open(cfg_path, "w") as f:
                yaml.dump(cfg, f)
            load_config(cfg_path)

            result = list_workflows()
            assert result.workflows == []


class TestDescribeWorkflow:
    def test_returns_full_metadata(self, loaded_workflow_config):
        result = describe_workflow("test_wf")
        assert result.id == "test_wf"
        assert result.title == "Test Workflow"
        assert result.summary == "A test."
        assert len(result.stages) == 2
        assert result.stages[0].name == "stage1"
        assert result.stages[0].script == "scripts/stage1.sh"
        assert result.stages[1].gated_by == "skip_stage2"
        assert result.stages[1].script is None
        assert len(result.parameter_groups) == 2
        assert result.caveats == ["Test caveat"]

    def test_parameter_groups_are_typed(self, loaded_workflow_config):
        """parameter_groups should be ParameterGroupInfo, not raw dicts."""
        result = describe_workflow("test_wf")
        grp = result.parameter_groups[0]
        assert grp.name == "inputs"
        assert grp.title == "Inputs"
        assert grp.parameters == ["input_dir", "threads"]

    def test_unknown_workflow_raises(self, loaded_workflow_config):
        with pytest.raises(ValueError, match="Unknown workflow"):
            describe_workflow("nonexistent")


class TestGetWorkflowParameterSchema:
    def test_returns_parameters(self, loaded_workflow_config):
        result = get_workflow_parameter_schema("test_wf")
        assert result.workflow_id == "test_wf"
        assert len(result.parameters) == 3
        names = [p.name for p in result.parameters]
        assert "input_dir" in names
        assert "threads" in names
        assert "dry_run" in names

    def test_parameter_metadata_fields(self, loaded_workflow_config):
        """category, stage, source, repeatable, path_must_exist are populated."""
        result = get_workflow_parameter_schema("test_wf")
        by_name = {p.name: p for p in result.parameters}

        inp = by_name["input_dir"]
        assert inp.category == "inputs"
        assert inp.stage == "top_level"
        assert inp.source == "workflow_wrapper"
        assert inp.repeatable is False
        assert inp.path_must_exist is True
        assert inp.type == "directory"
        assert inp.required is True

        thr = by_name["threads"]
        assert thr.category == "inputs"
        assert thr.default == 4
        assert thr.repeatable is False

        dr = by_name["dry_run"]
        assert dr.category == "mode"
        assert dr.type == "bool"
        assert dr.default is False

    def test_returns_groups(self, loaded_workflow_config):
        result = get_workflow_parameter_schema("test_wf")
        assert len(result.parameter_groups) == 2
        assert result.parameter_groups[0].name == "inputs"
        assert result.parameter_groups[0].title == "Inputs"
        assert result.parameter_groups[0].parameters == ["input_dir", "threads"]

    def test_returns_constraints(self, loaded_workflow_config):
        result = get_workflow_parameter_schema("test_wf")
        assert len(result.constraints) == 1
        c = result.constraints[0]
        assert c.kind == "positive"
        assert c.params == ["threads"]
        assert c.message == "Must be positive."

    def test_unknown_workflow_raises(self, loaded_workflow_config):
        with pytest.raises(ValueError, match="Unknown workflow"):
            get_workflow_parameter_schema("missing")
