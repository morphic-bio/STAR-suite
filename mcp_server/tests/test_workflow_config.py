"""Tests for workflow config loading and schema resolution."""

import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import get_workflow_schema, get_workflow_schemas, load_config


MINIMAL_WORKFLOW_SCHEMA = {
    "id": "test_workflow",
    "title": "Test Workflow",
    "summary": "A test workflow.",
    "entry_script": "scripts/test.sh",
    "parameters": [],
}


def _make_config(temp_dir, workflows=None, workflow_schemas=None):
    """Build a config dict + write schema files, return config path."""
    cfg = {
        "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
            "temp_root": str(temp_dir / "tmp"),
        },
        "trusted_roots": [str(temp_dir), "/tmp"],
    }

    if workflows:
        cfg["workflows"] = workflows

    # Write schema files
    if workflow_schemas:
        for rel_path, schema_data in workflow_schemas.items():
            full_path = temp_dir / rel_path
            full_path.parent.mkdir(parents=True, exist_ok=True)
            with open(full_path, "w") as f:
                yaml.dump(schema_data, f)

    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(cfg, f)
    return config_path


@pytest.fixture(autouse=True)
def _reset_config():
    """Reset global config state after each test."""
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


class TestWorkflowConfigLoading:
    def test_config_loads_without_workflows(self):
        """Config loads fine when no workflows section is present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cfg_path = _make_config(tmp)
            config = load_config(cfg_path)
            assert config.workflows == []
            assert get_workflow_schemas() == {}

    def test_config_loads_with_workflow(self):
        """Config loads and resolves a workflow schema file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            schema_rel = "workflows/test.yaml"
            cfg_path = _make_config(
                tmp,
                workflows=[
                    {
                        "id": "test_workflow",
                        "title": "Test",
                        "summary": "Test wf",
                        "entry_script": "scripts/test.sh",
                        "kind": "shell_workflow",
                        "schema_file": schema_rel,
                    }
                ],
                workflow_schemas={schema_rel: MINIMAL_WORKFLOW_SCHEMA},
            )
            config = load_config(cfg_path)

            assert len(config.workflows) == 1
            assert config.workflows[0].id == "test_workflow"

            schemas = get_workflow_schemas()
            assert "test_workflow" in schemas
            assert schemas["test_workflow"].title == "Test Workflow"

    def test_get_workflow_schema_by_id(self):
        """get_workflow_schema returns the right schema."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            schema_rel = "workflows/test.yaml"
            cfg_path = _make_config(
                tmp,
                workflows=[
                    {
                        "id": "my_wf",
                        "title": "My WF",
                        "entry_script": "scripts/x.sh",
                        "schema_file": schema_rel,
                    }
                ],
                workflow_schemas={
                    schema_rel: {**MINIMAL_WORKFLOW_SCHEMA, "id": "my_wf"}
                },
            )
            load_config(cfg_path)

            schema = get_workflow_schema("my_wf")
            assert schema is not None
            assert schema.id == "my_wf"

            assert get_workflow_schema("nonexistent") is None

    def test_missing_schema_file_raises(self):
        """Config loading fails if schema_file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cfg_path = _make_config(
                tmp,
                workflows=[
                    {
                        "id": "bad",
                        "title": "Bad",
                        "entry_script": "x.sh",
                        "schema_file": "workflows/nonexistent.yaml",
                    }
                ],
            )
            with pytest.raises(FileNotFoundError, match="nonexistent.yaml"):
                load_config(cfg_path)

    def test_invalid_workflow_config_rejected(self):
        """WorkflowConfig requires id, entry_script, schema_file."""
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
                "workflows": [{"id": "incomplete"}],  # missing required fields
            }
            cfg_path = tmp / "config.yaml"
            with open(cfg_path, "w") as f:
                yaml.dump(cfg, f)

            with pytest.raises(Exception):
                load_config(cfg_path)

    def test_get_workflow_on_mcpconfig(self):
        """MCPConfig.get_workflow() helper works."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            schema_rel = "workflows/test.yaml"
            cfg_path = _make_config(
                tmp,
                workflows=[
                    {
                        "id": "wf1",
                        "title": "WF1",
                        "entry_script": "s.sh",
                        "schema_file": schema_rel,
                    }
                ],
                workflow_schemas={
                    schema_rel: {**MINIMAL_WORKFLOW_SCHEMA, "id": "wf1"}
                },
            )
            config = load_config(cfg_path)
            assert config.get_workflow("wf1") is not None
            assert config.get_workflow("wf1").title == "WF1"
            assert config.get_workflow("missing") is None
