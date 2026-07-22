"""Integration tests: validate -> render -> preflight composition."""

import os
import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.schemas.run_config import RunConfig
from mcp_server.tools.preflight import run_preflight
from mcp_server.tools.workflows import (
    validate_workflow_parameters,
    render_workflow_command,
)


INTEGRATION_SCHEMA = {
    "id": "int_wf",
    "title": "Integration Test Workflow",
    "summary": "Test the validate -> render flow.",
    "entry_script": "scripts/entry.sh",
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
            "name": "threads",
            "cli_flag": "--threads",
            "type": "int",
            "default": 4,
            "category": "performance",
        },
        {
            "name": "dry_run",
            "cli_flag": "--dry-run",
            "type": "bool",
            "default": False,
            "category": "mode",
        },
    ],
    "constraints": [
        {"kind": "positive", "params": ["threads"], "message": "Positive."},
    ],
    "rendering": {
        "flag_order": ["input_dir", "threads", "dry_run"],
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
def integration_env():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Create real directories/files the workflow needs
        input_dir = tmp / "data"
        input_dir.mkdir()
        scripts_dir = tmp / "scripts"
        scripts_dir.mkdir()
        entry_script = scripts_dir / "entry.sh"
        entry_script.write_text("#!/bin/bash\necho test\n")
        os.chmod(entry_script, 0o755)

        schema_rel = "workflows/int_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(INTEGRATION_SCHEMA, f)

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
                    "id": "int_wf",
                    "title": "Integration",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": schema_rel,
                }
            ],
            "scripts": [
                {
                    "name": "entry_script",
                    "path": "scripts/entry.sh",
                    "module": "test",
                    "description": "Test entry",
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
            "entry_script": str(entry_script),
        }


class TestValidateThenRender:
    def test_valid_params_render_successfully(self, integration_env):
        """validate -> render pipeline produces correct argv."""
        params = {
            "input_dir": integration_env["input_dir"],
            "threads": 8,
            "dry_run": True,
        }

        # Step 1: validate
        val_result = validate_workflow_parameters("int_wf", params)
        assert val_result.valid
        assert val_result.errors == []

        # Step 2: render using normalized params
        render_result = render_workflow_command("int_wf", val_result.normalized_params)
        assert render_result.entry_script == integration_env["entry_script"]
        assert "--input-dir" in render_result.argv
        assert integration_env["input_dir"] in render_result.argv
        assert "--threads" in render_result.argv
        assert "8" in render_result.argv
        assert "--dry-run" in render_result.argv

    def test_invalid_params_do_not_render(self, integration_env):
        """Invalid params should fail validation before rendering."""
        params = {
            "input_dir": integration_env["input_dir"],
            "threads": -1,
        }

        val_result = validate_workflow_parameters("int_wf", params, check_paths=False)
        assert not val_result.valid
        assert any("positive" in e for e in val_result.errors)

    def test_defaults_flow_through(self, integration_env):
        """Defaults are applied in validation and used in rendering."""
        params = {"input_dir": integration_env["input_dir"]}

        val_result = validate_workflow_parameters("int_wf", params)
        assert val_result.valid
        assert val_result.normalized_params["threads"] == 4
        assert val_result.normalized_params["dry_run"] is False

        render_result = render_workflow_command("int_wf", val_result.normalized_params)
        assert "--threads" in render_result.argv
        assert "4" in render_result.argv
        assert "--dry-run" not in render_result.argv

    def test_rendered_script_is_preflight_compatible(self, integration_env):
        """The rendered entry script path matches a configured script."""
        params = {"input_dir": integration_env["input_dir"], "dry_run": True}
        val_result = validate_workflow_parameters("int_wf", params)
        render_result = render_workflow_command("int_wf", val_result.normalized_params)

        # The rendered entry script should be the resolved path
        script_path = Path(render_result.entry_script)
        assert script_path.exists()
        assert script_path.name == "entry.sh"

    def test_validate_render_preflight_composition(self, integration_env):
        """Full validate -> render -> preflight pipeline passes."""
        params = {
            "input_dir": integration_env["input_dir"],
            "threads": 8,
            "dry_run": True,
        }

        # Step 1: validate
        val_result = validate_workflow_parameters("int_wf", params)
        assert val_result.valid, f"Validation errors: {val_result.errors}"

        # Step 2: render
        render_result = render_workflow_command("int_wf", val_result.normalized_params)
        assert render_result.entry_script == integration_env["entry_script"]

        # Step 3: preflight — map rendered output into RunConfig
        run_config = RunConfig(
            script="entry_script",
            args=render_result.argv[1:],
            out_dir=str(Path(integration_env["tmp"]) / "output"),
            env_overrides=render_result.env_overrides or None,
        )
        # Create the output dir so preflight can check it
        Path(run_config.out_dir).mkdir(parents=True, exist_ok=True)

        preflight_result = run_preflight(run_config)
        assert preflight_result.valid, (
            f"Preflight errors: {preflight_result.errors}\n"
            f"Checks: {[(c.name, c.passed, c.message) for c in preflight_result.checks]}"
        )
