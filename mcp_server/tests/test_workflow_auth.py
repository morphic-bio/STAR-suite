"""Auth/HTTP-wrapper tests for workflow and scaffold tools in app.py.

Mirrors test_auth.py coverage: verifies that each app.py tool wrapper
enforces auth correctly and returns structured responses.
"""

import os
import shutil
import tempfile
from pathlib import Path
from textwrap import dedent

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.app import (
    list_workflows as list_workflows_tool,
    describe_workflow as describe_workflow_tool,
    get_workflow_scripts as get_workflow_scripts_tool,
    get_workflow_parameter_schema as get_workflow_parameter_schema_tool,
    validate_workflow_parameters as validate_workflow_parameters_tool,
    render_workflow_command as render_workflow_command_tool,
    scaffold_workflow_schema as scaffold_workflow_schema_tool,
    validate_draft_workflow_schema as validate_draft_workflow_schema_tool,
)


# Access underlying functions from the FunctionTool wrappers
def list_workflows(**kw):
    return list_workflows_tool.fn(**kw)

def describe_workflow(**kw):
    return describe_workflow_tool.fn(**kw)

def get_workflow_parameter_schema(**kw):
    return get_workflow_parameter_schema_tool.fn(**kw)

def validate_workflow_parameters(**kw):
    return validate_workflow_parameters_tool.fn(**kw)

def render_workflow_command(**kw):
    return render_workflow_command_tool.fn(**kw)

def get_workflow_scripts(**kw):
    return get_workflow_scripts_tool.fn(**kw)

def scaffold_workflow_schema(**kw):
    return scaffold_workflow_schema_tool.fn(**kw)

def validate_draft_workflow_schema(**kw):
    return validate_draft_workflow_schema_tool.fn(**kw)


# Minimal workflow schema for testing
MINI_SCHEMA = {
    "id": "auth_wf",
    "title": "Auth Test Workflow",
    "summary": "Workflow for auth tests.",
    "entry_script": "scripts/entry.sh",
    "parameters": [
        {
            "name": "input_dir",
            "cli_flag": "--input-dir",
            "type": "directory",
            "required": True,
        },
        {
            "name": "dry_run",
            "cli_flag": "--dry-run",
            "type": "bool",
            "default": False,
        },
    ],
    "rendering": {"flag_order": ["input_dir", "dry_run"]},
}

PRIVATE_SCHEMA = {
    "id": "private_wf",
    "title": "Private Workflow",
    "summary": "Only visible with auth.",
    "entry_script": "scripts/entry.sh",
    "parameters": [
        {
            "name": "input_dir",
            "cli_flag": "--input-dir",
            "type": "directory",
            "required": True,
        },
    ],
    "rendering": {"flag_order": ["input_dir"]},
}


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def auth_env():
    """Environment with auth enabled and a workflow configured."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Create entry script
        scripts_dir = tmp / "scripts"
        scripts_dir.mkdir()
        entry = scripts_dir / "entry.sh"
        entry.write_text("#!/bin/bash\necho stub\n")
        os.chmod(entry, 0o755)

        # Create a scaffoldable script
        scaffold_script = scripts_dir / "run_scaffold_target.sh"
        scaffold_script.write_text(dedent("""\
            #!/bin/bash
            INPUT=""
            usage() { cat <<EOF
            Options:
              --input DIR    Input directory
            EOF
            }
            while [[ $# -gt 0 ]]; do
              case "$1" in
                --input) INPUT="$2"; shift 2 ;;
                *) exit 1 ;;
              esac
            done
        """))
        os.chmod(scaffold_script, 0o755)

        # Create input dir
        input_dir = tmp / "data"
        input_dir.mkdir()

        # Write workflow schemas
        schema_dir = tmp / "workflows"
        schema_dir.mkdir()
        with open(schema_dir / "auth_wf.yaml", "w") as f:
            yaml.dump(MINI_SCHEMA, f)
        with open(schema_dir / "private_wf.yaml", "w") as f:
            yaml.dump(PRIVATE_SCHEMA, f)

        cfg = {
            "server": {
                "auth_token": "secret-token",
                "public_discovery": True,
            },
            "paths": {
                "repo_root": str(tmp),
                "artifact_log_root": str(tmp / "artifacts"),
                "temp_root": str(tmp / "tmp"),
            },
            "trusted_roots": [str(tmp), "/tmp"],
            "workflows": [
                {
                    "id": "auth_wf",
                    "title": "Auth Test Workflow",
                    "summary": "For auth tests.",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": "workflows/auth_wf.yaml",
                },
                {
                    "id": "private_wf",
                    "title": "Private Workflow",
                    "summary": "Only visible with auth.",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": "workflows/private_wf.yaml",
                    "visibility": "private",
                },
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield {
            "tmp": tmp,
            "input_dir": str(input_dir),
            "scaffold_script": "scripts/run_scaffold_target.sh",
        }


# ---------------------------------------------------------------------------
# Discovery tools: public by default (is_discovery=True)
# ---------------------------------------------------------------------------


class TestDiscoveryToolsPublicAccess:
    """Discovery tools should work without auth when public_discovery=True."""

    def test_list_workflows_no_auth(self, auth_env):
        result = list_workflows()
        assert "error" not in result or result.get("error") is not True
        assert "workflows" in result

    def test_describe_workflow_no_auth(self, auth_env):
        result = describe_workflow(workflow_id="auth_wf")
        assert "error" not in result or result.get("error") is not True
        assert result["id"] == "auth_wf"

    def test_get_workflow_scripts_no_auth(self, auth_env):
        result = get_workflow_scripts(workflow_id="auth_wf")
        assert "error" not in result or result.get("error") is not True
        assert result["workflow_id"] == "auth_wf"

    def test_get_workflow_parameter_schema_no_auth(self, auth_env):
        result = get_workflow_parameter_schema(workflow_id="auth_wf")
        assert "error" not in result or result.get("error") is not True
        assert result["workflow_id"] == "auth_wf"
        assert len(result["parameters"]) == 2

    def test_scaffold_workflow_schema_no_auth(self, auth_env):
        result = scaffold_workflow_schema(
            script_path=auth_env["scaffold_script"]
        )
        assert "error" not in result or result.get("error") is not True
        assert "draft_yaml" in result

    def test_validate_draft_no_auth(self, auth_env):
        draft = 'id: "x"\ntitle: "X"\nsummary: "X."\nentry_script: "x.sh"\n'
        result = validate_draft_workflow_schema(draft_yaml=draft)
        assert "error" not in result or result.get("error") is not True
        assert "valid" in result


# ---------------------------------------------------------------------------
# Mutating/sensitive tools: require auth
# ---------------------------------------------------------------------------


class TestAuthRequiredTools:
    """validate_workflow_parameters and render_workflow_command require auth."""

    def test_validate_requires_auth(self, auth_env):
        result = validate_workflow_parameters(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"]},
        )
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"

    def test_validate_with_correct_token(self, auth_env):
        result = validate_workflow_parameters(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"]},
            auth_token="secret-token",
        )
        assert result.get("error") is not True
        assert "valid" in result

    def test_validate_with_wrong_token(self, auth_env):
        result = validate_workflow_parameters(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"]},
            auth_token="wrong-token",
        )
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"

    def test_render_requires_auth(self, auth_env):
        result = render_workflow_command(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"], "dry_run": True},
        )
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"

    def test_render_with_correct_token(self, auth_env):
        result = render_workflow_command(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"], "dry_run": True},
            auth_token="secret-token",
        )
        assert result.get("error") is not True
        assert "argv" in result
        assert "--dry-run" in result["argv"]

    def test_render_with_wrong_token(self, auth_env):
        result = render_workflow_command(
            workflow_id="auth_wf",
            params={"input_dir": auth_env["input_dir"]},
            auth_token="wrong-token",
        )
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"


# ---------------------------------------------------------------------------
# Discovery tools: auth required when public_discovery=False
# ---------------------------------------------------------------------------


class TestDiscoveryToolsWithAuthRequired:
    """When public_discovery=False, discovery tools also require auth."""

    @pytest.fixture(autouse=True)
    def _disable_public(self, auth_env):
        config_module._config.server.public_discovery = False
        yield
        # auth_env fixture reset handles cleanup

    def test_list_workflows_requires_auth(self, auth_env):
        result = list_workflows()
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"

    def test_list_workflows_with_token(self, auth_env):
        result = list_workflows(auth_token="secret-token")
        assert "error" not in result or result.get("error") is not True
        assert "workflows" in result

    def test_describe_requires_auth(self, auth_env):
        result = describe_workflow(workflow_id="auth_wf")
        assert result.get("error") is True

    def test_describe_with_token(self, auth_env):
        result = describe_workflow(workflow_id="auth_wf", auth_token="secret-token")
        assert result.get("error") is not True
        assert result["id"] == "auth_wf"

    def test_schema_requires_auth(self, auth_env):
        result = get_workflow_parameter_schema(workflow_id="auth_wf")
        assert result.get("error") is True

    def test_scaffold_requires_auth(self, auth_env):
        result = scaffold_workflow_schema(
            script_path=auth_env["scaffold_script"]
        )
        assert result.get("error") is True

    def test_validate_draft_requires_auth(self, auth_env):
        result = validate_draft_workflow_schema(draft_yaml="id: x\n")
        assert result.get("error") is True


# ---------------------------------------------------------------------------
# Error handling: tools return ErrorResponse on bad input
# ---------------------------------------------------------------------------


class TestToolErrorHandling:
    """Tools should return structured ErrorResponse on errors, not crash."""

    def test_describe_unknown_workflow(self, auth_env):
        result = describe_workflow(workflow_id="nonexistent")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_scripts_unknown_workflow(self, auth_env):
        result = get_workflow_scripts(workflow_id="nonexistent")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_schema_unknown_workflow(self, auth_env):
        result = get_workflow_parameter_schema(workflow_id="nonexistent")
        assert result.get("error") is True

    def test_validate_unknown_workflow(self, auth_env):
        result = validate_workflow_parameters(
            workflow_id="nonexistent",
            params={},
            auth_token="secret-token",
        )
        assert result.get("error") is True

    def test_render_unknown_workflow(self, auth_env):
        result = render_workflow_command(
            workflow_id="nonexistent",
            params={},
            auth_token="secret-token",
        )
        assert result.get("error") is True

    def test_scaffold_missing_script(self, auth_env):
        result = scaffold_workflow_schema(script_path="scripts/nonexistent.sh")
        assert result.get("error") is True
        assert result["code"] == "SCAFFOLD_FAILED"

    def test_scaffold_outside_trusted_roots(self, auth_env):
        result = scaffold_workflow_schema(script_path="/etc/passwd")
        assert result.get("error") is True
        assert result["code"] == "SCAFFOLD_FAILED"


# ---------------------------------------------------------------------------
# Visibility: content-level filtering by workflow visibility + auth
# ---------------------------------------------------------------------------


class TestWorkflowVisibility:
    """Private workflows are hidden from unauthenticated callers."""

    def test_list_omits_private_without_auth(self, auth_env):
        result = list_workflows()
        ids = [w["id"] for w in result["workflows"]]
        assert "auth_wf" in ids
        assert "private_wf" not in ids

    def test_list_includes_private_with_auth(self, auth_env):
        result = list_workflows(auth_token="secret-token")
        ids = [w["id"] for w in result["workflows"]]
        assert "auth_wf" in ids
        assert "private_wf" in ids

    def test_describe_private_without_auth_returns_error(self, auth_env):
        result = describe_workflow(workflow_id="private_wf")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_describe_private_with_auth_succeeds(self, auth_env):
        result = describe_workflow(workflow_id="private_wf", auth_token="secret-token")
        assert result.get("error") is not True
        assert result["id"] == "private_wf"

    def test_scripts_private_without_auth_returns_error(self, auth_env):
        result = get_workflow_scripts(workflow_id="private_wf")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_scripts_private_with_auth_succeeds(self, auth_env):
        result = get_workflow_scripts(
            workflow_id="private_wf", auth_token="secret-token"
        )
        assert result.get("error") is not True
        assert result["workflow_id"] == "private_wf"

    def test_scripts_public_redacts_host_metadata(self, auth_env):
        """Public discovery of scripts must not expose host-specific metadata."""
        result = get_workflow_scripts(workflow_id="auth_wf")
        assert result.get("error") is not True
        # absolute_path must be None for public callers
        for s in result["scripts"]:
            assert s["absolute_path"] is None
        # provenance must not contain repo_root, git_commit, git_remote
        prov = result["provenance"]
        assert "repo_root" not in prov
        assert "git_commit" not in prov
        assert "git_remote" not in prov

    def test_scripts_authenticated_includes_host_metadata(self, auth_env):
        """Authenticated callers get full provenance and absolute paths."""
        result = get_workflow_scripts(
            workflow_id="auth_wf", auth_token="secret-token"
        )
        assert result.get("error") is not True
        for s in result["scripts"]:
            assert s["absolute_path"] is not None
        assert "repo_root" in result["provenance"]

    def test_schema_private_without_auth_returns_error(self, auth_env):
        result = get_workflow_parameter_schema(workflow_id="private_wf")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_schema_private_with_auth_succeeds(self, auth_env):
        result = get_workflow_parameter_schema(
            workflow_id="private_wf", auth_token="secret-token"
        )
        assert result.get("error") is not True
        assert result["workflow_id"] == "private_wf"

    def test_public_workflow_visible_without_auth(self, auth_env):
        """Public workflows remain visible without auth (no regression)."""
        result = describe_workflow(workflow_id="auth_wf")
        assert result.get("error") is not True
        assert result["id"] == "auth_wf"
