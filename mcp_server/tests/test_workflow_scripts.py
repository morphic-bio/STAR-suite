"""Tests for get_workflow_scripts tool."""

import os
import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import get_workflow_scripts
from mcp_server.app import (
    get_workflow_scripts as get_workflow_scripts_tool,
)


# Access the underlying function from the FunctionTool wrapper
def get_workflow_scripts_app(**kw):
    return get_workflow_scripts_tool.fn(**kw)


SCHEMA_WITH_SCRIPTS = {
    "id": "wf_scripts",
    "title": "Script Test Workflow",
    "summary": "Workflow with scripts section.",
    "entry_script": "scripts/entry.sh",
    "scripts": [
        {
            "role": "entry",
            "path": "scripts/entry.sh",
            "description": "Top-level orchestrator.",
            "language": "bash",
        },
        {
            "role": "batch_runner",
            "path": "scripts/batch.sh",
            "description": "Per-sample batch runner.",
            "language": "bash",
        },
        {
            "role": "helper",
            "path": "scripts/helper.py",
            "description": "Python helper.",
            "language": "python",
        },
    ],
    "parameters": [],
    "rendering": {"flag_order": []},
}

SCHEMA_WITHOUT_SCRIPTS = {
    "id": "wf_no_scripts",
    "title": "No Scripts Workflow",
    "summary": "Workflow without explicit scripts section.",
    "entry_script": "scripts/entry.sh",
    "parameters": [],
    "rendering": {"flag_order": []},
}

PRIVATE_SCHEMA_WITH_SCRIPTS = {
    "id": "wf_private_scripts",
    "title": "Private Scripts Workflow",
    "summary": "Private workflow with scripts.",
    "entry_script": "scripts/entry.sh",
    "scripts": [
        {
            "role": "entry",
            "path": "scripts/entry.sh",
            "description": "Entry.",
            "language": "bash",
        },
    ],
    "parameters": [],
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
def scripts_env():
    """Environment with workflows that have scripts sections."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Create script files so exists=True
        scripts_dir = tmp / "scripts"
        scripts_dir.mkdir()
        (scripts_dir / "entry.sh").write_text("#!/bin/bash\necho entry\n")
        os.chmod(scripts_dir / "entry.sh", 0o755)
        (scripts_dir / "batch.sh").write_text("#!/bin/bash\necho batch\n")
        os.chmod(scripts_dir / "batch.sh", 0o755)
        (scripts_dir / "helper.py").write_text("#!/usr/bin/env python3\nprint('helper')\n")
        os.chmod(scripts_dir / "helper.py", 0o755)

        # Write workflow schemas
        wf_dir = tmp / "workflows"
        wf_dir.mkdir()
        with open(wf_dir / "wf_scripts.yaml", "w") as f:
            yaml.dump(SCHEMA_WITH_SCRIPTS, f)
        with open(wf_dir / "wf_no_scripts.yaml", "w") as f:
            yaml.dump(SCHEMA_WITHOUT_SCRIPTS, f)
        with open(wf_dir / "wf_private_scripts.yaml", "w") as f:
            yaml.dump(PRIVATE_SCHEMA_WITH_SCRIPTS, f)

        cfg = {
            "server": {
                "auth_token": "test-token",
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
                    "id": "wf_scripts",
                    "title": "Script Test Workflow",
                    "summary": "With scripts.",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": "workflows/wf_scripts.yaml",
                },
                {
                    "id": "wf_no_scripts",
                    "title": "No Scripts Workflow",
                    "summary": "Without scripts.",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": "workflows/wf_no_scripts.yaml",
                },
                {
                    "id": "wf_private_scripts",
                    "title": "Private Scripts Workflow",
                    "summary": "Private with scripts.",
                    "entry_script": "scripts/entry.sh",
                    "schema_file": "workflows/wf_private_scripts.yaml",
                    "visibility": "private",
                },
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield tmp


class TestGetWorkflowScripts:
    """Tests for the get_workflow_scripts tool function."""

    def test_returns_all_scripts(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        assert result.workflow_id == "wf_scripts"
        assert result.title == "Script Test Workflow"
        assert result.entry_script == "scripts/entry.sh"
        assert len(result.scripts) == 3

    def test_script_roles(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        roles = [s.role for s in result.scripts]
        assert "entry" in roles
        assert "batch_runner" in roles
        assert "helper" in roles

    def test_script_paths_and_absolute_paths(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        entry = [s for s in result.scripts if s.role == "entry"][0]
        assert entry.path == "scripts/entry.sh"
        assert entry.absolute_path == str(scripts_env / "scripts" / "entry.sh")

    def test_script_languages(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        by_role = {s.role: s for s in result.scripts}
        assert by_role["entry"].language == "bash"
        assert by_role["helper"].language == "python"

    def test_script_exists_flag(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        for s in result.scripts:
            assert s.exists is True

    def test_missing_script_exists_false(self, scripts_env):
        """If a declared script doesn't exist on disk, exists=False."""
        # Remove one script file
        (scripts_env / "scripts" / "batch.sh").unlink()
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        by_role = {s.role: s for s in result.scripts}
        assert by_role["batch_runner"].exists is False
        assert by_role["entry"].exists is True

    def test_fallback_entry_script_when_no_scripts_section(self, scripts_env):
        """Workflows without a scripts section still return the entry script."""
        result = get_workflow_scripts("wf_no_scripts", authenticated=True)
        assert len(result.scripts) == 1
        assert result.scripts[0].role == "entry"
        assert result.scripts[0].path == "scripts/entry.sh"

    def test_provenance_includes_repo_root(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        assert result.provenance["repo_root"] == str(scripts_env)

    def test_provenance_includes_schema_file(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        assert result.provenance["workflow_schema"] == "workflows/wf_scripts.yaml"

    def test_unknown_workflow_raises(self, scripts_env):
        with pytest.raises(ValueError, match="Unknown workflow"):
            get_workflow_scripts("nonexistent", authenticated=True)

    def test_private_workflow_hidden_without_auth(self, scripts_env):
        with pytest.raises(ValueError, match="Unknown workflow"):
            get_workflow_scripts("wf_private_scripts", authenticated=False)

    def test_private_workflow_visible_with_auth(self, scripts_env):
        result = get_workflow_scripts("wf_private_scripts", authenticated=True)
        assert result.workflow_id == "wf_private_scripts"
        assert len(result.scripts) == 1


class TestPublicRedaction:
    """Public/unauthenticated calls must not expose host-specific metadata."""

    def test_public_omits_absolute_path(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        for s in result.scripts:
            assert s.absolute_path is None, (
                f"absolute_path should be None for public callers, "
                f"got {s.absolute_path!r} on role={s.role}"
            )

    def test_public_omits_repo_root(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        assert "repo_root" not in result.provenance

    def test_public_omits_git_commit(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        assert "git_commit" not in result.provenance

    def test_public_omits_git_remote(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        assert "git_remote" not in result.provenance

    def test_public_retains_structural_fields(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        assert result.workflow_id == "wf_scripts"
        assert result.title == "Script Test Workflow"
        assert result.entry_script == "scripts/entry.sh"
        assert len(result.scripts) == 3
        for s in result.scripts:
            assert s.role
            assert s.path
            assert s.language

    def test_public_retains_workflow_schema_in_provenance(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        assert result.provenance.get("workflow_schema") == "workflows/wf_scripts.yaml"

    def test_public_retains_exists_flag(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=False)
        for s in result.scripts:
            assert s.exists is True

    def test_authenticated_includes_absolute_path(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        for s in result.scripts:
            assert s.absolute_path is not None
            assert str(scripts_env) in s.absolute_path

    def test_authenticated_includes_repo_root(self, scripts_env):
        result = get_workflow_scripts("wf_scripts", authenticated=True)
        assert result.provenance["repo_root"] == str(scripts_env)

    def test_fallback_entry_public_omits_absolute_path(self, scripts_env):
        """Workflows without a scripts section: public still omits absolute_path."""
        result = get_workflow_scripts("wf_no_scripts", authenticated=False)
        assert len(result.scripts) == 1
        assert result.scripts[0].absolute_path is None

    def test_fallback_entry_authenticated_includes_absolute_path(self, scripts_env):
        result = get_workflow_scripts("wf_no_scripts", authenticated=True)
        assert result.scripts[0].absolute_path is not None


class TestGetWorkflowScriptsAppWrapper:
    """Tests for the app.py wrapper (auth enforcement)."""

    def test_public_discovery_no_auth_needed(self, scripts_env):
        result = get_workflow_scripts_app(workflow_id="wf_scripts")
        assert result.get("error") is not True
        assert result["workflow_id"] == "wf_scripts"
        assert len(result["scripts"]) == 3

    def test_public_discovery_redacts_absolute_path(self, scripts_env):
        """App-layer public call must redact absolute_path."""
        result = get_workflow_scripts_app(workflow_id="wf_scripts")
        assert result.get("error") is not True
        for s in result["scripts"]:
            assert s["absolute_path"] is None

    def test_public_discovery_redacts_provenance(self, scripts_env):
        """App-layer public call must not include repo_root/git_commit/git_remote."""
        result = get_workflow_scripts_app(workflow_id="wf_scripts")
        prov = result["provenance"]
        assert "repo_root" not in prov
        assert "git_commit" not in prov
        assert "git_remote" not in prov

    def test_authenticated_discovery_includes_absolute_path(self, scripts_env):
        result = get_workflow_scripts_app(
            workflow_id="wf_scripts", auth_token="test-token"
        )
        assert result.get("error") is not True
        for s in result["scripts"]:
            assert s["absolute_path"] is not None

    def test_authenticated_discovery_includes_provenance(self, scripts_env):
        result = get_workflow_scripts_app(
            workflow_id="wf_scripts", auth_token="test-token"
        )
        prov = result["provenance"]
        assert "repo_root" in prov

    def test_private_workflow_hidden_without_token(self, scripts_env):
        result = get_workflow_scripts_app(workflow_id="wf_private_scripts")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_private_workflow_visible_with_token(self, scripts_env):
        result = get_workflow_scripts_app(
            workflow_id="wf_private_scripts", auth_token="test-token"
        )
        assert result.get("error") is not True
        assert result["workflow_id"] == "wf_private_scripts"

    def test_unknown_workflow_returns_error(self, scripts_env):
        result = get_workflow_scripts_app(workflow_id="nonexistent")
        assert result.get("error") is True
        assert result["code"] == "WORKFLOW_ERROR"

    def test_auth_required_when_public_discovery_false(self, scripts_env):
        config_module._config.server.public_discovery = False
        result = get_workflow_scripts_app(workflow_id="wf_scripts")
        assert result.get("error") is True
        assert result["code"] == "AUTH_FAILED"

    def test_auth_passes_when_public_discovery_false(self, scripts_env):
        config_module._config.server.public_discovery = False
        result = get_workflow_scripts_app(
            workflow_id="wf_scripts", auth_token="test-token"
        )
        assert result.get("error") is not True
        assert result["workflow_id"] == "wf_scripts"
