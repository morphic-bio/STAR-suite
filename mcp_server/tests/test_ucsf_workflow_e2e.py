"""E2E dry-run contract test for the UCSF STAR-suite production workflow.

Uses the real workflow schema file (mcp_server/workflows/ucsf_star_suite_production.yaml)
with temporary directories and stub executables. Never executes the real workflow.
"""

import os
import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import get_workflow_schema, load_config
from mcp_server.schemas.run_config import RunConfig
from mcp_server.tools.preflight import run_preflight
from mcp_server.tools.workflows import (
    describe_workflow,
    get_workflow_parameter_schema,
    list_workflows,
    render_workflow_command,
    validate_workflow_parameters,
)


# Path to the real UCSF workflow schema (relative to repo root)
UCSF_SCHEMA_REL = "mcp_server/workflows/ucsf_star_suite_production.yaml"
UCSF_ENTRY_SCRIPT = "scripts/paper/run_ucsf_corrected_production_workflow.sh"
UCSF_WORKFLOW_ID = "ucsf_star_suite_production"


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def ucsf_e2e_env():
    """Set up an environment with temp paths for UCSF workflow E2E testing.

    Uses the real workflow schema file but with temporary fixture paths.
    """
    # Find the real repo root (we need the real schema file)
    this_file = Path(__file__).resolve()
    repo_root = this_file.parent.parent.parent  # mcp_server/tests/ -> mcp_server/ -> repo

    real_schema_path = repo_root / UCSF_SCHEMA_REL
    if not real_schema_path.exists():
        pytest.skip(f"UCSF schema not found at {real_schema_path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Create temp directories for required path params
        dataset_root = tmp / "ucsf-corrected"
        dataset_root.mkdir()
        genome_dir = tmp / "genome"
        genome_dir.mkdir()
        out_root = tmp / "output"
        out_root.mkdir()

        # Create temp files for required file params
        feature_ref = tmp / "feature_ref.csv"
        feature_ref.write_text("header\n")
        solo_cb_whitelist = tmp / "cb_whitelist.txt"
        solo_cb_whitelist.write_text("AAAA\n")
        cr_whitelist = tmp / "cr_whitelist.txt"
        cr_whitelist.write_text("AAAA\n")

        # Create a stub STAR binary
        star_bin = tmp / "STAR.release"
        star_bin.write_text("#!/bin/bash\necho stub\n")
        os.chmod(star_bin, 0o755)

        # Create entry script stub (for preflight-compatible path)
        scripts_dir = tmp / "scripts" / "paper"
        scripts_dir.mkdir(parents=True)
        entry_script = scripts_dir / "run_ucsf_corrected_production_workflow.sh"
        entry_script.write_text("#!/bin/bash\necho stub\n")
        os.chmod(entry_script, 0o755)

        # Copy the real schema file into our temp tree so config can find it
        wf_dir = tmp / "mcp_server" / "workflows"
        wf_dir.mkdir(parents=True)
        import shutil
        shutil.copy2(real_schema_path, wf_dir / "ucsf_star_suite_production.yaml")

        # Build config
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
                    "id": UCSF_WORKFLOW_ID,
                    "title": "UCSF STAR-suite Production Workflow",
                    "summary": "Corrected UCSF GEX+guide workflow.",
                    "entry_script": UCSF_ENTRY_SCRIPT,
                    "kind": "shell_workflow",
                    "schema_file": UCSF_SCHEMA_REL,
                }
            ],
            "scripts": [
                {
                    "name": "ucsf_corrected_production",
                    "path": UCSF_ENTRY_SCRIPT,
                    "module": "workflow",
                    "description": "UCSF corrected production workflow",
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)

        yield {
            "tmp": tmp,
            "dataset_root": str(dataset_root),
            "genome_dir": str(genome_dir),
            "out_root": str(out_root),
            "feature_ref": str(feature_ref),
            "solo_cb_whitelist": str(solo_cb_whitelist),
            "cr_whitelist": str(cr_whitelist),
            "star_bin": str(star_bin),
            "entry_script": str(entry_script),
        }


class TestUCSFWorkflowDiscovery:
    def test_workflow_is_discoverable(self, ucsf_e2e_env):
        result = list_workflows()
        ids = [w.id for w in result.workflows]
        assert UCSF_WORKFLOW_ID in ids

    def test_describe_returns_stages(self, ucsf_e2e_env):
        result = describe_workflow(UCSF_WORKFLOW_ID)
        assert result.id == UCSF_WORKFLOW_ID
        stage_names = [s.name for s in result.stages]
        assert "star_alignment" in stage_names
        assert "downstream" in stage_names

        # Verify stage metadata fields
        by_name = {s.name: s for s in result.stages}
        star_stage = by_name["star_alignment"]
        assert star_stage.script == "scripts/run_ucsf_perturb_yremove_batch.sh"
        assert star_stage.gated_by is None

        downstream = by_name["downstream"]
        assert downstream.gated_by == "star_only"
        assert downstream.script == "scripts/run_scrna_downstream_gene_full_velocyto.sh"

    def test_describe_returns_typed_parameter_groups(self, ucsf_e2e_env):
        result = describe_workflow(UCSF_WORKFLOW_ID)
        grp_names = [g.name for g in result.parameter_groups]
        assert "selection" in grp_names
        assert "inputs" in grp_names
        # Verify ParameterGroupInfo structure
        selection = next(g for g in result.parameter_groups if g.name == "selection")
        assert selection.title == "Sample Selection"
        assert "samples" in selection.parameters
        assert "all_samples" in selection.parameters

    def test_parameter_schema_has_groups(self, ucsf_e2e_env):
        result = get_workflow_parameter_schema(UCSF_WORKFLOW_ID)
        assert result.workflow_id == UCSF_WORKFLOW_ID
        group_names = [g.name for g in result.parameter_groups]
        assert "selection" in group_names
        assert "inputs" in group_names
        assert "performance" in group_names
        assert "execution_mode" in group_names
        assert "transfer" in group_names
        assert "output" in group_names

    def test_parameter_metadata_fields(self, ucsf_e2e_env):
        """Verify category, source, path_must_exist on UCSF params."""
        result = get_workflow_parameter_schema(UCSF_WORKFLOW_ID)
        by_name = {p.name: p for p in result.parameters}

        dr = by_name["dataset_root"]
        assert dr.category == "inputs"
        assert dr.source == "workflow_wrapper"
        assert dr.path_must_exist is True
        assert dr.type == "directory"

        thr = by_name["threads"]
        assert thr.category == "performance"
        assert thr.default == 24
        assert thr.type == "int"

        so = by_name["star_only"]
        assert so.category == "execution_mode"
        assert so.type == "bool"
        assert so.default is False

        samp = by_name["samples"]
        assert samp.type == "string_list"
        assert samp.category == "selection"
        assert samp.repeatable is False

    def test_declarative_rendering_hints(self, ucsf_e2e_env):
        """Verify env_var, skip_when_default, is_output_root on UCSF params."""
        result = get_workflow_parameter_schema(UCSF_WORKFLOW_ID)
        by_name = {p.name: p for p in result.parameters}

        ds = by_name["downsample_seed"]
        assert ds.env_var == "DOWNSAMPLE_SEED"

        dr = by_name["downsample_reads"]
        assert dr.skip_when_default is True

        out = by_name["out_root"]
        assert out.is_output_root is True

        # Non-annotated params should have defaults
        thr = by_name["threads"]
        assert thr.env_var is None
        assert thr.skip_when_default is False
        assert thr.is_output_root is False

    def test_parameter_schema_has_constraints(self, ucsf_e2e_env):
        result = get_workflow_parameter_schema(UCSF_WORKFLOW_ID)
        kinds = [c.kind for c in result.constraints]
        assert "mutual_exclusion" in kinds
        assert "group_required" in kinds
        assert "positive" in kinds

        # Verify constraint details
        me = next(c for c in result.constraints if c.kind == "mutual_exclusion")
        assert "samples" in me.params
        assert "all_samples" in me.params
        assert me.message != ""

        gr = next(c for c in result.constraints if c.kind == "group_required")
        assert "globus_src_endpoint" in gr.params
        assert "globus_dst_endpoint" in gr.params
        assert "globus_dst_root" in gr.params

        pos = next(c for c in result.constraints if c.kind == "positive")
        assert "threads" in pos.params


class TestUCSFValidation:
    def test_minimal_valid_request(self, ucsf_e2e_env):
        """Minimal valid params: just override paths to temp fixtures."""
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],
            "out_root": ucsf_e2e_env["out_root"],
            "all_samples": True,
            "dry_run": True,
        }
        result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params)
        assert result.valid, f"Validation errors: {result.errors}"
        assert result.errors == []

    def test_mutual_exclusion_enforced(self, ucsf_e2e_env):
        params = {
            "samples": "s1,s2",
            "all_samples": True,
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],
            "dry_run": True,
        }
        result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params, check_paths=True)
        assert not result.valid
        assert any("Mutually exclusive" in e for e in result.errors)


class TestUCSFRendering:
    def test_dry_run_renders_flag(self, ucsf_e2e_env):
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],
            "out_root": ucsf_e2e_env["out_root"],
            "all_samples": True,
            "dry_run": True,
        }
        result = render_workflow_command(UCSF_WORKFLOW_ID, params)
        assert "--dry-run" in result.argv
        assert "--all-samples" in result.argv

    def test_entry_script_is_ucsf_wrapper(self, ucsf_e2e_env):
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "all_samples": True,
            "dry_run": True,
        }
        result = render_workflow_command(UCSF_WORKFLOW_ID, params)
        assert result.entry_script.endswith(
            "run_ucsf_corrected_production_workflow.sh"
        )

    def test_star_only_rendered(self, ucsf_e2e_env):
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "star_only": True,
            "dry_run": True,
        }
        result = render_workflow_command(UCSF_WORKFLOW_ID, params)
        assert "--star-only" in result.argv

    def test_threads_rendered(self, ucsf_e2e_env):
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "threads": 16,
            "dry_run": True,
        }
        result = render_workflow_command(UCSF_WORKFLOW_ID, params)
        idx = result.argv.index("--threads")
        assert result.argv[idx + 1] == "16"

    def test_downsample_seed_env_override(self, ucsf_e2e_env):
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "downsample_seed": 42,
            "dry_run": True,
        }
        result = render_workflow_command(UCSF_WORKFLOW_ID, params)
        assert result.env_overrides.get("DOWNSAMPLE_SEED") == "42"


class TestUCSFE2EContractFlow:
    def test_schema_validate_render_pipeline(self, ucsf_e2e_env):
        """Full E2E: schema discovery -> validate -> render -> check argv."""
        # 1. Schema is discoverable
        schema = get_workflow_schema(UCSF_WORKFLOW_ID)
        assert schema is not None
        assert schema.id == UCSF_WORKFLOW_ID

        # 2. Validate minimal valid request
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],
            "out_root": ucsf_e2e_env["out_root"],
            "all_samples": True,
            "dry_run": True,
        }
        val_result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params)
        assert val_result.valid, f"Validation errors: {val_result.errors}"

        # 3. Render command
        render_result = render_workflow_command(
            UCSF_WORKFLOW_ID, val_result.normalized_params
        )

        # 4. Assert expected argv content
        assert render_result.entry_script.endswith(
            "run_ucsf_corrected_production_workflow.sh"
        )
        assert "--dry-run" in render_result.argv
        assert "--all-samples" in render_result.argv
        assert "--dataset-root" in render_result.argv
        assert ucsf_e2e_env["dataset_root"] in render_result.argv
        assert "--genome-dir" in render_result.argv
        assert "--star-bin" in render_result.argv
        assert render_result.output_root == ucsf_e2e_env["out_root"]

        # 5. Verify the rendered entry script path exists
        assert Path(render_result.entry_script).exists()

    def test_validate_render_preflight_composition(self, ucsf_e2e_env):
        """Full E2E: validate -> render -> preflight passes end-to-end."""
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],
            "out_root": ucsf_e2e_env["out_root"],
            "all_samples": True,
            "dry_run": True,
        }

        # Step 1: validate
        val_result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params)
        assert val_result.valid, f"Validation errors: {val_result.errors}"

        # Step 2: render
        render_result = render_workflow_command(
            UCSF_WORKFLOW_ID, val_result.normalized_params
        )

        # Step 3: map into RunConfig and run preflight
        run_config = RunConfig(
            script="ucsf_corrected_production",
            args=render_result.argv[1:],
            out_dir=render_result.output_root,
            env_overrides=render_result.env_overrides or None,
        )

        preflight_result = run_preflight(run_config)
        assert preflight_result.valid, (
            f"Preflight errors: {preflight_result.errors}\n"
            f"Checks: {[(c.name, c.passed, c.message) for c in preflight_result.checks]}"
        )

        # Verify key preflight checks passed
        check_names = {c.name: c for c in preflight_result.checks}
        assert check_names["script_allowed"].passed
        assert check_names["script_path"].passed
        assert check_names["output_directory"].passed


class TestUCSFExecutableValidation:
    def test_non_executable_star_bin_fails_validation(self, ucsf_e2e_env):
        """A non-executable star_bin should produce a validation error."""
        # Create a non-executable file
        tmp = ucsf_e2e_env["tmp"]
        bad_star = Path(tmp) / "STAR.noexec"
        bad_star.write_text("#!/bin/bash\necho stub\n")
        os.chmod(bad_star, 0o644)  # readable but NOT executable

        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": str(bad_star),
            "out_root": ucsf_e2e_env["out_root"],
            "dry_run": True,
        }
        result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params)
        assert not result.valid
        assert any("not executable" in e for e in result.errors)

    def test_executable_star_bin_passes_validation(self, ucsf_e2e_env):
        """A properly executable star_bin should pass validation."""
        params = {
            "dataset_root": ucsf_e2e_env["dataset_root"],
            "genome_dir": ucsf_e2e_env["genome_dir"],
            "feature_ref": ucsf_e2e_env["feature_ref"],
            "solo_cb_whitelist": ucsf_e2e_env["solo_cb_whitelist"],
            "cr_whitelist": ucsf_e2e_env["cr_whitelist"],
            "star_bin": ucsf_e2e_env["star_bin"],  # fixture creates it executable
            "out_root": ucsf_e2e_env["out_root"],
            "dry_run": True,
        }
        result = validate_workflow_parameters(UCSF_WORKFLOW_ID, params)
        assert result.valid, f"Validation errors: {result.errors}"
        assert not any("executable" in e for e in result.errors)
