"""Stress tests for scaffold_workflow_schema against real scripts.

Exercises the scaffold tool against every while/case-based script in the
repository — smoke tests, production benchmarks, and E2E pipelines.  Then
validates each scaffolded schema through validate_draft_workflow_schema and
(for production scripts) through the full validate->render pipeline.
"""

import os
import shutil
import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.scaffold import (
    scaffold_workflow_schema,
    validate_draft_workflow_schema,
)

REPO_ROOT = Path("/mnt/pikachu/STAR-suite")


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def real_repo_env():
    """Config fixture pointing at the real repo root so scaffold can resolve
    relative paths and trusted-root checks pass.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        cfg = {
            "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
            "paths": {
                "repo_root": str(REPO_ROOT),
                "artifact_log_root": str(tmp / "artifacts"),
                "temp_root": str(tmp / "tmp"),
            },
            "trusted_roots": [str(REPO_ROOT), "/tmp"],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)
        load_config(cfg_path)
        yield


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _skip_if_missing(path: Path):
    if not path.exists():
        pytest.skip(f"Script not available: {path}")


def _scaffold_and_validate(script_rel: str, min_params: int):
    """Scaffold a script and run it through validate_draft.

    Returns (scaffold_result, validate_result, parsed_yaml_dict).
    """
    script = REPO_ROOT / script_rel
    _skip_if_missing(script)

    result = scaffold_workflow_schema(script_rel)

    # Basic scaffold sanity
    assert result.parameter_count >= min_params, (
        f"{script_rel}: expected >= {min_params} params, got {result.parameter_count}"
    )
    assert result.draft_yaml
    assert result.source_script

    # Parse the YAML
    parsed = yaml.safe_load(result.draft_yaml)
    assert isinstance(parsed, dict)
    assert "parameters" in parsed
    assert "rendering" in parsed

    # Validate through the draft validator
    val = validate_draft_workflow_schema(result.draft_yaml)
    assert val.valid, f"{script_rel} validation errors: {val.errors}"
    assert val.parameter_count == result.parameter_count

    return result, val, parsed


# ---------------------------------------------------------------------------
# Smoke test scripts (tests/run_*.sh with while/case)
# ---------------------------------------------------------------------------


class TestScaffoldSmokeScripts:
    """Scaffold each smoke test that has a while/case flag block."""

    # (script_relative_path, minimum_expected_params)
    SMOKE_SCRIPTS = [
        ("tests/run_a375_public_smoke.sh", 5),
        ("tests/run_gex_binary_smoke.sh", 3),
        ("tests/run_flex_tiny_public_binary_smoke.sh", 3),
        ("tests/run_adapter_clip_synthetic_test.sh", 2),
        ("tests/run_binary_test_matrix.sh", 3),
        ("tests/run_genome_generate_validation.sh", 2),
        ("tests/run_slam_end_to_end.sh", 2),
        ("tests/run_slam_end_to_end_fixture.sh", 3),
        ("tests/run_slam_end_to_end_vb_commonmask.sh", 4),
        ("tests/run_a375_mex_comparison.sh", 2),
        ("tests/run_default_groups_smoke.sh", 2),
        ("tests/run_downsample10M_alignment_and_gedi_goalset.sh", 3),
        ("tests/run_parms_tests.sh", 2),
    ]

    @pytest.mark.parametrize("script_rel,min_params", SMOKE_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in SMOKE_SCRIPTS])
    def test_scaffold_and_validate(self, real_repo_env, script_rel, min_params):
        """Each script scaffolds without error and produces a valid schema."""
        _scaffold_and_validate(script_rel, min_params)

    @pytest.mark.parametrize("script_rel,min_params", SMOKE_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in SMOKE_SCRIPTS])
    def test_no_duplicate_params(self, real_repo_env, script_rel, min_params):
        """Scaffolded schemas should not have duplicate parameter names."""
        _, _, parsed = _scaffold_and_validate(script_rel, min_params)
        names = [p["name"] for p in parsed["parameters"]]
        assert len(names) == len(set(names)), (
            f"{script_rel}: duplicate params: "
            f"{[n for n in names if names.count(n) > 1]}"
        )

    @pytest.mark.parametrize("script_rel,min_params", SMOKE_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in SMOKE_SCRIPTS])
    def test_flag_order_matches_params(self, real_repo_env, script_rel, min_params):
        """flag_order should list exactly the scaffolded parameter names."""
        _, _, parsed = _scaffold_and_validate(script_rel, min_params)
        param_names = [p["name"] for p in parsed["parameters"]]
        flag_order = parsed["rendering"]["flag_order"]
        assert flag_order == param_names, (
            f"{script_rel}: flag_order mismatch"
        )

    @pytest.mark.parametrize("script_rel,min_params", SMOKE_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in SMOKE_SCRIPTS])
    def test_no_shell_variable_defaults(self, real_repo_env, script_rel, min_params):
        """No parameter default should contain raw ${VAR} or $(...) references."""
        _, _, parsed = _scaffold_and_validate(script_rel, min_params)
        for p in parsed["parameters"]:
            default = p.get("default")
            if default is not None:
                default_str = str(default)
                assert "${" not in default_str, (
                    f"{script_rel}: param '{p['name']}' has raw shell var "
                    f"in default: {default_str}"
                )
                assert "$(" not in default_str, (
                    f"{script_rel}: param '{p['name']}' has command sub "
                    f"in default: {default_str}"
                )


# ---------------------------------------------------------------------------
# E2E production scripts (scripts/paper/run_*.sh + flex)
# ---------------------------------------------------------------------------


class TestScaffoldProductionScripts:
    """Scaffold production/benchmark scripts and verify quality."""

    PROD_SCRIPTS = [
        ("scripts/paper/run_msk_30polyko_benchmark.sh", 4),
        ("scripts/paper/run_pe_bulk_feature_benchmark.sh", 15),
        ("scripts/paper/run_a375_gexonly_legacy_benchmark.sh", 3),  # mode flags are unrecognized pattern
        ("scripts/paper/run_msk_gexonly_legacy_benchmark.sh", 3),  # mode flags are unrecognized pattern
        ("scripts/paper/run_ucsf_corrected_production_workflow.sh", 14),
        ("scripts/paper/run_ucsf_ebs2_2_benchmark.sh", 4),
        ("scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh", 3),  # mode flags are unrecognized pattern
        ("scripts/run_flex_cr_config.sh", 8),
    ]

    @pytest.mark.parametrize("script_rel,min_params", PROD_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in PROD_SCRIPTS])
    def test_scaffold_and_validate(self, real_repo_env, script_rel, min_params):
        """Production scripts scaffold without error and produce valid schemas."""
        _scaffold_and_validate(script_rel, min_params)

    @pytest.mark.parametrize("script_rel,min_params", PROD_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in PROD_SCRIPTS])
    def test_no_duplicate_params(self, real_repo_env, script_rel, min_params):
        _, _, parsed = _scaffold_and_validate(script_rel, min_params)
        names = [p["name"] for p in parsed["parameters"]]
        assert len(names) == len(set(names))

    @pytest.mark.parametrize("script_rel,min_params", PROD_SCRIPTS,
                             ids=[s[0].split("/")[-1] for s in PROD_SCRIPTS])
    def test_no_shell_variable_defaults(self, real_repo_env, script_rel, min_params):
        _, _, parsed = _scaffold_and_validate(script_rel, min_params)
        for p in parsed["parameters"]:
            default = p.get("default")
            if default is not None:
                default_str = str(default)
                assert "${" not in default_str, (
                    f"{script_rel}: '{p['name']}' raw shell var: {default_str}"
                )
                assert "$(" not in default_str, (
                    f"{script_rel}: '{p['name']}' command sub: {default_str}"
                )


class TestProductionDescriptionQuality:
    """Production scripts with usage() heredocs should have real descriptions."""

    # Scripts with verified usage() heredoc
    HEREDOC_SCRIPTS = [
        "scripts/paper/run_msk_30polyko_benchmark.sh",
        "scripts/paper/run_pe_bulk_feature_benchmark.sh",
        "scripts/paper/run_a375_gexonly_legacy_benchmark.sh",
        "scripts/paper/run_msk_gexonly_legacy_benchmark.sh",
        "scripts/paper/run_ucsf_corrected_production_workflow.sh",
        "scripts/paper/run_ucsf_ebs2_2_benchmark.sh",
        "scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh",
        "scripts/run_flex_cr_config.sh",
    ]

    @pytest.mark.parametrize("script_rel", HEREDOC_SCRIPTS,
                             ids=[s.split("/")[-1] for s in HEREDOC_SCRIPTS])
    def test_most_descriptions_populated(self, real_repo_env, script_rel):
        """At least 50% of params should have real (non-TODO) descriptions."""
        result = scaffold_workflow_schema(script_rel)
        parsed = yaml.safe_load(result.draft_yaml)
        params = parsed["parameters"]
        todo_count = sum(1 for p in params if "TODO" in p.get("description", ""))
        real_count = len(params) - todo_count
        assert real_count >= len(params) / 2, (
            f"{script_rel}: only {real_count}/{len(params)} params have real "
            f"descriptions (TODO count: {todo_count})"
        )

    @pytest.mark.parametrize("script_rel", HEREDOC_SCRIPTS,
                             ids=[s.split("/")[-1] for s in HEREDOC_SCRIPTS])
    def test_descriptions_not_garbled(self, real_repo_env, script_rel):
        """No description should start with '--' (cross-contamination sign)."""
        result = scaffold_workflow_schema(script_rel)
        parsed = yaml.safe_load(result.draft_yaml)
        for p in parsed["parameters"]:
            desc = p.get("description", "")
            assert not desc.startswith("--"), (
                f"{script_rel}: param '{p['name']}' description garbled: {desc!r}"
            )


# ---------------------------------------------------------------------------
# Cross-check: Flex, MSK, PE specific parameters
# ---------------------------------------------------------------------------


class TestFlexSchemaAccuracy:
    """Verify key Flex parameters are captured correctly."""

    def test_flex_cr_config_key_params(self, real_repo_env):
        script = "scripts/run_flex_cr_config.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 8)
        names = {p["name"] for p in parsed["parameters"]}

        # Core Flex config params
        expected_core = {"cr_config", "threads", "genome_dir", "cb_whitelist"}
        missing = expected_core - names
        assert not missing, f"Missing expected Flex params: {missing}"

    def test_flex_smoke_key_params(self, real_repo_env):
        script = "tests/run_flex_tiny_public_binary_smoke.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 3)
        names = {p["name"] for p in parsed["parameters"]}

        # Core params: outdir, threads, read_limit
        expected_core = {"outdir", "threads"}
        missing = expected_core - names
        assert not missing, f"Missing expected Flex smoke params: {missing}"


class TestMSKSchemaAccuracy:
    """Verify key MSK benchmark parameters."""

    def test_msk_30polyko_key_params(self, real_repo_env):
        script = "scripts/paper/run_msk_30polyko_benchmark.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 4)
        names = {p["name"] for p in parsed["parameters"]}

        expected = {"outdir", "threads", "star_bin"}
        missing = expected - names
        assert not missing, f"Missing MSK params: {missing}"

    def test_msk_gexonly_legacy_key_params(self, real_repo_env):
        script = "scripts/paper/run_msk_gexonly_legacy_benchmark.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 3)
        names = {p["name"] for p in parsed["parameters"]}

        expected = {"outdir", "threads", "star_bin"}
        missing = expected - names
        assert not missing, f"Missing MSK legacy params: {missing}"


class TestPESchemaAccuracy:
    """Verify key PE (paired-end) benchmark parameters."""

    def test_pe_bulk_feature_key_params(self, real_repo_env):
        script = "scripts/paper/run_pe_bulk_feature_benchmark.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 15)
        names = {p["name"] for p in parsed["parameters"]}

        expected_core = {
            "out_root", "threads", "star_index", "gtf", "dry_run",
        }
        missing = expected_core - names
        assert not missing, f"Missing PE params: {missing}"

    def test_pe_bool_flags_typed_correctly(self, real_repo_env):
        """Boolean flags like --dry-run, --yremove should be type=bool."""
        script = "scripts/paper/run_pe_bulk_feature_benchmark.sh"
        _skip_if_missing(REPO_ROOT / script)
        _, _, parsed = _scaffold_and_validate(script, 15)
        by_name = {p["name"]: p for p in parsed["parameters"]}

        if "dry_run" in by_name:
            assert by_name["dry_run"]["type"] == "bool", (
                f"dry_run should be bool, got {by_name['dry_run']['type']}"
            )


# ---------------------------------------------------------------------------
# UCSF production schema: compare scaffold output to hand-curated schema
# ---------------------------------------------------------------------------


class TestUCSFSchemaConsistency:
    """Compare scaffolded UCSF schema against the hand-curated production YAML."""

    @pytest.fixture
    def curated_schema(self):
        schema_path = REPO_ROOT / "mcp_server" / "workflows" / "ucsf_star_suite_production.yaml"
        if not schema_path.exists():
            pytest.skip("Curated UCSF schema not available")
        with open(schema_path) as f:
            return yaml.safe_load(f)

    def test_scaffolded_covers_curated_params(self, real_repo_env, curated_schema):
        """Every param in the hand-curated schema should appear in the scaffold."""
        script = "scripts/paper/run_ucsf_corrected_production_workflow.sh"
        _skip_if_missing(REPO_ROOT / script)
        result = scaffold_workflow_schema(script)
        parsed = yaml.safe_load(result.draft_yaml)
        scaffolded_names = {p["name"] for p in parsed["parameters"]}
        curated_names = {p["name"] for p in curated_schema["parameters"]}

        missing = curated_names - scaffolded_names
        assert not missing, (
            f"Scaffold missed params from curated schema: {missing}"
        )

    def test_scaffolded_flags_match_curated(self, real_repo_env, curated_schema):
        """CLI flags should match between scaffold and curated schema."""
        script = "scripts/paper/run_ucsf_corrected_production_workflow.sh"
        _skip_if_missing(REPO_ROOT / script)
        result = scaffold_workflow_schema(script)
        parsed = yaml.safe_load(result.draft_yaml)

        scaffolded_flags = {p["name"]: p["cli_flag"] for p in parsed["parameters"]}
        curated_flags = {p["name"]: p["cli_flag"] for p in curated_schema["parameters"]}

        # Only check params that exist in both
        common = set(scaffolded_flags) & set(curated_flags)
        mismatches = {}
        for name in common:
            if scaffolded_flags[name] != curated_flags[name]:
                mismatches[name] = (scaffolded_flags[name], curated_flags[name])

        assert not mismatches, (
            f"CLI flag mismatches between scaffold and curated: {mismatches}"
        )

    def test_curated_schema_validates(self, real_repo_env, curated_schema):
        """The hand-curated UCSF schema should also pass draft validation."""
        draft_yaml = yaml.dump(curated_schema, default_flow_style=False, sort_keys=False)
        val = validate_draft_workflow_schema(draft_yaml)
        assert val.valid, f"Curated schema validation errors: {val.errors}"


# ---------------------------------------------------------------------------
# Full pipeline: scaffold -> validate -> render (production scripts)
# ---------------------------------------------------------------------------


class TestScaffoldToRenderPipeline:
    """For production scripts, test the full scaffold->validate->render path.

    This registers the scaffolded schema as a workflow and exercises the
    validate_workflow_parameters and render_workflow_command tools.
    """

    @pytest.fixture
    def pipeline_env(self):
        """Config with real repo root and workflow registration capability."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cfg = {
                "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
                "paths": {
                    "repo_root": str(REPO_ROOT),
                    "artifact_log_root": str(tmp / "artifacts"),
                    "temp_root": str(tmp / "tmp"),
                },
                "trusted_roots": [
                    str(REPO_ROOT), "/tmp", str(tmp),
                    "/storage", "/mnt/pikachu", "/home",
                ],
                "scripts": [
                    {
                        "name": "ucsf_corrected_production",
                        "path": "scripts/paper/run_ucsf_corrected_production_workflow.sh",
                        "module": "workflow",
                        "description": "UCSF corrected production workflow",
                    },
                    {
                        "name": "msk_30polyko",
                        "path": "scripts/paper/run_msk_30polyko_benchmark.sh",
                        "module": "workflow",
                        "description": "MSK 30-polyko benchmark",
                    },
                    {
                        "name": "pe_bulk_feature",
                        "path": "scripts/paper/run_pe_bulk_feature_benchmark.sh",
                        "module": "workflow",
                        "description": "PE bulk feature benchmark",
                    },
                    {
                        "name": "flex_cr_config",
                        "path": "scripts/run_flex_cr_config.sh",
                        "module": "workflow",
                        "description": "Flex CR config workflow",
                    },
                ],
            }
            cfg_path = tmp / "config.yaml"
            with open(cfg_path, "w") as f:
                yaml.dump(cfg, f)
            load_config(cfg_path)
            yield tmp

    def _scaffold_register_and_test(self, pipeline_env, script_rel, workflow_id,
                                     test_params):
        """Scaffold a script, register as workflow, validate params, render."""
        from mcp_server.tools.workflows import (
            validate_workflow_parameters,
            render_workflow_command,
        )

        _skip_if_missing(REPO_ROOT / script_rel)

        # Step 1: Scaffold
        result = scaffold_workflow_schema(script_rel)
        assert result.parameter_count > 0

        # Step 2: Validate draft
        val = validate_draft_workflow_schema(result.draft_yaml)
        assert val.valid, f"Draft validation failed: {val.errors}"

        # Step 3: Parse and register as a workflow schema
        parsed = yaml.safe_load(result.draft_yaml)
        parsed["id"] = workflow_id

        from mcp_server.schemas.workflow import WorkflowSchema
        schema = WorkflowSchema(**parsed)

        # Register in the config module's schema cache
        config_module._workflow_schemas[workflow_id] = schema

        # Also register in config's workflows list if not there
        config = config_module._config
        workflow_entry_exists = any(
            w.id == workflow_id for w in (config.workflows or [])
        )
        if not workflow_entry_exists:
            from mcp_server.schemas.config import WorkflowConfig
            config.workflows = list(config.workflows or [])
            config.workflows.append(WorkflowConfig(
                id=workflow_id,
                title=schema.title,
                summary=schema.summary,
                entry_script=schema.entry_script,
                kind=schema.kind,
                schema_file="<in-memory>",  # not loaded from disk
            ))

        # Step 4: Validate workflow parameters
        val_result = validate_workflow_parameters(workflow_id, test_params)
        # Result may be a Pydantic model or a dict depending on the code path
        if hasattr(val_result, "valid"):
            assert val_result.valid is True, (
                f"Parameter validation failed: {val_result.errors}"
            )
        else:
            assert val_result.get("valid") is True, (
                f"Parameter validation failed: {val_result.get('errors', [])}"
            )

        # Step 5: Render command
        render_result = render_workflow_command(workflow_id, test_params)
        if hasattr(render_result, "argv"):
            assert len(render_result.argv) > 0
            # Convert to dict for consistent downstream checks
            render_result = render_result.model_dump()
        else:
            assert "argv" in render_result
            assert len(render_result["argv"]) > 0

        return render_result

    def test_msk_30polyko_pipeline(self, pipeline_env):
        """MSK benchmark: scaffold -> validate -> render."""
        result = self._scaffold_register_and_test(
            pipeline_env,
            "scripts/paper/run_msk_30polyko_benchmark.sh",
            "msk_30polyko_scaffold",
            {"outdir": "/tmp/msk_out", "threads": 8},
        )
        assert any("/tmp/msk_out" in a for a in result["argv"])

    def test_pe_bulk_feature_pipeline(self, pipeline_env):
        """PE benchmark: scaffold -> validate -> render."""
        result = self._scaffold_register_and_test(
            pipeline_env,
            "scripts/paper/run_pe_bulk_feature_benchmark.sh",
            "pe_bulk_scaffold",
            {"out_root": "/tmp/pe_out", "threads": 4},
        )
        assert any("/tmp/pe_out" in a for a in result["argv"])

    def test_flex_cr_config_pipeline(self, pipeline_env):
        """Flex config: scaffold -> validate -> render."""
        result = self._scaffold_register_and_test(
            pipeline_env,
            "scripts/run_flex_cr_config.sh",
            "flex_cr_scaffold",
            {"threads": 4},
        )
        argv_str = " ".join(result["argv"])
        assert "--threads" in argv_str

    def test_ucsf_pipeline(self, pipeline_env):
        """UCSF production: scaffold -> validate -> render."""
        result = self._scaffold_register_and_test(
            pipeline_env,
            "scripts/paper/run_ucsf_corrected_production_workflow.sh",
            "ucsf_scaffold",
            {"out_root": "/tmp/ucsf_out", "threads": 16, "all_samples": True},
        )
        argv_str = " ".join(result["argv"])
        assert "--out-root" in argv_str


# ---------------------------------------------------------------------------
# Aggregate statistics: report for manual review
# ---------------------------------------------------------------------------


class TestScaffoldInventoryReport:
    """Generate a summary report of all scaffolded scripts.

    This test always passes but prints a diagnostic table useful for
    manual review of the scaffold tool's coverage.
    """

    ALL_SCRIPTS = (
        TestScaffoldSmokeScripts.SMOKE_SCRIPTS
        + TestScaffoldProductionScripts.PROD_SCRIPTS
    )

    def test_print_scaffold_summary(self, real_repo_env, capsys):
        """Print a summary table of all scaffolded scripts."""
        rows = []
        for script_rel, _ in self.ALL_SCRIPTS:
            script = REPO_ROOT / script_rel
            if not script.exists():
                rows.append((script_rel, "MISSING", 0, 0, ""))
                continue
            try:
                result = scaffold_workflow_schema(script_rel)
                val = validate_draft_workflow_schema(result.draft_yaml)
                parsed = yaml.safe_load(result.draft_yaml)
                todo_count = sum(
                    1 for p in parsed["parameters"]
                    if "TODO" in p.get("description", "")
                )
                rows.append((
                    script_rel,
                    "VALID" if val.valid else "INVALID",
                    result.parameter_count,
                    todo_count,
                    "; ".join(val.errors) if val.errors else "",
                ))
            except Exception as e:
                rows.append((script_rel, f"ERROR: {e}", 0, 0, ""))

        # Print report
        print("\n" + "=" * 100)
        print("SCAFFOLD STRESS TEST REPORT")
        print("=" * 100)
        print(f"{'Script':<60} {'Status':<10} {'Params':>6} {'TODOs':>5} Errors")
        print("-" * 100)
        for script, status, params, todos, errors in rows:
            name = script.split("/")[-1]
            print(f"{name:<60} {status:<10} {params:>6} {todos:>5} {errors}")
        print("=" * 100)

        # All should at least scaffold without error
        error_scripts = [r[0] for r in rows if r[1].startswith("ERROR")]
        assert not error_scripts, f"Scripts that failed to scaffold: {error_scripts}"
