"""Tests for workflow schema scaffolding and draft validation tools."""

import os
import tempfile
from pathlib import Path
from textwrap import dedent

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.scaffold import (
    scaffold_workflow_schema,
    validate_draft_workflow_schema,
    _clean_default,
    _extract_defaults,
    _extract_parameters,
    _extract_usage_block,
    _extract_usage_descriptions,
    _extract_exports,
)


# ---------------------------------------------------------------------------
# Fixture: minimal config with a temp repo root
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def scaffold_env():
    """Set up a temp repo root with config loaded."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        cfg = {
            "server": {"host": "127.0.0.1", "port": 9999, "transport": "http"},
            "paths": {
                "repo_root": str(tmp),
                "artifact_log_root": str(tmp / "artifacts"),
                "temp_root": str(tmp / "tmp"),
            },
            "trusted_roots": [str(tmp), "/tmp"],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)
        yield tmp


# ---------------------------------------------------------------------------
# Sample script content for testing
# ---------------------------------------------------------------------------

SAMPLE_SCRIPT = dedent("""\
    #!/bin/bash
    set -euo pipefail

    REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
    STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
    DATASET_ROOT="${UCSF_DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected}"
    THREADS="${UCSF_THREADS:-24}"
    CELLBENDER_GPU=0
    TRIM_QC=0
    SAMPLES_CSV=""
    ALL_SAMPLES=0
    DRY_RUN=0
    OUT_ROOT=""
    DOWNSAMPLE_READS=0
    DOWNSAMPLE_SEED=1
    GENOME_DIR="/storage/autoindex_110_44/bulk_index"
    FEATURE_REF="/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref.csv"

    usage() {
      cat <<EOF
    Usage: $(basename "$0") [options]

    Options:
      --samples CSV            Comma-separated sample IDs
      --all-samples            Run every corrected UCSF sample
      --dataset-root DIR       Root of corrected UCSF dataset
      --genome-dir DIR         STAR genomeDir
      --feature-ref FILE       Feature reference CSV
      --star-bin PATH          STAR binary path
      --threads N              STAR threads
      --cellbender-gpu         Enable GPU CellBender
      --trim-qc                Emit read-level trim-QC reports
      --star-only              Run STAR alignment only
      --dry-run                Prepare manifests only, do not execute
      --out-root DIR           Output root directory
      --downsample-reads N     Downsample each library to N reads
      --downsample-seed N      Downsampling seed
      -h, --help               Show help
    EOF
    }

    export DOWNSAMPLE_SEED

    while [[ $# -gt 0 ]]; do
      case "$1" in
        --samples) SAMPLES_CSV="$2"; shift 2 ;;
        --all-samples) ALL_SAMPLES=1; shift ;;
        --dataset-root) DATASET_ROOT="$2"; shift 2 ;;
        --genome-dir) GENOME_DIR="$2"; shift 2 ;;
        --feature-ref) FEATURE_REF="$2"; shift 2 ;;
        --star-bin) STAR_BIN="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --cellbender-gpu) CELLBENDER_GPU=1; shift ;;
        --trim-qc) TRIM_QC=1; shift ;;
        --star-only) STAR_ONLY=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --downsample-reads) DOWNSAMPLE_READS="$2"; shift 2 ;;
        --downsample-seed) DOWNSAMPLE_SEED="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
      esac
    done

    echo "Running with THREADS=$THREADS"
""")


MINIMAL_SCRIPT = dedent("""\
    #!/bin/bash
    INPUT=""
    VERBOSE=0

    while [[ $# -gt 0 ]]; do
      case "$1" in
        --input) INPUT="$2"; shift 2 ;;
        --verbose) VERBOSE=1; shift ;;
        -h|--help) echo "help"; exit 0 ;;
        *) exit 1 ;;
      esac
    done
""")


# ---------------------------------------------------------------------------
# Unit tests: internal parsers
# ---------------------------------------------------------------------------


class TestExtractDefaults:
    def test_env_fallback(self):
        text = 'THREADS="${UCSF_THREADS:-24}"'
        defaults = _extract_defaults(text)
        assert defaults["THREADS"] == "24"

    def test_env_fallback_path(self):
        text = 'DATASET_ROOT="${UCSF_DATASET_ROOT:-/mnt/pikachu/data}"'
        defaults = _extract_defaults(text)
        assert defaults["DATASET_ROOT"] == "/mnt/pikachu/data"

    def test_numeric_default(self):
        text = "DOWNSAMPLE_READS=0"
        defaults = _extract_defaults(text)
        assert defaults["DOWNSAMPLE_READS"] == "0"

    def test_empty_string_default(self):
        text = 'SAMPLES_CSV=""'
        defaults = _extract_defaults(text)
        assert defaults["SAMPLES_CSV"] == ""

    def test_quoted_literal(self):
        text = 'GENOME_DIR="/storage/autoindex_110_44/bulk_index"'
        defaults = _extract_defaults(text)
        assert defaults["GENOME_DIR"] == "/storage/autoindex_110_44/bulk_index"

    def test_multiple_defaults(self):
        defaults = _extract_defaults(SAMPLE_SCRIPT)
        assert "THREADS" in defaults
        assert "CELLBENDER_GPU" in defaults
        assert "SAMPLES_CSV" in defaults


class TestExtractUsageBlock:
    def test_extracts_heredoc(self):
        block = _extract_usage_block(SAMPLE_SCRIPT)
        assert "--samples" in block
        assert "Comma-separated sample IDs" in block

    def test_no_usage_function(self):
        block = _extract_usage_block("#!/bin/bash\necho hello\n")
        assert block == ""


class TestExtractUsageDescriptions:
    def test_value_param_description(self):
        descs = _extract_usage_descriptions(SAMPLE_SCRIPT)
        assert descs["--samples"] == "Comma-separated sample IDs"

    def test_bool_param_description(self):
        descs = _extract_usage_descriptions(SAMPLE_SCRIPT)
        assert descs["--cellbender-gpu"] == "Enable GPU CellBender"

    def test_dir_param_description(self):
        descs = _extract_usage_descriptions(SAMPLE_SCRIPT)
        assert descs["--dataset-root"] == "Root of corrected UCSF dataset"

    def test_strips_arg_hint(self):
        """Arg hint like N, CSV, DIR should be stripped from description."""
        descs = _extract_usage_descriptions(SAMPLE_SCRIPT)
        # --threads N  STAR threads → should be "STAR threads", not "N STAR threads"
        assert descs["--threads"].startswith("STAR threads")
        # --samples CSV → "Comma-separated..."
        assert descs["--samples"].startswith("Comma")

    def test_no_cross_contamination(self):
        """Descriptions should only come from usage() block, not ARGS arrays."""
        # Build a script with an ARGS block that looks like usage lines
        script = dedent("""\
            #!/bin/bash
            VAR=""
            usage() { cat <<EOF
            Options:
              --input DIR    Input directory
            EOF
            }
            ARGS=(
              --input "${INPUT}"
              --output "${OUTPUT}"
            )
            while [[ $# -gt 0 ]]; do
              case "$1" in
                --input) VAR="$2"; shift 2 ;;
                *) exit 1 ;;
              esac
            done
        """)
        descs = _extract_usage_descriptions(script)
        assert descs.get("--input") == "Input directory"
        # --output from ARGS block should NOT appear
        assert "--output" not in descs

    def test_all_params_extracted(self):
        descs = _extract_usage_descriptions(SAMPLE_SCRIPT)
        expected_flags = [
            "--samples", "--all-samples", "--dataset-root", "--genome-dir",
            "--feature-ref", "--star-bin", "--threads", "--cellbender-gpu",
            "--trim-qc", "--star-only", "--dry-run", "--out-root",
            "--downsample-reads", "--downsample-seed",
        ]
        for flag in expected_flags:
            assert flag in descs, f"Missing usage description for {flag}"


class TestCleanDefault:
    def test_literal_passthrough(self):
        assert _clean_default("/mnt/pikachu/data") == "/mnt/pikachu/data"

    def test_empty_string(self):
        assert _clean_default("") == ""

    def test_none(self):
        assert _clean_default(None) is None

    def test_pure_variable_ref(self):
        assert _clean_default("${THREADS}") is None

    def test_command_substitution(self):
        assert _clean_default("/mnt/out_$(date +%Y%m%d)") is None

    def test_repo_root_prefix_stripped(self):
        assert _clean_default("${REPO_ROOT}/core/legacy/source/STAR.release") == "core/legacy/source/STAR.release"

    def test_embedded_variable(self):
        assert _clean_default("/path/${USER}/data") is None

    def test_numeric_string(self):
        assert _clean_default("24") == "24"


class TestExtractExports:
    def test_finds_export(self):
        exports = _extract_exports(SAMPLE_SCRIPT)
        assert "DOWNSAMPLE_SEED" in exports

    def test_no_exports(self):
        exports = _extract_exports(MINIMAL_SCRIPT)
        assert len(exports) == 0


class TestExtractParameters:
    def test_value_params_detected(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert "threads" in by_name
        assert by_name["threads"]["type"] == "int"

    def test_bool_params_detected(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert "cellbender_gpu" in by_name
        assert by_name["cellbender_gpu"]["type"] == "bool"

    def test_help_excluded(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        names = [p["name"] for p in params]
        assert "help" not in names

    def test_descriptions_populated(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["threads"]["description"].startswith("STAR threads")
        assert by_name["dry_run"]["description"] == "Prepare manifests only, do not execute"

    def test_defaults_populated(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["threads"]["default"] == 24  # coerced to int

    def test_variable_defaults_cleaned(self):
        """Shell variable references in defaults should be cleaned."""
        script = dedent("""\
            #!/bin/bash
            REPO_ROOT="/repo"
            STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/STAR}"
            CPU_CORES="${CPU_CORES:-${THREADS}}"
            OUT="${OUT:-/tmp/out_$(date +%s)}"
            usage() { cat <<EOF
            Options:
              --star-bin PATH     STAR binary
              --cpu-cores N       CPU cores
              --out DIR           Output
            EOF
            }
            while [[ $# -gt 0 ]]; do
              case "$1" in
                --star-bin) STAR_BIN="$2"; shift 2 ;;
                --cpu-cores) CPU_CORES="$2"; shift 2 ;;
                --out) OUT="$2"; shift 2 ;;
                *) exit 1 ;;
              esac
            done
        """)
        params = _extract_parameters(script)
        by_name = {p["name"]: p for p in params}
        assert by_name["star_bin"]["default"] == "core/STAR"  # stripped ${REPO_ROOT}/
        assert by_name["cpu_cores"]["default"] is None  # pure ${THREADS} ref
        assert by_name["out"]["default"] is None  # command substitution

    def test_env_var_detected(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["downsample_seed"].get("env_var") == "DOWNSAMPLE_SEED"

    def test_dir_type_inferred(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["dataset_root"]["type"] == "directory"
        assert by_name["genome_dir"]["type"] == "directory"
        assert by_name["out_root"]["type"] == "directory"

    def test_file_type_inferred(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["feature_ref"]["type"] == "file"
        assert by_name["star_bin"]["type"] == "file"

    def test_executable_hint_on_bin(self):
        params = _extract_parameters(SAMPLE_SCRIPT)
        by_name = {p["name"]: p for p in params}
        assert by_name["star_bin"].get("must_be_executable") is True

    def test_minimal_script(self):
        params = _extract_parameters(MINIMAL_SCRIPT)
        assert len(params) == 2
        by_name = {p["name"]: p for p in params}
        assert "input" in by_name
        assert "verbose" in by_name
        assert by_name["verbose"]["type"] == "bool"


# ---------------------------------------------------------------------------
# Integration tests: scaffold_workflow_schema
# ---------------------------------------------------------------------------


class TestScaffoldWorkflowSchema:
    def test_scaffold_produces_valid_yaml(self, scaffold_env):
        script = scaffold_env / "scripts" / "run_test_workflow.sh"
        script.parent.mkdir(parents=True)
        script.write_text(SAMPLE_SCRIPT)
        os.chmod(script, 0o755)

        result = scaffold_workflow_schema(str(script))
        assert result.parameter_count > 0
        assert result.source_script == str(script)

        # The YAML should be parseable
        parsed = yaml.safe_load(result.draft_yaml)
        assert parsed["id"] == "test_workflow"
        assert parsed["entry_script"] == "scripts/run_test_workflow.sh"
        assert len(parsed["parameters"]) == result.parameter_count

    def test_scaffold_relative_path(self, scaffold_env):
        script = scaffold_env / "scripts" / "run_my_script.sh"
        script.parent.mkdir(parents=True)
        script.write_text(MINIMAL_SCRIPT)
        os.chmod(script, 0o755)

        result = scaffold_workflow_schema("scripts/run_my_script.sh")
        assert result.parameter_count == 2

    def test_scaffold_outside_trusted_roots(self, scaffold_env):
        with pytest.raises(ValueError, match="outside trusted roots"):
            scaffold_workflow_schema("/etc/passwd")

    def test_scaffold_missing_script(self, scaffold_env):
        with pytest.raises(FileNotFoundError):
            scaffold_workflow_schema("scripts/nonexistent.sh")

    def test_scaffold_has_review_notes(self, scaffold_env):
        script = scaffold_env / "scripts" / "run_example.sh"
        script.parent.mkdir(parents=True)
        script.write_text(MINIMAL_SCRIPT)
        os.chmod(script, 0o755)

        result = scaffold_workflow_schema(str(script))
        assert len(result.notes) > 0
        assert any("best-effort" in n.lower() for n in result.notes)

    def test_scaffold_strips_run_prefix(self, scaffold_env):
        script = scaffold_env / "scripts" / "run_my_pipeline.sh"
        script.parent.mkdir(parents=True)
        script.write_text(MINIMAL_SCRIPT)
        os.chmod(script, 0o755)

        result = scaffold_workflow_schema(str(script))
        parsed = yaml.safe_load(result.draft_yaml)
        assert parsed["id"] == "my_pipeline"

    def test_scaffold_includes_flag_order(self, scaffold_env):
        script = scaffold_env / "scripts" / "run_ordered.sh"
        script.parent.mkdir(parents=True)
        script.write_text(SAMPLE_SCRIPT)
        os.chmod(script, 0o755)

        result = scaffold_workflow_schema(str(script))
        parsed = yaml.safe_load(result.draft_yaml)
        assert "flag_order" in parsed["rendering"]
        assert len(parsed["rendering"]["flag_order"]) == result.parameter_count

    def test_scaffold_yaml_roundtrips_through_validate(self, scaffold_env):
        """Scaffolded YAML should pass validate_draft with only warnings, no errors."""
        script = scaffold_env / "scripts" / "run_pipeline.sh"
        script.parent.mkdir(parents=True)
        script.write_text(SAMPLE_SCRIPT)
        os.chmod(script, 0o755)

        scaffold_result = scaffold_workflow_schema(str(script))
        validate_result = validate_draft_workflow_schema(scaffold_result.draft_yaml)
        assert validate_result.valid, f"Validation errors: {validate_result.errors}"
        assert validate_result.parameter_count == scaffold_result.parameter_count


# ---------------------------------------------------------------------------
# Tests: validate_draft_workflow_schema
# ---------------------------------------------------------------------------


VALID_DRAFT = dedent("""\
    id: "test_wf"
    title: "Test Workflow"
    summary: "A test workflow."
    entry_script: "scripts/test.sh"
    parameters:
      - name: "input_dir"
        cli_flag: "--input-dir"
        type: "directory"
        required: true
        description: "Input directory."
      - name: "threads"
        cli_flag: "--threads"
        type: "int"
        default: 4
        description: "Thread count."
      - name: "dry_run"
        cli_flag: "--dry-run"
        type: "bool"
        default: false
        description: "Dry run mode."
    rendering:
      flag_order: ["input_dir", "threads", "dry_run"]
""")


class TestValidateDraftWorkflowSchema:
    def test_valid_draft_passes(self):
        result = validate_draft_workflow_schema(VALID_DRAFT)
        assert result.valid
        assert result.errors == []
        assert result.parameter_count == 3

    def test_invalid_yaml_syntax(self):
        result = validate_draft_workflow_schema("{{bad yaml")
        assert not result.valid
        assert any("YAML parse" in e for e in result.errors)

    def test_non_dict_yaml(self):
        result = validate_draft_workflow_schema("- a list\n- not a dict")
        assert not result.valid
        assert any("mapping" in e for e in result.errors)

    def test_missing_required_fields(self):
        result = validate_draft_workflow_schema("id: test\n")
        assert not result.valid
        # entry_script is required
        assert any("entry_script" in e for e in result.errors)

    def test_invalid_parameter_type(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "invalid_type"
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid

    def test_todo_warnings(self):
        draft = dedent("""\
            id: "test"
            title: "TODO: fix this"
            summary: "TODO: add summary"
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
                description: "TODO"
        """)
        result = validate_draft_workflow_schema(draft)
        assert result.valid  # TODOs are warnings, not errors
        assert any("title contains TODO" in w for w in result.warnings)
        assert any("summary contains TODO" in w for w in result.warnings)
        assert any("'x' description contains TODO" in w for w in result.warnings)

    def test_flag_order_references_unknown_param(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
            rendering:
              flag_order: ["x", "nonexistent"]
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("nonexistent" in e for e in result.errors)

    def test_constraint_references_unknown_param(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
            constraints:
              - kind: "mutual_exclusion"
                params: ["x", "ghost"]
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("ghost" in e for e in result.errors)

    def test_duplicate_param_names(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
              - name: "x"
                cli_flag: "--y"
                type: "int"
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("duplicate parameter name" in e for e in result.errors)

    def test_duplicate_cli_flags(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--same"
                type: "string"
              - name: "y"
                cli_flag: "--same"
                type: "int"
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("duplicate cli_flag" in e for e in result.errors)

    def test_enum_without_choices(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "mode"
                cli_flag: "--mode"
                type: "enum"
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("enum" in e and "choices" in e for e in result.errors)

    def test_path_must_exist_wrong_type_warns(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
                path_must_exist: true
        """)
        result = validate_draft_workflow_schema(draft)
        assert result.valid  # warning, not error
        assert any("path_must_exist" in w for w in result.warnings)

    def test_must_be_executable_wrong_type_warns(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "directory"
                must_be_executable: true
        """)
        result = validate_draft_workflow_schema(draft)
        assert result.valid  # warning, not error
        assert any("must_be_executable" in w for w in result.warnings)

    def test_no_flag_order_warns(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
        """)
        result = validate_draft_workflow_schema(draft)
        assert result.valid
        assert any("flag_order" in w for w in result.warnings)

    def test_group_references_unknown_param(self):
        draft = dedent("""\
            id: "test"
            title: "Test"
            summary: "Test."
            entry_script: "test.sh"
            parameters:
              - name: "x"
                cli_flag: "--x"
                type: "string"
            parameter_groups:
              - name: "grp"
                title: "Group"
                parameters: ["x", "missing"]
        """)
        result = validate_draft_workflow_schema(draft)
        assert not result.valid
        assert any("missing" in e for e in result.errors)


class TestScaffoldThenValidateRealScript:
    """Test scaffolding the actual UCSF production script if available."""

    def test_scaffold_real_ucsf_script(self, scaffold_env):
        """Scaffold the actual UCSF script and validate the output."""
        real_script = Path("/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh")
        if not real_script.exists():
            pytest.skip("Real UCSF script not available")

        # Copy into scaffold_env trusted root
        import shutil
        dest = scaffold_env / "scripts" / "paper" / real_script.name
        dest.parent.mkdir(parents=True)
        shutil.copy2(real_script, dest)

        result = scaffold_workflow_schema(str(dest))
        assert result.parameter_count > 10  # UCSF script has ~14 params

        # Validate the scaffolded output
        val_result = validate_draft_workflow_schema(result.draft_yaml)
        assert val_result.valid, f"Errors: {val_result.errors}"

    def test_real_ucsf_descriptions_not_garbled(self, scaffold_env):
        """Descriptions should come from usage(), not from ARGS or case blocks."""
        real_script = Path("/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh")
        if not real_script.exists():
            pytest.skip("Real UCSF script not available")

        import shutil
        dest = scaffold_env / "scripts" / "paper" / real_script.name
        dest.parent.mkdir(parents=True)
        shutil.copy2(real_script, dest)

        result = scaffold_workflow_schema(str(dest))
        parsed = yaml.safe_load(result.draft_yaml)
        by_name = {p["name"]: p for p in parsed["parameters"]}

        # These descriptions were garbled before the parser fix
        assert not by_name["star_bin"]["description"].startswith("--")
        assert not by_name["out_root"]["description"].startswith("--")
        assert not by_name["feature_ref"]["description"].startswith("--")
        assert not by_name["solo_cb_whitelist"]["description"].startswith("--")

        # Verify a few known-good descriptions
        assert "Comma-separated" in by_name["samples"]["description"]
        assert by_name["dry_run"]["description"] == "Prepare manifests only"

    def test_real_ucsf_defaults_cleaned(self, scaffold_env):
        """Shell variable defaults should be cleaned, not passed raw."""
        real_script = Path("/mnt/pikachu/STAR-suite/scripts/paper/run_ucsf_corrected_production_workflow.sh")
        if not real_script.exists():
            pytest.skip("Real UCSF script not available")

        import shutil
        dest = scaffold_env / "scripts" / "paper" / real_script.name
        dest.parent.mkdir(parents=True)
        shutil.copy2(real_script, dest)

        result = scaffold_workflow_schema(str(dest))
        parsed = yaml.safe_load(result.draft_yaml)
        by_name = {p["name"]: p for p in parsed["parameters"]}

        # star_bin default should be relative path, not ${REPO_ROOT}/...
        star_default = by_name["star_bin"].get("default")
        if star_default is not None:
            assert "${" not in str(star_default)

        # cellbender_cpu_cores default was ${THREADS} — should be null
        cpu_default = by_name.get("cellbender_cpu_cores", {}).get("default")
        assert cpu_default is None

        # out_root default had $(date ...) — should be null
        out_default = by_name["out_root"].get("default")
        assert out_default is None
