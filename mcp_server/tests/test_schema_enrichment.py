"""Tests for enriched schema metadata on star_genome_generate and star_bulk_pe_batch.

Verifies that:
- advisory fields (label, help, example, widget_hint, aliases, advanced,
  display_order, min_value, max_value) are present where defined
- canonical fields (name, type, required, default, description) remain stable
- alias resolution works in validate_workflow_parameters
- structured field_errors are populated alongside flat errors/warnings
- parameter group descriptions are surfaced
- backward compatibility: old consumers see the same core shape
"""

import tempfile
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import (
    get_workflow_parameter_schema,
    validate_workflow_parameters,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _load_real_workflow(workflow_id: str, tmp_path: Path):
    """Load a real workflow YAML from the repo into the test config."""
    repo_root = Path(__file__).resolve().parent.parent.parent  # STAR-suite root
    yaml_path = repo_root / "mcp_server" / "workflows" / f"{workflow_id}.yaml"
    assert yaml_path.exists(), f"Workflow YAML not found: {yaml_path}"

    # Read title and entry_script from the YAML schema itself
    with open(yaml_path) as f:
        wf_data = yaml.safe_load(f)

    config_dict = {
        "server": {
            "host": "127.0.0.1",
            "port": 9999,
            "auth_token": "",
            "transport": "http",
        },
        "paths": {
            "repo_root": str(repo_root),
            "artifact_log_root": str(tmp_path / "artifacts"),
            "temp_root": str(tmp_path / "tmp"),
        },
        "trusted_roots": [str(repo_root), str(tmp_path), "/tmp", "/storage"],
        "datasets": [],
        "scripts": [],
        "test_suites": [],
        "required_binaries": [],
        "workflows": [
            {
                "id": workflow_id,
                "title": wf_data.get("title", workflow_id),
                "entry_script": wf_data.get("entry_script", ""),
                "schema_file": str(yaml_path),
                "visibility": "public",
            }
        ],
    }

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)

    load_config(config_path)


@pytest.fixture
def genome_generate_config(tmp_path):
    _load_real_workflow("star_genome_generate", tmp_path)
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def bulk_pe_batch_config(tmp_path):
    _load_real_workflow("star_bulk_pe_batch", tmp_path)
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


# ---------------------------------------------------------------------------
# star_genome_generate enrichment tests
# ---------------------------------------------------------------------------


class TestGenomeGenerateEnrichment:
    """Advisory metadata on star_genome_generate."""

    def test_canonical_fields_stable(self, genome_generate_config):
        """Core fields (name, type, required, default, description) remain present."""
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert "run_mode" in params
        assert params["run_mode"].type == "string"
        assert params["run_mode"].required is True
        assert params["run_mode"].default == "genomeGenerate"

        assert "genome_dir" in params
        assert params["genome_dir"].type == "directory"
        assert params["genome_dir"].required is True

        assert "genome_fasta" in params
        assert params["genome_fasta"].type == "file"
        assert params["genome_fasta"].required is True
        assert params["genome_fasta"].path_must_exist is True

        assert "run_thread_n" in params
        assert params["run_thread_n"].type == "int"
        assert params["run_thread_n"].default == 8

    def test_label_present(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].label == "Run mode"
        assert params["genome_dir"].label == "Genome directory"
        assert params["genome_fasta"].label == "Reference FASTA"
        assert params["run_thread_n"].label == "Threads"

    def test_help_present(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["genome_dir"].help is not None
        assert len(params["genome_dir"].help) > len(params["genome_dir"].description)

    def test_example_present(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["genome_dir"].example == "/tmp/star-demo-index"
        assert params["genome_fasta"].example is not None
        assert ".fa" in params["genome_fasta"].example

    def test_widget_hint_present(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].widget_hint == "readonly"
        assert params["genome_dir"].widget_hint == "path"
        assert params["genome_fasta"].widget_hint == "path"
        assert params["run_thread_n"].widget_hint == "number"

    def test_aliases_present(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].aliases == ["runMode"]
        assert params["genome_dir"].aliases == ["genomeDir"]
        assert params["genome_fasta"].aliases == ["genomeFastaFiles"]
        assert params["run_thread_n"].aliases == ["runThreadN"]

    def test_advanced_flag(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].advanced is False
        assert params["run_thread_n"].advanced is True

    def test_display_order(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].display_order == 0
        assert params["genome_dir"].display_order == 1
        assert params["genome_fasta"].display_order == 2
        assert params["run_thread_n"].display_order == 3

    def test_min_max_value(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        params = {p.name: p for p in result.parameters}

        assert params["run_thread_n"].min_value == 1
        assert params["run_thread_n"].max_value == 128

    def test_group_description(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        groups = {g.name: g for g in result.parameter_groups}

        assert "core" in groups
        assert groups["core"].description != ""
        assert "index" in groups["core"].description.lower()


# ---------------------------------------------------------------------------
# star_bulk_pe_batch enrichment tests
# ---------------------------------------------------------------------------


class TestBulkPeBatchEnrichment:
    """Advisory metadata on star_bulk_pe_batch."""

    def test_canonical_fields_stable(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert "genome_dir" in params
        assert params["genome_dir"].type == "directory"
        assert params["genome_dir"].required is True
        assert params["genome_dir"].path_must_exist is True

        assert "read_files_mate1" in params
        assert params["read_files_mate1"].operand_group == "readfiles_pe"

        assert "run_thread_n" in params
        assert params["run_thread_n"].default == 16

    def test_label_present(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert params["genome_dir"].label == "Genome directory"
        assert params["read_files_mate1"].label == "Mate 1 FASTQs"
        assert params["read_files_mate2"].label == "Mate 2 FASTQs"
        assert params["out_file_name_prefix"].label == "Output prefix"

    def test_help_present(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert params["batch_mode"].help is not None
        assert "genome" in params["batch_mode"].help.lower()

    def test_example_present(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert params["read_files_mate1"].example is not None
        assert "R1" in params["read_files_mate1"].example

    def test_widget_hint_present(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert params["run_mode"].widget_hint == "readonly"
        assert params["genome_dir"].widget_hint == "path"
        assert params["run_thread_n"].widget_hint == "number"

    def test_aliases_on_bulk(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        assert params["genome_dir"].aliases == ["genomeDir"]
        assert params["run_thread_n"].aliases == ["runThreadN"]
        assert params["batch_mode"].aliases == ["batchMode"]
        assert params["out_file_name_prefix"].aliases == ["outFileNamePrefix"]

    def test_no_alias_on_multi_operand_params(self, bulk_pe_batch_config):
        """Multi-operand params (readFilesIn, outSAMtype) must not carry aliases
        because the alias maps one key to one param, losing the other operand."""
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        # readFilesIn is shared by mate1 and mate2 via operand_group — no alias
        assert params["read_files_mate1"].aliases is None
        assert params["read_files_mate2"].aliases is None

        # outSAMtype is shared by out_sam_kind and out_sam_sort — no alias
        assert params["out_sam_kind"].aliases is None
        assert params["out_sam_sort"].aliases is None

    def test_grouped_sections(self, bulk_pe_batch_config):
        """Groups have been reorganized into inputs / output / performance."""
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        group_names = [g.name for g in result.parameter_groups]

        assert "inputs" in group_names
        assert "output" in group_names
        assert "performance" in group_names

    def test_group_descriptions(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        groups = {g.name: g for g in result.parameter_groups}

        for name in ["inputs", "output", "performance"]:
            assert groups[name].description != "", f"Group '{name}' missing description"

    def test_display_order_sequential(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        orders = [p.display_order for p in result.parameters if p.display_order is not None]
        assert orders == sorted(orders), "display_order values should be sequential"
        assert len(set(orders)) == len(orders), "display_order values should be unique"

    def test_advanced_flags(self, bulk_pe_batch_config):
        result = get_workflow_parameter_schema("star_bulk_pe_batch")
        params = {p.name: p for p in result.parameters}

        # Core user-facing params are not advanced
        assert params["genome_dir"].advanced is False
        assert params["read_files_mate1"].advanced is False

        # Tuning params are advanced
        assert params["run_thread_n"].advanced is True
        assert params["batch_mode"].advanced is True
        assert params["out_sam_sort"].advanced is True


# ---------------------------------------------------------------------------
# Alias normalization tests
# ---------------------------------------------------------------------------


class TestAliasNormalization:
    """Alias resolution in validate_workflow_parameters."""

    def test_alias_resolves_to_canonical(self, genome_generate_config):
        """Passing a camelCase alias should be normalized to the canonical name."""
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "runMode": "genomeGenerate",
                "genomeDir": "/tmp/test-idx",
                "genomeFastaFiles": "/tmp/fake.fa",
                "runThreadN": 4,
            },
            check_paths=False,
        )
        # Aliases should be resolved to canonical names in normalized_params
        assert "run_mode" in result.normalized_params
        assert "genome_dir" in result.normalized_params
        assert "genome_fasta" in result.normalized_params
        assert "run_thread_n" in result.normalized_params

        # Alias keys should not appear
        assert "runMode" not in result.normalized_params
        assert "genomeDir" not in result.normalized_params

    def test_alias_and_canonical_collision_is_error(self, genome_generate_config):
        """Supplying both alias and canonical for the same param is a collision error."""
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "runMode": "genomeGenerate",
                "genomeDir": "/tmp/alias-idx",
                "genome_dir": "/tmp/canonical-idx",
                "genomeFastaFiles": "/tmp/fake.fa",
                "runThreadN": 4,
            },
            check_paths=False,
        )
        assert result.valid is False
        collision_msgs = [e for e in result.errors if "Alias collision" in e]
        assert len(collision_msgs) == 1
        assert "genome_dir" in collision_msgs[0]


# ---------------------------------------------------------------------------
# Structured field_errors tests
# ---------------------------------------------------------------------------


class TestFieldErrors:
    """Structured field_errors alongside flat errors/warnings."""

    def test_missing_required_has_field_error(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {"run_mode": "genomeGenerate"},
            check_paths=False,
        )
        assert result.valid is False
        assert len(result.errors) > 0

        # field_errors should contain entries for missing required params
        field_names = [fe.field for fe in result.field_errors]
        assert "genome_dir" in field_names
        assert "genome_fasta" in field_names

    def test_type_error_has_field_error(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "run_mode": "genomeGenerate",
                "genome_dir": "/tmp/test",
                "genome_fasta": "/tmp/fake.fa",
                "run_thread_n": "not_a_number",
            },
            check_paths=False,
        )
        assert result.valid is False

        field_names = [fe.field for fe in result.field_errors]
        assert "run_thread_n" in field_names

        # Check the kind is "error"
        for fe in result.field_errors:
            if fe.field == "run_thread_n":
                assert fe.kind == "error"

    def test_positive_constraint_field_error(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "run_mode": "genomeGenerate",
                "genome_dir": "/tmp/test",
                "genome_fasta": "/tmp/fake.fa",
                "run_thread_n": -1,
            },
            check_paths=False,
        )
        assert result.valid is False

        field_names = [fe.field for fe in result.field_errors]
        assert "run_thread_n" in field_names

    def test_valid_params_no_field_errors(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "run_mode": "genomeGenerate",
                "genome_dir": "/tmp/test",
                "genome_fasta": "/tmp/fake.fa",
                "run_thread_n": 4,
            },
            check_paths=False,
        )
        assert result.valid is True
        assert len(result.field_errors) == 0
        assert len(result.errors) == 0

    def test_unknown_param_field_error(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {
                "run_mode": "genomeGenerate",
                "genome_dir": "/tmp/test",
                "genome_fasta": "/tmp/fake.fa",
                "bogus_param": "value",
            },
            check_paths=False,
        )
        assert result.valid is False
        field_names = [fe.field for fe in result.field_errors]
        assert "bogus_param" in field_names


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------


class TestBackwardCompatibility:
    """Old consumers expecting name, type, required, default, description still work."""

    def test_old_shape_preserved(self, genome_generate_config):
        result = get_workflow_parameter_schema("star_genome_generate")
        data = result.model_dump()

        for p in data["parameters"]:
            assert "name" in p
            assert "type" in p
            assert "required" in p
            assert "default" in p
            assert "description" in p
            assert "cli_flag" in p

    def test_validation_response_old_fields(self, genome_generate_config):
        result = validate_workflow_parameters(
            "star_genome_generate",
            {"run_mode": "genomeGenerate"},
            check_paths=False,
        )
        data = result.model_dump()

        assert "valid" in data
        assert "normalized_params" in data
        assert "warnings" in data
        assert "errors" in data
        # field_errors is additive, not replacing errors
        assert isinstance(data["field_errors"], list)

    def test_advisory_fields_are_optional_in_model_dump(self, genome_generate_config):
        """Advisory fields default to None/False so old consumers can ignore them."""
        result = get_workflow_parameter_schema("star_genome_generate")
        data = result.model_dump()

        # All advisory fields should be present in the dict (Pydantic includes defaults)
        for p in data["parameters"]:
            # These are the advisory keys -- they should exist but can be null
            for key in ["label", "help", "example", "widget_hint", "aliases",
                        "advanced", "display_order", "min_value", "max_value"]:
                assert key in p, f"Advisory field '{key}' missing from parameter '{p['name']}'"
