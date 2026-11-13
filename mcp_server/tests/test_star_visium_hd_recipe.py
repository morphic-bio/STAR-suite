"""Contract tests for the star_visium_hd recipe (scripts/run_visium_hd.sh).

Schema load -> validate -> render -> run the rendered argv with --dry-run
against stub inputs and a stub STAR binary.
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

import mcp_server.config as config_module
from mcp_server.config import load_config
from mcp_server.tools.workflows import (
    render_workflow_command,
    validate_workflow_parameters,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
WORKFLOW_ID = "star_visium_hd"
CONTRACT_FILES = (
    "barcode_coords.tsv",
    "bc1_full_oligos.txt",
    "bc2_full_oligos.txt",
    "bc1_whitelist.txt",
    "bc2_whitelist.txt",
)


@pytest.fixture
def recipe(tmp_path):
    yaml_path = REPO_ROOT / "mcp_server" / "workflows" / f"{WORKFLOW_ID}.yaml"
    with open(yaml_path) as f:
        wf_data = yaml.safe_load(f)
    config_dict = {
        "server": {"host": "127.0.0.1", "port": 9999, "auth_token": "", "transport": "http"},
        "paths": {
            "repo_root": str(REPO_ROOT),
            "artifact_log_root": str(tmp_path / "artifacts"),
            "temp_root": str(tmp_path / "tmp"),
        },
        "trusted_roots": [str(REPO_ROOT), str(tmp_path), "/tmp"],
        "datasets": [],
        "scripts": [],
        "test_suites": [],
        "required_binaries": [],
        "workflows": [
            {
                "id": WORKFLOW_ID,
                "title": wf_data["title"],
                "entry_script": wf_data["entry_script"],
                "kind": wf_data["kind"],
                "schema_file": str(yaml_path),
                "visibility": "public",
            }
        ],
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)
    load_config(config_path)
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def inputs(tmp_path):
    fastqs = tmp_path / "fastqs"
    fastqs.mkdir()
    for lane in ("L001", "L002"):
        for read in ("R1", "R2"):
            with gzip.open(fastqs / f"Slide_S1_{lane}_{read}_001.fastq.gz", "wb") as f:
                f.write(b"@r\nACGT\n+\nIIII\n")
    contract = tmp_path / "contract"
    contract.mkdir()
    for name in CONTRACT_FILES:
        (contract / name).write_text("x\n")
    genome = tmp_path / "genome"
    genome.mkdir()
    (genome / "SA").write_text("")
    star = tmp_path / "STAR"
    star.write_text("#!/bin/sh\nexit 0\n")
    star.chmod(0o755)
    probe_set = tmp_path / "probe_set.csv"
    probe_set.write_text("gene_id,probe_seq,probe_id\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("")
    return {
        "fastq_dir": str(fastqs),
        "contract": str(contract),
        "genome_dir": str(genome),
        "star": str(star),
        "probe_set": str(probe_set),
        "gtf": str(gtf),
        "out_dir": str(tmp_path / "out"),
    }


def gex_params(inputs):
    return {
        "mode": "gex",
        "fastq_dir": inputs["fastq_dir"],
        "contract": inputs["contract"],
        "out_dir": inputs["out_dir"],
        "genome_dir": inputs["genome_dir"],
        "star": inputs["star"],
        "dry_run": True,
    }


def test_gex_renders_the_wrapper_argv(recipe, inputs):
    params = gex_params(inputs)
    assert validate_workflow_parameters(WORKFLOW_ID, params).valid
    rendered = render_workflow_command(WORKFLOW_ID, params)
    assert rendered.argv == [
        str(REPO_ROOT / "scripts" / "run_visium_hd.sh"),
        "--mode", "gex",
        "--fastq-dir", inputs["fastq_dir"],
        "--contract", inputs["contract"],
        "--out-dir", inputs["out_dir"],
        "--genome-dir", inputs["genome_dir"],
        "--threads", "16",
        "--star", inputs["star"],
        "--dry-run",
    ]
    assert rendered.output_root == inputs["out_dir"]


def test_slide_or_contract_is_required(recipe, inputs):
    params = gex_params(inputs)
    del params["contract"]
    result = validate_workflow_parameters(WORKFLOW_ID, params)
    assert not result.valid
    assert any("slide" in error and "contract" in error for error in result.errors)


def test_gex_dry_run_prints_the_star_recipe(recipe, inputs):
    rendered = render_workflow_command(WORKFLOW_ID, gex_params(inputs))
    result = subprocess.run(rendered.argv, capture_output=True, text=True, check=True)
    assert "--soloSpatialGexIntegrated yes" in result.stdout
    assert "Slide_S1_L001_R2_001.fastq.gz" in result.stdout
    assert not Path(inputs["out_dir"]).exists()


@pytest.mark.skipif(shutil.which("g++") is None, reason="building a Flex reference needs g++")
def test_flex_dry_run_builds_a_missing_reference(recipe, inputs, tmp_path):
    params = gex_params(inputs)
    del params["genome_dir"]
    params.update(
        mode="flex",
        flex_reference=str(tmp_path / "flex_reference"),
        probe_set=inputs["probe_set"],
        gtf=inputs["gtf"],
        reads=1000,
    )
    assert validate_workflow_parameters(WORKFLOW_ID, params).valid
    rendered = render_workflow_command(WORKFLOW_ID, params)
    result = subprocess.run(rendered.argv, capture_output=True, text=True, check=True)
    assert "--runMode hashCacheGenerate" in result.stdout
    assert "--soloSpatialFlexIntegrated yes" in result.stdout
    assert "--soloSpatialExpectedReads 1000 --soloSpatialExpectedCandidates 1500" in result.stdout
    assert not (tmp_path / "flex_reference").exists()


def test_missing_genome_index_is_rejected(recipe, inputs):
    params = gex_params(inputs)
    del params["genome_dir"]
    rendered = render_workflow_command(WORKFLOW_ID, params)
    result = subprocess.run(rendered.argv, capture_output=True, text=True)
    assert result.returncode == 2
    assert "--genome-dir" in result.stderr
