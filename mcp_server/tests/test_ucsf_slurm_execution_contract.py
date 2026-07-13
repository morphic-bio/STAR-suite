from __future__ import annotations

import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import yaml
from scipy import sparse

from mcp_server.schemas.workflow import WorkflowSchema


REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "mcp_server" / "workflows" / "ucsf_star_suite_production.yaml"


def test_ucsf_workflow_publishes_per_sample_cellbender_decomposition():
    workflow = WorkflowSchema.model_validate(yaml.safe_load(SCHEMA_PATH.read_text()))

    assert workflow.execution is not None
    assert workflow.execution.decomposition_id == "ucsf_star_suite_per_sample_cellbender/v1"
    assert workflow.execution.scatter.parameter == "samples"
    stages = {stage.name: stage for stage in workflow.execution.stages}
    assert set(stages) == {
        "star_alignment",
        "downstream_prepare",
        "cellbender",
        "downstream_finalize",
        "gather",
    }
    assert stages["star_alignment"].resource_class == "star_cpu_highmem"
    assert stages["cellbender"].resource_class == "cellbender_gpu"
    assert stages["cellbender"].foreach == "sample"
    assert workflow.execution.gather.stage == "gather"


def test_ucsf_slurm_stage_entrypoints_have_reviewable_help():
    scripts = [
        REPO_ROOT / "scripts" / "slurm" / "run_ucsf_downstream_prepare_stage.sh",
        REPO_ROOT / "scripts" / "slurm" / "run_ucsf_cellbender_stage.sh",
        REPO_ROOT / "scripts" / "slurm" / "run_ucsf_downstream_finalize_stage.sh",
        REPO_ROOT / "scripts" / "slurm" / "gather_ucsf_tutorial.py",
        REPO_ROOT / "scripts" / "slurm" / "combine_ucsf_filters.py",
    ]
    for script in scripts:
        result = subprocess.run([str(script), "--help"], text=True, capture_output=True, check=False)
        assert result.returncode == 0, (script, result.stderr)
        assert "usage" in result.stdout.lower()


def test_ucsf_prepare_stage_has_no_hidden_image_helper_dependency():
    source = (REPO_ROOT / "scripts" / "slurm" / "run_ucsf_downstream_prepare_stage.sh").read_text()

    assert "/usr/local/bin/combineFilters.py" not in source
    assert "${TOOLS_DIR}/combine_ucsf_filters.py" in source


def test_ucsf_prepare_stage_applies_adaptive_max_before_baseline_masks():
    source = (REPO_ROOT / "scripts" / "slurm" / "run_ucsf_downstream_prepare_stage.sh").read_text()

    threshold_call = source.index('python3 "${COMPUTE_ADAPTIVE_QC}"')
    combine_call = source.index('python3 "${COMBINE_FILTERS}"')
    assert threshold_call < combine_call
    assert '--max-genes "${EFFECTIVE_MAX_GENES}"' in source
    assert '--doublet-scores "${OUTPUT_DIR}/filtered_barcodes_with_scores.txt"' in source


def test_cellbender_stage_relocates_temporary_files_to_the_output_volume():
    source = (REPO_ROOT / "scripts" / "slurm" / "run_ucsf_cellbender_stage.sh").read_text()

    assert 'export TMPDIR="${OUTPUT_DIR}/.tmp"' in source


def test_doublet_stage_uses_a_stable_configurable_seed():
    source = (REPO_ROOT / "scripts" / "run_star_cell_doublets.R").read_text()

    assert 'Sys.getenv("SCDBLFINDER_SEED", "1")' in source
    assert "set.seed(scdblfinder_seed)" in source
    assert "marking all STAR cells as singlets" not in source


def test_feature_integration_defers_annotations_for_cellbender_python():
    source = (REPO_ROOT / "scripts" / "integrate_feature_library.py").read_text()

    assert source.splitlines()[1] == "from __future__ import annotations"


def test_portable_filter_computes_mito_percentage_from_total_counts(tmp_path: Path):
    counts = ad.AnnData(
        X=sparse.csr_matrix(np.array([[2, 8]], dtype=np.int32)),
        obs={"sample": ["fixture"]},
        var={"gene": ["mt", "nuclear"]},
    )
    counts.obs_names = ["barcode-1"]
    counts.var_names = ["ENSG00000198888", "ENSG00000000001"]
    input_path = tmp_path / "counts.h5ad"
    non_empty = tmp_path / "non_empty.txt"
    doublets = tmp_path / "doublets.txt"
    doublet_scores = tmp_path / "doublet_scores.txt"
    output_dir = tmp_path / "output"
    counts.write_h5ad(input_path)
    non_empty.write_text("barcode-1\n")
    doublets.write_text("")
    doublet_scores.write_text(
        "Barcode\tClassification\tScore\nbarcode-1\tsinglet\t0.125\n"
    )

    subprocess.run(
        [
            str(REPO_ROOT / "scripts" / "slurm" / "combine_ucsf_filters.py"),
            "--input-file",
            str(input_path),
            "--non-empty-barcodes",
            str(non_empty),
            "--doublet-barcodes",
            str(doublets),
            "--doublet-scores",
            str(doublet_scores),
            "--output-dir",
            str(output_dir),
            "--min-genes",
            "1",
            "--max-genes",
            "3",
            "--mt-pct-cutoff",
            "50",
        ],
        check=True,
    )

    result = ad.read_h5ad(output_dir / "unfiltered_counts.h5ad")
    assert result.obs.loc["barcode-1", "total_counts"] == 10
    assert result.obs.loc["barcode-1", "mt_pct"] == 20
    assert result.obs.loc["barcode-1", "doublet_scores"] == 0.125
    assert bool(result.obs.loc["barcode-1", "filter"])


def test_portable_filter_rejects_incomplete_doublet_scores(tmp_path: Path):
    counts = ad.AnnData(
        X=sparse.csr_matrix([[2]], dtype=np.int32),
        obs={"sample": ["fixture"]},
        var={"gene": ["gene-1"]},
    )
    counts.obs_names = ["barcode-1"]
    counts.var_names = ["gene-1"]
    input_path = tmp_path / "counts.h5ad"
    non_empty = tmp_path / "non_empty.txt"
    doublets = tmp_path / "doublets.txt"
    doublet_scores = tmp_path / "doublet_scores.txt"
    output_dir = tmp_path / "output"
    counts.write_h5ad(input_path)
    non_empty.write_text("barcode-1\n")
    doublets.write_text("")
    doublet_scores.write_text("Barcode\tClassification\tScore\n")

    result = subprocess.run(
        [
            str(REPO_ROOT / "scripts" / "slurm" / "combine_ucsf_filters.py"),
            "--input-file",
            str(input_path),
            "--non-empty-barcodes",
            str(non_empty),
            "--doublet-barcodes",
            str(doublets),
            "--doublet-scores",
            str(doublet_scores),
            "--output-dir",
            str(output_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "cover exactly the non-empty barcode set" in result.stderr
