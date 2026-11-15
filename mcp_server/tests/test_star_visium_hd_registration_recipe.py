"""Contract tests for the star_visium_hd_registration recipe
(scripts/run_visium_hd_registration.sh, a wrapper for scripts/run_registration.sh
of BioDepot/stargrid).

Schema load -> validate -> render -> run the rendered argv with --dry-run
against a stub stargrid checkout, stub tools and stub inputs.
"""

import os
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
WORKFLOW_ID = "star_visium_hd_registration"
WRAPPER = REPO_ROOT / "scripts" / "run_visium_hd_registration.sh"
TOOLS = (
    "visium_hd_tiff_decimate",
    "visium_hd_capture_grid",
    "visium_hd_register",
    "visium_hd_count_register",
)
REAL_SCRIPT = Path(os.environ.get("STARGRID_REPO", "/nonexistent")) / "scripts" / "run_registration.sh"


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
    """A stub stargrid checkout whose run_registration.sh is the real one when
    available, else a stub that echoes its arguments; stub tools; stub inputs."""
    stargrid = tmp_path / "stargrid"
    (stargrid / "scripts").mkdir(parents=True)
    script = stargrid / "scripts" / "run_registration.sh"
    script.write_text("#!/usr/bin/env bash\necho \"ARGS: $*\"\n")
    script.chmod(0o755)
    bins = tmp_path / "bin"
    bins.mkdir()
    for name in TOOLS:
        tool = bins / name
        tool.write_text("#!/bin/sh\nexit 0\n")
        tool.chmod(0o755)
    he = tmp_path / "he.btf"
    he.write_text("x")
    cytassist = tmp_path / "cytassist.tif"
    cytassist.write_text("x")
    layout = tmp_path / "H1-TEST.vlf"
    layout.write_text('{"areas": ["D1"]}\n')
    counts = tmp_path / "SpatialGex.out"
    for size in ("002", "008", "016"):
        d = counts / "hard" / f"square_{size}um"
        d.mkdir(parents=True)
        (d / "matrix.mtx").write_text("")
    return {
        "stargrid": str(stargrid),
        "bin_dir": str(bins),
        "he": str(he),
        "cytassist": str(cytassist),
        "layout": str(layout),
        "counts": str(counts),
        "out_dir": str(tmp_path / "out"),
    }


def params_for(inputs):
    return {
        "stargrid": inputs["stargrid"],
        "he": inputs["he"],
        "cytassist": inputs["cytassist"],
        "layout": inputs["layout"],
        "area": "D1",
        "counts": inputs["counts"],
        "out_dir": inputs["out_dir"],
        "bin_dir": inputs["bin_dir"],
        "dry_run": True,
    }


def test_schema_loads_with_the_script_options(recipe):
    from mcp_server.tools.workflows import get_workflow_parameter_schema

    schema = get_workflow_parameter_schema(WORKFLOW_ID)
    names = {p.name for p in schema.parameters}
    assert names == {
        "stargrid", "he", "cytassist", "layout", "area", "counts", "out_dir",
        "segmentation", "pixel_size_um", "cellpose_model", "python", "container",
        "cpu", "threads", "bin_dir", "dry_run",
    }
    required = {p.name for p in schema.parameters if p.required}
    assert required == {"stargrid", "he", "cytassist", "layout", "area", "counts", "out_dir"}


def test_complete_parameters_validate_and_render(recipe, inputs):
    params = params_for(inputs)
    assert validate_workflow_parameters(WORKFLOW_ID, params).valid
    rendered = render_workflow_command(WORKFLOW_ID, params)
    assert rendered.argv == [
        str(WRAPPER),
        "--stargrid", inputs["stargrid"],
        "--he", inputs["he"],
        "--cytassist", inputs["cytassist"],
        "--layout", inputs["layout"],
        "--area", "D1",
        "--counts", inputs["counts"],
        "--out-dir", inputs["out_dir"],
        "--threads", "16",
        "--bin-dir", inputs["bin_dir"],
        "--dry-run",
    ]
    assert rendered.output_root == inputs["out_dir"]


def test_nuclei_options_render_in_order(recipe, inputs, tmp_path):
    seg = tmp_path / "seg"
    seg.mkdir()
    model = tmp_path / "cpsam"
    model.write_text("x")
    params = params_for(inputs)
    params.update(segmentation=str(seg), pixel_size_um=0.2738, cellpose_model=str(model),
                  python="/usr/bin/python3", container="img:tag", cpu=True, threads=12)
    assert validate_workflow_parameters(WORKFLOW_ID, params).valid
    argv = render_workflow_command(WORKFLOW_ID, params).argv
    tail = argv[argv.index("--out-dir") + 2:]
    assert tail == [
        "--segmentation", str(seg), "--pixel-size-um", "0.2738",
        "--cellpose-model", str(model), "--python", "/usr/bin/python3",
        "--container", "img:tag", "--cpu", "--threads", "12",
        "--bin-dir", inputs["bin_dir"], "--dry-run",
    ]


@pytest.mark.parametrize("missing", ["stargrid", "he", "cytassist", "layout", "area", "counts", "out_dir"])
def test_missing_required_parameter_is_refused(recipe, inputs, missing):
    params = params_for(inputs)
    del params[missing]
    result = validate_workflow_parameters(WORKFLOW_ID, params)
    assert not result.valid
    assert any(missing in error for error in result.errors)


@pytest.mark.parametrize("bad", [{"threads": 0}, {"threads": -4}, {"pixel_size_um": 0},
                                 {"pixel_size_um": -0.27}, {"threads": "many"}])
def test_bad_values_are_refused(recipe, inputs, bad):
    params = params_for(inputs)
    params.update(bad)
    result = validate_workflow_parameters(WORKFLOW_ID, params)
    assert not result.valid
    assert any(next(iter(bad)) in error for error in result.errors)


def test_missing_input_files_are_refused(recipe, inputs, tmp_path):
    params = params_for(inputs)
    params["he"] = str(tmp_path / "absent.btf")
    params["stargrid"] = str(tmp_path / "no_such_checkout")
    result = validate_workflow_parameters(WORKFLOW_ID, params)
    assert not result.valid
    assert any("he" in e and "does not exist" in e for e in result.errors)
    assert any("stargrid" in e and "does not exist" in e for e in result.errors)


def test_rendered_argv_reaches_run_registration_unchanged(recipe, inputs):
    params = params_for(inputs)
    rendered = render_workflow_command(WORKFLOW_ID, params)
    result = subprocess.run(rendered.argv, capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "ARGS: " + " ".join(rendered.argv[3:])
    assert "--stargrid" not in result.stdout


def test_stargrid_comes_from_the_environment(inputs):
    env = {"PATH": "/usr/bin:/bin", "STARGRID_REPO": inputs["stargrid"]}
    result = subprocess.run([str(WRAPPER), "--he", "a", "--dry-run"], capture_output=True,
                            text=True, check=True, env=env)
    assert result.stdout.strip() == "ARGS: --he a --dry-run"


def test_wrapper_fails_clearly_without_a_checkout(tmp_path):
    env = {"PATH": "/usr/bin:/bin"}
    result = subprocess.run([str(WRAPPER), "--he", "a"], capture_output=True, text=True, env=env)
    assert result.returncode == 2 and "STARGRID_REPO" in result.stderr
    result = subprocess.run([str(WRAPPER), "--stargrid", str(tmp_path / "none")],
                            capture_output=True, text=True, env=env)
    assert result.returncode == 2 and "not a folder" in result.stderr
    (tmp_path / "empty").mkdir()
    result = subprocess.run([str(WRAPPER), "--stargrid", str(tmp_path / "empty")],
                            capture_output=True, text=True, env=env)
    assert result.returncode == 2 and "run_registration.sh not found" in result.stderr


def test_wrapper_prints_usage(inputs):
    result = subprocess.run([str(WRAPPER), "--help"], capture_output=True, text=True,
                            check=True, env={"PATH": "/usr/bin:/bin"})
    assert "--stargrid DIR" in result.stdout


@pytest.mark.skipif(not REAL_SCRIPT.is_file(), reason="set STARGRID_REPO to a checkout of BioDepot/stargrid")
def test_dry_run_with_the_real_script_prints_the_stages(recipe, inputs):
    params = params_for(inputs)
    params["stargrid"] = str(REAL_SCRIPT.parent.parent)
    params["container"] = "some-image:tag"
    params["pixel_size_um"] = 0.27
    cellpose_model = Path(inputs["he"]).parent / "cpsam"
    cellpose_model.write_text("x")
    params["cellpose_model"] = str(cellpose_model)
    rendered = render_workflow_command(WORKFLOW_ID, params)
    result = subprocess.run(rendered.argv, capture_output=True, text=True, check=True)
    for tool in TOOLS:
        assert tool in result.stdout
    assert "segment_hd_he_cellpose.py" in result.stdout
    assert "# in container some-image:tag" in result.stdout
    assert not Path(inputs["out_dir"]).exists()
    params["area"] = "Z9"
    bad = subprocess.run(render_workflow_command(WORKFLOW_ID, params).argv,
                         capture_output=True, text=True)
    assert bad.returncode == 2 and "capture area" in bad.stderr
