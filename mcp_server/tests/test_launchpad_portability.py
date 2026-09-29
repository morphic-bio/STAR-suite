"""Portable Launchpad contracts, using tiny processes rather than biological runs."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import socket
import sys
import time
from urllib.request import urlopen

import pytest
import yaml
from starlette.testclient import TestClient

from mcp_server import config as config_module
from mcp_server.app import build_http_app
from mcp_server.config import get_workflow_schema, load_config
from mcp_server.launchpad.multiome import FASTQ_GROUPS, REFERENCES, runtime_check
from mcp_server.tools.workflows import render_workflow_command, validate_workflow_parameters

ROOT = Path(__file__).resolve().parents[2]
API = "/launchpad/api/workflows/morphic_multiome"


@pytest.fixture
def portable(tmp_path, monkeypatch):
    root = tmp_path / "relocated suite"
    (root / "mcp_server").mkdir(parents=True)
    shutil.copytree(ROOT / "mcp_server/workflows", root / "mcp_server/workflows")
    shutil.copytree(ROOT / "share", root / "share")
    raw = yaml.safe_load((ROOT / "mcp_server/config.yaml").read_text())
    raw["paths"]["artifact_log_root"] = str(tmp_path / "records")
    path = root / "mcp_server/config.yaml"
    path.write_text(yaml.safe_dump(raw))
    for key in ("_config", "_config_path", "_config_loaded_at", "_workflow_schemas",
                "_workflow_configs", "_workflow_origins", "_recipe_catalogs"):
        monkeypatch.setattr(config_module, key, getattr(config_module, key))
    monkeypatch.chdir(tmp_path)  # Deliberately outside both checkout and config directory.
    load_config(path)
    return root


@pytest.fixture
def inputs(tmp_path):
    params = {name: str(tmp_path / name) for group in FASTQ_GROUPS for name in group}
    params.update({name: str(tmp_path / name) for name in REFERENCES})
    for name, value in params.items():
        if name == "genome_dir":
            Path(value).mkdir()
        else:
            Path(value).write_text("fixture\n")
    params["out_dir"] = str(tmp_path / "new run")
    return params


@pytest.fixture
def fake_runtime(portable, inputs, tmp_path):
    star = tmp_path / "STAR"
    star.write_text('#!/usr/bin/env python3\nprint(\'{"suite_version":"1.9.5.b","chromap_atac":true}\')\n')
    star.chmod(0o755)
    helper = tmp_path / "star_multiome_atac_peak_mex"
    helper.write_text("#!/bin/sh\nexit 0\n")
    helper.chmod(0o755)
    return {**inputs, "star_bin": str(star), "atac_peak_mex_bin": str(helper)}


def test_defaults_relocate_and_have_no_workstation_paths(portable):
    cfg = config_module.get_config()
    assert cfg.paths.repo_root == portable
    assert not cfg.datasets
    assert str(portable) in cfg.trusted_roots
    paths = [ROOT / "mcp_server/config.yaml", ROOT / "mcp_server/schemas/config.py",
             *list((ROOT / "mcp_server/workflows").glob("*.yaml")),
             *list((ROOT / "mcp_server/launchpad/static").glob("*"))]
    for path in paths:
        if not path.is_file():
            continue
        text = path.read_text()
        assert all(token not in text for token in ("/mnt/pikachu", "/storage/", "/home/lhhung", "/Users/")), path


def test_data_roots_are_explicit_site_configuration(portable, tmp_path, monkeypatch):
    data = tmp_path / "site data"
    monkeypatch.setenv("STAR_SUITE_DATA_ROOTS", str(data))
    cfg = load_config(portable / "mcp_server/config.yaml")
    assert str(data) in cfg.trusted_roots


def test_labels_required_inputs_and_namespaced_route(portable, inputs):
    with TestClient(build_http_app(), base_url="http://127.0.0.1") as client:
        for workflow in ("morphic_multiome", "starsuite.official%2Fmultiome"):
            response = client.get(f"/launchpad/api/workflows/{workflow}/schema")
            assert response.status_code == 200, response.text
            fields = {p["name"]: p for p in response.json()["parameters"]}
            assert "R2" in fields["atac_barcode"]["label"]
            assert "R3" in fields["atac_r2"]["label"]
            assert not fields["dry_run"]["default"]
            assert all(fields[name]["required"] and fields[name]["default"] is None for name in REFERENCES)
        assert client.get("/launchpad/").status_code == 200
        assert client.post(API + "/render", json={"params": {}}).status_code == 400
    rendered = render_workflow_command("morphic_multiome", inputs)
    assert str(portable / "share/star-suite/catalogs/official/scripts") in rendered.entry_script
    assert "--dry-run" not in rendered.argv
    assert "--stop-after-local-mex" in rendered.argv


@pytest.mark.parametrize("change", [
    {"gex_r1": ""}, {"gex_r1": "/missing"}, {"gex_r1": "/one,/two"},
    {"gex_r1": "/tmp/$(bad)"}, {"genome_dir": None}, {"atac_whitelist": ["wrong type"]},
    {"chromap_threads": 0}, {"skip_build": False}, {"stop_after_local_mex": False},
])
def test_reject_bad_inputs_before_execution(portable, inputs, change):
    result = validate_workflow_parameters("morphic_multiome", {**inputs, **change}, check_paths=True)
    assert not result.valid


def test_lane_normalization_and_real_recipe_dry_run(portable, inputs, fake_runtime):
    # No engine is invoked. The bundled recipe produces its actual command file.
    inputs.update({key: fake_runtime[key] for key in ("star_bin", "atac_peak_mex_bin")})
    extra = Path(inputs["gex_r1"]).parent / "lane two"
    extra.touch()
    inputs["gex_r1"] += "\n" + str(extra)
    inputs["gex_r2"] += ", " + str(extra)
    result = validate_workflow_parameters("morphic_multiome", inputs, check_paths=True)
    assert result.valid, result.errors
    assert "\n" not in result.normalized_params["gex_r1"]
    params = {**result.normalized_params, "dry_run": True}
    rendered = render_workflow_command("morphic_multiome", params)
    env = {k: v for k, v in os.environ.items() if not k.startswith("STAR_MULTIOME_")}
    env.update(rendered.env_overrides)
    assert "star_bin" not in rendered.argv and "atac_peak_mex_bin" not in rendered.argv
    run = subprocess.run(rendered.argv, env=env, capture_output=True, text=True, timeout=20)
    assert run.returncode == 0, run.stderr
    command = (Path(inputs["out_dir"]) / "RUN_STAR_MULTIOME.sh").read_text()
    assert f'--readFilesIn "{params["gex_r2"]}" "{params["gex_r1"]}"' in command
    assert '--chromapAtacReadFormat "bc:8:23:-"' in command
    assert '--clipAdapterType CellRanger4' in command
    assert not (Path(inputs["out_dir"]) / "LOCAL_MEX_READY.txt").exists()


def test_runtime_distinguishes_portable_and_chromap_builds(portable, fake_runtime):
    identity = runtime_check(fake_runtime, portable)
    assert len(identity["binary_sha256"]) == 64
    assert len(identity["atac_helper_sha256"]) == 64
    binary = Path(fake_runtime["star_bin"])
    binary.write_text(binary.read_text().replace('true', 'false'))
    with pytest.raises(ValueError, match="without Chromap"):
        runtime_check(fake_runtime, portable)
    with TestClient(build_http_app(), base_url="http://127.0.0.1") as client:
        response = client.post(API + "/validate", json={"params": fake_runtime, "check_paths": True})
        assert not response.json()["valid"]
        assert client.post(API + "/launch", json={"params": fake_runtime}).status_code == 400
    assert not Path(fake_runtime["out_dir"]).exists()


def wait_job(client, job_id):
    deadline = time.monotonic() + 8
    while time.monotonic() < deadline:
        job = client.get(f"/launchpad/api/jobs/{job_id}").json()
        if job["status"] not in {"running", "cancelling"}:
            return job
        time.sleep(.025)
    pytest.fail("fixture did not finish")


def test_job_result_logs_record_and_fresh_output(portable, fake_runtime):
    entry = portable / get_workflow_schema("morphic_multiome").entry_script
    entry.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
out = Path(sys.argv[sys.argv.index('--out-dir') + 1])
(out / 'LOCAL_MEX_READY.txt').write_text('fixture only')
print('fixture completed', flush=True)
''')
    with TestClient(build_http_app(), base_url="http://127.0.0.1") as client:
        response = client.post(API + "/launch", json={"params": fake_runtime})
        assert response.status_code == 202, response.text
        job = wait_job(client, response.json()["id"])
        assert job["status"] == "succeeded" and "fixture completed" in job["log_tail"]
        record = json.loads(Path(job["record_path"]).read_text())
        assert record["status"] == "succeeded" and record["runtime"]["features"]["chromap_atac"]
        assert client.post(API + "/launch", json={"params": fake_runtime}).status_code == 409


def test_serial_jobs_cancel_and_access_guards(portable, fake_runtime):
    entry = portable / get_workflow_schema("morphic_multiome").entry_script
    entry.write_text("#!/usr/bin/env python3\nimport time\ntime.sleep(60)\n")
    with TestClient(build_http_app(), base_url="http://127.0.0.1") as client:
        for headers in ({"Host": "evil.example"}, {"Origin": "http://evil.example"}):
            assert client.post(API + "/launch", json={"params": fake_runtime}, headers=headers).status_code == 403
        response = client.post(API + "/launch", json={"params": fake_runtime})
        assert response.status_code == 202, response.text
        job_id = response.json()["id"]
        competing = {**fake_runtime, "out_dir": fake_runtime["out_dir"] + "2"}
        assert client.post(API + "/launch", json={"params": competing}).status_code == 409
        assert not Path(competing["out_dir"]).exists()
        assert client.post(f"/launchpad/api/jobs/{job_id}/cancel", json={}).status_code == 200
        assert wait_job(client, job_id)["status"] == "cancelled"


@pytest.fixture
def installed(portable, tmp_path):
    # Exercise the same payload staging code used by both tarballs and .deb.
    prefix = tmp_path / "installed suite"
    subprocess.run([sys.executable, str(ROOT / "scripts/release/stage_launchpad.py"),
                    "--stage-root", str(prefix)], check=True, capture_output=True)
    cfg_path = prefix / "share/star-suite/launchpad/mcp_server/config.yaml"
    raw = yaml.safe_load(cfg_path.read_text())
    raw["paths"]["artifact_log_root"] = str(tmp_path / "installed records")
    cfg_path.write_text(yaml.safe_dump(raw))
    cfg = load_config(cfg_path)
    return prefix, cfg


def test_installed_payload_resolves_without_checkout(installed):
    prefix, cfg = installed
    root = prefix / "share/star-suite/launchpad"
    assert cfg.paths.repo_root == root
    assert str(prefix / "bin") in cfg.trusted_roots
    assert not cfg.scripts and not cfg.datasets
    for workflow in ("morphic_multiome", "starsuite.official/multiome"):
        schema = get_workflow_schema(workflow)
        origin = config_module.get_workflow_root(workflow)
        assert (origin / schema.entry_script).is_file()
    star = get_workflow_schema("star_genome_generate")
    assert (root / star.entry_script).resolve() == prefix / "bin/STAR"
    result = subprocess.run([str(prefix / "bin/star-suite-launchpad"), "--help"],
                            capture_output=True, text=True, check=True)
    assert "--setup" in result.stdout
    assert (root / "mcp_server/launchpad/static/vendor/alpine-3.14.3.min.js").is_file()


def test_installed_browser_form_run_and_reload(installed, fake_runtime, tmp_path):
    playwright = pytest.importorskip("playwright.sync_api")
    prefix, _ = installed
    entry = prefix / "share/star-suite/catalogs/official/scripts/run_star_multiome_lane_smoke.sh"
    entry.write_text('''#!/usr/bin/env python3
import sys, time
from pathlib import Path
out = Path(sys.argv[sys.argv.index('--out-dir') + 1])
print('fixture recipe started', flush=True)
time.sleep(1)
(out / 'LOCAL_MEX_READY.txt').write_text('fixture only')
''')
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    env = {**os.environ, "STAR_SUITE_LAUNCHPAD_PYTHON": sys.executable}
    log = tmp_path / "installed-server.log"
    with log.open("w") as stream:
        process = subprocess.Popen([str(prefix / "bin/star-suite-launchpad"), "--port", str(port)],
                                   cwd=tmp_path, env=env, stdout=stream, stderr=stream)
    base = f"http://127.0.0.1:{port}"
    try:
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            assert process.poll() is None, log.read_text()
            try:
                with urlopen(base + "/launchpad/", timeout=.5):
                    break
            except OSError:
                time.sleep(.05)
        else:
            pytest.fail(log.read_text())
        with playwright.sync_playwright() as pw:
            browser = pw.chromium.launch(headless=True, args=["--no-sandbox"])
            page = browser.new_page(viewport={"width": 1280, "height": 960})
            errors = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            # No CDN or other external browser resource is needed.
            page.route("**/*", lambda route: route.continue_() if route.request.url.startswith(base + "/") else route.abort())
            page.goto(base + "/launchpad/")
            page.locator('#recipe-select option[value="morphic_multiome"]').wait_for(state="attached")
            page.locator("#recipe-select").select_option("morphic_multiome")
            field = lambda name: page.locator(f'[data-param="{name}"]')
            field("atac_r2").wait_for()
            assert "ATAC R3" in field("atac_r2").inner_text()
            assert not field("gex_cbq").is_visible()
            field("input_format").locator("select").select_option("cbq")
            playwright.expect(field("gex_cbq")).to_be_visible()
            assert not field("atac_r1").is_visible()
            field("input_format").locator("select").select_option("fastq")
            for name, value in fake_runtime.items():
                field(name).locator('textarea:visible, input[type="text"]:visible').fill(value)
            page.get_by_role("button", name="Generate command", exact=True).click()
            playwright.expect(page.locator("#cmd-block")).to_contain_text("--atac-r2")
            assert "--dry-run" not in page.locator("#cmd-block").inner_text()
            assert str(entry) in page.locator("#cmd-block").inner_text()
            page.get_by_role("button", name="Run sample", exact=True).click()
            playwright.expect(page.locator("#multiome-job-status")).to_contain_text("succeeded", timeout=15000)
            assert "fixture recipe started" in page.locator("#multiome-job-log").inner_text()
            page.reload()
            playwright.expect(page.locator("#multiome-job-status")).to_contain_text("succeeded", timeout=10000)
            playwright.expect(page.locator("#recipe-select")).to_have_value("morphic_multiome")
            playwright.expect(field("gex_r1").locator("textarea:visible")).to_have_value(fake_runtime["gex_r1"])
            evidence = ROOT / "tests/launchpad_portability_output"
            evidence.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(evidence / "installed-browser.png"), full_page=True)
            assert not errors
            browser.close()
    finally:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=5)
