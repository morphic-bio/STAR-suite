"""HTTP tests for STAR Launchpad API (TestClient)."""

import tempfile
from pathlib import Path
from unittest import mock

import pytest
import yaml

import mcp_server.config as config_module
import mcp_server.launchpad.api as launchpad_api
from mcp_server.app import build_http_app
from mcp_server.config import load_config
from mcp_server.tools.workflows import get_workflow_parameter_schema
from starlette.requests import Request

from .test_workflow_discovery import SAMPLE_WORKFLOW_SCHEMA

BWB_NEXTFLOW_UTILS_ROOT = Path(__file__).resolve().parents[2].parent / "bwb-nextflow-utils"

LAUNCH_TEST_SCHEMA = {
    "id": "launch_test_wf",
    "title": "Launch test",
    "summary": "No-op binary for launch endpoint tests.",
    "kind": "shell_workflow",
    "entry_script": "/bin/true",
    "parameters": [],
    "constraints": [],
    "rendering": {"flag_order": [], "omit_absent_optionals": True},
}

PRIVATE_MINI_SCHEMA = {
    "id": "private_wf",
    "title": "Private Workflow",
    "summary": "Internal only unless Launchpad trusts loopback.",
    "kind": "shell_workflow",
    "entry_script": "scripts/private.sh",
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


def _request_with_client(client):
    return Request(
        {
            "type": "http",
            "asgi": {"version": "3.0", "spec_version": "2.3"},
            "http_version": "1.1",
            "method": "GET",
            "path": "/launchpad/api/workflows",
            "raw_path": b"/launchpad/api/workflows",
            "root_path": "",
            "scheme": "http",
            "query_string": b"",
            "headers": [],
            "client": client,
            "server": ("127.0.0.1", 8765),
        }
    )


@pytest.fixture(autouse=True)
def _reset():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}


@pytest.fixture
def launchpad_client():
    """Config with ``test_wf`` and Starlette TestClient for Launchpad routes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        schema_rel = "workflows/test_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(SAMPLE_WORKFLOW_SCHEMA, f)

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
                    "id": "test_wf",
                    "title": "Test Workflow",
                    "summary": "A test.",
                    "entry_script": "scripts/test.sh",
                    "kind": "shell_workflow",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)

        (tmp / "scripts").mkdir()
        (tmp / "scripts" / "test.sh").write_text("#!/bin/bash\necho ok\n")

        app = build_http_app()
        from starlette.testclient import TestClient

        with TestClient(app) as client:
            yield client


@pytest.fixture
def launchpad_client_launch_noop():
    """Minimal workflow whose entry is ``/bin/true`` for safe ``launch`` tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        schema_rel = "workflows/launch_test_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / schema_rel, "w") as f:
            yaml.dump(LAUNCH_TEST_SCHEMA, f)

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
                    "id": "launch_test_wf",
                    "title": "Launch test",
                    "summary": "No-op.",
                    "entry_script": "/bin/true",
                    "kind": "shell_workflow",
                    "schema_file": schema_rel,
                }
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)

        app = build_http_app()
        from starlette.testclient import TestClient

        with TestClient(app) as client:
            yield client


@pytest.fixture
def launchpad_client_public_private():
    """One public and one private workflow (for trust / visibility tests)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        pub_rel = "workflows/test_wf.yaml"
        priv_rel = "workflows/private_wf.yaml"
        (tmp / "workflows").mkdir()
        with open(tmp / pub_rel, "w") as f:
            yaml.dump(SAMPLE_WORKFLOW_SCHEMA, f)
        with open(tmp / priv_rel, "w") as f:
            yaml.dump(PRIVATE_MINI_SCHEMA, f)

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
                    "id": "test_wf",
                    "title": "Test Workflow",
                    "summary": "A test.",
                    "entry_script": "scripts/test.sh",
                    "kind": "shell_workflow",
                    "schema_file": pub_rel,
                },
                {
                    "id": "private_wf",
                    "title": "Private Workflow",
                    "summary": "Hidden from remote Launchpad.",
                    "entry_script": "scripts/private.sh",
                    "kind": "shell_workflow",
                    "schema_file": priv_rel,
                    "visibility": "private",
                },
            ],
        }
        cfg_path = tmp / "config.yaml"
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        load_config(cfg_path)

        (tmp / "scripts").mkdir()
        (tmp / "scripts" / "test.sh").write_text("#!/bin/bash\necho ok\n")
        (tmp / "scripts" / "private.sh").write_text("#!/bin/bash\necho ok\n")

        app = build_http_app()
        from starlette.testclient import TestClient

        with TestClient(app) as client:
            yield client


class TestLaunchpadTrustedLocal:
    def test_loopback_hosts_trusted(self):
        assert launchpad_api.launchpad_request_trusted_local(
            _request_with_client(("127.0.0.1", 12345))
        )
        assert launchpad_api.launchpad_request_trusted_local(
            _request_with_client(("::1", 12345))
        )
        assert launchpad_api.launchpad_request_trusted_local(
            _request_with_client(("::ffff:127.0.0.1", 12345))
        )

    def test_testclient_trusted(self):
        assert launchpad_api.launchpad_request_trusted_local(
            _request_with_client(("testclient", 50000))
        )

    def test_lan_not_trusted(self):
        assert not launchpad_api.launchpad_request_trusted_local(
            _request_with_client(("192.0.2.10", 12345))
        )

    def test_missing_client_not_trusted(self):
        scope = {
            "type": "http",
            "asgi": {"version": "3.0", "spec_version": "2.3"},
            "http_version": "1.1",
            "method": "GET",
            "path": "/",
            "raw_path": b"/",
            "root_path": "",
            "scheme": "http",
            "query_string": b"",
            "headers": [],
            "server": ("127.0.0.1", 8765),
        }
        assert not launchpad_api.launchpad_request_trusted_local(Request(scope))


class TestLaunchpadApi:
    def test_list_workflows(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/workflows")
        assert r.status_code == 200
        data = r.json()
        assert "workflows" in data
        assert len(data["workflows"]) == 1
        assert data["workflows"][0]["id"] == "test_wf"

    def test_schema_ok(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/workflows/test_wf/schema")
        assert r.status_code == 200
        data = r.json()
        assert data["workflow_id"] == "test_wf"
        assert "required_files" in data
        names = [x["name"] for x in data["required_files"]]
        assert "input_dir" in names
        assert data["required_files"][0]["type"] == "directory"

    def test_schema_matches_mcp_contract(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/workflows/test_wf/schema")
        assert r.status_code == 200

        expected = get_workflow_parameter_schema("test_wf").model_dump()
        assert r.json() == expected

    def test_schema_unknown(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/workflows/missing_wf/schema")
        assert r.status_code == 404
        assert r.json()["code"] == "NOT_FOUND"

    def test_describe_ok(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/workflows/test_wf/describe")
        assert r.status_code == 200
        data = r.json()
        assert data["id"] == "test_wf"
        assert data["caveats"] == ["Test caveat"]

    def test_validate_no_path_checks(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/validate",
            json={"params": {"input_dir": "/nonexistent/path/that/does/not/exist"}},
        )
        assert r.status_code == 200
        data = r.json()
        assert data["valid"] is True
        assert data["errors"] == []

    def test_validate_checks_paths_when_requested(self, launchpad_client):
        missing = "/tmp/launchpad-test-missing-input-dir"
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/validate",
            json={
                "params": {"input_dir": missing},
                "check_paths": True,
            },
        )
        assert r.status_code == 200
        data = r.json()
        assert data["valid"] is False
        assert any("does not exist" in msg for msg in data["errors"])
        assert any(fe["field"] == "input_dir" for fe in data["field_errors"])

    def test_validate_rejects_path_checks_off_loopback(self, launchpad_client):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.post(
                "/launchpad/api/workflows/test_wf/validate",
                json={"params": {"input_dir": "/tmp"}, "check_paths": True},
            )
            assert r.status_code == 403
            assert r.json()["code"] == "FORBIDDEN"

    def test_validate_bad_body(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/validate",
            json={},
        )
        assert r.status_code == 400

    def test_validate_bad_check_paths_type(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/validate",
            json={"params": {"input_dir": "/tmp"}, "check_paths": "yes"},
        )
        assert r.status_code == 400
        assert r.json()["code"] == "BAD_REQUEST"

    def test_render(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/render",
            json={"params": {"input_dir": "/tmp"}},
        )
        assert r.status_code == 200
        data = r.json()
        assert data["workflow_id"] == "test_wf"
        assert "argv" in data
        # Launchpad should not leak host repo_root via absolute entry_script paths.
        assert data["entry_script"] == "scripts/test.sh"
        assert data["argv"][0] == "scripts/test.sh"
        assert data["shell_preview"]
        assert "--input-dir" in data["shell_preview"]

    def test_static_index(self, launchpad_client):
        r = launchpad_client.get("/launchpad/")
        assert r.status_code == 200
        assert b"STAR Launchpad" in r.content

    def test_list_all_workflows_when_trusted_local_mock(
        self, launchpad_client_public_private
    ):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=True
        ):
            r = launchpad_client_public_private.get("/launchpad/api/workflows")
            assert r.status_code == 200
            ids = {w["id"] for w in r.json()["workflows"]}
            assert ids == {"test_wf", "private_wf"}

    def test_list_public_only_when_not_trusted_local_mock(
        self, launchpad_client_public_private
    ):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client_public_private.get("/launchpad/api/workflows")
            assert r.status_code == 200
            ids = {w["id"] for w in r.json()["workflows"]}
            assert ids == {"test_wf"}

    def test_capabilities_loopback(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/capabilities")
        assert r.status_code == 200
        data = r.json()
        assert data["launch_supported"] is True
        assert data["server_path_check_supported"] is True
        assert "script_lane_translate_supported" in data
        assert "script_lane_utils_ready" in data
        assert "script_lane_viability_inputs_supported" in data
        assert "script_lane_execute_supported" in data

    def test_capabilities_not_trusted_when_mocked(self, launchpad_client):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.get("/launchpad/api/capabilities")
            assert r.status_code == 200
            data = r.json()
            assert data["launch_supported"] is False
            assert data["server_path_check_supported"] is False
            assert data["script_lane_viability_inputs_supported"] is False
            assert data["script_lane_execute_supported"] is False

    def test_launch_starts_process(self, launchpad_client_launch_noop):
        r = launchpad_client_launch_noop.post(
            "/launchpad/api/workflows/launch_test_wf/launch",
            json={"params": {}},
        )
        assert r.status_code == 200
        data = r.json()
        assert data.get("ok") is True
        assert isinstance(data.get("pid"), int)
        assert data.get("message")

    def test_launch_forbidden_off_loopback(self, launchpad_client_launch_noop):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client_launch_noop.post(
                "/launchpad/api/workflows/launch_test_wf/launch",
                json={"params": {}},
            )
            assert r.status_code == 403
            assert r.json().get("code") == "FORBIDDEN"

    def test_launch_checks_paths_before_start(self, launchpad_client):
        missing = "/tmp/launchpad-test-missing-input-dir"
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/launch",
            json={"params": {"input_dir": missing}},
        )
        assert r.status_code == 400
        data = r.json()
        assert data.get("code") == "VALIDATION_FAILED"
        assert data["validation"]["valid"] is False
        assert any("input_dir" == fe["field"] for fe in data["validation"]["field_errors"])

    def test_quit_requires_confirm(self, launchpad_client):
        r = launchpad_client.post("/launchpad/api/quit", json={})
        assert r.status_code == 400
        assert r.json().get("code") == "BAD_REQUEST"

    def test_quit_ok_trusted_local(self, launchpad_client):
        with mock.patch.object(launchpad_api, "_schedule_self_terminate") as sched:
            r = launchpad_client.post("/launchpad/api/quit", json={"confirm": True})
            assert r.status_code == 200
            data = r.json()
            assert data.get("ok") is True
            assert isinstance(data.get("pid"), int)
            assert data.get("message")
            assert sched.called

    def test_quit_forbidden_off_loopback(self, launchpad_client):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.post("/launchpad/api/quit", json={"confirm": True})
            assert r.status_code == 403
            assert r.json().get("code") == "FORBIDDEN"


def _minimal_script_lane_bash() -> str:
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "#@workflow id=launchpad_sl_test",
            "#@input genome file required",
            "#@output bam file",
            "#@step id=copy tool=cp",
            "#@uses genome",
            "#@produces bam",
            'cp "${genome}" "${bam}"',
        ]
    )


@pytest.fixture
def script_lane_bwb(monkeypatch):
    if not BWB_NEXTFLOW_UTILS_ROOT.is_dir():
        pytest.skip("bwb-nextflow-utils not found (set BWB_NEXTFLOW_UTILS_ROOT or sibling checkout)")
    monkeypatch.setenv("BWB_NEXTFLOW_UTILS_ROOT", str(BWB_NEXTFLOW_UTILS_ROOT))
    import mcp_server.launchpad.script_lane_bridge as slb

    slb._svc = None
    yield
    slb._svc = None


class TestLaunchpadBrowseApi:
    def test_capabilities_exposes_browse_roots(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/capabilities")
        assert r.status_code == 200
        data = r.json()
        assert data["trusted_local"] is True
        assert data["browse_supported"] is True
        assert data["upload_supported"] is True
        assert isinstance(data["browse_roots"], list)
        assert data["browse_roots"], "expected at least one trusted root"

    def test_browse_without_path_lists_roots(self, launchpad_client):
        r = launchpad_client.get("/launchpad/api/browse")
        assert r.status_code == 200
        data = r.json()
        assert data["is_root_list"] is True
        assert data["parent"] is None
        names = {e["name"] for e in data["entries"]}
        assert names, "expected trusted roots to be listed"
        for e in data["entries"]:
            assert e["kind"] == "dir"

    def test_browse_directory_lists_entries(self, launchpad_client, tmp_path):
        # tmp_path is under pytest's /tmp, which is in trusted_roots.
        sub = tmp_path / "browse_fixture"
        sub.mkdir()
        (sub / "a.txt").write_text("hello")
        (sub / "nested").mkdir()
        r = launchpad_client.get(
            "/launchpad/api/browse", params={"path": str(sub)}
        )
        assert r.status_code == 200, r.text
        data = r.json()
        assert data["is_root_list"] is False
        assert data["path"] == str(sub)
        kinds = {e["name"]: e["kind"] for e in data["entries"]}
        assert kinds["a.txt"] == "file"
        assert kinds["nested"] == "dir"
        # Directories must sort before files.
        order = [e["name"] for e in data["entries"]]
        assert order.index("nested") < order.index("a.txt")

    def test_browse_rejects_path_outside_trusted_roots(self, launchpad_client):
        r = launchpad_client.get(
            "/launchpad/api/browse", params={"path": "/etc"}
        )
        assert r.status_code == 403
        assert r.json()["code"] == "FORBIDDEN"

    def test_browse_rejects_traversal(self, launchpad_client):
        # "/tmp/../etc" resolves to "/etc" (outside trusted_roots for this fixture).
        r = launchpad_client.get(
            "/launchpad/api/browse", params={"path": "/tmp/../etc"}
        )
        assert r.status_code == 403

    def test_browse_rejects_symlink_escape(self, launchpad_client, tmp_path):
        # A symlink that points outside trusted_roots must be rejected even though
        # its nominal path lives under a trusted root.
        target = Path("/etc")
        if not target.is_dir():
            pytest.skip("no /etc to symlink")
        link = tmp_path / "escape"
        link.symlink_to(target)
        r = launchpad_client.get(
            "/launchpad/api/browse", params={"path": str(link)}
        )
        assert r.status_code == 403

    def test_browse_file_path_rejected(self, launchpad_client, tmp_path):
        f = tmp_path / "solo.txt"
        f.write_text("x")
        r = launchpad_client.get(
            "/launchpad/api/browse", params={"path": str(f)}
        )
        assert r.status_code == 400
        assert r.json()["code"] == "BAD_REQUEST"

    def test_browse_missing_path(self, launchpad_client, tmp_path):
        r = launchpad_client.get(
            "/launchpad/api/browse",
            params={"path": str(tmp_path / "does_not_exist")},
        )
        assert r.status_code == 404

    def test_upload_writes_file_to_temp_root(self, launchpad_client):
        content = b"uploaded bytes\n"
        r = launchpad_client.post(
            "/launchpad/api/upload",
            files={"file": ("hello.bin", content, "application/octet-stream")},
        )
        assert r.status_code == 200, r.text
        data = r.json()
        assert data["ok"] is True
        assert data["filename"] == "hello.bin"
        assert data["size"] == len(content)
        written = Path(data["path"])
        assert written.exists()
        assert written.read_bytes() == content
        # Must land under configured temp_root/launchpad_uploads/<uuid>/
        cfg = config_module.get_config()
        temp_root = Path(str(cfg.paths.temp_root)).resolve()
        assert str(written).startswith(str(temp_root))
        assert "launchpad_uploads" in written.parts

    def test_upload_sanitizes_filename(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/upload",
            files={"file": ("../weird name!.tar.gz", b"x", "application/octet-stream")},
        )
        assert r.status_code == 200, r.text
        data = r.json()
        # No path separators; hostile chars replaced with underscores.
        assert "/" not in data["filename"]
        assert "\\" not in data["filename"]
        assert data["filename"].startswith("weird") or data["filename"].startswith("_weird")

    def test_upload_rejects_non_multipart(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/upload",
            json={"file": "not a real upload"},
        )
        assert r.status_code == 400
        assert r.json()["code"] == "BAD_REQUEST"

    def test_upload_requires_file_field(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/upload",
            files={"other": ("nope.txt", b"x", "text/plain")},
        )
        assert r.status_code == 400


@pytest.mark.usefixtures("script_lane_bwb")
class TestLaunchpadScriptLaneApi:
    def test_translate_success(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        assert r.status_code == 200
        data = r.json()
        assert data.get("ok") is True
        assert data.get("ir", {}).get("metadata", {}).get("id") == "launchpad_sl_test"

    def test_translate_refusal(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": "not a bash workflow script"},
        )
        assert r.status_code == 200
        data = r.json()
        assert data.get("ok") is False
        assert data.get("fallback_to_temporal") is True
        assert data.get("reasons")

    def test_viability_loopback_with_inputs(self, launchpad_client, tmp_path):
        tr = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        assert tr.status_code == 200
        ir = tr.json()["ir"]
        src = tmp_path / "g.bin"
        src.write_text("x")
        r = launchpad_client.post(
            "/launchpad/api/script-lane/viability",
            json={"ir": ir, "input_values": {"genome": str(src)}},
        )
        assert r.status_code == 200
        assert r.json().get("ok") is True

    def test_viability_empty_input_values_dict_matches_omitted(self, launchpad_client):
        """``input_values: {}`` must behave like omitted (no required-input probe)."""
        tr = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        ir = tr.json()["ir"]
        r_omit = launchpad_client.post(
            "/launchpad/api/script-lane/viability",
            json={"ir": ir},
        )
        r_empty = launchpad_client.post(
            "/launchpad/api/script-lane/viability",
            json={"ir": ir, "input_values": {}},
        )
        assert r_omit.status_code == 200
        assert r_empty.status_code == 200
        assert r_omit.json() == r_empty.json()
        assert r_omit.json().get("ok") is True

    def test_viability_remote_forbidden_with_nonempty_inputs(self, launchpad_client, tmp_path):
        tr = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        ir = tr.json()["ir"]
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.post(
                "/launchpad/api/script-lane/viability",
                json={"ir": ir, "input_values": {"genome": str(tmp_path / "nope")}},
            )
            assert r.status_code == 403
            assert r.json().get("code") == "FORBIDDEN"

    def test_execute_loopback(self, launchpad_client, tmp_path):
        tr = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        ir = tr.json()["ir"]
        src = tmp_path / "in.txt"
        src.write_text("hello")
        wd = tmp_path / "sl_run"
        wd.mkdir()
        r = launchpad_client.post(
            "/launchpad/api/script-lane/execute",
            json={"ir": ir, "input_values": {"genome": str(src)}, "workdir": str(wd)},
        )
        assert r.status_code == 200
        data = r.json()
        assert data.get("ok") is True
        assert data.get("steps")

    def test_execute_remote_forbidden(self, launchpad_client, tmp_path):
        tr = launchpad_client.post(
            "/launchpad/api/script-lane/translate",
            json={"script_text": _minimal_script_lane_bash()},
        )
        ir = tr.json()["ir"]
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.post(
                "/launchpad/api/script-lane/execute",
                json={
                    "ir": ir,
                    "input_values": {"genome": str(tmp_path / "x")},
                    "workdir": str(tmp_path / "w"),
                },
            )
            assert r.status_code == 403
            assert r.json().get("code") == "FORBIDDEN"
