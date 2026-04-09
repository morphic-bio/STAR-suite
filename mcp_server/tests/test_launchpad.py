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
from starlette.requests import Request

from .test_workflow_discovery import SAMPLE_WORKFLOW_SCHEMA

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

    def test_validate_bad_body(self, launchpad_client):
        r = launchpad_client.post(
            "/launchpad/api/workflows/test_wf/validate",
            json={},
        )
        assert r.status_code == 400

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
        assert r.json() == {"launch_supported": True}

    def test_capabilities_not_trusted_when_mocked(self, launchpad_client):
        with mock.patch.object(
            launchpad_api, "launchpad_request_trusted_local", return_value=False
        ):
            r = launchpad_client.get("/launchpad/api/capabilities")
            assert r.status_code == 200
            assert r.json() == {"launch_supported": False}

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
