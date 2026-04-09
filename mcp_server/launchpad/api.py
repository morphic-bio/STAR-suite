"""HTTP JSON API for STAR Launchpad (thin wrappers over shared workflow core).

Loopback clients see all workflows (including private), like authenticated MCP
discovery. Non-loopback clients see public workflows only.

Rendered commands still use schema-relative entry paths (no host path leak).
Validation uses check_paths=False for path planning on the user's machine.

``POST .../launch`` starts the rendered argv on the server host with
``subprocess.Popen`` (detached, no shell). It is allowed only for loopback
clients so arbitrary remote browsers cannot trigger execution.
"""

import os
import signal
import shlex
import subprocess
import threading
import time
from pathlib import Path

from starlette.requests import Request
from starlette.responses import JSONResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles

from ..tools.workflows import (
    describe_workflow,
    get_workflow_parameter_schema,
    list_workflows,
    render_workflow_command,
    validate_workflow_parameters,
)

_LAUNCHPAD_STATIC = Path(__file__).resolve().parent / "static"


def launchpad_request_trusted_local(request: Request) -> bool:
    """True when the HTTP peer is loopback (same machine as the server).

    Used only for Launchpad JSON API visibility (public vs all workflows).
    Do not trust ``X-Forwarded-For`` here: it can be spoofed from off-host clients.
    """
    client = request.client
    if client is None:
        return False
    host = (client.host or "").strip().lower()
    if not host:
        return False
    if host in ("127.0.0.1", "::1", "localhost"):
        return True
    # Starlette/httpx TestClient uses a synthetic peer (not a real remote IP).
    if host == "testclient":
        return True
    # IPv4-mapped IPv6 loopback
    if host.startswith("::ffff:"):
        tail = host.split("::ffff:", 1)[-1]
        if tail == "127.0.0.1":
            return True
    return False


def _json_error(code: str, message: str, status_code: int = 400) -> JSONResponse:
    return JSONResponse(
        {"error": True, "code": code, "message": message},
        status_code=status_code,
    )


async def lp_capabilities(request: Request) -> JSONResponse:
    """Whether this browser session may use server-side Run in shell (loopback only)."""
    return JSONResponse(
        {"launch_supported": launchpad_request_trusted_local(request)}
    )


def _schedule_self_terminate(
    pid: int,
    *,
    delay_seconds: float = 0.25,
    sig: int = signal.SIGTERM,
) -> None:
    """Terminate the current server process after returning a response.

    This lets the HTTP handler acknowledge the request before stopping the
    underlying Uvicorn process.
    """

    def _worker() -> None:
        time.sleep(delay_seconds)
        try:
            os.kill(pid, sig)
        except Exception:
            # Best-effort shutdown: if the process already exited, ignore.
            pass

    threading.Thread(target=_worker, daemon=True).start()


async def lp_list_workflows(request: Request) -> JSONResponse:
    try:
        authed = launchpad_request_trusted_local(request)
        result = list_workflows(authenticated=authed)
        return JSONResponse(result.model_dump())
    except Exception as e:
        return _json_error("INTERNAL_ERROR", str(e), 500)


async def lp_get_schema(request: Request) -> JSONResponse:
    wf_id = request.path_params["workflow_id"]
    try:
        authed = launchpad_request_trusted_local(request)
        result = get_workflow_parameter_schema(wf_id, authenticated=authed)
        return JSONResponse(result.model_dump())
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)


async def lp_describe(request: Request) -> JSONResponse:
    wf_id = request.path_params["workflow_id"]
    try:
        authed = launchpad_request_trusted_local(request)
        result = describe_workflow(wf_id, authenticated=authed)
        return JSONResponse(result.model_dump())
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)


async def lp_validate(request: Request) -> JSONResponse:
    wf_id = request.path_params["workflow_id"]
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    params = body.get("params")
    if params is None or not isinstance(params, dict):
        return _json_error("BAD_REQUEST", "Body must contain an object 'params'", 400)
    try:
        result = validate_workflow_parameters(
            wf_id, params, check_paths=False
        )
        return JSONResponse(result.model_dump())
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)


async def lp_render(request: Request) -> JSONResponse:
    wf_id = request.path_params["workflow_id"]
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    params = body.get("params")
    if params is None or not isinstance(params, dict):
        return _json_error("BAD_REQUEST", "Body must contain an object 'params'", 400)
    try:
        result = render_workflow_command(wf_id, params)
        # Launchpad is intentionally unauthenticated and should not leak host
        # filesystem layout (e.g. repo_root). Rewrite the entry script from the
        # schema path and rebuild argv/shell_preview accordingly.
        authed = launchpad_request_trusted_local(request)
        public_entry = describe_workflow(wf_id, authenticated=authed).entry_script
        payload = result.model_dump()
        argv = list(payload.get("argv") or [])
        if argv:
            argv[0] = public_entry
        payload["entry_script"] = public_entry
        payload["argv"] = argv
        payload["shell_preview"] = " ".join(shlex.quote(a) for a in argv)
        return JSONResponse(payload)
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)


async def lp_launch(request: Request) -> JSONResponse:
    """Validate, render, and start argv on the server (loopback peers only)."""
    if not launchpad_request_trusted_local(request):
        return _json_error(
            "FORBIDDEN",
            "Run in shell is only available when Launchpad is opened from the same "
            "machine as the server (e.g. http://127.0.0.1:…/launchpad/).",
            403,
        )
    wf_id = request.path_params["workflow_id"]
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    params = body.get("params")
    if params is None or not isinstance(params, dict):
        return _json_error("BAD_REQUEST", "Body must contain an object 'params'", 400)

    try:
        val = validate_workflow_parameters(wf_id, params, check_paths=False)
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)

    if not val.valid:
        return JSONResponse(
            {
                "error": True,
                "code": "VALIDATION_FAILED",
                "validation": val.model_dump(),
            },
            status_code=400,
        )

    try:
        result = render_workflow_command(wf_id, params)
    except ValueError as e:
        return _json_error("NOT_FOUND", str(e), 404)

    argv = list(result.argv or [])
    if not argv:
        return _json_error("INTERNAL_ERROR", "Empty argv after render", 500)

    env = os.environ.copy()
    for k, v in (result.env_overrides or {}).items():
        env[str(k)] = str(v)

    popen_kw: dict = {
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "env": env,
    }
    if os.name != "nt":
        popen_kw["start_new_session"] = True

    try:
        proc = subprocess.Popen(argv, **popen_kw)
    except OSError as e:
        return _json_error("EXEC_FAILED", str(e), 500)

    return JSONResponse(
        {
            "ok": True,
            "pid": proc.pid,
            "message": (
                f"Started process {proc.pid} on the server host (detached). "
                "STAR writes logs under your chosen output prefix; stdout/stderr are not captured here."
            ),
        }
    )


async def lp_quit(request: Request) -> JSONResponse:
    """Gracefully stop the STAR Server (Launchpad + MCP) on this host.

    Loopback-only: prevents remote browsers from shutting down the service.
    """
    if not launchpad_request_trusted_local(request):
        return _json_error(
            "FORBIDDEN",
            "Quit server is only available when Launchpad is opened from the same "
            "machine as the server (e.g. http://127.0.0.1:…/launchpad/).",
            403,
        )
    try:
        body = await request.json()
    except Exception:
        body = {}
    if not isinstance(body, dict) or body.get("confirm") is not True:
        return _json_error(
            "BAD_REQUEST",
            "Set JSON body to {\"confirm\": true} to stop the server.",
            400,
        )

    pid = os.getpid()
    _schedule_self_terminate(pid)
    return JSONResponse(
        {
            "ok": True,
            "pid": pid,
            "message": "Server is shutting down (SIGTERM).",
        }
    )


def get_launchpad_routes() -> list:
    """Routes for Launchpad API (must be registered before static /launchpad mount)."""
    return [
        Route(
            "/launchpad/api/capabilities",
            endpoint=lp_capabilities,
            methods=["GET"],
        ),
        Route(
            "/launchpad/api/quit",
            endpoint=lp_quit,
            methods=["POST"],
        ),
        Route("/launchpad/api/workflows", endpoint=lp_list_workflows, methods=["GET"]),
        Route(
            "/launchpad/api/workflows/{workflow_id}/schema",
            endpoint=lp_get_schema,
            methods=["GET"],
        ),
        Route(
            "/launchpad/api/workflows/{workflow_id}/describe",
            endpoint=lp_describe,
            methods=["GET"],
        ),
        Route(
            "/launchpad/api/workflows/{workflow_id}/validate",
            endpoint=lp_validate,
            methods=["POST"],
        ),
        Route(
            "/launchpad/api/workflows/{workflow_id}/render",
            endpoint=lp_render,
            methods=["POST"],
        ),
        Route(
            "/launchpad/api/workflows/{workflow_id}/launch",
            endpoint=lp_launch,
            methods=["POST"],
        ),
        Mount(
            "/launchpad",
            app=StaticFiles(directory=str(_LAUNCHPAD_STATIC), html=True),
        ),
    ]
