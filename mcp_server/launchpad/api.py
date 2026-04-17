"""HTTP JSON API for STAR Launchpad (thin wrappers over shared workflow core).

Loopback clients see all workflows (including private), like authenticated MCP
discovery. Non-loopback clients see public workflows only.

Rendered commands still use schema-relative entry paths (no host path leak).
Validation defaults to ``check_paths=False`` for browser-side planning, but
trusted-local callers may opt into server-host path checks.

``POST .../launch`` starts the rendered argv on the server host with
``subprocess.Popen`` (detached, no shell). It is allowed only for loopback
clients so arbitrary remote browsers cannot trigger execution.

Script Lane endpoints (``/launchpad/api/script-lane/*``) delegate annotated Bash,
simple-local viability, and simple-local execution to ``bwb-nextflow-utils``
via ``script_lane_bridge``; they do not reuse ``/launch`` or workflow render.
"""

import os
import re
import signal
import shlex
import subprocess
import threading
import time
import uuid
from pathlib import Path

from starlette.requests import Request
from starlette.responses import FileResponse, JSONResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles

from ..config import get_config
from ..tools.workflows import (
    describe_workflow,
    get_workflow_parameter_schema,
    list_workflows,
    render_workflow_command,
    validate_workflow_parameters,
)

from .script_lane_bridge import (
    ScriptLaneBackendUnavailable,
    backend_status,
    check_simple_local_ir_viability,
    execute_simple_local_ir,
    translate_simple_script_to_ir,
)

_LAUNCHPAD_STATIC = Path(__file__).resolve().parent / "static"
_BRIDGE_DIST = Path(__file__).resolve().parents[2] / "extension" / "dist"


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


# --- Filesystem browse / upload helpers --------------------------------------
# Browsing is restricted to the config's ``trusted_roots``. Paths are resolved
# (symlinks followed) and re-checked so clients can't escape via ``..`` or
# symlinked directories.

_UPLOAD_SUBDIR = "launchpad_uploads"
_MAX_BROWSE_ENTRIES = 5000
_MAX_UPLOAD_BYTES = 4 * 1024 * 1024 * 1024  # 4 GiB hard ceiling per file
# Filesystems are byte-oriented, but we keep upload filenames ASCII-safe to
# avoid path surprises when the same string later appears in shell commands.
_FILENAME_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def _trusted_roots_resolved() -> list[Path]:
    """Resolved trusted-root Paths (config must already be loaded)."""
    cfg = get_config()
    roots: list[Path] = []
    for raw in cfg.trusted_roots or []:
        try:
            p = Path(str(raw)).expanduser().resolve()
        except OSError:
            continue
        if p.is_dir():
            roots.append(p)
    return roots


def _resolve_under_trusted(path_str: str) -> Path | None:
    """Resolve ``path_str`` and return it only if it lives under a trusted root."""
    try:
        resolved = Path(path_str).expanduser().resolve()
    except OSError:
        return None
    for root in _trusted_roots_resolved():
        if resolved == root:
            return resolved
        try:
            resolved.relative_to(root)
        except ValueError:
            continue
        return resolved
    return None


def _entry_payload(name: str, child: Path) -> dict:
    """Stat one directory entry; tolerate permission / race errors."""
    try:
        st = child.lstat()
    except OSError:
        return {"name": name, "kind": "other", "size": None, "mtime": None}
    kind: str
    size: int | None
    if child.is_dir():
        kind = "dir"
        size = None
    elif child.is_file():
        kind = "file"
        size = int(st.st_size)
    elif child.is_symlink():
        # Broken or non-resolving symlink.
        kind = "other"
        size = None
    else:
        kind = "other"
        size = None
    return {
        "name": name,
        "kind": kind,
        "size": size,
        "mtime": float(st.st_mtime),
    }


def _sanitize_upload_filename(raw: str) -> str:
    """Strip directory separators and risky characters from an uploaded filename."""
    base = os.path.basename(raw or "").strip()
    if not base or base in (".", ".."):
        return "upload.bin"
    cleaned = _FILENAME_SAFE.sub("_", base)
    # Avoid leading dots producing hidden files.
    cleaned = cleaned.lstrip(".") or "upload.bin"
    return cleaned[:200]


def _upload_dir() -> Path:
    """Create (if needed) and return the base directory for uploaded files."""
    cfg = get_config()
    base = Path(str(cfg.paths.temp_root)).expanduser().resolve() / _UPLOAD_SUBDIR
    base.mkdir(parents=True, exist_ok=True)
    return base


def _pick_bridge_vsix() -> Path | None:
    """Return the newest ``.vsix`` in extension/dist/, or None if none built."""
    if not _BRIDGE_DIST.is_dir():
        return None
    candidates = sorted(
        _BRIDGE_DIST.glob("*.vsix"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return candidates[0] if candidates else None


async def lp_bridge_extension_info(_request: Request) -> JSONResponse:
    """Return metadata so the Launchpad can show a helpful install box."""
    vsix = _pick_bridge_vsix()
    if vsix is None:
        return JSONResponse(
            {
                "available": False,
                "message": (
                    "No pre-built bridge extension is shipped with this STAR "
                    "Server install. Build it with `cd extension && npm install "
                    "&& npm run package` to populate extension/dist/."
                ),
            }
        )
    stat = vsix.stat()
    return JSONResponse(
        {
            "available": True,
            "filename": vsix.name,
            "size": int(stat.st_size),
            "download_url": "/launchpad/api/bridge/extension",
            "install_command": f"code --install-extension {vsix.name}",
        }
    )


async def lp_bridge_extension_download(_request: Request) -> FileResponse | JSONResponse:
    """Stream the packaged bridge .vsix to the browser."""
    vsix = _pick_bridge_vsix()
    if vsix is None:
        return _json_error(
            "NOT_FOUND",
            "No bridge extension .vsix has been built on this server.",
            404,
        )
    return FileResponse(
        path=str(vsix),
        filename=vsix.name,
        media_type="application/octet-stream",
    )


async def lp_capabilities(request: Request) -> JSONResponse:
    """Whether this browser session may use server-side Run in shell (loopback only)."""
    trusted_local = launchpad_request_trusted_local(request)
    bs = backend_status()
    ready = bool(bs.get("ready"))
    try:
        roots = [str(p) for p in _trusted_roots_resolved()]
    except Exception:
        roots = []
    return JSONResponse(
        {
            "launch_supported": trusted_local,
            "server_path_check_supported": trusted_local,
            "script_lane_utils_ready": ready,
            "script_lane_translate_supported": ready,
            "script_lane_viability_inputs_supported": bool(ready and trusted_local),
            "script_lane_execute_supported": bool(ready and trusted_local),
            "trusted_local": trusted_local,
            "browse_supported": bool(roots),
            "upload_supported": True,
            "browse_roots": roots,
            "bridge_extension_available": _pick_bridge_vsix() is not None,
        }
    )


async def lp_browse(request: Request) -> JSONResponse:
    """List directory entries under a ``trusted_roots`` path.

    With no ``path`` query param, returns the configured roots as virtual entries
    so clients can render a "pick a root" first step.
    """
    raw_path = (request.query_params.get("path") or "").strip()
    roots = _trusted_roots_resolved()
    if not roots:
        return _json_error(
            "NO_TRUSTED_ROOTS",
            "Server has no trusted_roots configured; filesystem browsing is disabled.",
            503,
        )

    if not raw_path:
        return JSONResponse(
            {
                "path": None,
                "parent": None,
                "is_root_list": True,
                "entries": [
                    {"name": str(p), "kind": "dir", "size": None, "mtime": None}
                    for p in roots
                ],
                "roots": [str(p) for p in roots],
            }
        )

    target = _resolve_under_trusted(raw_path)
    if target is None:
        return _json_error(
            "FORBIDDEN",
            "Path is outside the server's trusted_roots.",
            403,
        )
    if not target.exists():
        return _json_error("NOT_FOUND", f"Path does not exist: {target}", 404)
    if not target.is_dir():
        return _json_error(
            "BAD_REQUEST",
            "Path is a file; browse only lists directories.",
            400,
        )

    try:
        names = sorted(os.listdir(target), key=lambda n: n.lower())
    except PermissionError:
        return _json_error("FORBIDDEN", f"Permission denied reading: {target}", 403)
    except OSError as e:
        return _json_error("INTERNAL_ERROR", f"Could not list directory: {e}", 500)

    entries: list[dict] = []
    truncated = False
    for name in names:
        if len(entries) >= _MAX_BROWSE_ENTRIES:
            truncated = True
            break
        entries.append(_entry_payload(name, target / name))
    # Directories before files, each group alphabetical (case-insensitive).
    entries.sort(key=lambda e: (0 if e["kind"] == "dir" else 1, e["name"].lower()))

    parent: str | None
    try:
        # Only surface a parent link when it is itself within trusted_roots.
        p_candidate = target.parent
        parent = str(p_candidate) if _resolve_under_trusted(str(p_candidate)) else None
    except Exception:
        parent = None

    return JSONResponse(
        {
            "path": str(target),
            "parent": parent,
            "is_root_list": False,
            "entries": entries,
            "truncated": truncated,
            "roots": [str(p) for p in roots],
        }
    )


async def lp_upload(request: Request) -> JSONResponse:
    """Accept a multipart file upload and stage it under ``temp_root``.

    Returns the absolute server path that workflow params can reference.
    """
    content_type = (request.headers.get("content-type") or "").lower()
    if not content_type.startswith("multipart/form-data"):
        return _json_error(
            "BAD_REQUEST",
            "Upload requires multipart/form-data with a 'file' field.",
            400,
        )
    declared = request.headers.get("content-length")
    if declared and declared.isdigit() and int(declared) > _MAX_UPLOAD_BYTES:
        return _json_error(
            "PAYLOAD_TOO_LARGE",
            f"Upload exceeds server limit of {_MAX_UPLOAD_BYTES} bytes.",
            413,
        )
    try:
        form = await request.form()
    except Exception as e:
        return _json_error("BAD_REQUEST", f"Could not parse multipart body: {e}", 400)

    upload = form.get("file")
    try:
        from starlette.datastructures import UploadFile
    except ImportError:  # pragma: no cover
        UploadFile = None  # type: ignore[assignment]

    if UploadFile is None or not isinstance(upload, UploadFile):
        return _json_error(
            "BAD_REQUEST",
            "Form must include a 'file' field with the uploaded content.",
            400,
        )

    filename = _sanitize_upload_filename(upload.filename or "upload.bin")
    try:
        base = _upload_dir()
    except Exception as e:
        return _json_error("INTERNAL_ERROR", f"Upload staging unavailable: {e}", 500)

    staging_dir = base / uuid.uuid4().hex
    staging_dir.mkdir(parents=True, exist_ok=False)
    dest = staging_dir / filename

    total = 0
    try:
        with open(dest, "wb") as sink:
            while True:
                chunk = await upload.read(1024 * 1024)
                if not chunk:
                    break
                total += len(chunk)
                if total > _MAX_UPLOAD_BYTES:
                    sink.close()
                    try:
                        dest.unlink(missing_ok=True)
                        staging_dir.rmdir()
                    except OSError:
                        pass
                    return _json_error(
                        "PAYLOAD_TOO_LARGE",
                        f"Upload exceeds server limit of {_MAX_UPLOAD_BYTES} bytes.",
                        413,
                    )
                sink.write(chunk)
    except Exception as e:
        try:
            dest.unlink(missing_ok=True)
            staging_dir.rmdir()
        except OSError:
            pass
        return _json_error("INTERNAL_ERROR", f"Could not write upload: {e}", 500)
    finally:
        try:
            await upload.close()
        except Exception:
            pass

    return JSONResponse(
        {
            "ok": True,
            "path": str(dest.resolve()),
            "filename": filename,
            "size": total,
        }
    )


def _script_lane_nonempty_input_values(body: dict) -> bool:
    """True when viability should treat ``input_values`` as host-sensitive."""
    raw = body.get("input_values")
    return isinstance(raw, dict) and bool(raw)


async def lp_script_lane_translate(request: Request) -> JSONResponse:
    """Annotated Bash -> IR (or refusal). No host path probes; any client may call."""
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    if not isinstance(body, dict):
        return _json_error("BAD_REQUEST", "Body must be a JSON object", 400)
    script_text = body.get("script_text")
    if not isinstance(script_text, str):
        return _json_error("BAD_REQUEST", "Body must contain string 'script_text'", 400)
    try:
        result = translate_simple_script_to_ir(script_text)
    except ScriptLaneBackendUnavailable as e:
        return JSONResponse(
            {"error": True, "code": "SERVICE_UNAVAILABLE", "message": str(e)},
            status_code=503,
        )
    except Exception as e:
        return _json_error("INTERNAL_ERROR", str(e), 500)
    return JSONResponse(result)


async def lp_script_lane_viability(request: Request) -> JSONResponse:
    """Simple-local viability; non-empty ``input_values`` is loopback-only (path checks)."""
    trusted_local = launchpad_request_trusted_local(request)
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    if not isinstance(body, dict):
        return _json_error("BAD_REQUEST", "Body must be a JSON object", 400)
    ir = body.get("ir")
    if not isinstance(ir, dict):
        return _json_error("BAD_REQUEST", "Body must contain object 'ir'", 400)
    if _script_lane_nonempty_input_values(body) and not trusted_local:
        return _json_error(
            "FORBIDDEN",
            "Viability with input_values (server-host path checks) is only available "
            "from the same machine as the server.",
            403,
        )
    raw_iv = body.get("input_values")
    input_values: dict | None
    if raw_iv is None:
        input_values = None
    elif isinstance(raw_iv, dict):
        # Treat {} like omitted: backend distinguishes None (skip value checks) from
        # a dict (run required-input / path validation even when keys are missing).
        input_values = dict(raw_iv) if raw_iv else None
    else:
        return _json_error("BAD_REQUEST", "'input_values' must be an object or omitted", 400)
    raw_vs = body.get("validate_schema", True)
    validate_schema = bool(raw_vs) if isinstance(raw_vs, bool) else True
    try:
        result = check_simple_local_ir_viability(
            ir,
            input_values=input_values,
            validate_schema=validate_schema,
        )
    except ScriptLaneBackendUnavailable as e:
        return JSONResponse(
            {"error": True, "code": "SERVICE_UNAVAILABLE", "message": str(e)},
            status_code=503,
        )
    except Exception as e:
        return _json_error("INTERNAL_ERROR", str(e), 500)
    return JSONResponse(result)


async def lp_script_lane_execute(request: Request) -> JSONResponse:
    """Run simple-local IR on the server host (loopback peers only)."""
    if not launchpad_request_trusted_local(request):
        return _json_error(
            "FORBIDDEN",
            "Script Lane local execute is only available when Launchpad is opened from "
            "the same machine as the server (e.g. http://127.0.0.1:…/launchpad/).",
            403,
        )
    try:
        body = await request.json()
    except Exception:
        return _json_error("BAD_REQUEST", "Expected JSON body", 400)
    if not isinstance(body, dict):
        return _json_error("BAD_REQUEST", "Body must be a JSON object", 400)
    ir = body.get("ir")
    if not isinstance(ir, dict):
        return _json_error("BAD_REQUEST", "Body must contain object 'ir'", 400)
    raw_iv = body.get("input_values")
    if not isinstance(raw_iv, dict):
        return _json_error("BAD_REQUEST", "Body must contain object 'input_values'", 400)
    workdir = body.get("workdir")
    if not isinstance(workdir, str) or not workdir.strip():
        return _json_error(
            "BAD_REQUEST",
            "Body must contain non-empty string 'workdir'",
            400,
        )
    vac = body.get("viability_already_checked", False)
    viability_already_checked = bool(vac) if isinstance(vac, bool) else False
    try:
        result = execute_simple_local_ir(
            ir,
            input_values=dict(raw_iv),
            workdir=workdir.strip(),
            viability_already_checked=viability_already_checked,
        )
    except ScriptLaneBackendUnavailable as e:
        return JSONResponse(
            {"error": True, "code": "SERVICE_UNAVAILABLE", "message": str(e)},
            status_code=503,
        )
    except Exception as e:
        return _json_error("INTERNAL_ERROR", str(e), 500)
    return JSONResponse(result)


def _extract_check_paths(body: dict, *, trusted_local: bool) -> tuple[bool, JSONResponse | None]:
    """Return requested path-check mode, rejecting non-local server path probes."""
    raw = body.get("check_paths", False)
    if not isinstance(raw, bool):
        return False, _json_error("BAD_REQUEST", "'check_paths' must be a boolean", 400)
    if raw and not trusted_local:
        return (
            False,
            _json_error(
                "FORBIDDEN",
                "Server-host path checks are only available from the same machine as the server.",
                403,
            ),
        )
    return raw, None


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
    trusted_local = launchpad_request_trusted_local(request)
    check_paths, error = _extract_check_paths(body, trusted_local=trusted_local)
    if error is not None:
        return error
    try:
        result = validate_workflow_parameters(
            wf_id, params, check_paths=check_paths
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
        val = validate_workflow_parameters(wf_id, params, check_paths=True)
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
            "/launchpad/api/browse",
            endpoint=lp_browse,
            methods=["GET"],
        ),
        Route(
            "/launchpad/api/upload",
            endpoint=lp_upload,
            methods=["POST"],
        ),
        Route(
            "/launchpad/api/bridge/extension/info",
            endpoint=lp_bridge_extension_info,
            methods=["GET"],
        ),
        Route(
            "/launchpad/api/bridge/extension",
            endpoint=lp_bridge_extension_download,
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
        Route(
            "/launchpad/api/script-lane/translate",
            endpoint=lp_script_lane_translate,
            methods=["POST"],
        ),
        Route(
            "/launchpad/api/script-lane/viability",
            endpoint=lp_script_lane_viability,
            methods=["POST"],
        ),
        Route(
            "/launchpad/api/script-lane/execute",
            endpoint=lp_script_lane_execute,
            methods=["POST"],
        ),
        Mount(
            "/launchpad",
            app=StaticFiles(directory=str(_LAUNCHPAD_STATIC), html=True),
        ),
    ]
