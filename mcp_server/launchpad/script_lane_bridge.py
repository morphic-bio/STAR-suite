"""Thin adapter from STAR Launchpad to bwb-nextflow-utils simple-script lane.

Canonical translation, viability, and execution live in bwb-nextflow-utils
(``tools/mcp_workflow_services.py``). This module only locates that checkout
and re-exports calls for the Launchpad HTTP layer.

Resolution order for the utils checkout:

1. Environment variable ``BWB_NEXTFLOW_UTILS_ROOT`` (absolute path to repo root).
2. Sibling directory ``../bwb-nextflow-utils`` next to the STAR-suite repo root
   (typical coordinated stack layout on one machine).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

_svc: ModuleType | None = None


class ScriptLaneBackendUnavailable(RuntimeError):
    """Raised when no usable bwb-nextflow-utils checkout is on ``sys.path``."""


def _candidate_roots() -> list[Path]:
    out: list[Path] = []
    env = os.environ.get("BWB_NEXTFLOW_UTILS_ROOT", "").strip()
    if env:
        out.append(Path(env).resolve())
    star_root = Path(__file__).resolve().parent.parent.parent
    out.append((star_root.parent / "bwb-nextflow-utils").resolve())
    # De-duplicate while preserving order
    seen: set[Path] = set()
    uniq: list[Path] = []
    for p in out:
        if p in seen:
            continue
        seen.add(p)
        uniq.append(p)
    return uniq


def _ensure_services() -> ModuleType:
    global _svc
    if _svc is not None:
        return _svc
    last_err: Exception | None = None
    for root in _candidate_roots():
        tools = root / "tools"
        if not (tools / "mcp_workflow_services.py").is_file():
            continue
        rs, ts = str(root), str(tools)
        inserted: list[str] = []
        for p in (ts, rs):
            if p not in sys.path:
                sys.path.insert(0, p)
                inserted.append(p)
        try:
            import mcp_workflow_services as svc  # type: ignore[import-untyped, import-not-found]
        except Exception as e:
            last_err = e
            for p in inserted:
                try:
                    sys.path.remove(p)
                except ValueError:
                    pass
            continue
        _svc = svc
        return svc
    msg = (
        "bwb-nextflow-utils not found or failed to import. Set BWB_NEXTFLOW_UTILS_ROOT "
        "to the repository root (directory containing tools/mcp_workflow_services.py), "
        "or place a checkout beside STAR-suite as ../bwb-nextflow-utils."
    )
    if last_err is not None:
        msg = f"{msg} Last import error: {last_err!r}"
    raise ScriptLaneBackendUnavailable(msg)


def backend_status() -> dict[str, Any]:
    """Whether the adapter can load ``mcp_workflow_services``."""
    try:
        s = _ensure_services()
        root = Path(s.__file__).resolve().parent.parent
        return {"ready": True, "utils_root": str(root)}
    except ScriptLaneBackendUnavailable as e:
        return {"ready": False, "message": str(e)}


def get_mcp_workflow_services() -> ModuleType:
    """Return the loaded ``mcp_workflow_services`` module (shared by Script Lane + Composition)."""
    return _ensure_services()


def translate_simple_script_to_ir(script_text: str) -> dict[str, Any]:
    return _ensure_services().translate_simple_script_to_ir(script_text)


def check_simple_local_ir_viability(
    ir_document: dict[str, Any],
    *,
    input_values: dict[str, Any] | None = None,
    validate_schema: bool = True,
) -> dict[str, Any]:
    return _ensure_services().check_simple_local_ir_viability(
        ir_document,
        input_values=input_values,
        validate_schema=validate_schema,
    )


def execute_simple_local_ir(
    ir_document: dict[str, Any],
    *,
    input_values: dict[str, Any],
    workdir: str,
    viability_already_checked: bool = False,
) -> dict[str, Any]:
    return _ensure_services().execute_simple_local_ir(
        ir_document,
        input_values=input_values,
        workdir=workdir,
        viability_already_checked=viability_already_checked,
    )
