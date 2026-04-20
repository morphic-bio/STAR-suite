"""Launchpad Composition tab: thin adapter to ``bwb-nextflow-utils`` headless constructor.

Resolution order for the utils checkout matches :mod:`script_lane_bridge` (env
``BWB_NEXTFLOW_UTILS_ROOT`` or sibling ``../bwb-nextflow-utils``).

Bulk STAR composition defaults load Workflow IR from **STAR-suite** JSON under
``launchpad/composition_ir/`` (maintained with ``mcp_server/workflows/star_*.yaml``),
not from ``bwb-nextflow-utils`` test fixtures—so Launchpad stays STAR-first for
bulk IR drift signals (missing files, constructor failures after YAML changes).
"""

from __future__ import annotations

import copy
import json
import re
from pathlib import Path
from typing import Any

from .script_lane_bridge import ScriptLaneBackendUnavailable, get_mcp_workflow_services

HEADLESS_COMPOSITION_ARTIFACTS_SCHEMA_ID = "biodepot.launchpad_composition_artifacts/v0"

_BULK_ARCHETYPE_ID = "star_index_bulk_pe_then_bulk_de"
_DEFAULT_BULK_COMPONENT_BINDING: dict[str, str] = {
    "mkindex": "star_genome_generate",
    "align": "star_bulk_pe_batch",
}
_BULK_DEFAULT_WORKFLOW_IDS = ("star_genome_generate", "star_bulk_pe_batch")

_RECIPE_INPUT_FROM = re.compile(r"^recipe\.inputs\.(?P<key>[A-Za-z0-9_.-]+)$")

_COMPOSITION_IR_DIR = Path(__file__).resolve().parent / "composition_ir"


def _load_star_bulk_composition_ir_map() -> dict[str, Any]:
    """Load STAR-suite bulk composition IR (see ``composition_ir/README.md``)."""
    out: dict[str, Any] = {}
    for wid in _BULK_DEFAULT_WORKFLOW_IDS:
        path = _COMPOSITION_IR_DIR / f"{wid}.json"
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing STAR-suite composition IR for {wid!r}: {path}",
            )
        doc = json.loads(path.read_text())
        meta = doc.get("metadata") if isinstance(doc.get("metadata"), dict) else {}
        meta_id = str(meta.get("id") or "")
        if meta_id != wid:
            raise ValueError(
                f"Composition IR metadata.id mismatch for {path}: "
                f"expected {wid!r}, got {meta_id!r}",
            )
        out[wid] = doc
    return out


def describe_composition_recipe_inputs(
    profile_id: str,
    workflow_ir_by_id: dict[str, Any] | None = None,
    component_binding: dict[str, str] | None = None,
) -> dict[str, Any]:
    """``describe_constructor_profile_recipe_inputs`` with merged IR + binding defaults."""
    _bundle, merged_ir, merged_binding = merge_composition_context(
        profile_id,
        workflow_ir_by_id,
        component_binding,
    )
    return get_mcp_workflow_services().describe_constructor_profile_recipe_inputs(
        profile_id,
        merged_ir,
        component_binding=merged_binding,
    )


def merge_composition_context(
    profile_id: str,
    workflow_ir_by_id: dict[str, Any] | None,
    component_binding: dict[str, str] | None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    """Resolve profile bundle and merge caller IR/binding over archetype defaults."""
    svc = get_mcp_workflow_services()
    bundle = svc.resolve_profile_to_constructor_bundle(profile_id)
    archetype_id = str(bundle.get("archetype_id") or "")

    default_ir: dict[str, Any] = {}
    default_binding: dict[str, str] = {}
    if archetype_id == _BULK_ARCHETYPE_ID:
        default_ir = _load_star_bulk_composition_ir_map()
        default_binding = dict(_DEFAULT_BULK_COMPONENT_BINDING)

    merged_ir = dict(default_ir)
    if workflow_ir_by_id:
        merged_ir.update(dict(workflow_ir_by_id))

    merged_binding = dict(default_binding)
    if component_binding:
        merged_binding.update({str(k): str(v) for k, v in component_binding.items()})

    return bundle, merged_ir, merged_binding


def _coerce_recipe_input_value(raw: Any, kind: str) -> Any:
    if raw is None:
        return None
    k = (kind or "string").strip().lower()
    if k == "boolean":
        if isinstance(raw, bool):
            return raw
        s = str(raw).strip().lower()
        if s in ("true", "1", "yes", "on"):
            return True
        if s in ("false", "0", "no", "off", ""):
            return False
        raise ValueError(f"not a boolean: {raw!r}")
    if k == "integer":
        if isinstance(raw, bool):
            raise ValueError("integer value cannot be bool")
        if isinstance(raw, int):
            return raw
        s = str(raw).strip()
        if not s:
            raise ValueError("empty integer")
        return int(s, 10)
    if k == "number":
        if isinstance(raw, bool):
            raise ValueError("number value cannot be bool")
        if isinstance(raw, (int, float)):
            return float(raw)
        s = str(raw).strip()
        if not s:
            raise ValueError("empty number")
        return float(s)
    if k in ("path_list", "string_list"):
        if isinstance(raw, list):
            return [str(x) for x in raw]
        s = str(raw).strip()
        if not s:
            return []
        if s.startswith("["):
            parsed = json.loads(s)
            if not isinstance(parsed, list):
                raise ValueError("JSON list required for array input")
            return [str(x) for x in parsed]
        return [p.strip() for p in s.split(",") if p.strip()]
    # path, string, default
    return str(raw)


def apply_recipe_input_literals(
    recipe: dict[str, Any],
    recipe_input_values: dict[str, Any],
) -> dict[str, Any]:
    """Replace ``from: recipe.inputs.*`` bindings with ``value`` literals where provided.

    Purely structural: does not interpret profile semantics beyond matching
    binding shapes already emitted by the headless constructor.
    """
    if not recipe_input_values:
        return recipe
    specs = recipe.get("inputs")
    if not isinstance(specs, dict):
        specs = {}

    out = copy.deepcopy(recipe)
    steps = out.get("steps")
    if not isinstance(steps, list):
        return out

    for step in steps:
        if not isinstance(step, dict):
            continue
        bindings = step.get("inputs")
        if not isinstance(bindings, dict):
            continue
        for port, binding in list(bindings.items()):
            if not isinstance(binding, dict):
                continue
            ref = binding.get("from")
            if not isinstance(ref, str):
                continue
            m = _RECIPE_INPUT_FROM.match(ref.strip())
            if not m:
                continue
            key = m.group("key")
            if key not in recipe_input_values:
                continue
            raw_val = recipe_input_values[key]
            if raw_val is None or (isinstance(raw_val, str) and not raw_val.strip()):
                continue
            kind = "string"
            sp = specs.get(key)
            if isinstance(sp, dict) and isinstance(sp.get("kind"), str):
                kind = str(sp["kind"])
            try:
                coerced = _coerce_recipe_input_value(raw_val, kind)
            except (TypeError, ValueError, json.JSONDecodeError) as e:
                raise ValueError(f"recipe_input_values[{key!r}]: {e}") from e
            bindings[port] = {"value": coerced}
    return out


def build_composition_draft(
    profile_id: str,
    *,
    workflow_ir_by_id: dict[str, Any] | None = None,
    component_binding: dict[str, str] | None = None,
    recipe_input_values: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Call ``build_headless_recipe_draft`` with merged context + optional literals."""
    svc = get_mcp_workflow_services()
    bundle, merged_ir, merged_binding = merge_composition_context(
        profile_id,
        workflow_ir_by_id,
        component_binding,
    )
    draft = svc.build_headless_recipe_draft(
        archetype_id=str(bundle["archetype_id"]),
        component_binding=merged_binding,
        workflow_ir_by_id=merged_ir,
        recipe_id=str(bundle.get("default_recipe_id") or "") or None,
        helper_ids_by_role=dict(bundle.get("helper_ids_by_role") or {}),
    )
    if recipe_input_values:
        draft = apply_recipe_input_literals(draft, dict(recipe_input_values))
    return {
        "profile_id": bundle["profile_id"],
        "archetype_id": bundle["archetype_id"],
        "helper_ids_by_role": dict(bundle.get("helper_ids_by_role") or {}),
        "component_binding": dict(merged_binding),
        "workflow_ir_by_id_keys": sorted(merged_ir.keys()),
        "recipe_draft": draft,
    }


def build_composition_artifacts(
    profile_id: str,
    *,
    workflow_ir_by_id: dict[str, Any] | None = None,
    component_binding: dict[str, str] | None = None,
    recipe_input_values: dict[str, Any] | None = None,
    strict: bool = True,
) -> dict[str, Any]:
    """Draft + validate + normalize + preview + script + manifest (side-effect free)."""
    svc = get_mcp_workflow_services()
    built = build_composition_draft(
        profile_id,
        workflow_ir_by_id=workflow_ir_by_id,
        component_binding=component_binding,
        recipe_input_values=recipe_input_values,
    )
    draft = built["recipe_draft"]
    _, merged_ir, _ = merge_composition_context(
        profile_id,
        workflow_ir_by_id,
        component_binding,
    )
    expanded = svc.ensure_workflow_ir_map_for_recipe(draft, merged_ir)
    validation = svc.validate_execution_recipe(
        draft,
        workflow_ir_by_id=expanded,
        strict=bool(strict),
    )
    normalized: dict[str, Any] | None = None
    preview: dict[str, Any] | None = None
    script: dict[str, Any] | None = None
    manifest: dict[str, Any] | None = None
    if validation.get("ok"):
        normalized = svc.normalize_execution_recipe(
            draft,
            workflow_ir_by_id=expanded,
            strict=bool(strict),
        )
        preview = svc.preview_execution_recipe(
            draft,
            workflow_ir_by_id=expanded,
            strict=bool(strict),
        )
        script = svc.generate_execution_recipe_script(
            recipe=draft,
            normalized_recipe=normalized,
            preview=preview,
            workflow_ir_by_id=expanded,
        )
        manifest = svc.generate_execution_recipe_manifest(
            recipe=draft,
            normalized_recipe=normalized,
            preview=preview,
            script=script,
        )
    return {
        "schema": HEADLESS_COMPOSITION_ARTIFACTS_SCHEMA_ID,
        "profile_id": built["profile_id"],
        "archetype_id": built["archetype_id"],
        "helper_ids_by_role": built["helper_ids_by_role"],
        "component_binding": built["component_binding"],
        "workflow_ir_by_id_keys": built["workflow_ir_by_id_keys"],
        "recipe_id": draft.get("id"),
        "recipe_draft": draft,
        "validation": validation,
        "normalized_recipe": normalized,
        "preview": preview,
        "generated_script": script,
        "generated_manifest": manifest,
    }


__all__ = [
    "HEADLESS_COMPOSITION_ARTIFACTS_SCHEMA_ID",
    "ScriptLaneBackendUnavailable",
    "apply_recipe_input_literals",
    "build_composition_artifacts",
    "build_composition_draft",
    "describe_composition_recipe_inputs",
    "merge_composition_context",
]
