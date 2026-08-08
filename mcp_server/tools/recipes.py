"""Recipe catalog discovery, source reconciliation, locks, and bundles."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional

from packaging.version import Version

from ..config import (
    LoadedRecipeCatalog,
    WorkflowOrigin,
    get_config,
    get_recipe_catalogs,
    get_workflow_config,
    get_workflow_configs,
    get_workflow_origin,
    get_workflow_schema,
)
from ..schemas.config import RecipeResolutionPolicy, WorkflowConfig
from .workflows import render_workflow_command, validate_workflow_parameters


RECIPE_LOCK_SCHEMA = "biodepot.recipe_lock/v1"
RESOLVED_WORKFLOW_SCHEMA = "biodepot.resolved_workflow/v1"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_digest(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _git_identity(root: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, command in (
        ("commit", ["git", "rev-parse", "HEAD"]),
        ("remote", ["git", "remote", "get-url", "origin"]),
    ):
        try:
            completed = subprocess.run(
                command,
                cwd=str(root),
                capture_output=True,
                text=True,
                timeout=5,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if completed.returncode == 0 and completed.stdout.strip():
            result[key] = completed.stdout.strip()
    return result


def _suite_identity() -> dict[str, Any]:
    """Return portable STAR Suite source/version identity for a lock."""
    root = get_config().paths.repo_root.resolve()
    identity: dict[str, Any] = {
        "source_revision": _git_identity(root).get("commit"),
    }
    version_file = root / "core" / "legacy" / "source" / "VERSION"
    if version_file.is_file():
        text = version_file.read_text(encoding="utf-8", errors="replace")
        match = re.search(r'STAR_SUITE_VERSION\s+"([^"]+)"', text)
        if match is not None:
            identity["version"] = match.group(1)
        identity["version_file_sha256"] = _sha256_file(version_file)
    return identity


def _visible(workflow: WorkflowConfig, authenticated: bool) -> bool:
    return authenticated or workflow.visibility != "private"


def _compatible(workflow: WorkflowConfig, application: Optional[str]) -> bool:
    return (
        application is None
        or not workflow.applications
        or application in workflow.applications
    )


def _logical_id(workflow: WorkflowConfig) -> str:
    return workflow.logical_id or workflow.id


def _source_ref(workflow_id: str, origin: WorkflowOrigin) -> str:
    return f"{origin.catalog.id}::{workflow_id}"


def _lineage(workflow_id: str) -> list[str]:
    lineage: list[str] = []
    current: Optional[str] = workflow_id
    while current is not None:
        lineage.append(current)
        workflow = get_workflow_config(current)
        current = workflow.extends if workflow is not None else None
    lineage.reverse()
    return lineage


def _candidate(workflow_id: str) -> dict[str, Any]:
    workflow = get_workflow_config(workflow_id)
    origin = get_workflow_origin(workflow_id)
    if workflow is None or origin is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")
    return {
        "workflow_id": workflow.id,
        "source_ref": _source_ref(workflow.id, origin),
        "logical_id": _logical_id(workflow),
        "version": workflow.version,
        "title": workflow.title,
        "summary": workflow.summary,
        "visibility": workflow.visibility,
        "applications": list(workflow.applications),
        "extends": workflow.extends,
        "replaces": list(workflow.replaces),
        "lineage": _lineage(workflow.id),
        "image": workflow.image,
        "catalog": {
            "id": origin.catalog.id,
            "namespace": origin.catalog.namespace,
            "version": origin.catalog.version,
            "order": origin.catalog.order,
            "built_in": origin.catalog.built_in,
            "trust": origin.catalog.trust,
        },
    }


def list_recipe_catalogs(authenticated: bool = False) -> dict[str, Any]:
    """Return loaded catalogs in deterministic discovery order."""
    workflows = get_workflow_configs()
    origins = {
        workflow_id: get_workflow_origin(workflow_id) for workflow_id in workflows
    }
    items: list[dict[str, Any]] = []
    for catalog in get_recipe_catalogs():
        visible_ids = [
            workflow_id
            for workflow_id, workflow in workflows.items()
            if origins[workflow_id] is not None
            and origins[workflow_id].catalog.id == catalog.id  # type: ignore[union-attr]
            and _visible(workflow, authenticated)
        ]
        item: dict[str, Any] = {
            "id": catalog.id,
            "namespace": catalog.namespace,
            "version": catalog.version,
            "order": catalog.order,
            "built_in": catalog.built_in,
            "workflow_count": len(visible_ids),
            "source_kind": "builtin" if catalog.built_in else "local",
            "trust": catalog.trust,
        }
        if catalog.manifest_path is not None:
            item["manifest"] = catalog.manifest_path.name
        if authenticated:
            item["root"] = str(catalog.root)
            if catalog.manifest_path is not None:
                item["manifest_path"] = str(catalog.manifest_path)
            item["git"] = _git_identity(catalog.root)
        items.append(item)
    return {"catalogs": items}


def describe_recipe_catalog(
    catalog_id: str, authenticated: bool = False
) -> dict[str, Any]:
    """Describe one catalog and its visible recipes."""
    catalog = next((item for item in get_recipe_catalogs() if item.id == catalog_id), None)
    if catalog is None:
        raise ValueError(f"Unknown recipe catalog: {catalog_id}")

    catalog_item = next(
        item
        for item in list_recipe_catalogs(authenticated=authenticated)["catalogs"]
        if item["id"] == catalog_id
    )
    candidates = [
        _candidate(workflow_id)
        for workflow_id, workflow in get_workflow_configs().items()
        if (origin := get_workflow_origin(workflow_id)) is not None
        and origin.catalog.id == catalog_id
        and _visible(workflow, authenticated)
    ]
    return {"catalog": catalog_item, "recipes": candidates}


def list_provenance_repositories(
    authenticated: bool = False,
) -> dict[str, Any]:
    """Describe the ordered provenance search path and explicit write target."""
    provenance = get_config().provenance
    search = []
    for order, repository in enumerate(provenance.search):
        item: dict[str, Any] = {"id": repository.id, "order": order}
        if authenticated:
            item["root"] = str(repository.root)
            item["exists"] = repository.root.is_dir()
        search.append(item)

    write: Optional[dict[str, Any]] = None
    if provenance.write is not None:
        write = {"id": provenance.write.id}
        if authenticated:
            write["root"] = str(provenance.write.root)
            write["exists"] = provenance.write.root.is_dir()
    return {"search": search, "write": write}


def list_recipe_candidates(
    logical_id: Optional[str] = None,
    application: Optional[str] = None,
    authenticated: bool = False,
) -> dict[str, Any]:
    """List source-specific candidates grouped by logical recipe identity."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for workflow_id, workflow in get_workflow_configs().items():
        candidate_logical_id = _logical_id(workflow)
        if logical_id is not None and candidate_logical_id != logical_id:
            continue
        if not _visible(workflow, authenticated) or not _compatible(
            workflow, application
        ):
            continue
        grouped.setdefault(candidate_logical_id, []).append(_candidate(workflow_id))

    groups = []
    for candidate_logical_id in sorted(grouped):
        candidates = sorted(
            grouped[candidate_logical_id],
            key=lambda item: (item["catalog"]["order"], item["workflow_id"]),
        )
        groups.append(
            {
                "logical_id": candidate_logical_id,
                "conflict": len(candidates) > 1,
                "candidates": candidates,
            }
        )
    return {"application": application, "groups": groups}


def list_recipe_conflicts(
    application: Optional[str] = None,
    authenticated: bool = False,
) -> dict[str, Any]:
    """Return only logical identities supplied by multiple visible sources."""
    listed = list_recipe_candidates(
        application=application, authenticated=authenticated
    )
    return {
        "application": application,
        "conflicts": [group for group in listed["groups"] if group["conflict"]],
    }


def _select_explicit(
    candidates: list[dict[str, Any]], selected_workflow_id: str
) -> Optional[dict[str, Any]]:
    return next(
        (
            item
            for item in candidates
            if selected_workflow_id in (item["workflow_id"], item["source_ref"])
        ),
        None,
    )


def resolve_recipe(
    reference: str,
    *,
    application: Optional[str] = None,
    policy: Optional[RecipeResolutionPolicy] = None,
    selected_workflow_id: Optional[str] = None,
    authenticated: bool = False,
) -> dict[str, Any]:
    """Resolve a logical id according to an explicit consumer policy.

    ``catalog-id::workflow-id`` is an unambiguous source reference. A plain
    fully qualified workflow id is used directly when it is not also a shared
    logical id. Shared logical identities are reconciled by policy.
    """
    configured_policy = get_config().recipe_resolution.policy_for(application)
    effective_policy = policy or configured_policy

    if "::" in reference:
        catalog_id, workflow_id = reference.split("::", 1)
        workflow = get_workflow_config(workflow_id)
        origin = get_workflow_origin(workflow_id)
        if (
            workflow is None
            or origin is None
            or origin.catalog.id != catalog_id
            or not _visible(workflow, authenticated)
            or not _compatible(workflow, application)
        ):
            raise ValueError(f"Unknown or incompatible recipe source: {reference}")
        candidate = _candidate(workflow_id)
        return {
            "status": "selected",
            "reference": reference,
            "application": application,
            "policy": effective_policy,
            "reason": "explicit_source",
            "selected": candidate,
            "candidates": [candidate],
        }

    logical_candidates = list_recipe_candidates(
        logical_id=reference,
        application=application,
        authenticated=authenticated,
    )["groups"]
    candidates = logical_candidates[0]["candidates"] if logical_candidates else []

    if not candidates:
        workflow = get_workflow_config(reference)
        if (
            workflow is None
            or not _visible(workflow, authenticated)
            or not _compatible(workflow, application)
        ):
            raise ValueError(f"Unknown or incompatible recipe: {reference}")
        candidates = [_candidate(reference)]

    if selected_workflow_id is not None:
        selected = _select_explicit(candidates, selected_workflow_id)
        if selected is None:
            raise ValueError(
                f"Selected workflow '{selected_workflow_id}' is not a candidate "
                f"for '{reference}'"
            )
        return {
            "status": "selected",
            "reference": reference,
            "application": application,
            "policy": effective_policy,
            "reason": "user_selection",
            "selected": selected,
            "candidates": candidates,
        }

    if len(candidates) == 1:
        return {
            "status": "selected",
            "reference": reference,
            "application": application,
            "policy": effective_policy,
            "reason": "only_candidate",
            "selected": candidates[0],
            "candidates": candidates,
        }

    if effective_policy == "prefer_newest":
        selected = max(
            candidates,
            key=lambda item: (
                Version(item["version"]),
                item["catalog"]["order"],
                item["workflow_id"],
            ),
        )
        return {
            "status": "selected",
            "reference": reference,
            "application": application,
            "policy": effective_policy,
            "reason": "highest_version_then_later_catalog",
            "selected": selected,
            "candidates": candidates,
        }

    status = "selection_required" if effective_policy == "prompt" else "separate"
    return {
        "status": status,
        "reference": reference,
        "application": application,
        "policy": effective_policy,
        "reason": (
            "consumer_must_prompt"
            if effective_policy == "prompt"
            else "use_a_source_ref_or_fully_qualified_workflow_id"
        ),
        "selected": None,
        "candidates": candidates,
    }


def _catalog_relative(path: Path, catalog: LoadedRecipeCatalog, declared: str) -> str:
    try:
        return path.resolve().relative_to(catalog.root).as_posix()
    except ValueError:
        return declared


def _workflow_sources(workflow_id: str, strict_files: bool) -> list[dict[str, Any]]:
    schema = get_workflow_schema(workflow_id)
    origin = get_workflow_origin(workflow_id)
    if schema is None or origin is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")

    declared_scripts = [
        (script.role, script.path, script.language) for script in schema.scripts
    ]
    declared_scripts.extend(
        (f"stage:{stage.name}", stage.script, "bash")
        for stage in schema.stages
        if stage.script is not None
    )
    if schema.execution is not None:
        declared_scripts.extend(
            (
                f"execution_stage:{stage.name}",
                stage.entry_script,
                "bash",
            )
            for stage in schema.execution.stages
        )
    if not declared_scripts:
        declared_scripts = [("entry", schema.entry_script, "bash")]

    sources: list[dict[str, Any]] = []
    seen: set[Path] = set()
    schema_path = origin.schema_path.resolve()
    for role, declared, language, path in [
        ("workflow_schema", origin.schema_file, "yaml", schema_path),
        *[
            (
                role,
                declared,
                language,
                (
                    Path(declared)
                    if Path(declared).is_absolute()
                    else origin.catalog.root / declared
                ).resolve(),
            )
            for role, declared, language in declared_scripts
        ],
    ]:
        if path in seen:
            continue
        seen.add(path)
        exists = path.is_file()
        if strict_files and not exists:
            raise FileNotFoundError(
                f"Locked source for workflow '{workflow_id}' does not exist: {path}"
            )
        entry: dict[str, Any] = {
            "role": role,
            "path": _catalog_relative(path, origin.catalog, declared),
            "declared_path": declared,
            "language": language,
            "exists": exists,
        }
        if exists:
            entry["sha256"] = _sha256_file(path)
            entry["size"] = path.stat().st_size
        sources.append(entry)
    return sources


def _locked_lineage(workflow_id: str, strict_files: bool) -> list[dict[str, Any]]:
    locked: list[dict[str, Any]] = []
    for lineage_id in _lineage(workflow_id):
        workflow = get_workflow_config(lineage_id)
        origin = get_workflow_origin(lineage_id)
        if workflow is None or origin is None:
            raise ValueError(f"Unknown workflow: {lineage_id}")
        locked.append(
            {
                "workflow_id": lineage_id,
                "logical_id": _logical_id(workflow),
                "version": workflow.version,
                "catalog_id": origin.catalog.id,
                "catalog_namespace": origin.catalog.namespace,
                "catalog_version": origin.catalog.version,
                "catalog_trust": origin.catalog.trust,
                "catalog_manifest_sha256": (
                    _sha256_file(origin.catalog.manifest_path)
                    if origin.catalog.manifest_path is not None
                    else None
                ),
                "source_revision": _git_identity(origin.catalog.root).get("commit"),
                "extends": workflow.extends,
                "sources": _workflow_sources(lineage_id, strict_files),
            }
        )
    return locked


def build_recipe_lock(
    reference: str,
    params: dict[str, Any],
    *,
    application: Optional[str] = None,
    policy: Optional[RecipeResolutionPolicy] = None,
    selected_workflow_id: Optional[str] = None,
    authenticated: bool = True,
    check_paths: bool = False,
    strict_files: bool = True,
) -> dict[str, Any]:
    """Build a deterministic lock without executing the recipe."""
    resolution = resolve_recipe(
        reference,
        application=application,
        policy=policy,
        selected_workflow_id=selected_workflow_id,
        authenticated=authenticated,
    )
    if resolution["selected"] is None:
        raise ValueError(
            f"Recipe '{reference}' is unresolved ({resolution['status']}); "
            "select a source explicitly"
        )
    workflow_id = resolution["selected"]["workflow_id"]
    workflow = get_workflow_config(workflow_id)
    origin = get_workflow_origin(workflow_id)
    if workflow is None or origin is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")

    validation = validate_workflow_parameters(
        workflow_id, params, check_paths=check_paths
    )
    if not validation.valid:
        raise ValueError("Invalid recipe parameters: " + "; ".join(validation.errors))
    rendered = render_workflow_command(workflow_id, validation.normalized_params)
    portable_argv = [workflow.entry_script, *rendered.argv[1:]]
    provenance = get_config().provenance

    manifest: dict[str, Any] = {
        "schema": RECIPE_LOCK_SCHEMA,
        "suite": _suite_identity(),
        "recipe": {
            "workflow_id": workflow_id,
            "logical_id": _logical_id(workflow),
            "version": workflow.version,
            "catalog_id": origin.catalog.id,
            "catalog_namespace": origin.catalog.namespace,
            "catalog_version": origin.catalog.version,
            "catalog_trust": origin.catalog.trust,
            "source_revision": _git_identity(origin.catalog.root).get("commit"),
            "image": workflow.image,
            "lineage": _locked_lineage(workflow_id, strict_files),
        },
        "resolution": {
            "reference": reference,
            "application": application,
            "policy": resolution["policy"],
            "reason": resolution["reason"],
            "selected_source_ref": resolution["selected"]["source_ref"],
        },
        "parameters": validation.normalized_params,
        "rendered": {
            "argv": portable_argv,
            "shell_preview": shlex.join(portable_argv),
            "env_overrides": rendered.env_overrides,
            "output_root": rendered.output_root,
        },
        "provenance": {
            "search": [repository.id for repository in provenance.search],
            "write": provenance.write.id if provenance.write is not None else None,
        },
    }
    if origin.catalog.manifest_path is not None:
        manifest["recipe"]["catalog_manifest_sha256"] = _sha256_file(
            origin.catalog.manifest_path
        )
    manifest["digest"] = _canonical_digest(manifest)
    return manifest


def _source_absolute_path(workflow_id: str, declared_path: str) -> Path:
    origin = get_workflow_origin(workflow_id)
    if origin is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")
    path = Path(declared_path)
    return path.resolve() if path.is_absolute() else (origin.catalog.root / path).resolve()


def create_recipe_bundle(
    target: Path,
    reference: str,
    params: dict[str, Any],
    *,
    run_id: str,
    application: Optional[str] = "workbench",
    policy: Optional[RecipeResolutionPolicy] = None,
    selected_workflow_id: Optional[str] = None,
    image: Optional[str] = None,
    resources: Optional[dict[str, Any]] = None,
    check_paths: bool = False,
) -> dict[str, Any]:
    """Atomically create a portable, immutable-input execution bundle."""
    if not run_id.strip():
        raise ValueError("run_id must be non-empty")
    target = target.resolve()
    if target.exists():
        raise FileExistsError(f"Bundle target already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)

    recipe_lock = build_recipe_lock(
        reference,
        params,
        application=application,
        policy=policy,
        selected_workflow_id=selected_workflow_id,
        authenticated=True,
        check_paths=check_paths,
        strict_files=True,
    )
    workflow_id = recipe_lock["recipe"]["workflow_id"]
    workflow = get_workflow_config(workflow_id)
    schema = get_workflow_schema(workflow_id)
    origin = get_workflow_origin(workflow_id)
    if workflow is None or schema is None or origin is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")
    resolved_image = image or workflow.image
    if not resolved_image:
        raise ValueError(
            f"Recipe '{workflow_id}' has no execution image; provide --image or "
            "declare image in its catalog entry"
        )

    temp_root = Path(
        tempfile.mkdtemp(prefix=f".{target.name}.tmp-", dir=str(target.parent))
    )
    copied: list[dict[str, Any]] = []
    external_dependencies: list[dict[str, Any]] = []
    try:
        for lineage_entry in recipe_lock["recipe"]["lineage"]:
            lineage_id = lineage_entry["workflow_id"]
            lineage_origin = get_workflow_origin(lineage_id)
            if lineage_origin is None:
                raise ValueError(f"Unknown workflow: {lineage_id}")
            for source in lineage_entry["sources"]:
                if source["language"] == "native":
                    external_dependencies.append(
                        {
                            "workflow_id": lineage_id,
                            "path": source["declared_path"],
                            "sha256": source.get("sha256"),
                        }
                    )
                    continue
                source_path = _source_absolute_path(
                    lineage_id, source["declared_path"]
                )
                relative = Path(source["path"])
                if relative.is_absolute() or ".." in relative.parts:
                    relative = Path(source_path.name)
                destination = (
                    temp_root
                    / "catalogs"
                    / lineage_origin.catalog.id
                    / relative
                )
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_path, destination)
                copied.append(
                    {
                        "workflow_id": lineage_id,
                        "path": destination.relative_to(temp_root).as_posix(),
                        "sha256": source["sha256"],
                    }
                )

        entry_language = next(
            (
                script.language
                for script in schema.scripts
                if script.path == schema.entry_script
            ),
            "bash",
        )
        command_argv = list(recipe_lock["rendered"]["argv"])
        if entry_language != "native":
            entry_relative = Path(schema.entry_script)
            if entry_relative.is_absolute() or ".." in entry_relative.parts:
                entry_relative = Path(entry_relative.name)
            command_argv[0] = (
                Path("catalogs") / origin.catalog.id / entry_relative
            ).as_posix()

        resolved_resources = {"cores": 1, "mem_mb": 2048, "gpus": 0}
        if resources:
            resolved_resources.update(resources)
        for field in ("cores", "mem_mb", "gpus"):
            value = resolved_resources.get(field)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(
                    f"Bundle resource '{field}' must be a non-negative integer"
                )
        node = {
            "id": "1",
            "title": schema.title,
            "description": schema.summary,
            "image_name": resolved_image,
            "launch": {
                "command": [shlex.join(command_argv)],
                "environment": recipe_lock["rendered"]["env_overrides"],
            },
            "inputs": {},
            "outputs": {},
            "resources": resolved_resources,
        }
        resolved_workflow = {
            "schema": RESOLVED_WORKFLOW_SCHEMA,
            "resolved_workflow": {
                "run_id": run_id,
                "recipe_lock_digest": recipe_lock["digest"],
                "nodes": {"1": node},
                "links": [],
            },
        }
        provenance = get_config().provenance
        provenance_template = {
            "schema": "biodepot.provenance_run/v1",
            "run_id": run_id,
            "recipe_lock_digest": recipe_lock["digest"],
            "search_repositories": [item.id for item in provenance.search],
            "write_repository": (
                {
                    "id": provenance.write.id,
                    "root": str(provenance.write.root),
                }
                if provenance.write is not None
                else None
            ),
            "status": "planned",
        }

        files = {
            "recipe.lock.json": recipe_lock,
            "manifests/resolved_workflow_v1.json": resolved_workflow,
            "commands/rendered-command.json": {
                "argv": command_argv,
                "shell_preview": shlex.join(command_argv),
                "env_overrides": recipe_lock["rendered"]["env_overrides"],
            },
            "provenance/run.template.json": provenance_template,
            "bundle.json": {
                "schema": "biodepot.recipe_bundle/v1",
                "recipe_lock_digest": recipe_lock["digest"],
                "run_id": run_id,
                "copied_sources": copied,
                "external_dependencies": external_dependencies,
            },
        }
        for relative, payload in files.items():
            output_path = temp_root / relative
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(
                json.dumps(payload, sort_keys=True, indent=2) + "\n",
                encoding="utf-8",
            )

        os.replace(temp_root, target)
    except Exception:
        shutil.rmtree(temp_root, ignore_errors=True)
        raise

    return {
        "bundle": str(target),
        "run_id": run_id,
        "workflow_id": workflow_id,
        "recipe_lock_digest": recipe_lock["digest"],
        "resolved_workflow": str(target / "manifests/resolved_workflow_v1.json"),
    }
