"""Workflow discovery, validation, and rendering tools."""

import os
import shlex
from pathlib import Path
from typing import Any, Optional

from ..config import get_config, get_workflow_schema, get_workflow_schemas
from ..schemas.config import WorkflowConfig
from ..schemas.responses import (
    ConstraintInfo,
    DescribeWorkflowResponse,
    GetWorkflowScriptsResponse,
    ListWorkflowsResponse,
    ParameterGroupInfo,
    ParameterInfo,
    RenderWorkflowResponse,
    RequiredFileParamInfo,
    ValidateWorkflowResponse,
    WorkflowInfo,
    WorkflowParameterSchemaResponse,
    WorkflowScriptDetail,
    WorkflowStageInfo,
)
from ..schemas.workflow import WorkflowParameterDef, WorkflowSchema
from .utils import is_path_allowed, validate_path


def _get_workflow_config(workflow_id: str) -> Optional[WorkflowConfig]:
    """Look up the WorkflowConfig entry for *workflow_id*."""
    config = get_config()
    return config.get_workflow(workflow_id)


def _is_visible(workflow_id: str, authenticated: bool) -> bool:
    """Return True if *workflow_id* should be exposed given the auth state."""
    wf_cfg = _get_workflow_config(workflow_id)
    if wf_cfg is None:
        return True
    if wf_cfg.visibility == "private" and not authenticated:
        return False
    return True


def list_workflows(authenticated: bool = False) -> ListWorkflowsResponse:
    """Return configured workflow templates visible to the caller."""
    schemas = get_workflow_schemas()
    items = []
    for schema in schemas.values():
        if not _is_visible(schema.id, authenticated):
            continue
        items.append(
            WorkflowInfo(
                id=schema.id,
                title=schema.title,
                summary=schema.summary,
                kind=schema.kind,
                entry_script=schema.entry_script,
                supported_modes=schema.supported_modes,
            )
        )
    return ListWorkflowsResponse(workflows=items)


def describe_workflow(
    workflow_id: str,
    authenticated: bool = False,
) -> DescribeWorkflowResponse:
    """Return full metadata for a workflow."""
    schema = get_workflow_schema(workflow_id)
    if schema is None or not _is_visible(workflow_id, authenticated):
        raise ValueError(f"Unknown workflow: {workflow_id}")

    stages = [
        WorkflowStageInfo(
            name=s.name,
            title=s.title,
            script=s.script,
            description=s.description,
            gated_by=s.gated_by,
        )
        for s in schema.stages
    ]

    groups = [
        ParameterGroupInfo(
            name=g.name,
            title=g.title,
            parameters=g.parameters,
            gated_by=g.gated_by,
        )
        for g in schema.parameter_groups
    ]

    return DescribeWorkflowResponse(
        id=schema.id,
        title=schema.title,
        summary=schema.summary,
        kind=schema.kind,
        entry_script=schema.entry_script,
        supported_modes=schema.supported_modes,
        caveats=schema.caveats,
        default_output_layout=schema.default_output_layout,
        stages=stages,
        parameter_groups=groups,
    )


def get_workflow_scripts(
    workflow_id: str,
    include_content: bool = False,
    authenticated: bool = False,
) -> GetWorkflowScriptsResponse:
    """Return the scripts composing a workflow with provenance metadata.

    When *authenticated* is False (public / unauthenticated callers),
    host-specific details are redacted:

    * ``absolute_path`` is omitted from each script entry.
    * ``repo_root``, ``git_commit``, and ``git_remote`` are omitted from
      provenance.

    Structural information (role, relative path, description, language,
    exists) is always returned so downstream agents can reason about
    the workflow's script composition without gaining host metadata.
    """
    schema = get_workflow_schema(workflow_id)
    if schema is None or not _is_visible(workflow_id, authenticated):
        raise ValueError(f"Unknown workflow: {workflow_id}")

    config = get_config()
    repo_root = Path(config.paths.repo_root)

    scripts: list[WorkflowScriptDetail] = []
    for s in schema.scripts:
        abs_path = Path(s.path)
        if not abs_path.is_absolute():
            abs_path = repo_root / s.path

        detail = WorkflowScriptDetail(
            role=s.role,
            path=s.path,
            absolute_path=str(abs_path) if authenticated else None,
            description=s.description,
            language=s.language,
            exists=abs_path.exists(),
        )
        scripts.append(detail)

    # If no scripts section in the schema, at minimum return the entry script
    if not scripts:
        entry_abs = repo_root / schema.entry_script
        scripts.append(
            WorkflowScriptDetail(
                role="entry",
                path=schema.entry_script,
                absolute_path=str(entry_abs) if authenticated else None,
                description="Entry script",
                language="bash",
                exists=entry_abs.exists(),
            )
        )

    # Build provenance — public callers get only the schema file reference;
    # authenticated callers additionally get repo_root, git_commit, git_remote.
    provenance: dict[str, Any] = {
        "workflow_schema": _get_workflow_config(workflow_id).schema_file
        if _get_workflow_config(workflow_id)
        else None,
    }

    if authenticated:
        provenance["repo_root"] = str(repo_root)

        # Best-effort git info (authenticated only)
        import subprocess

        try:
            commit = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                cwd=str(repo_root),
                timeout=5,
            )
            if commit.returncode == 0:
                provenance["git_commit"] = commit.stdout.strip()
        except Exception:
            pass

        try:
            remote = subprocess.run(
                ["git", "remote", "get-url", "origin"],
                capture_output=True,
                text=True,
                cwd=str(repo_root),
                timeout=5,
            )
            if remote.returncode == 0:
                provenance["git_remote"] = remote.stdout.strip()
        except Exception:
            pass

    return GetWorkflowScriptsResponse(
        workflow_id=schema.id,
        title=schema.title,
        entry_script=schema.entry_script,
        scripts=scripts,
        provenance=provenance,
    )


def _required_file_params(schema: WorkflowSchema) -> list[RequiredFileParamInfo]:
    """Derive a required-file checklist from schema metadata.

    This is intended for UX (Launchpad / agents): parameters where path existence
    matters when check_paths=True.
    """
    out: list[RequiredFileParamInfo] = []
    for p in schema.parameters:
        if p.type in ("file", "directory") and p.path_must_exist:
            out.append(
                RequiredFileParamInfo(
                    name=p.name,
                    cli_flag=p.cli_flag,
                    type=p.type,
                    description=p.description,
                )
            )
    return out


def get_workflow_parameter_schema(
    workflow_id: str,
    authenticated: bool = False,
) -> WorkflowParameterSchemaResponse:
    """Return the machine-readable parameter schema for a workflow."""
    schema = get_workflow_schema(workflow_id)
    if schema is None or not _is_visible(workflow_id, authenticated):
        raise ValueError(f"Unknown workflow: {workflow_id}")

    params = [
        ParameterInfo(
            name=p.name,
            cli_flag=p.cli_flag,
            type=p.type,
            required=p.required,
            default=p.default,
            description=p.description,
            choices=p.choices,
            repeatable=p.repeatable,
            operand_group=p.operand_group,
            path_must_exist=p.path_must_exist,
            must_be_executable=p.must_be_executable,
            category=p.category,
            stage=p.stage,
            source=p.source,
            env_var=p.env_var,
            skip_when_default=p.skip_when_default,
            is_output_root=p.is_output_root,
            ui_gated_by=p.ui_gated_by,
        )
        for p in schema.parameters
    ]

    groups = [
        ParameterGroupInfo(
            name=g.name,
            title=g.title,
            parameters=g.parameters,
            gated_by=g.gated_by,
        )
        for g in schema.parameter_groups
    ]

    constraints = [
        ConstraintInfo(kind=c.kind, params=c.params, message=c.message)
        for c in schema.constraints
    ]

    return WorkflowParameterSchemaResponse(
        workflow_id=schema.id,
        parameters=params,
        parameter_groups=groups,
        constraints=constraints,
        required_files=_required_file_params(schema),
    )


def _normalize_params(
    schema: WorkflowSchema, params: dict[str, Any]
) -> dict[str, Any]:
    """Normalize parameters: apply defaults and coerce types."""
    normalized: dict[str, Any] = {}

    for p_def in schema.parameters:
        if p_def.name in params:
            val = params[p_def.name]
            # Coerce types (skip on failure — validation catches type errors)
            if p_def.type == "int" and val is not None:
                try:
                    val = int(val)
                except (TypeError, ValueError):
                    pass
            elif p_def.type == "float" and val is not None:
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    pass
            elif p_def.type == "bool":
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                else:
                    val = bool(val)
            normalized[p_def.name] = val
        elif p_def.default is not None:
            normalized[p_def.name] = p_def.default
        # Omit params with no value and no default

    return normalized


def validate_workflow_parameters(
    workflow_id: str,
    params: dict[str, Any],
    check_paths: bool = True,
) -> ValidateWorkflowResponse:
    """Validate a structured params dict against a workflow schema."""
    schema = get_workflow_schema(workflow_id)
    if schema is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")

    errors: list[str] = []
    warnings: list[str] = []

    # Check for unknown params
    known_names = {p.name for p in schema.parameters}
    for key in params:
        if key not in known_names:
            errors.append(f"Unknown parameter: {key}")

    # Check required params
    for p_def in schema.parameters:
        if p_def.required and p_def.name not in params:
            errors.append(f"Missing required parameter: {p_def.name}")

    # Type checks
    for p_def in schema.parameters:
        if p_def.name not in params:
            continue
        val = params[p_def.name]
        if val is None:
            continue

        if p_def.type == "int":
            try:
                int(val)
            except (TypeError, ValueError):
                errors.append(f"Parameter '{p_def.name}' must be an integer, got: {val!r}")
        elif p_def.type == "float":
            try:
                float(val)
            except (TypeError, ValueError):
                errors.append(f"Parameter '{p_def.name}' must be a float, got: {val!r}")
        elif p_def.type == "bool":
            if not isinstance(val, (bool, int, str)):
                errors.append(f"Parameter '{p_def.name}' must be a boolean, got: {val!r}")
        elif p_def.type == "enum":
            if p_def.choices and str(val) not in p_def.choices:
                errors.append(
                    f"Parameter '{p_def.name}' must be one of {p_def.choices}, got: {val!r}"
                )

    # Normalize for constraint checking
    normalized = _normalize_params(schema, params)

    # Constraint validation
    for constraint in schema.constraints:
        if constraint.kind == "mutual_exclusion":
            present = [
                p for p in constraint.params
                if p in params and params[p] is not None and params[p] is not False
                and params[p] != "" and params[p] != 0
            ]
            if len(present) > 1:
                errors.append(
                    f"Mutually exclusive parameters both provided: "
                    f"{', '.join(present)}. {constraint.message}"
                )

        elif constraint.kind == "group_required":
            present = [
                p for p in constraint.params
                if p in params and params[p] is not None
                and params[p] != ""
            ]
            if 0 < len(present) < len(constraint.params):
                missing = [p for p in constraint.params if p not in present]
                warnings.append(
                    f"Incomplete parameter group: provided {present}, "
                    f"missing {missing}. {constraint.message}"
                )

        elif constraint.kind == "dependency":
            # For the cellbender_gpu / star_only dependency:
            # warn if cellbender_gpu is set but star_only is also set
            if len(constraint.params) >= 2:
                dep_param = constraint.params[0]
                gate_param = constraint.params[1]
                dep_val = normalized.get(dep_param)
                gate_val = normalized.get(gate_param)
                if dep_val and gate_val:
                    warnings.append(
                        f"Parameter '{dep_param}' has no effect when "
                        f"'{gate_param}' is set. {constraint.message}"
                    )

        elif constraint.kind == "positive":
            for p_name in constraint.params:
                val = normalized.get(p_name)
                if val is not None:
                    try:
                        if int(val) <= 0:
                            errors.append(
                                f"Parameter '{p_name}' must be positive, got: {val}"
                            )
                    except (TypeError, ValueError):
                        pass

        elif constraint.kind == "non_negative":
            for p_name in constraint.params:
                val = normalized.get(p_name)
                if val is not None:
                    try:
                        if int(val) < 0:
                            errors.append(
                                f"Parameter '{p_name}' must be non-negative, got: {val}"
                            )
                    except (TypeError, ValueError):
                        pass

    # Path validation (file/directory params with path_must_exist)
    if check_paths:
        config = get_config()
        repo_root = Path(config.paths.repo_root)
        for p_def in schema.parameters:
            if p_def.name not in normalized:
                continue
            val = normalized[p_def.name]
            if val is None or val == "":
                continue

            if p_def.type in ("file", "directory") and p_def.path_must_exist:
                path = Path(str(val))
                # Resolve relative paths against repo root
                if not path.is_absolute():
                    path = repo_root / path
                    normalized[p_def.name] = str(path)
                if not is_path_allowed(path):
                    errors.append(
                        f"Parameter '{p_def.name}' path outside trusted roots: {path}"
                    )
                elif not path.exists():
                    errors.append(
                        f"Parameter '{p_def.name}' path does not exist: {path}"
                    )
                elif p_def.type == "file" and not path.is_file():
                    errors.append(
                        f"Parameter '{p_def.name}' is not a file: {path}"
                    )
                elif p_def.type == "directory" and not path.is_dir():
                    errors.append(
                        f"Parameter '{p_def.name}' is not a directory: {path}"
                    )
                elif p_def.must_be_executable and not os.access(path, os.X_OK):
                    errors.append(
                        f"Parameter '{p_def.name}' is not executable: {path}"
                    )

    return ValidateWorkflowResponse(
        valid=len(errors) == 0,
        normalized_params=normalized,
        warnings=warnings,
        errors=errors,
    )


def render_workflow_command(
    workflow_id: str,
    params: dict[str, Any],
) -> RenderWorkflowResponse:
    """Render validated params into a command invocation.

    This does not execute anything.
    """
    schema = get_workflow_schema(workflow_id)
    if schema is None:
        raise ValueError(f"Unknown workflow: {workflow_id}")

    config = get_config()
    repo_root = Path(config.paths.repo_root)

    # Resolve entry script path
    entry_script = schema.entry_script
    entry_path = Path(entry_script)
    if not entry_path.is_absolute():
        entry_path = repo_root / entry_script

    # Normalize params
    normalized = _normalize_params(schema, params)

    # Build argv
    argv: list[str] = [str(entry_path)]
    env_overrides: dict[str, str] = {}
    output_root: Optional[str] = None

    # Build a name -> param_def lookup
    param_lookup = {p.name: p for p in schema.parameters}

    # Determine flag order
    flag_order = schema.rendering.flag_order
    if not flag_order:
        flag_order = [p.name for p in schema.parameters]

    # Read rendering configuration
    omit_absent = schema.rendering.omit_absent_optionals
    csv_style = schema.rendering.csv_style

    # Collect env_var overrides from schema (declarative, not hardcoded)
    for p_def in schema.parameters:
        if p_def.env_var and p_def.name in normalized and normalized[p_def.name] is not None:
            val = normalized[p_def.name]
            # Respect skip_when_default for env_var params too
            if p_def.skip_when_default and p_def.default is not None:
                try:
                    if p_def.type == "int" and int(val) == int(p_def.default):
                        continue
                    elif p_def.type == "float" and float(val) == float(p_def.default):
                        continue
                    elif p_def.type not in ("int", "float") and val == p_def.default:
                        continue
                except (TypeError, ValueError):
                    pass
            env_overrides[p_def.env_var] = str(val)

    def _operand_token(q: WorkflowParameterDef, q_name: str) -> Optional[str]:
        """Single argv token for STAR-style multi-operand flags (--readFilesIn a b)."""
        has_v = q_name in normalized
        v = normalized.get(q_name)
        if q.env_var and q.source != "workflow_wrapper":
            return None
        if omit_absent and not has_v and not q.required:
            return None
        if v is None:
            return None
        if q.skip_when_default and q.default is not None:
            try:
                if q.type == "int" and int(v) == int(q.default):
                    return None
                if q.type == "float" and float(v) == float(q.default):
                    return None
                if q.type not in ("int", "float") and v == q.default:
                    return None
            except (TypeError, ValueError):
                pass
        if q.type == "string_list":
            if not v:
                return None
            if isinstance(v, list):
                return ",".join(str(x) for x in v)
            return str(v)
        if q.type == "bool":
            return None
        if q.type == "string":
            if not v or str(v) == "":
                return None
            return str(v)
        if q.type in ("int", "float"):
            return str(v)
        if q.type in ("file", "directory"):
            return str(v)
        return str(v)

    n_flags = len(flag_order)
    fi = 0
    while fi < n_flags:
        p_name = flag_order[fi]
        if p_name not in param_lookup:
            fi += 1
            continue

        p_def = param_lookup[p_name]

        # Params delivered via env_var from sub-scripts (source != workflow_wrapper)
        # are already in env_overrides — do not emit as top-level CLI flags.
        if p_def.env_var and p_def.source != "workflow_wrapper":
            fi += 1
            continue

        if p_def.operand_group:
            og = p_def.operand_group
            cflag = p_def.cli_flag
            group_names: list[str] = []
            fj = fi
            while fj < n_flags:
                pn2 = flag_order[fj]
                q = param_lookup.get(pn2)
                if not q or q.operand_group != og or q.cli_flag != cflag:
                    break
                group_names.append(pn2)
                fj += 1

            pieces: list[str] = []
            skip_group = False
            for pn2 in group_names:
                q = param_lookup[pn2]
                val = normalized.get(pn2)
                if q.is_output_root and val:
                    output_root = str(val)
                tok = _operand_token(q, pn2)
                if tok is None:
                    skip_group = True
                    break
                pieces.append(tok)
            if not skip_group and pieces:
                argv.append(cflag)
                argv.extend(pieces)
            fi = fj
            continue

        has_value = p_name in normalized
        val = normalized.get(p_name)

        # When omit_absent_optionals is true, skip absent optional params
        if omit_absent and not has_value and not p_def.required:
            fi += 1
            continue

        # Skip None/empty values for non-required params
        if val is None:
            fi += 1
            continue

        # Track output root (declarative via is_output_root)
        if p_def.is_output_root and val:
            output_root = str(val)

        # skip_when_default: omit flag when value equals the schema default
        if p_def.skip_when_default and p_def.default is not None:
            try:
                if p_def.type == "int" and int(val) == int(p_def.default):
                    fi += 1
                    continue
                elif p_def.type == "float" and float(val) == float(p_def.default):
                    fi += 1
                    continue
                elif p_def.type not in ("int", "float") and val == p_def.default:
                    fi += 1
                    continue
            except (TypeError, ValueError):
                pass

        if p_def.type == "bool":
            if val:
                argv.append(p_def.cli_flag)
        elif p_def.type == "string_list":
            if val:
                argv.append(p_def.cli_flag)
                if isinstance(val, list):
                    joined = ",".join(str(v) for v in val)
                else:
                    joined = str(val)
                if csv_style == "quoted":
                    argv.append(joined)
                else:
                    argv.append(joined)
        elif p_def.type == "int":
            argv.append(p_def.cli_flag)
            argv.append(str(val))
        elif p_def.type == "string":
            if val and str(val) != "":
                argv.append(p_def.cli_flag)
                argv.append(str(val))
        else:
            # file, directory, float
            if val is not None:
                argv.append(p_def.cli_flag)
                argv.append(str(val))

        fi += 1

    # Build shell preview
    shell_preview = " ".join(shlex.quote(a) for a in argv)

    return RenderWorkflowResponse(
        workflow_id=schema.id,
        entry_script=str(entry_path),
        argv=argv,
        shell_preview=shell_preview,
        env_overrides=env_overrides,
        output_root=output_root,
    )
