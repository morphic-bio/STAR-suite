"""Workflow schema scaffolding and draft validation tools.

These tools propose schemas — they never modify the running config.
The draft must be committed to the repo and loaded via reload_config
to become active.
"""

import re
from pathlib import Path
from typing import Any, Optional

import yaml
from pydantic import ValidationError

from ..config import get_config
from ..schemas.responses import (
    ScaffoldWorkflowResponse,
    ValidateDraftWorkflowResponse,
)
from ..schemas.workflow import WorkflowSchema
from .utils import is_path_allowed


# ---------------------------------------------------------------------------
# Shell script parser helpers
# ---------------------------------------------------------------------------

# Matches case arms: --flag-name) body ;;
# Also handles -h|--help) style combined short+long
# DOTALL lets .*? span multi-line case bodies (body on separate lines from ;;)
_CASE_ARM_RE = re.compile(
    r"^\s*(?:-\w[\s,|]+)?(--[\w][\w-]*)\)\s*(.*?);;",
    re.MULTILINE | re.DOTALL,
)

# In a case body: VAR="$2" or VAR=${2} or VAR=$2
_VALUE_BODY_RE = re.compile(r"(\w+)=[\"']?\$\{?2\}?[\"']?")

# In a case body: VAR=1 or VAR=0 or VAR=true or VAR=false
# \b ensures we don't partially match multi-digit numbers (e.g. VAR=10)
_BOOL_BODY_RE = re.compile(r"(\w+)=(0|1|true|false)\b")

# Default with env fallback: VAR="${ENV_VAR:-default}"
_DEFAULT_ENV_FALLBACK_RE = re.compile(
    r'^(\w+)="\$\{\w+:-(.*?)\}"', re.MULTILINE
)

# Default with simple quoted string: VAR="literal"
_DEFAULT_QUOTED_RE = re.compile(
    r'^(\w+)="([^$"][^"]*)"', re.MULTILINE
)

# Default with simple number: VAR=0 or VAR=24
_DEFAULT_NUMBER_RE = re.compile(r"^(\w+)=(\d+)\s*$", re.MULTILINE)

# Default with empty string: VAR=""
_DEFAULT_EMPTY_RE = re.compile(r'^(\w+)=""\s*$', re.MULTILINE)

# Exported variables: export VAR
_EXPORT_RE = re.compile(r"^\s*export\s+(\w+)", re.MULTILINE)

# Inside usage block: a line starting with --flag, optional arg hint, then description.
# The hint is a single short uppercase token (N, CSV, DIR, FILE, PATH, ID, P).
_USAGE_FLAG_LINE_RE = re.compile(
    r"^\s+(--[\w-]+)\s+(.*?)\s*$", re.MULTILINE
)

# Matches a leading arg-hint token: 1-5 all-uppercase letters
_ARG_HINT_RE = re.compile(r"^([A-Z]{1,5})\s+")


def _extract_usage_block(text: str) -> str:
    """Extract the heredoc body from a usage() function.

    Handles: usage() { cat <<EOF ... EOF }
             usage() { cat <<USAGE ... USAGE }
    """
    m = re.search(
        r"usage\s*\(\)\s*\{[^}]*?cat\s+<<\s*[\"']?(\w+)[\"']?\s*\n(.*?)\n\s*\1",
        text,
        re.DOTALL,
    )
    if m:
        return m.group(2)
    return ""


def _extract_defaults(text: str) -> dict[str, str]:
    """Extract variable defaults from script initialization section."""
    defaults: dict[str, str] = {}

    # Env fallback: VAR="${ENV:-default}"
    for m in _DEFAULT_ENV_FALLBACK_RE.finditer(text):
        defaults[m.group(1)] = m.group(2)

    # Quoted literal: VAR="value" (skip lines that reference $)
    for m in _DEFAULT_QUOTED_RE.finditer(text):
        var = m.group(1)
        if var not in defaults:
            defaults[var] = m.group(2)

    # Numeric: VAR=24
    for m in _DEFAULT_NUMBER_RE.finditer(text):
        var = m.group(1)
        if var not in defaults:
            defaults[var] = m.group(2)

    # Empty string: VAR=""
    for m in _DEFAULT_EMPTY_RE.finditer(text):
        var = m.group(1)
        if var not in defaults:
            defaults[var] = ""

    return defaults


def _extract_usage_descriptions(text: str) -> dict[str, str]:
    """Extract flag descriptions from the usage() heredoc only.

    Restricts search to the usage() function body to avoid matching
    flag-like lines elsewhere in the script (e.g., ARGS arrays).
    Strips a single leading uppercase arg-hint token (N, CSV, DIR, etc.)
    to handle lines where the flag+hint consume the alignment space.
    """
    usage_text = _extract_usage_block(text)
    if not usage_text:
        # Fallback: no usage block found — return empty
        return {}

    descs: dict[str, str] = {}
    for m in _USAGE_FLAG_LINE_RE.finditer(usage_text):
        flag = m.group(1)
        rest = m.group(2).strip()

        # Skip -h/--help
        if flag == "--help":
            continue

        # Strip a single leading uppercase arg-hint token (CSV, N, DIR, etc.)
        hint_m = _ARG_HINT_RE.match(rest)
        if hint_m:
            rest = rest[hint_m.end():].strip()

        if rest:
            descs[flag] = rest

    return descs


def _extract_exports(text: str) -> set[str]:
    """Extract exported variable names."""
    return {m.group(1) for m in _EXPORT_RE.finditer(text)}


def _flag_to_param_name(flag: str) -> str:
    """Convert --flag-name to flag_name."""
    return flag.lstrip("-").replace("-", "_")


def _clean_default(raw: Optional[str]) -> Any:
    """Clean a raw default value extracted from a shell script.

    Handles shell variable references and command substitutions:
    - "${REPO_ROOT}/path" → "path" (strip ${REPO_ROOT}/ prefix, keep relative)
    - "${VAR}" (pure variable ref) → None
    - "$(date ...)" (command substitution) → None
    - Literals and empty strings pass through unchanged.

    Returns the cleaned default, or None if it can't be resolved.
    """
    if raw is None:
        return None

    # Pure variable reference: ${VAR} or $VAR
    if re.match(r"^\$\{?\w+\}?$", raw):
        return None

    # Contains command substitution: $(...)
    if "$(" in raw:
        return None

    # Starts with ${REPO_ROOT}/ — strip to relative path
    m = re.match(r"^\$\{(?:REPO_ROOT|SCRIPT_DIR)\}/(.*)", raw)
    if m:
        return m.group(1)

    # Contains other ${VAR} references embedded in the value
    if "${" in raw:
        return None

    return raw


def _infer_type(param: dict[str, Any]) -> None:
    """Refine parameter type based on default value and name heuristics."""
    if param["type"] == "bool":
        return

    default = param.get("default")
    name = param["name"]

    # Integer if default is numeric
    if default is not None and default != "":
        try:
            int(default)
            param["type"] = "int"
            param["default"] = int(default)
            return
        except (ValueError, TypeError):
            pass

    # Path heuristics from the name
    dir_hints = ("root", "dir", "directory", "genome_dir", "out_root")
    file_hints = ("ref", "whitelist", "bin", "csv", "file", "path")

    if any(hint in name for hint in dir_hints):
        param["type"] = "directory"
        if default and default.startswith("/"):
            param["path_must_exist"] = True
        return

    if any(hint in name for hint in file_hints):
        param["type"] = "file"
        if default and default.startswith("/"):
            param["path_must_exist"] = True
        if "bin" in name:
            param["must_be_executable"] = True
        return

    # Comma-separated list hint
    if default is not None and "," in str(default):
        param["type"] = "string_list"
        return

    # CSV in the name or usage hint suggests string_list
    if "csv" in name or param.get("_usage_hint", "").upper() in ("CSV", "LIST"):
        param["type"] = "string_list"
        return


def _extract_parameters(text: str) -> list[dict[str, Any]]:
    """Extract parameters from case arms and correlate with defaults."""
    defaults = _extract_defaults(text)
    descs = _extract_usage_descriptions(text)
    exports = _extract_exports(text)

    params: list[dict[str, Any]] = []
    seen_flags: set[str] = set()

    for m in _CASE_ARM_RE.finditer(text):
        flag = m.group(1)
        body = m.group(2)

        # Skip help and duplicates
        if flag in ("--help",) or flag in seen_flags:
            continue
        seen_flags.add(flag)

        param_name = _flag_to_param_name(flag)

        # Detect value vs bool from body
        val_m = _VALUE_BODY_RE.search(body)
        bool_m = _BOOL_BODY_RE.search(body)

        if val_m:
            var_name = val_m.group(1)
            raw_default = defaults.get(var_name)
            default_val = _clean_default(raw_default)
            param = {
                "name": param_name,
                "cli_flag": flag,
                "type": "string",
                "required": False,
                "default": default_val,
                "description": descs.get(flag, "TODO"),
                "_var": var_name,
            }
            # Check if this var is exported
            if var_name in exports:
                param["env_var"] = var_name
        elif bool_m:
            var_name = bool_m.group(1)
            param = {
                "name": param_name,
                "cli_flag": flag,
                "type": "bool",
                "required": False,
                "default": False,
                "description": descs.get(flag, "TODO"),
                "_var": var_name,
            }
        else:
            # Unrecognized body pattern — include as string with a note
            param = {
                "name": param_name,
                "cli_flag": flag,
                "type": "string",
                "required": False,
                "default": None,
                "description": descs.get(flag, "TODO (unrecognized pattern)"),
            }

        params.append(param)

    # Infer types and clean up internal fields
    for p in params:
        _infer_type(p)
        p.pop("_var", None)

    return params


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def scaffold_workflow_schema(script_path: str) -> ScaffoldWorkflowResponse:
    """Parse a shell script and return a draft workflow schema YAML.

    Best-effort extraction from while/case flag-parsing patterns.
    The output is a starting point for manual refinement.

    Args:
        script_path: Path to the shell script (absolute or relative to repo root).

    Returns:
        ScaffoldWorkflowResponse with draft YAML and extraction notes.
    """
    config = get_config()
    path = Path(script_path)

    # Resolve relative paths
    if not path.is_absolute():
        path = Path(config.paths.repo_root) / path

    # Security check
    if not is_path_allowed(path):
        raise ValueError(f"Script path outside trusted roots: {path}")

    if not path.exists():
        raise FileNotFoundError(f"Script not found: {path}")

    if not path.is_file():
        raise ValueError(f"Not a file: {path}")

    text = path.read_text()
    params = _extract_parameters(text)

    # Compute relative entry_script path
    repo_root = Path(config.paths.repo_root)
    try:
        entry_script = str(path.relative_to(repo_root))
    except ValueError:
        entry_script = str(path)

    # Derive a workflow ID from the script name
    stem = path.stem
    # Strip common prefixes
    for prefix in ("run_", "run-"):
        if stem.startswith(prefix):
            stem = stem[len(prefix) :]
    workflow_id = stem.replace("-", "_")

    # Build the draft schema dict
    schema: dict[str, Any] = {
        "id": workflow_id,
        "title": f"TODO: Human-readable title for {path.name}",
        "summary": "TODO: one-line summary",
        "kind": "shell_workflow",
        "entry_script": entry_script,
        "supported_modes": ["local"],
        "parameters": _clean_params_for_yaml(params),
        "rendering": {
            "flag_order": [p["name"] for p in params],
        },
    }

    draft_yaml = yaml.dump(schema, default_flow_style=False, sort_keys=False)

    notes = [
        "This is a best-effort extraction. Review all types, defaults, and descriptions.",
        "Parameters with 'TODO' descriptions need manual documentation.",
        "Add path_must_exist: true to file/directory params that must exist at runtime.",
        "Add stages, parameter_groups, and constraints sections as needed.",
        "See mcp_server/workflows/AUTHORING.md for the full field reference.",
    ]

    # Flag params that couldn't get descriptions
    todo_count = sum(1 for p in params if "TODO" in p.get("description", ""))
    if todo_count:
        notes.append(
            f"{todo_count} parameter(s) have TODO descriptions — "
            f"check the script's usage() text or add manually."
        )

    return ScaffoldWorkflowResponse(
        draft_yaml=draft_yaml,
        parameter_count=len(params),
        source_script=str(path),
        notes=notes,
    )


def _clean_params_for_yaml(params: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Prepare parameter dicts for YAML output — remove None values."""
    cleaned = []
    for p in params:
        out: dict[str, Any] = {}
        for key in (
            "name",
            "cli_flag",
            "type",
            "required",
            "default",
            "description",
            "path_must_exist",
            "must_be_executable",
            "env_var",
        ):
            if key in p and p[key] is not None:
                out[key] = p[key]
        cleaned.append(out)
    return cleaned


def validate_draft_workflow_schema(
    draft_yaml: str,
) -> ValidateDraftWorkflowResponse:
    """Validate a draft YAML string against the WorkflowSchema model.

    Does NOT load the schema into the running server config.

    Args:
        draft_yaml: YAML string to validate.

    Returns:
        ValidateDraftWorkflowResponse with validation results.
    """
    errors: list[str] = []
    warnings: list[str] = []
    parameter_count: Optional[int] = None

    # Step 1: Parse YAML
    try:
        data = yaml.safe_load(draft_yaml)
    except yaml.YAMLError as e:
        return ValidateDraftWorkflowResponse(
            valid=False,
            errors=[f"YAML parse error: {e}"],
        )

    if not isinstance(data, dict):
        return ValidateDraftWorkflowResponse(
            valid=False,
            errors=["YAML must be a mapping (dict), not a scalar or list"],
        )

    # Step 2: Validate against Pydantic model
    try:
        schema = WorkflowSchema(**data)
    except ValidationError as e:
        for err in e.errors():
            loc = " -> ".join(str(x) for x in err["loc"])
            errors.append(f"{loc}: {err['msg']}")
        return ValidateDraftWorkflowResponse(
            valid=False,
            errors=errors,
        )

    parameter_count = len(schema.parameters)

    # Step 3: Semantic checks (warnings, not hard errors)

    # Check for TODO placeholders
    if "TODO" in schema.title:
        warnings.append("title contains TODO placeholder")
    if "TODO" in schema.summary:
        warnings.append("summary contains TODO placeholder")
    for p in schema.parameters:
        if "TODO" in p.description:
            warnings.append(f"parameter '{p.name}' description contains TODO")

    # Check that all params in flag_order exist
    param_names = {p.name for p in schema.parameters}
    for flag_name in schema.rendering.flag_order:
        if flag_name not in param_names:
            errors.append(
                f"rendering.flag_order references unknown parameter: {flag_name}"
            )

    # Check that all params in groups exist
    for group in schema.parameter_groups:
        for gp_name in group.parameters:
            if gp_name not in param_names:
                errors.append(
                    f"parameter_group '{group.name}' references unknown parameter: {gp_name}"
                )

    # Check that all params in constraints exist
    for constraint in schema.constraints:
        for cp_name in constraint.params:
            if cp_name not in param_names:
                errors.append(
                    f"constraint ({constraint.kind}) references unknown parameter: {cp_name}"
                )

    # Check for duplicate parameter names
    seen_names: set[str] = set()
    for p in schema.parameters:
        if p.name in seen_names:
            errors.append(f"duplicate parameter name: {p.name}")
        seen_names.add(p.name)

    # Check for duplicate cli_flags
    seen_flags: set[str] = set()
    for p in schema.parameters:
        if p.cli_flag in seen_flags:
            errors.append(f"duplicate cli_flag: {p.cli_flag}")
        seen_flags.add(p.cli_flag)

    # Check that parameters with path_must_exist are file or directory type
    for p in schema.parameters:
        if p.path_must_exist and p.type not in ("file", "directory"):
            warnings.append(
                f"parameter '{p.name}' has path_must_exist but type is '{p.type}' "
                f"(expected 'file' or 'directory')"
            )

    # Check that must_be_executable is only on file params
    for p in schema.parameters:
        if p.must_be_executable and p.type != "file":
            warnings.append(
                f"parameter '{p.name}' has must_be_executable but type is '{p.type}' "
                f"(expected 'file')"
            )

    # Check that enum params have choices
    for p in schema.parameters:
        if p.type == "enum" and not p.choices:
            errors.append(
                f"parameter '{p.name}' is type 'enum' but has no choices defined"
            )

    # Warn if no flag_order is set but there are parameters
    if schema.parameters and not schema.rendering.flag_order:
        warnings.append(
            "No rendering.flag_order defined — flag emission order will be "
            "non-deterministic. Add flag_order for reproducible commands."
        )

    return ValidateDraftWorkflowResponse(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        parameter_count=parameter_count,
    )
