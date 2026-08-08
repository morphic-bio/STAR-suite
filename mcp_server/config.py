"""Configuration loading and management for MCP server."""

import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import yaml
from packaging.version import InvalidVersion, Version

from .schemas.config import (
    MCPConfig,
    ProvenanceHierarchyConfig,
    RecipeCatalogManifest,
    WorkflowConfig,
)
from .schemas.workflow import WorkflowSchema


# Global config instance
_config: Optional[MCPConfig] = None
_config_path: Optional[Path] = None
_config_loaded_at: Optional[datetime] = None
_workflow_schemas: dict[str, WorkflowSchema] = {}
_workflow_configs: dict[str, WorkflowConfig] = {}
_workflow_origins: dict[str, "WorkflowOrigin"] = {}
_recipe_catalogs: tuple["LoadedRecipeCatalog", ...] = ()


BUILTIN_CATALOG_ID = "starsuite-builtin"
BUILTIN_CATALOG_NAMESPACE = "starsuite.core"


@dataclass(frozen=True)
class LoadedRecipeCatalog:
    """Resolved identity and filesystem root for one loaded recipe catalog."""

    id: str
    namespace: str
    version: str
    root: Path
    manifest_path: Optional[Path]
    built_in: bool = False
    order: int = 0
    trust: str = "untrusted"


@dataclass(frozen=True)
class WorkflowOrigin:
    """Catalog ownership and schema location for one workflow."""

    catalog: LoadedRecipeCatalog
    schema_file: str
    schema_path: Path


def expand_env_vars(value: Any) -> Any:
    """Recursively expand environment variables in config values.

    Supports ${VAR_NAME} syntax. Missing vars expand to empty string.
    """
    if isinstance(value, str):
        # Match ${VAR_NAME} pattern
        pattern = re.compile(r"\$\{([^}]+)\}")

        def replacer(match: re.Match) -> str:
            var_name = match.group(1)
            return os.environ.get(var_name, "")

        return pattern.sub(replacer, value)
    elif isinstance(value, dict):
        return {k: expand_env_vars(v) for k, v in value.items()}
    elif isinstance(value, list):
        return [expand_env_vars(item) for item in value]
    else:
        return value


def _resolved_path(base: Path, value: Path | str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _confined_catalog_path(root: Path, value: str, label: str) -> Path:
    raw_path = Path(value)
    if raw_path.is_absolute():
        raise ValueError(f"{label} must be relative to its catalog root: {value}")

    resolved = (root / raw_path).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{label} escapes its catalog root: {value}") from exc
    return resolved


def _read_yaml(path: Path, label: str) -> Any:
    with open(path) as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"{label} must contain a YAML object: {path}")
    return expand_env_vars(loaded)


def _load_workflow_schema(path: Path, workflow_id: str) -> WorkflowSchema:
    if not path.exists():
        raise FileNotFoundError(
            f"Workflow schema file not found for '{workflow_id}': {path}"
        )
    return WorkflowSchema(**_read_yaml(path, "workflow schema"))


def _validate_external_workflow_paths(
    schema: WorkflowSchema,
    workflow: WorkflowConfig,
    root: Path,
) -> None:
    if schema.id != workflow.id:
        raise ValueError(
            f"catalog workflow id '{workflow.id}' does not match schema id "
            f"'{schema.id}'"
        )
    if schema.entry_script != workflow.entry_script:
        raise ValueError(
            f"catalog workflow '{workflow.id}' entry_script does not match its "
            "workflow schema"
        )

    _confined_catalog_path(
        root, workflow.entry_script, f"catalog workflow '{workflow.id}' entry_script"
    )
    for script in schema.scripts:
        _confined_catalog_path(
            root,
            script.path,
            f"catalog workflow '{workflow.id}' script '{script.role}'",
        )
    for stage in schema.stages:
        if stage.script:
            _confined_catalog_path(
                root,
                stage.script,
                f"catalog workflow '{workflow.id}' stage '{stage.name}' script",
            )
    if schema.execution is not None:
        for stage in schema.execution.stages:
            _confined_catalog_path(
                root,
                stage.entry_script,
                (
                    f"catalog workflow '{workflow.id}' execution stage "
                    f"'{stage.name}' entry_script"
                ),
            )


def _load_workflow_schemas(
    config: MCPConfig,
    repo_root: Path,
    config_dir: Path,
) -> tuple[
    dict[str, WorkflowSchema],
    dict[str, WorkflowConfig],
    dict[str, WorkflowOrigin],
    tuple[LoadedRecipeCatalog, ...],
]:
    """Load the built-in compatibility catalog and ordered external catalogs."""

    schemas: dict[str, WorkflowSchema] = {}
    workflow_configs: dict[str, WorkflowConfig] = {}
    origins: dict[str, WorkflowOrigin] = {}

    builtin = LoadedRecipeCatalog(
        id=BUILTIN_CATALOG_ID,
        namespace=BUILTIN_CATALOG_NAMESPACE,
        version="configured",
        root=repo_root.resolve(),
        manifest_path=None,
        built_in=True,
        trust="trusted",
    )
    catalogs: list[LoadedRecipeCatalog] = [builtin]
    catalog_ids = {builtin.id}
    catalog_namespaces = {builtin.namespace}

    def register(
        workflow: WorkflowConfig,
        schema: WorkflowSchema,
        origin: WorkflowOrigin,
    ) -> None:
        if workflow.id in schemas:
            previous = origins[workflow.id].catalog.id
            raise ValueError(
                f"duplicate workflow id '{workflow.id}' in catalogs "
                f"'{previous}' and '{origin.catalog.id}'"
            )
        try:
            Version(workflow.version)
        except InvalidVersion as exc:
            raise ValueError(
                f"workflow '{workflow.id}' has invalid PEP 440 version: "
                f"{workflow.version}"
            ) from exc
        schemas[workflow.id] = schema
        workflow_configs[workflow.id] = workflow
        origins[workflow.id] = origin

    # Backward-compatible inline workflows retain existing absolute-path behavior.
    for workflow in config.workflows:
        schema_path = Path(workflow.schema_file)
        if not schema_path.is_absolute():
            schema_path = builtin.root / schema_path
        schema_path = schema_path.resolve()
        register(
            workflow,
            _load_workflow_schema(schema_path, workflow.id),
            WorkflowOrigin(
                catalog=builtin,
                schema_file=workflow.schema_file,
                schema_path=schema_path,
            ),
        )

    for source in config.recipe_catalogs:
        manifest_path = _resolved_path(config_dir, source.manifest)
        if not manifest_path.exists():
            raise FileNotFoundError(f"Recipe catalog manifest not found: {manifest_path}")

        manifest = RecipeCatalogManifest(**_read_yaml(manifest_path, "recipe catalog"))
        if manifest.id in catalog_ids:
            raise ValueError(f"duplicate recipe catalog id: {manifest.id}")
        if manifest.namespace in catalog_namespaces:
            raise ValueError(
                f"duplicate recipe catalog namespace: {manifest.namespace}"
            )

        catalog = LoadedRecipeCatalog(
            id=manifest.id,
            namespace=manifest.namespace,
            version=manifest.version,
            root=manifest_path.parent.resolve(),
            manifest_path=manifest_path,
            order=len(catalogs),
            trust=source.trust,
        )
        catalogs.append(catalog)
        catalog_ids.add(catalog.id)
        catalog_namespaces.add(catalog.namespace)

        workflow_prefix = f"{catalog.namespace}/"
        for workflow in manifest.workflows:
            if not workflow.id.startswith(workflow_prefix):
                raise ValueError(
                    f"catalog '{catalog.id}' workflow id '{workflow.id}' must use "
                    f"namespace prefix '{workflow_prefix}'"
                )
            schema_path = _confined_catalog_path(
                catalog.root,
                workflow.schema_file,
                f"catalog workflow '{workflow.id}' schema_file",
            )
            schema = _load_workflow_schema(schema_path, workflow.id)
            _validate_external_workflow_paths(schema, workflow, catalog.root)
            register(
                workflow,
                schema,
                WorkflowOrigin(
                    catalog=catalog,
                    schema_file=workflow.schema_file,
                    schema_path=schema_path,
                ),
            )

    _validate_workflow_relationships(workflow_configs)

    return schemas, workflow_configs, origins, tuple(catalogs)


def _validate_workflow_relationships(
    workflows: dict[str, WorkflowConfig],
) -> None:
    """Validate explicit lineage after the complete catalog stack is loaded."""
    for workflow in workflows.values():
        references = ([workflow.extends] if workflow.extends is not None else []) + list(
            workflow.replaces
        )
        for reference in references:
            if reference not in workflows:
                raise ValueError(
                    f"workflow '{workflow.id}' references unknown workflow "
                    f"'{reference}'"
                )
            if (
                workflow.visibility == "public"
                and workflows[reference].visibility == "private"
            ):
                raise ValueError(
                    f"public workflow '{workflow.id}' cannot reference private "
                    f"workflow '{reference}'"
                )

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(workflow_id: str, chain: list[str]) -> None:
        if workflow_id in visited:
            return
        if workflow_id in visiting:
            cycle_start = chain.index(workflow_id)
            cycle = chain[cycle_start:] + [workflow_id]
            raise ValueError(
                "recipe inheritance cycle: " + " -> ".join(cycle)
            )

        visiting.add(workflow_id)
        workflow = workflows[workflow_id]
        if workflow.extends is not None:
            visit(workflow.extends, chain + [workflow_id])
        visiting.remove(workflow_id)
        visited.add(workflow_id)

    for workflow_id in workflows:
        visit(workflow_id, [])


def _resolve_provenance_paths(
    provenance: ProvenanceHierarchyConfig,
    config_dir: Path,
) -> None:
    for repository in provenance.search:
        repository.root = _resolved_path(config_dir, repository.root)
    if provenance.write is not None:
        provenance.write.root = _resolved_path(config_dir, provenance.write.root)


def load_config(config_path: Optional[Path] = None) -> MCPConfig:
    """Load configuration from YAML file.

    Args:
        config_path: Path to config file. If None, uses default location.

    Returns:
        MCPConfig instance.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        ValueError: If config is invalid.
    """
    global _config, _config_path, _config_loaded_at
    global _workflow_schemas, _workflow_configs, _workflow_origins, _recipe_catalogs

    if config_path is None:
        # Default to config.yaml in the same directory as this file
        config_path = Path(__file__).parent / "config.yaml"

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        raw_config = yaml.safe_load(f)

    # Expand environment variables
    expanded_config = expand_env_vars(raw_config)

    # Parse into Pydantic model
    config = MCPConfig(**expanded_config)
    config_dir = config_path.parent.resolve()
    _resolve_provenance_paths(config.provenance, config_dir)

    # Load workflow schemas
    repo_root = config.paths.repo_root.resolve()
    schemas, workflow_configs, origins, catalogs = _load_workflow_schemas(
        config, repo_root, config_dir
    )

    # Store globally
    _config = config
    _config_path = config_path
    _config_loaded_at = datetime.now(timezone.utc)
    _workflow_schemas = schemas
    _workflow_configs = workflow_configs
    _workflow_origins = origins
    _recipe_catalogs = catalogs

    return config


def get_config() -> MCPConfig:
    """Get the current configuration.

    Returns:
        MCPConfig instance.

    Raises:
        RuntimeError: If config hasn't been loaded yet.
    """
    if _config is None:
        raise RuntimeError("Configuration not loaded. Call load_config() first.")
    return _config


def get_workflow_schemas() -> dict[str, WorkflowSchema]:
    """Get all loaded workflow schemas."""
    return _workflow_schemas


def get_workflow_schema(workflow_id: str) -> Optional[WorkflowSchema]:
    """Get a workflow schema by ID."""
    return _workflow_schemas.get(workflow_id)


def get_workflow_config(workflow_id: str) -> Optional[WorkflowConfig]:
    """Get the configuration entry for a built-in or external workflow."""
    return _workflow_configs.get(workflow_id)


def get_workflow_configs() -> dict[str, WorkflowConfig]:
    """Get all built-in and external workflow configuration entries."""
    return _workflow_configs


def get_workflow_origin(workflow_id: str) -> Optional[WorkflowOrigin]:
    """Get catalog ownership and schema location for a workflow."""
    return _workflow_origins.get(workflow_id)


def get_workflow_root(workflow_id: str) -> Optional[Path]:
    """Get the filesystem root against which workflow-relative paths resolve."""
    origin = get_workflow_origin(workflow_id)
    if origin is not None:
        return origin.catalog.root

    # Draft/scaffold tests and propose-only tooling may temporarily register a
    # validated schema directly in the in-memory registry before it is added to
    # a catalog. Preserve the historical repo-root resolution for that narrow
    # compatibility path.
    if workflow_id in _workflow_schemas and _config is not None:
        return _config.paths.repo_root.resolve()
    return None


def get_recipe_catalogs() -> tuple[LoadedRecipeCatalog, ...]:
    """Return loaded recipe catalogs in discovery order."""
    return _recipe_catalogs


def get_config_path() -> Optional[Path]:
    """Get the path to the loaded config file."""
    return _config_path


def get_config_loaded_at() -> Optional[datetime]:
    """Get when the config was last loaded."""
    return _config_loaded_at


def reload_config() -> tuple[MCPConfig, datetime]:
    """Reload configuration from the same file.

    Returns:
        Tuple of (new config, load time).

    Raises:
        RuntimeError: If no config was previously loaded.
        ValueError: If new config is invalid (old config is retained).
    """
    global _config, _config_loaded_at
    global _workflow_schemas, _workflow_configs, _workflow_origins, _recipe_catalogs

    if _config_path is None:
        raise RuntimeError("No config previously loaded. Call load_config() first.")

    # Load into a new variable first to validate
    old_config = _config
    old_schemas = _workflow_schemas
    old_workflow_configs = _workflow_configs
    old_workflow_origins = _workflow_origins
    old_recipe_catalogs = _recipe_catalogs
    try:
        new_config = load_config(_config_path)
        return new_config, _config_loaded_at  # type: ignore
    except Exception:
        # Restore old config on failure
        _config = old_config
        _workflow_schemas = old_schemas
        _workflow_configs = old_workflow_configs
        _workflow_origins = old_workflow_origins
        _recipe_catalogs = old_recipe_catalogs
        raise
