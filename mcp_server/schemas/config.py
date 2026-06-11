"""Pydantic models for MCP server configuration."""

from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


class ServerConfig(BaseModel):
    """Server connection settings."""

    host: str = "0.0.0.0"
    port: int = 8765
    auth_token: str = Field(default="", description="Bearer token for auth")
    transport: str = Field(default="http", pattern="^(http|stdio)$")
    public_discovery: bool = Field(
        default=True,
        description="If true, discovery tools (list_datasets, list_test_suites, find_docs, find_tests) don't require auth",
    )


class PathsConfig(BaseModel):
    """Path configuration for the server."""

    repo_root: Path = Field(description="Root of the STAR-suite repository")
    artifact_log_root: Path = Field(description="Where to store run logs")
    temp_root: Path = Field(default=Path("/tmp"), description="Temp directory root")


class DatasetConfig(BaseModel):
    """Configuration for a dataset."""

    id: str = Field(description="Unique identifier for the dataset")
    path: Path = Field(description="Absolute path to the dataset")
    description: str = Field(default="", description="Human-readable description")
    metadata_file: Optional[str] = Field(
        default=None, description="Optional metadata file relative to path"
    )


class ExecutionConfig(BaseModel):
    """Execution constraints."""

    default_timeout_seconds: int = Field(default=7200, description="Default 2 hours")
    max_concurrent_jobs: int = Field(default=1, description="Max concurrent jobs")
    queue_max_size: int = Field(default=10, description="Max queued jobs")


class DiskSpaceConfig(BaseModel):
    """Disk space thresholds."""

    default_min_gb: int = Field(default=20, description="Default minimum GB")
    full_depth_min_gb: int = Field(default=100, description="Full-depth minimum GB")


class ScriptConfig(BaseModel):
    """Configuration for an allowlisted script."""

    name: str = Field(description="Unique script identifier")
    path: str = Field(description="Path to script (relative to repo root or absolute)")
    module: str = Field(description="Module this script belongs to")
    description: str = Field(default="", description="Human-readable description")
    timeout_seconds: Optional[int] = Field(
        default=None, description="Override default timeout"
    )
    args: Optional[list[str]] = Field(
        default=None, description="Default arguments for the script"
    )
    working_dir: Optional[str] = Field(
        default=None, description="Working directory (relative to repo root)"
    )
    fixtures: Optional[list[str]] = Field(
        default=None, description="Required fixture paths"
    )
    runnable: bool = Field(
        default=True, description="Whether this script is currently runnable"
    )


class TestSuiteConfig(BaseModel):
    """Configuration for a test suite."""

    module: str = Field(description="Module identifier")
    description: str = Field(default="", description="Human-readable description")
    scripts: list[str] = Field(description="List of script names in this suite")
    disk_space_min_gb: Optional[int] = Field(
        default=None, description="Override disk space requirement"
    )
    runnable: bool = Field(default=True, description="Whether suite is runnable")


class BinaryConfig(BaseModel):
    """Configuration for a required binary."""

    name: str = Field(description="Binary name (e.g., 'STAR')")
    paths: list[str] = Field(description="Paths to search for the binary")


class WorkflowConfig(BaseModel):
    """Configuration entry for a structured workflow."""

    id: str = Field(description="Unique workflow identifier")
    title: str = Field(description="Human-readable title")
    summary: str = Field(default="", description="One-line summary")
    entry_script: str = Field(description="Entry script path (relative to repo root)")
    kind: str = Field(default="shell_workflow", description="Workflow kind")
    schema_file: str = Field(description="Path to workflow schema YAML (relative to repo root)")
    visibility: str = Field(
        default="public",
        pattern="^(public|private)$",
        description="public = visible without auth; private = requires valid auth_token",
    )


# Site-specific paths that fill the ${...} placeholders in DEFAULT_AGENT_PROTOCOL.
# These are morphic's values, mirrored from the canonical
# morphic-recipes/docs/authoring/localization.example.yaml. Other sites override
# via `locations:` in config.yaml (the protocol text itself stays unchanged).
DEFAULT_LOCATIONS = {
    "PROVENANCE_ROOT": "/mnt/pikachu/morphic-provenance",
    "RECIPES_ROOT": "/mnt/pikachu/morphic-recipes",
}

# Lab-agnostic protocol TEMPLATE. ${PROVENANCE_ROOT}/${RECIPES_ROOT} are filled at
# load time from `locations` (see MCPConfig.rendered_agent_protocol). Canonical
# copy: morphic-recipes/docs/authoring/agent_protocol.txt.
DEFAULT_AGENT_PROTOCOL = (
    "PROVENANCE-FIRST EXECUTION. Recipes and workflows are STARTING POINTS, not "
    "turnkey commands. Resource/scale parameters (thread counts, memory, "
    "--*-low-mem flags, start mode) MUST come from a known-good run, never invented. "
    "Before running at non-trivial scale: (1) consult your provenance record "
    "(${PROVENANCE_ROOT}) runs/<project>/<run_id>/{run.json,commands/} "
    "for the exact parameters that actually worked; (2) reproduce them, adapting only "
    "what the new input/machine requires and noting the deviation; (3) if none exists "
    "at your scale, start from the closest run, scale conservatively, and record a new "
    "run. Recipes live in your recipe repo (${RECIPES_ROOT}); this "
    "server exposes suite workflows only — cross-reference recipes + provenance before "
    "executing. Inventing thread/memory params and running blind is a known OOM "
    "failure mode. "
    "COMPOSE-UP EXECUTION. A recipe often emits a SUPERSET of outputs; do not run the "
    "maximal set blindly. (1) Start from the MINIMAL functional core (the tested floor); "
    "(2) ADD only the output layers your target workflow actually consumes — check the "
    "recipe's COMPOSITION block for what is optional, how to add/drop it, and which "
    "provenance run is the parameter oracle for each layer; (3) preview with --dry-run "
    "(hand to a human if useful) before the real run. Emitting layers the target never "
    "uses wastes compute and distorts benchmark comparisons against tools that emitted "
    "less (e.g. RNA-velocity / extra BAMs vs Cell Ranger ARC --no-bam). Step back and "
    "match outputs to the task instead of running a recipe verbatim."
)


class MCPConfig(BaseModel):
    """Root configuration model."""

    server: ServerConfig = Field(default_factory=ServerConfig)
    paths: PathsConfig
    trusted_roots: list[str] = Field(description="Trusted path roots")
    datasets: list[DatasetConfig] = Field(default_factory=list)
    execution: ExecutionConfig = Field(default_factory=ExecutionConfig)
    disk_space: DiskSpaceConfig = Field(default_factory=DiskSpaceConfig)
    scripts: list[ScriptConfig] = Field(default_factory=list)
    test_suites: list[TestSuiteConfig] = Field(default_factory=list)
    required_binaries: list[BinaryConfig] = Field(default_factory=list)
    workflows: list[WorkflowConfig] = Field(default_factory=list)
    agent_protocol: str = Field(
        default=DEFAULT_AGENT_PROTOCOL,
        description="Provenance-first + compose-up protocol TEMPLATE surfaced to agents on workflow discovery.",
    )
    locations: dict[str, str] = Field(
        default_factory=lambda: dict(DEFAULT_LOCATIONS),
        description="Site-specific paths that fill the ${...} placeholders in agent_protocol.",
    )

    @property
    def rendered_agent_protocol(self) -> str:
        """agent_protocol with ${PROVENANCE_ROOT}/${RECIPES_ROOT} filled from locations."""
        from string import Template
        return Template(self.agent_protocol).safe_substitute(self.locations)

    def get_script(self, name: str) -> Optional[ScriptConfig]:
        """Get a script by name."""
        for script in self.scripts:
            if script.name == name:
                return script
        return None

    def get_test_suite(self, module: str) -> Optional[TestSuiteConfig]:
        """Get a test suite by module name."""
        for suite in self.test_suites:
            if suite.module == module:
                return suite
        return None

    def get_dataset(self, dataset_id: str) -> Optional[DatasetConfig]:
        """Get a dataset by ID."""
        for dataset in self.datasets:
            if dataset.id == dataset_id:
                return dataset
        return None

    def get_workflow(self, workflow_id: str) -> Optional["WorkflowConfig"]:
        """Get a workflow config by ID."""
        for wf in self.workflows:
            if wf.id == workflow_id:
                return wf
        return None
