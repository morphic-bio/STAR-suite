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
