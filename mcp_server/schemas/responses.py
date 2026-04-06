"""Pydantic models for MCP tool responses."""

from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel, Field


# --- list_datasets ---


class DatasetInfo(BaseModel):
    """Information about a dataset."""

    id: str
    path: str
    description: str
    size_bytes: Optional[int] = None
    mtime: Optional[datetime] = None
    metadata: Optional[dict[str, Any]] = None


class ListDatasetsResponse(BaseModel):
    """Response from list_datasets."""

    datasets: list[DatasetInfo]


# --- list_test_suites ---


class ScriptInfo(BaseModel):
    """Information about a script within a test suite."""

    name: str
    path: str
    description: str
    runnable: bool
    missing_fixtures: list[str] = Field(default_factory=list)


class TestSuiteInfo(BaseModel):
    """Information about a test suite."""

    module: str
    description: str
    disk_space_min_gb: int
    runnable: bool
    scripts: list[ScriptInfo]


class ListTestSuitesResponse(BaseModel):
    """Response from list_test_suites."""

    suites: list[TestSuiteInfo]


# --- find_docs ---


class DocMatch(BaseModel):
    """A documentation match."""

    path: str
    title: Optional[str] = None
    snippet: Optional[str] = None


class FindDocsResponse(BaseModel):
    """Response from find_docs."""

    matches: list[DocMatch]


# --- find_tests ---


class TestMatch(BaseModel):
    """A test script match."""

    name: str
    path: str
    module: str
    description: str


class FindTestsResponse(BaseModel):
    """Response from find_tests."""

    matches: list[TestMatch]


# --- preflight ---


class PreflightCheck(BaseModel):
    """Result of a single preflight check."""

    name: str
    passed: bool
    details: Optional[dict[str, Any]] = None
    message: Optional[str] = None


class PreflightResponse(BaseModel):
    """Response from preflight validation."""

    valid: bool
    checks: list[PreflightCheck]
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


# --- run_script ---


class RunScriptResponse(BaseModel):
    """Immediate response from run_script."""

    run_id: str
    status: str  # "running", "queued"
    log_path: str
    stream_url: Optional[str] = None


class RunCompletionResponse(BaseModel):
    """Response when a run completes."""

    run_id: str
    status: str  # "completed", "failed", "timeout"
    exit_code: Optional[int] = None
    duration_seconds: Optional[float] = None
    log_path: str
    log_tail: Optional[str] = None
    outputs: list[str] = Field(default_factory=list)


# --- collect_outputs ---


class OutputFile(BaseModel):
    """Information about an output file."""

    path: str
    size_bytes: int
    mtime: datetime


class CollectOutputsResponse(BaseModel):
    """Response from collect_outputs."""

    run_id: str
    status: str
    exit_code: Optional[int] = None
    duration_seconds: Optional[float] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    script: str
    args: list[str] = Field(default_factory=list)
    log_files: dict[str, str] = Field(default_factory=dict)
    outputs: list[OutputFile] = Field(default_factory=list)


# --- reload_config ---


class ReloadConfigResponse(BaseModel):
    """Response from reload_config."""

    reloaded: bool
    config_path: str
    loaded_at: datetime
    error: Optional[str] = None


# --- errors ---


# --- workflows ---


class WorkflowInfo(BaseModel):
    """Summary information about a workflow."""

    id: str
    title: str
    summary: str
    kind: str
    entry_script: str
    supported_modes: list[str] = Field(default_factory=list)


class ListWorkflowsResponse(BaseModel):
    """Response from list_workflows."""

    workflows: list[WorkflowInfo]


class WorkflowStageInfo(BaseModel):
    """Stage info returned by describe_workflow."""

    name: str
    title: str
    script: Optional[str] = None
    description: str = ""
    gated_by: Optional[str] = None


class DescribeWorkflowResponse(BaseModel):
    """Response from describe_workflow."""

    id: str
    title: str
    summary: str
    kind: str
    entry_script: str
    supported_modes: list[str] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)
    default_output_layout: str = ""
    stages: list[WorkflowStageInfo] = Field(default_factory=list)
    parameter_groups: list["ParameterGroupInfo"] = Field(default_factory=list)


class ParameterInfo(BaseModel):
    """Machine-readable parameter definition."""

    name: str
    cli_flag: str
    type: str
    required: bool
    default: Any = None
    description: str = ""
    choices: Optional[list[str]] = None
    repeatable: bool = False
    path_must_exist: bool = False
    must_be_executable: bool = False
    category: str = "general"
    stage: str = "top_level"
    source: str = "workflow_wrapper"
    env_var: Optional[str] = None
    skip_when_default: bool = False
    is_output_root: bool = False


class ConstraintInfo(BaseModel):
    """Constraint between parameters."""

    kind: str
    params: list[str]
    message: str = ""


class ParameterGroupInfo(BaseModel):
    """Ordered parameter group."""

    name: str
    title: str
    parameters: list[str]


class WorkflowParameterSchemaResponse(BaseModel):
    """Response from get_workflow_parameter_schema."""

    workflow_id: str
    parameters: list[ParameterInfo]
    parameter_groups: list[ParameterGroupInfo] = Field(default_factory=list)
    constraints: list[ConstraintInfo] = Field(default_factory=list)


class ValidateWorkflowResponse(BaseModel):
    """Response from validate_workflow_parameters."""

    valid: bool
    normalized_params: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class RenderWorkflowResponse(BaseModel):
    """Response from render_workflow_command."""

    workflow_id: str
    entry_script: str
    argv: list[str]
    shell_preview: str
    env_overrides: dict[str, str] = Field(default_factory=dict)
    output_root: Optional[str] = None


# --- errors ---


class ScaffoldWorkflowResponse(BaseModel):
    """Response from scaffold_workflow_schema."""

    draft_yaml: str
    parameter_count: int
    source_script: str
    notes: list[str] = Field(default_factory=list)


class ValidateDraftWorkflowResponse(BaseModel):
    """Response from validate_draft_workflow_schema."""

    valid: bool
    errors: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    parameter_count: Optional[int] = None


class ErrorResponse(BaseModel):
    """Structured error response."""

    error: bool = True
    code: str
    message: str
    details: Optional[dict[str, Any]] = None
