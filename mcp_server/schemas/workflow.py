"""Pydantic models for workflow metadata and parameter schemas."""

from typing import Any, Optional

from pydantic import BaseModel, Field


class WorkflowParameterDef(BaseModel):
    """Definition of a single workflow parameter."""

    name: str = Field(description="Parameter identifier (snake_case)")
    cli_flag: str = Field(description="CLI flag as used by the entry script (e.g. --samples)")
    type: str = Field(
        description="Parameter type",
        pattern="^(string|int|float|bool|enum|file|directory|string_list)$",
    )
    required: bool = Field(default=False, description="Whether this parameter is required")
    default: Any = Field(default=None, description="Default value (null if none)")
    description: str = Field(default="", description="Human-readable description")
    choices: Optional[list[str]] = Field(
        default=None, description="Allowed values for enum type"
    )
    repeatable: bool = Field(default=False, description="Whether the flag can appear multiple times")
    operand_group: Optional[str] = Field(
        default=None,
        description=(
            "If set, consecutive params in flag_order with the same operand_group and "
            "cli_flag render as one flag followed by each value (e.g. STAR --readFilesIn mate1 mate2)."
        ),
    )
    path_must_exist: bool = Field(
        default=False, description="If true, path must exist at validation time"
    )
    must_be_executable: bool = Field(
        default=False, description="If true, file must have execute permission"
    )
    category: str = Field(default="general", description="Parameter group/category")
    stage: str = Field(
        default="top_level", description="Which stage of the workflow owns this param"
    )
    source: str = Field(
        default="workflow_wrapper",
        description="Which script layer defines this param",
    )

    # --- Rendering hints (declarative, not hardcoded in the renderer) ---
    env_var: Optional[str] = Field(
        default=None,
        description="If set, export this value as the named env var instead of (or in addition to) a CLI flag",
    )
    skip_when_default: bool = Field(
        default=False,
        description="If true, omit the flag when the value equals the schema default",
    )
    is_output_root: bool = Field(
        default=False,
        description="If true, this parameter's value is reported as the workflow output_root",
    )
    ui_gated_by: Optional[str] = Field(
        default=None,
        description=(
            "Optional boolean parameter name: human-facing UIs may hide this field "
            "while that parameter is true (no effect on server validation/rendering)."
        ),
    )

    # --- Planner/form metadata and optional numeric validation hints ---
    label: Optional[str] = Field(
        default=None,
        description="Stable human-readable label for form/UI display (advisory).",
    )
    help: Optional[str] = Field(
        default=None,
        description="Extended help text beyond description; may contain usage guidance (advisory).",
    )
    example: Optional[str] = Field(
        default=None,
        description="Example value for documentation and placeholder hints (advisory).",
    )
    widget_hint: Optional[str] = Field(
        default=None,
        description=(
            "Hint for UI rendering: text, number, checkbox, select, path, textarea, "
            "readonly, hidden. Clients may ignore (advisory)."
        ),
    )
    aliases: Optional[list[str]] = Field(
        default=None,
        description=(
            "Canonical alternate names for this parameter (e.g. camelCase STAR flag names). "
            "Used for alias normalization by planners and agents (advisory)."
        ),
    )
    advanced: bool = Field(
        default=False,
        description="If true, UIs may collapse this field into an advanced section (advisory).",
    )
    display_order: Optional[int] = Field(
        default=None,
        description="Explicit sort key within its group for form rendering (advisory).",
    )
    min_value: Optional[float] = Field(
        default=None,
        description="Minimum allowed numeric value, if applicable.",
    )
    max_value: Optional[float] = Field(
        default=None,
        description="Maximum allowed numeric value, if applicable.",
    )


class WorkflowConstraint(BaseModel):
    """A constraint between parameters."""

    kind: str = Field(
        description="Constraint type",
        pattern="^(mutual_exclusion|at_least_one|dependency|group_required|positive|non_negative)$",
    )
    params: list[str] = Field(description="Parameters involved in the constraint")
    message: str = Field(default="", description="Human-readable constraint description")


class WorkflowParameterGroup(BaseModel):
    """A named group of parameters for display ordering."""

    name: str = Field(description="Group identifier")
    title: str = Field(description="Human-readable group title")
    description: str = Field(default="", description="Help text for this group (advisory).")
    parameters: list[str] = Field(description="Parameter names in this group (display order)")
    gated_by: Optional[str] = Field(
        default=None,
        description=(
            "When set, UIs may hide this entire group while the named boolean "
            "parameter is true (e.g. star_only skips downstream-related fields)."
        ),
    )


class WorkflowScriptDef(BaseModel):
    """A script that participates in the workflow."""

    role: str = Field(description="Script role (e.g. entry, batch_runner, downstream, helper)")
    path: str = Field(description="Path relative to repo root")
    description: str = Field(default="", description="What this script does")
    language: str = Field(default="bash", description="Script language (bash, python)")


class WorkflowStage(BaseModel):
    """A semantic stage within the workflow."""

    name: str = Field(description="Stage identifier")
    title: str = Field(description="Human-readable title")
    script: Optional[str] = Field(default=None, description="Script used by this stage")
    description: str = Field(default="")
    gated_by: Optional[str] = Field(
        default=None,
        description="Parameter that gates this stage (e.g. star_only=false enables downstream)",
    )


class WorkflowRenderingRule(BaseModel):
    """Rules for how parameters are rendered to CLI flags."""

    bool_style: str = Field(
        default="flag_only",
        description="How booleans render: 'flag_only' = emit flag when true, omit when false",
    )
    flag_order: list[str] = Field(
        default_factory=list,
        description="Deterministic flag emission order (by parameter name)",
    )
    omit_absent_optionals: bool = Field(
        default=True, description="Skip flags for unset optional params"
    )
    csv_style: str = Field(
        default="quoted", description="How string_list params are emitted"
    )


class WorkflowSchema(BaseModel):
    """Complete machine-readable schema for a workflow.

    Loaded from a per-workflow YAML file (e.g. workflows/ucsf_star_suite_production.yaml).
    """

    id: str = Field(description="Workflow identifier")
    title: str = Field(description="Display name")
    summary: str = Field(description="One-line summary")
    kind: str = Field(default="shell_workflow", description="Workflow kind")
    entry_script: str = Field(description="Path to entry script relative to repo root")
    supported_modes: list[str] = Field(
        default_factory=lambda: ["local", "dry-run-capable"],
        description="Execution modes this workflow supports",
    )
    caveats: list[str] = Field(default_factory=list, description="Known caveats")
    default_output_layout: str = Field(
        default="",
        description="Description of default output directory structure",
    )

    scripts: list[WorkflowScriptDef] = Field(default_factory=list)
    stages: list[WorkflowStage] = Field(default_factory=list)
    parameter_groups: list[WorkflowParameterGroup] = Field(default_factory=list)
    parameters: list[WorkflowParameterDef] = Field(default_factory=list)
    constraints: list[WorkflowConstraint] = Field(default_factory=list)
    rendering: WorkflowRenderingRule = Field(default_factory=WorkflowRenderingRule)

    def get_parameter(self, name: str) -> Optional[WorkflowParameterDef]:
        """Get a parameter definition by name."""
        for p in self.parameters:
            if p.name == name:
                return p
        return None

    def get_required_parameters(self) -> list[WorkflowParameterDef]:
        """Get all required parameters."""
        return [p for p in self.parameters if p.required]

    def get_parameters_by_category(self, category: str) -> list[WorkflowParameterDef]:
        """Get parameters in a given category."""
        return [p for p in self.parameters if p.category == category]
