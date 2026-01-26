"""Pydantic models for run configuration."""

from typing import Optional

from pydantic import BaseModel, Field


class RunConfig(BaseModel):
    """Configuration for a script run request."""

    script: str = Field(description="Name of the script to run")
    args: Optional[list[str]] = Field(
        default=None, description="Additional arguments to pass"
    )
    module: Optional[str] = Field(
        default=None, description="Module context (for validation)"
    )
    dataset_id: Optional[str] = Field(
        default=None, description="Dataset to use (if applicable)"
    )
    out_dir: Optional[str] = Field(
        default=None, description="Output directory (defaults to temp)"
    )
    env_overrides: Optional[dict[str, str]] = Field(
        default=None, description="Environment variable overrides"
    )
    # Note: stream_logs was removed - streaming requires SSE infrastructure
    # that will be added in a future phase. Logs are always written to files
    # and can be retrieved via collect_outputs() or get_log_tail().
