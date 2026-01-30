"""MCP tools for STAR-suite."""

from .discovery import list_datasets, list_test_suites, find_docs, find_tests
from .executor import run_script, collect_outputs, get_log_tail
from .preflight import run_preflight
from .reload import reload_config
from .build import build_target, needs_rebuild, ensure_clean_build, smart_build

__all__ = [
    "list_datasets",
    "list_test_suites",
    "find_docs",
    "find_tests",
    "run_preflight",
    "run_script",
    "collect_outputs",
    "get_log_tail",
    "reload_config",
    "build_target",
    "needs_rebuild",
    "ensure_clean_build",
    "smart_build",
]
