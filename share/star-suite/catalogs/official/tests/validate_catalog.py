#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
STAR_SUITE_ROOT = Path(
    os.environ.get("STAR_SUITE_ROOT", ROOT.parent / "STAR-suite")
).resolve()
if not (STAR_SUITE_ROOT / "mcp_server/schemas/config.py").is_file():
    raise SystemExit(
        "ERROR: set STAR_SUITE_ROOT to a STAR Suite checkout with catalog support"
    )
sys.path.insert(0, str(STAR_SUITE_ROOT))

from mcp_server.config import _load_workflow_schemas  # noqa: E402
from mcp_server.schemas.config import MCPConfig, RecipeCatalogManifest  # noqa: E402
from mcp_server.schemas.workflow import WorkflowSchema  # noqa: E402


manifest = RecipeCatalogManifest(**yaml.safe_load((ROOT / "catalog.yaml").read_text()))
if manifest.id != "starsuite-official" or manifest.namespace != "starsuite.official":
    raise SystemExit("ERROR: unexpected official catalog identity")
if not manifest.workflows:
    raise SystemExit("ERROR: catalog is empty")

for workflow in manifest.workflows:
    if not workflow.id.startswith(f"{manifest.namespace}/"):
        raise SystemExit(f"ERROR: workflow outside namespace: {workflow.id}")
    for value in (workflow.entry_script, workflow.schema_file):
        path = Path(value)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"ERROR: non-confined path in {workflow.id}: {value}")
        if not (ROOT / path).is_file():
            raise SystemExit(f"ERROR: missing catalog file: {value}")
    schema = WorkflowSchema(**yaml.safe_load((ROOT / workflow.schema_file).read_text()))
    if schema.id != workflow.id or schema.entry_script != workflow.entry_script:
        raise SystemExit(f"ERROR: manifest/schema mismatch: {workflow.id}")

config = MCPConfig(
    paths={
        "repo_root": str(STAR_SUITE_ROOT),
        "artifact_log_root": str(ROOT / "run-output"),
    },
    trusted_roots=[str(ROOT), str(STAR_SUITE_ROOT)],
    recipe_catalogs=[{"manifest": str(ROOT / "catalog.yaml"), "trust": "trusted"}],
)
schemas, _, _, catalogs = _load_workflow_schemas(
    config, STAR_SUITE_ROOT, ROOT / "tests"
)
if len(schemas) != len(manifest.workflows):
    raise SystemExit("ERROR: external catalog loader returned the wrong workflow count")
if catalogs[-1].id != manifest.id:
    raise SystemExit("ERROR: external catalog identity was not retained")
print(f"PASS: loaded {len(schemas)} official workflows")

