"""Release integration tests for the pinned official recipe/evidence snapshots."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml

from mcp_server.config import _load_workflow_schemas
from mcp_server.schemas.config import MCPConfig


REPO_ROOT = Path(__file__).resolve().parents[2]
SHARE_ROOT = REPO_ROOT / "share/star-suite"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_official_snapshot_metadata_matches_vendored_files():
    metadata = json.loads((SHARE_ROOT / "SNAPSHOTS.json").read_text())
    recipes = metadata["recipes"]
    provenance = metadata["provenance"]

    catalog = SHARE_ROOT / recipes["catalog"]
    record_schema = SHARE_ROOT / provenance["record_schema"]
    assert _sha256(catalog) == recipes["catalog_sha256"]
    assert _sha256(record_schema) == provenance["record_schema_sha256"]
    assert len(list((SHARE_ROOT / "evidence/official/records").rglob("*.json"))) == 4


def test_default_config_loads_pinned_official_catalog_and_evidence():
    config_path = REPO_ROOT / "mcp_server/config.yaml"
    config = MCPConfig(**yaml.safe_load(config_path.read_text()))
    schemas, _, _, catalogs = _load_workflow_schemas(
        config, REPO_ROOT, config_path.parent
    )

    assert catalogs[-1].id == "starsuite-official"
    assert catalogs[-1].trust == "trusted"
    assert "starsuite.official/partitioned-bulk-pe-alignment" in schemas
    assert len([name for name in schemas if name.startswith("starsuite.official/")]) == 11
    assert config.provenance.search[0].id == "official-evidence"
    assert (config_path.parent / config.provenance.search[0].root).resolve() == (
        SHARE_ROOT / "evidence/official"
    ).resolve()
