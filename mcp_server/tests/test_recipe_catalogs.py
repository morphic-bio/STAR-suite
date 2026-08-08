"""Tests for ordered external recipe catalogs and provenance roots."""

from pathlib import Path

import pytest
import yaml
import json

import mcp_server.config as config_module
from mcp_server.config import (
    BUILTIN_CATALOG_ID,
    get_recipe_catalogs,
    get_workflow_config,
    get_workflow_origin,
    get_workflow_root,
    get_workflow_schema,
    load_config,
    reload_config,
)
from mcp_server.tools.workflows import (
    get_workflow_scripts,
    list_workflows,
    render_workflow_command,
)
from mcp_server.tools.recipes import (
    build_recipe_lock,
    create_recipe_bundle,
    describe_recipe_catalog,
    list_recipe_candidates,
    list_recipe_catalogs,
    list_recipe_conflicts,
    list_provenance_repositories,
    resolve_recipe,
)
from mcp_server.recipe_cli import main as recipe_cli_main


@pytest.fixture(autouse=True)
def _reset_config_state():
    yield
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None
    config_module._workflow_schemas = {}
    config_module._workflow_configs = {}
    config_module._workflow_origins = {}
    config_module._recipe_catalogs = ()


def _workflow_schema(workflow_id: str, entry_script: str = "scripts/run.sh") -> dict:
    return {
        "id": workflow_id,
        "title": f"Workflow {workflow_id}",
        "summary": "External catalog workflow.",
        "entry_script": entry_script,
        "scripts": [
            {
                "role": "entry",
                "path": entry_script,
                "description": "Catalog entry point.",
                "language": "bash",
            },
            {
                "role": "helper",
                "path": "helpers/helper.py",
                "description": "Catalog-local helper.",
                "language": "python",
            },
        ],
        "parameters": [],
        "rendering": {"flag_order": []},
    }


def _write_catalog(
    root: Path,
    *,
    catalog_id: str,
    namespace: str,
    workflow_local_name: str = "demo",
    entry_script: str = "scripts/run.sh",
    logical_id: str | None = None,
    recipe_version: str = "0",
    applications: list[str] | None = None,
    extends: str | None = None,
    replaces: list[str] | None = None,
    image: str | None = None,
) -> tuple[Path, str]:
    workflow_id = f"{namespace}/{workflow_local_name}"
    (root / "workflows").mkdir(parents=True)
    (root / "scripts").mkdir(parents=True)
    (root / "helpers").mkdir(parents=True)
    (root / "scripts" / "run.sh").write_text("#!/bin/sh\necho catalog\n")
    (root / "helpers" / "helper.py").write_text("print('helper')\n")
    schema_file = f"workflows/{workflow_local_name}.yaml"
    (root / schema_file).write_text(
        yaml.safe_dump(_workflow_schema(workflow_id, entry_script))
    )
    workflow = {
        "id": workflow_id,
        "title": f"Workflow {workflow_id}",
        "summary": "External catalog workflow.",
        "entry_script": entry_script,
        "kind": "shell_workflow",
        "schema_file": schema_file,
        "version": recipe_version,
    }
    if logical_id is not None:
        workflow["logical_id"] = logical_id
    if applications is not None:
        workflow["applications"] = applications
    if extends is not None:
        workflow["extends"] = extends
    if replaces is not None:
        workflow["replaces"] = replaces
    if image is not None:
        workflow["image"] = image
    manifest = {
        "schema": "biodepot.recipe_catalog/v1",
        "id": catalog_id,
        "namespace": namespace,
        "version": "1",
        "workflows": [workflow],
    }
    manifest_path = root / "catalog.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return manifest_path, workflow_id


def _write_config(
    root: Path,
    *,
    recipe_catalogs: list[dict] | None = None,
    workflows: list[dict] | None = None,
    provenance: dict | None = None,
    recipe_resolution: dict | None = None,
    auth_token: str = "",
) -> Path:
    repo_root = root / "suite"
    repo_root.mkdir(exist_ok=True)
    config = {
        "server": {
            "host": "127.0.0.1",
            "port": 9999,
            "transport": "http",
            "auth_token": auth_token,
            "public_discovery": True,
        },
        "paths": {
            "repo_root": str(repo_root),
            "artifact_log_root": str(root / "artifacts"),
            "temp_root": str(root / "tmp"),
        },
        "trusted_roots": [str(root), "/tmp"],
    }
    if recipe_catalogs is not None:
        config["recipe_catalogs"] = recipe_catalogs
    if workflows is not None:
        config["workflows"] = workflows
    if provenance is not None:
        config["provenance"] = provenance
    if recipe_resolution is not None:
        config["recipe_resolution"] = recipe_resolution
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path


def test_legacy_workflows_are_the_builtin_compatibility_catalog(tmp_path: Path):
    repo_root = tmp_path / "suite"
    (repo_root / "workflows").mkdir(parents=True)
    schema = _workflow_schema("legacy", "scripts/legacy.sh")
    (repo_root / "workflows" / "legacy.yaml").write_text(yaml.safe_dump(schema))
    config_path = _write_config(
        tmp_path,
        workflows=[
            {
                "id": "legacy",
                "title": "Legacy",
                "entry_script": "scripts/legacy.sh",
                "schema_file": "workflows/legacy.yaml",
            }
        ],
    )

    load_config(config_path)

    origin = get_workflow_origin("legacy")
    assert origin is not None
    assert origin.catalog.id == BUILTIN_CATALOG_ID
    assert origin.catalog.built_in is True
    assert get_workflow_root("legacy") == repo_root.resolve()
    assert [catalog.id for catalog in get_recipe_catalogs()] == [BUILTIN_CATALOG_ID]


def test_external_catalogs_load_in_declared_order(tmp_path: Path):
    first_manifest, first_id = _write_catalog(
        tmp_path / "first", catalog_id="first", namespace="example.first"
    )
    second_manifest, second_id = _write_catalog(
        tmp_path / "second", catalog_id="second", namespace="example.second"
    )
    config_path = _write_config(
        tmp_path,
        recipe_catalogs=[
            {"manifest": str(first_manifest.relative_to(tmp_path))},
            {"manifest": str(second_manifest.relative_to(tmp_path))},
        ],
    )

    load_config(config_path)

    assert [catalog.id for catalog in get_recipe_catalogs()] == [
        BUILTIN_CATALOG_ID,
        "first",
        "second",
    ]
    assert [item.id for item in list_workflows().workflows] == [first_id, second_id]
    assert get_workflow_config(first_id) is not None
    assert get_workflow_schema(second_id) is not None


def test_external_paths_and_provenance_use_owning_catalog_root(tmp_path: Path):
    manifest, workflow_id = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    load_config(
        _write_config(
            tmp_path,
            recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
        )
    )

    rendered = render_workflow_command(workflow_id, {})
    assert rendered.argv[0] == str(tmp_path / "external" / "scripts" / "run.sh")

    public = get_workflow_scripts(workflow_id, authenticated=False)
    assert public.provenance["catalog_id"] == "external"
    assert public.provenance["catalog_namespace"] == "example.external"
    assert public.provenance["catalog_trust"] == "untrusted"
    assert public.provenance["catalog_manifest"] == "catalog.yaml"
    assert "catalog_root" not in public.provenance
    assert "repo_root" not in public.provenance
    assert all(script.absolute_path is None for script in public.scripts)

    authenticated = get_workflow_scripts(workflow_id, authenticated=True)
    assert authenticated.provenance["catalog_root"] == str(
        (tmp_path / "external").resolve()
    )
    assert authenticated.provenance["repo_root"] == str(
        (tmp_path / "external").resolve()
    )
    assert all(
        script.absolute_path.startswith(str(tmp_path / "external"))
        for script in authenticated.scripts
    )

    summary = list_workflows().workflows[0]
    assert summary.catalog_id == "external"
    assert summary.catalog_namespace == "example.external"


@pytest.mark.parametrize("bad_path", ["../outside.yaml", "/tmp/outside.yaml"])
def test_external_schema_path_must_stay_in_catalog(
    tmp_path: Path, bad_path: str
):
    manifest, _ = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    raw = yaml.safe_load(manifest.read_text())
    raw["workflows"][0]["schema_file"] = bad_path
    manifest.write_text(yaml.safe_dump(raw, sort_keys=False))

    with pytest.raises(ValueError, match="schema_file.*catalog root"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


@pytest.mark.parametrize("bad_entry", ["../run.sh", "/usr/bin/env"])
def test_external_entry_script_must_stay_in_catalog(
    tmp_path: Path, bad_entry: str
):
    manifest, workflow_id = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    raw_manifest = yaml.safe_load(manifest.read_text())
    raw_manifest["workflows"][0]["entry_script"] = bad_entry
    manifest.write_text(yaml.safe_dump(raw_manifest, sort_keys=False))
    schema_path = tmp_path / "external" / "workflows" / "demo.yaml"
    raw_schema = yaml.safe_load(schema_path.read_text())
    raw_schema["entry_script"] = bad_entry
    raw_schema["scripts"][0]["path"] = bad_entry
    schema_path.write_text(yaml.safe_dump(raw_schema, sort_keys=False))

    with pytest.raises(ValueError, match=f"workflow '{workflow_id}'.*catalog root"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


def test_external_execution_stage_must_stay_in_catalog(tmp_path: Path):
    manifest, workflow_id = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    schema_path = tmp_path / "external" / "workflows" / "demo.yaml"
    raw_schema = yaml.safe_load(schema_path.read_text())
    raw_schema["execution"] = {
        "decomposition_id": "demo/v1",
        "scatter": {"parameter": "items", "item_name": "item"},
        "stages": [
            {
                "name": "scatter",
                "title": "Scatter",
                "entry_script": "../outside.sh",
                "resource_class": "small",
                "image": "example.invalid/image@sha256:test",
            }
        ],
        "gather": {"stage": "scatter"},
    }
    schema_path.write_text(yaml.safe_dump(raw_schema, sort_keys=False))

    with pytest.raises(ValueError, match=f"workflow '{workflow_id}'.*execution stage"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


def test_duplicate_catalog_id_fails_closed(tmp_path: Path):
    first, _ = _write_catalog(
        tmp_path / "first", catalog_id="duplicate", namespace="example.first"
    )
    second, _ = _write_catalog(
        tmp_path / "second", catalog_id="duplicate", namespace="example.second"
    )
    with pytest.raises(ValueError, match="duplicate recipe catalog id"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[
                    {"manifest": str(first.relative_to(tmp_path))},
                    {"manifest": str(second.relative_to(tmp_path))},
                ],
            )
        )


def test_duplicate_catalog_namespace_fails_closed(tmp_path: Path):
    first, _ = _write_catalog(
        tmp_path / "first", catalog_id="first", namespace="example.shared"
    )
    second, _ = _write_catalog(
        tmp_path / "second", catalog_id="second", namespace="example.shared"
    )
    with pytest.raises(ValueError, match="duplicate recipe catalog namespace"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[
                    {"manifest": str(first.relative_to(tmp_path))},
                    {"manifest": str(second.relative_to(tmp_path))},
                ],
            )
        )


def test_duplicate_workflow_id_within_catalog_fails_closed(tmp_path: Path):
    manifest, _ = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    raw = yaml.safe_load(manifest.read_text())
    raw["workflows"].append(dict(raw["workflows"][0]))
    manifest.write_text(yaml.safe_dump(raw, sort_keys=False))

    with pytest.raises(ValueError, match="duplicate workflow id"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


def test_workflow_id_must_use_catalog_namespace(tmp_path: Path):
    manifest, _ = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    raw = yaml.safe_load(manifest.read_text())
    raw["workflows"][0]["id"] = "wrong/demo"
    manifest.write_text(yaml.safe_dump(raw, sort_keys=False))

    with pytest.raises(ValueError, match="must use namespace prefix"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


def test_missing_catalog_manifest_fails_with_resolved_source(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="Recipe catalog manifest not found"):
        load_config(
            _write_config(
                tmp_path, recipe_catalogs=[{"manifest": "missing/catalog.yaml"}]
            )
        )


def test_failed_reload_restores_catalog_registry_atomically(tmp_path: Path):
    manifest, workflow_id = _write_catalog(
        tmp_path / "external", catalog_id="external", namespace="example.external"
    )
    load_config(
        _write_config(
            tmp_path,
            recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
        )
    )
    prior_origin = get_workflow_origin(workflow_id)
    prior_catalogs = get_recipe_catalogs()

    raw = yaml.safe_load(manifest.read_text())
    raw["workflows"][0]["schema_file"] = "../escape.yaml"
    manifest.write_text(yaml.safe_dump(raw, sort_keys=False))

    with pytest.raises(ValueError, match="escapes its catalog root"):
        reload_config()

    assert get_workflow_schema(workflow_id) is not None
    assert get_workflow_origin(workflow_id) == prior_origin
    assert get_recipe_catalogs() == prior_catalogs


def test_provenance_roots_are_ordered_and_relative_to_config(tmp_path: Path):
    config = load_config(
        _write_config(
            tmp_path,
            provenance={
                "search": [
                    {"id": "site", "root": "provenance/site"},
                    {"id": "project", "root": "provenance/project"},
                ],
                "write": {"id": "project", "root": "provenance/project"},
            },
        )
    )

    assert [repository.id for repository in config.provenance.search] == [
        "site",
        "project",
    ]
    assert config.provenance.search[0].root == (tmp_path / "provenance/site").resolve()
    assert config.provenance.write is not None
    assert config.provenance.write.root == (
        tmp_path / "provenance/project"
    ).resolve()

    public = list_provenance_repositories(authenticated=False)
    assert public == {
        "search": [{"id": "site", "order": 0}, {"id": "project", "order": 1}],
        "write": {"id": "project"},
    }
    authenticated = list_provenance_repositories(authenticated=True)
    assert authenticated["search"][0]["root"] == str(
        (tmp_path / "provenance/site").resolve()
    )
    assert authenticated["write"]["root"] == str(
        (tmp_path / "provenance/project").resolve()
    )


def test_duplicate_provenance_search_id_is_rejected(tmp_path: Path):
    with pytest.raises(ValueError, match="duplicate provenance search repository id"):
        load_config(
            _write_config(
                tmp_path,
                provenance={
                    "search": [
                        {"id": "site", "root": "one"},
                        {"id": "site", "root": "two"},
                    ]
                },
            )
        )


def test_provenance_write_id_cannot_point_to_different_root(tmp_path: Path):
    with pytest.raises(ValueError, match="different root"):
        load_config(
            _write_config(
                tmp_path,
                provenance={
                    "search": [{"id": "project", "root": "read-root"}],
                    "write": {"id": "project", "root": "write-root"},
                },
            )
        )


def _load_two_source_logical_recipe(
    tmp_path: Path,
    *,
    first_version: str = "1.0",
    second_version: str = "2.0",
    first_applications: list[str] | None = None,
    second_applications: list[str] | None = None,
    recipe_resolution: dict | None = None,
) -> tuple[str, str, str]:
    logical_id = "bulk-rna"
    first, first_id = _write_catalog(
        tmp_path / "first",
        catalog_id="first",
        namespace="example.first",
        logical_id=logical_id,
        recipe_version=first_version,
        applications=first_applications,
        image="example/star-suite@sha256:first",
    )
    second, second_id = _write_catalog(
        tmp_path / "second",
        catalog_id="second",
        namespace="example.second",
        logical_id=logical_id,
        recipe_version=second_version,
        applications=second_applications,
        image="example/star-suite@sha256:second",
    )
    load_config(
        _write_config(
            tmp_path,
            recipe_catalogs=[
                {"manifest": str(first.relative_to(tmp_path))},
                {"manifest": str(second.relative_to(tmp_path))},
            ],
            recipe_resolution=recipe_resolution,
        )
    )
    return logical_id, first_id, second_id


def test_catalog_discovery_redacts_roots_and_reports_conflicts(tmp_path: Path):
    logical_id, _, _ = _load_two_source_logical_recipe(tmp_path)

    public = list_recipe_catalogs(authenticated=False)
    external = public["catalogs"][1]
    assert external["id"] == "first"
    assert public["catalogs"][0]["trust"] == "trusted"
    assert external["trust"] == "untrusted"
    assert "root" not in external
    assert "git" not in external

    authenticated = describe_recipe_catalog("first", authenticated=True)
    assert authenticated["catalog"]["root"] == str((tmp_path / "first").resolve())
    assert authenticated["recipes"][0]["logical_id"] == logical_id

    conflicts = list_recipe_conflicts(authenticated=True)
    assert [item["logical_id"] for item in conflicts["conflicts"]] == [logical_id]


def test_resolution_supports_separate_prompt_newest_and_user_selection(
    tmp_path: Path,
):
    logical_id, first_id, second_id = _load_two_source_logical_recipe(
        tmp_path,
        recipe_resolution={
            "default_policy": "keep_separate",
            "applications": {
                "workbench": "prompt",
                "batch": "prefer_newest",
            },
        },
    )

    separate = resolve_recipe(logical_id, authenticated=True)
    assert separate["status"] == "separate"
    assert separate["selected"] is None

    prompt = resolve_recipe(
        logical_id, application="workbench", authenticated=True
    )
    assert prompt["status"] == "selection_required"

    newest = resolve_recipe(logical_id, application="batch", authenticated=True)
    assert newest["status"] == "selected"
    assert newest["selected"]["workflow_id"] == second_id

    chosen = resolve_recipe(
        logical_id,
        application="workbench",
        selected_workflow_id=first_id,
        authenticated=True,
    )
    assert chosen["reason"] == "user_selection"
    assert chosen["selected"]["workflow_id"] == first_id

    explicit = resolve_recipe(f"second::{second_id}", authenticated=True)
    assert explicit["reason"] == "explicit_source"


def test_prefer_newest_tie_uses_later_catalog_order(tmp_path: Path):
    logical_id, _, second_id = _load_two_source_logical_recipe(
        tmp_path, first_version="2.0", second_version="2.0"
    )
    resolved = resolve_recipe(
        logical_id, policy="prefer_newest", authenticated=True
    )
    assert resolved["selected"]["workflow_id"] == second_id
    assert resolved["reason"] == "highest_version_then_later_catalog"


def test_application_compatibility_keeps_variants_separate(tmp_path: Path):
    logical_id, first_id, second_id = _load_two_source_logical_recipe(
        tmp_path,
        first_applications=["workbench"],
        second_applications=["temporal"],
    )

    workbench = list_recipe_candidates(
        logical_id=logical_id, application="workbench", authenticated=True
    )
    temporal = list_recipe_candidates(
        logical_id=logical_id, application="temporal", authenticated=True
    )
    assert [item["workflow_id"] for item in workbench["groups"][0]["candidates"]] == [
        first_id
    ]
    assert [item["workflow_id"] for item in temporal["groups"][0]["candidates"]] == [
        second_id
    ]


def test_unknown_relationship_and_inheritance_cycle_fail_closed(tmp_path: Path):
    missing, _ = _write_catalog(
        tmp_path / "missing",
        catalog_id="missing",
        namespace="example.missing",
        extends="example.unknown/demo",
    )
    with pytest.raises(ValueError, match="references unknown workflow"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(missing.relative_to(tmp_path))}],
            )
        )

    first, first_id = _write_catalog(
        tmp_path / "cycle-first",
        catalog_id="cycle-first",
        namespace="example.cycle_first",
        extends="example.cycle_second/demo",
    )
    second, second_id = _write_catalog(
        tmp_path / "cycle-second",
        catalog_id="cycle-second",
        namespace="example.cycle_second",
        extends=first_id,
    )
    assert second_id == "example.cycle_second/demo"
    with pytest.raises(ValueError, match="recipe inheritance cycle"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[
                    {"manifest": str(first.relative_to(tmp_path))},
                    {"manifest": str(second.relative_to(tmp_path))},
                ],
            )
        )


def test_invalid_recipe_version_fails_during_catalog_load(tmp_path: Path):
    manifest, workflow_id = _write_catalog(
        tmp_path / "bad-version",
        catalog_id="bad-version",
        namespace="example.bad_version",
        recipe_version="not a version",
    )
    with pytest.raises(ValueError, match=f"workflow '{workflow_id}'.*PEP 440"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            )
        )


def test_public_recipe_cannot_expose_private_lineage(tmp_path: Path):
    parent, parent_id = _write_catalog(
        tmp_path / "private-parent",
        catalog_id="private-parent",
        namespace="example.private_parent",
    )
    raw_parent = yaml.safe_load(parent.read_text())
    raw_parent["workflows"][0]["visibility"] = "private"
    parent.write_text(yaml.safe_dump(raw_parent, sort_keys=False))
    child, _ = _write_catalog(
        tmp_path / "public-child",
        catalog_id="public-child",
        namespace="example.public_child",
        extends=parent_id,
    )

    with pytest.raises(ValueError, match="public workflow.*private workflow"):
        load_config(
            _write_config(
                tmp_path,
                recipe_catalogs=[
                    {"manifest": str(parent.relative_to(tmp_path))},
                    {"manifest": str(child.relative_to(tmp_path))},
                ],
            )
        )


def test_recipe_lock_is_deterministic_and_hashes_helpers(tmp_path: Path):
    logical_id, _, second_id = _load_two_source_logical_recipe(tmp_path)

    first_lock = build_recipe_lock(
        logical_id,
        {},
        policy="prefer_newest",
        application="batch",
    )
    second_lock = build_recipe_lock(
        logical_id,
        {},
        policy="prefer_newest",
        application="batch",
    )
    assert first_lock == second_lock
    assert first_lock["schema"] == "biodepot.recipe_lock/v1"
    assert "source_revision" in first_lock["suite"]
    assert first_lock["recipe"]["workflow_id"] == second_id
    roles = {
        source["role"]
        for source in first_lock["recipe"]["lineage"][-1]["sources"]
    }
    assert roles == {"workflow_schema", "entry", "helper"}
    assert first_lock["rendered"]["argv"][0] == "scripts/run.sh"

    (tmp_path / "second" / "helpers" / "helper.py").write_text(
        "print('changed')\n"
    )
    changed_lock = build_recipe_lock(
        logical_id,
        {},
        policy="prefer_newest",
        application="batch",
    )
    assert changed_lock["digest"] != first_lock["digest"]


def test_recipe_bundle_is_atomic_and_emits_workbench_temporal_packet(
    tmp_path: Path,
):
    logical_id, _, second_id = _load_two_source_logical_recipe(tmp_path)
    bundle = tmp_path / "bundle"
    result = create_recipe_bundle(
        bundle,
        logical_id,
        {},
        run_id="paper-run-001",
        application="workbench",
        policy="prefer_newest",
    )

    assert result["workflow_id"] == second_id
    assert (bundle / "recipe.lock.json").is_file()
    assert (bundle / "catalogs/second/scripts/run.sh").is_file()
    assert (bundle / "catalogs/second/helpers/helper.py").is_file()
    packet = json.loads(
        (bundle / "manifests/resolved_workflow_v1.json").read_text()
    )
    assert packet["schema"] == "biodepot.resolved_workflow/v1"
    assert packet["resolved_workflow"]["run_id"] == "paper-run-001"
    node = packet["resolved_workflow"]["nodes"]["1"]
    assert node["image_name"] == "example/star-suite@sha256:second"
    assert node["launch"]["command"][0].startswith(
        "catalogs/second/scripts/run.sh"
    )
    assert packet["resolved_workflow"]["links"] == []

    with pytest.raises(FileExistsError, match="already exists"):
        create_recipe_bundle(
            bundle,
            logical_id,
            {},
            run_id="paper-run-002",
            policy="prefer_newest",
        )


def test_recipe_cli_lists_and_resolves_catalogs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    logical_id, _, second_id = _load_two_source_logical_recipe(tmp_path)
    config_path = config_module.get_config_path()
    assert config_path is not None

    assert recipe_cli_main(
        ["--config", str(config_path), "catalog", "list"]
    ) == 0
    catalog_payload = json.loads(capsys.readouterr().out)
    assert [item["id"] for item in catalog_payload["catalogs"]] == [
        BUILTIN_CATALOG_ID,
        "first",
        "second",
    ]

    assert recipe_cli_main(
        ["--config", str(config_path), "provenance", "list"]
    ) == 0
    provenance_payload = json.loads(capsys.readouterr().out)
    assert provenance_payload == {"search": [], "write": None}

    assert recipe_cli_main(
        [
            "--config",
            str(config_path),
            "recipe",
            "resolve",
            logical_id,
            "--policy",
            "prefer_newest",
        ]
    ) == 0
    resolution = json.loads(capsys.readouterr().out)
    assert resolution["selected"]["workflow_id"] == second_id

    bundle = tmp_path / "cli-bundle"
    assert recipe_cli_main(
        [
            "--config",
            str(config_path),
            "recipe",
            "bundle",
            logical_id,
            "--policy",
            "prefer_newest",
            "--run-id",
            "cli-contract",
            "--output-dir",
            str(bundle),
        ]
    ) == 0
    bundle_result = json.loads(capsys.readouterr().out)
    assert bundle_result["workflow_id"] == second_id
    assert (bundle / "manifests/resolved_workflow_v1.json").is_file()


def test_recipe_mcp_discovery_redacts_and_lock_requires_auth(tmp_path: Path):
    from mcp_server.app import build_recipe_lock as lock_tool
    from mcp_server.app import list_recipe_catalogs as catalog_tool

    manifest, workflow_id = _write_catalog(
        tmp_path / "external",
        catalog_id="external",
        namespace="example.external",
        image="example/star-suite@sha256:external",
    )
    load_config(
        _write_config(
            tmp_path,
            recipe_catalogs=[{"manifest": str(manifest.relative_to(tmp_path))}],
            auth_token="secret",
        )
    )

    public = catalog_tool.fn()
    assert "root" not in public["catalogs"][1]
    authenticated = catalog_tool.fn(auth_token="secret")
    assert authenticated["catalogs"][1]["root"] == str(
        (tmp_path / "external").resolve()
    )

    denied = lock_tool.fn(reference=workflow_id, params={})
    assert denied["code"] == "AUTH_FAILED"
    locked = lock_tool.fn(
        reference=workflow_id, params={}, auth_token="secret"
    )
    assert locked["schema"] == "biodepot.recipe_lock/v1"
