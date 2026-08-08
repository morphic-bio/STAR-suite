"""Command-line discovery and locking for STAR Suite recipe catalogs."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

from .config import load_config
from .tools.recipes import (
    build_recipe_lock,
    create_recipe_bundle,
    describe_recipe_catalog,
    list_recipe_candidates,
    list_recipe_catalogs,
    list_recipe_conflicts,
    list_provenance_repositories,
    resolve_recipe,
)
from .tools.workflows import (
    describe_workflow,
    get_workflow_parameter_schema,
    render_workflow_command,
    validate_workflow_parameters,
)


def _parameters(args: argparse.Namespace) -> dict[str, Any]:
    if getattr(args, "params_json", None) and getattr(args, "params_file", None):
        raise ValueError("Use only one of --params-json and --params-file")
    if getattr(args, "params_file", None):
        value = json.loads(Path(args.params_file).read_text(encoding="utf-8"))
    elif getattr(args, "params_json", None):
        value = json.loads(args.params_json)
    else:
        value = {}
    if not isinstance(value, dict):
        raise ValueError("Recipe parameters must be a JSON object")
    return value


def _write_json(path: Path, payload: Any) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.tmp-", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as output:
            json.dump(payload, output, sort_keys=True, indent=2)
            output.write("\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _resolution_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--application", default="cli")
    parser.add_argument(
        "--policy", choices=("keep_separate", "prompt", "prefer_newest")
    )
    parser.add_argument(
        "--select",
        dest="selected_workflow_id",
        help="Candidate workflow id or catalog-id::workflow-id source reference",
    )


def _parameter_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--params-json", help="Parameters as a JSON object")
    parser.add_argument("--params-file", help="Path to a JSON parameter object")
    parser.add_argument(
        "--check-paths", action="store_true", help="Validate input path existence"
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="star-suite-recipes")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("config.yaml"),
        help="MCP configuration containing the catalog stack",
    )
    commands = parser.add_subparsers(dest="scope", required=True)

    catalogs = commands.add_parser("catalog", help="Inspect recipe catalogs")
    catalog_commands = catalogs.add_subparsers(dest="action", required=True)
    catalog_commands.add_parser("list")
    catalog_describe = catalog_commands.add_parser("describe")
    catalog_describe.add_argument("catalog_id")

    provenance = commands.add_parser(
        "provenance", help="Inspect provenance repository hierarchy"
    )
    provenance_commands = provenance.add_subparsers(dest="action", required=True)
    provenance_commands.add_parser("list")

    recipes = commands.add_parser("recipe", help="Resolve and lock recipes")
    recipe_commands = recipes.add_subparsers(dest="action", required=True)
    recipe_list = recipe_commands.add_parser("list")
    recipe_list.add_argument("--logical-id")
    recipe_list.add_argument("--application")

    recipe_conflicts = recipe_commands.add_parser("conflicts")
    recipe_conflicts.add_argument("--application")

    for action in ("resolve", "show", "validate", "render", "lock", "bundle"):
        subparser = recipe_commands.add_parser(action)
        subparser.add_argument("reference")
        _resolution_options(subparser)
        if action in ("validate", "render", "lock", "bundle"):
            _parameter_options(subparser)
        if action == "lock":
            subparser.add_argument("--output", type=Path)
        if action == "bundle":
            subparser.add_argument("--output-dir", type=Path, required=True)
            subparser.add_argument("--run-id", required=True)
            subparser.add_argument("--image")
            subparser.add_argument("--cores", type=int, default=1)
            subparser.add_argument("--mem-mb", type=int, default=2048)
            subparser.add_argument("--gpus", type=int, default=0)

    return parser


def _selected_resolution(args: argparse.Namespace) -> dict[str, Any]:
    return resolve_recipe(
        args.reference,
        application=args.application,
        policy=args.policy,
        selected_workflow_id=args.selected_workflow_id,
        authenticated=True,
    )


def _dispatch(args: argparse.Namespace) -> dict[str, Any]:
    if args.scope == "catalog":
        if args.action == "list":
            return list_recipe_catalogs(authenticated=True)
        return describe_recipe_catalog(args.catalog_id, authenticated=True)

    if args.scope == "provenance":
        return list_provenance_repositories(authenticated=True)

    if args.action == "list":
        return list_recipe_candidates(
            logical_id=args.logical_id,
            application=args.application,
            authenticated=True,
        )
    if args.action == "conflicts":
        return list_recipe_conflicts(
            application=args.application, authenticated=True
        )

    resolution = _selected_resolution(args)
    if args.action == "resolve":
        return resolution
    if args.action == "show":
        if resolution["selected"] is None:
            return resolution
        workflow_id = resolution["selected"]["workflow_id"]
        return {
            "resolution": resolution,
            "workflow": describe_workflow(workflow_id, authenticated=True).model_dump(),
            "parameters": get_workflow_parameter_schema(
                workflow_id, authenticated=True
            ).model_dump(),
        }

    if resolution["selected"] is None:
        raise ValueError(
            f"Recipe '{args.reference}' is unresolved ({resolution['status']}); "
            "use --select"
        )
    workflow_id = resolution["selected"]["workflow_id"]
    params = _parameters(args)
    if args.action == "validate":
        return validate_workflow_parameters(
            workflow_id, params, check_paths=args.check_paths
        ).model_dump()
    if args.action == "render":
        validation = validate_workflow_parameters(
            workflow_id, params, check_paths=args.check_paths
        )
        if not validation.valid:
            return validation.model_dump()
        return render_workflow_command(
            workflow_id, validation.normalized_params
        ).model_dump()
    if args.action == "lock":
        result = build_recipe_lock(
            args.reference,
            params,
            application=args.application,
            policy=args.policy,
            selected_workflow_id=args.selected_workflow_id,
            authenticated=True,
            check_paths=args.check_paths,
        )
        if args.output:
            _write_json(args.output, result)
            return {
                "output": str(args.output.resolve()),
                "digest": result["digest"],
                "workflow_id": result["recipe"]["workflow_id"],
            }
        return result
    if args.action == "bundle":
        return create_recipe_bundle(
            args.output_dir,
            args.reference,
            params,
            run_id=args.run_id,
            application=args.application,
            policy=args.policy,
            selected_workflow_id=args.selected_workflow_id,
            image=args.image,
            resources={
                "cores": args.cores,
                "mem_mb": args.mem_mb,
                "gpus": args.gpus,
            },
            check_paths=args.check_paths,
        )
    raise ValueError(f"Unsupported command: {args.scope} {args.action}")


def main(argv: Optional[list[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        load_config(args.config)
        result = _dispatch(args)
    except Exception as exc:
        json.dump(
            {"error": {"type": type(exc).__name__, "message": str(exc)}},
            sys.stderr,
            sort_keys=True,
        )
        sys.stderr.write("\n")
        return 2
    json.dump(result, sys.stdout, sort_keys=True, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
