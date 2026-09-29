#!/usr/bin/env python3
"""Stage a relocatable Launchpad payload into a tarball or Debian install prefix."""
import argparse
from pathlib import Path
import shutil

import yaml

ROOT = Path(__file__).resolve().parents[2]


def stage(prefix):
    prefix = Path(prefix).resolve()
    data = prefix / "share/star-suite"
    payload = data / "launchpad"
    shutil.copytree(ROOT / "share/star-suite", data, dirs_exist_ok=True)
    shutil.copytree(ROOT / "mcp_server", payload / "mcp_server", dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache", "tests"))
    config_path = payload / "mcp_server/config.yaml"
    config = yaml.safe_load(config_path.read_text())
    # Installed packages expose executable workflows, not checkout-only tests.
    kinds = {w["id"]: yaml.safe_load((payload / w["schema_file"]).read_text()).get("kind")
             for w in config["workflows"]}
    config["workflows"] = [w for w in config["workflows"] if kinds[w["id"]] == "star_cli" or w["id"] == "morphic_multiome"]
    config["scripts"] = []
    config["test_suites"] = []
    config["required_binaries"] = []
    config["recipe_catalogs"] = [{"manifest": "../../catalogs/official/catalog.yaml", "trust": "trusted"}]
    config["provenance"] = {"search": [{"id": "official-evidence", "root": "../../evidence/official"}]}
    config["trusted_roots"].append("../../../../bin")
    for workflow in config["workflows"]:
        if kinds[workflow["id"]] == "star_cli":
            workflow["entry_script"] = "../../../bin/STAR"
        else:
            workflow["entry_script"] = "../catalogs/official/scripts/run_star_multiome_lane_smoke.sh"
        schema_path = payload / workflow["schema_file"]
        text = schema_path.read_text().replace("core/legacy/source/STAR", "../../../bin/STAR")
        text = text.replace("share/star-suite/catalogs/official/", "../catalogs/official/")
        schema_path.write_text(text)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    bindir = prefix / "bin"
    bindir.mkdir(parents=True, exist_ok=True)
    launcher = bindir / "star-suite-launchpad"
    shutil.copy2(ROOT / "scripts/launchpad_cli.py", launcher)
    launcher.chmod(0o755)
    print(f"Staged Launchpad: {payload}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-root", type=Path, required=True)
    stage(parser.parse_args().stage_root)


if __name__ == "__main__":
    main()
