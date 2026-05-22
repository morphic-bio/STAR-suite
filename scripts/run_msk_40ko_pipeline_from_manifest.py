#!/usr/bin/env python3
"""Compatibility launcher for the canonical Morphic recipe copy."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> int:
    os.environ.setdefault("STAR_SUITE_ROOT", str(Path(__file__).resolve().parents[1]))
    recipe_root = Path(os.environ.get("MORPHIC_RECIPES_ROOT", "/mnt/pikachu/morphic-recipes"))
    target = recipe_root / "scripts" / Path(__file__).name
    if not target.is_file():
        sys.stderr.write(f"ERROR: canonical Morphic recipe not found: {target}\n")
        sys.stderr.write("Set MORPHIC_RECIPES_ROOT=/path/to/morphic-recipes or run the recipe there directly.\n")
        return 1
    os.execv(sys.executable, [sys.executable, str(target), *sys.argv[1:]])
    return 127


if __name__ == "__main__":
    raise SystemExit(main())
