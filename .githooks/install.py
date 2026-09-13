#!/usr/bin/env python3
"""Install attribution hooks in the shared Git hooks directory, preserving existing hooks."""
from pathlib import Path
import shutil
import subprocess

source = Path(__file__).resolve().parent
hooks = Path(subprocess.check_output(
    ["git", "rev-parse", "--path-format=absolute", "--git-path", "hooks"], text=True
).strip())
hooks.mkdir(parents=True, exist_ok=True)
names = ["check-attribution.py", "commit-msg", "pre-push"]
for name in names:
    target = hooks / name
    if target.exists() and target.read_bytes() != (source / name).read_bytes():
        raise SystemExit("Existing hook preserved; integrate the attribution check manually: " + str(target))
for name in names:
    shutil.copy2(source / name, hooks / name)
    (hooks / name).chmod(0o755)
print("Installed attribution checks in " + str(hooks))
