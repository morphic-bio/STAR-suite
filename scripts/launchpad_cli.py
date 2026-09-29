#!/usr/bin/env python3
"""Start Launchpad from a source checkout or an installed STAR Suite package."""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Launch the STAR Suite browser UI.")
    parser.add_argument("--setup", action="store_true", help="Create a user-owned Python environment and install UI dependencies, then exit.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--config", type=Path, help="Optional site configuration.")
    parser.add_argument("--env-dir", type=Path, default=Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local/share")) / "star-suite/launchpad-venv")
    args = parser.parse_args()
    script = Path(__file__).resolve()
    source = script.parent.name == "scripts" and (script.parents[1] / "mcp_server").is_dir()
    prefix = script.parents[1]
    root = prefix if source else prefix / "share/star-suite/launchpad"
    if not (root / "mcp_server/config.yaml").is_file():
        parser.error(f"Launchpad assets are missing from {root}")
    venv = args.env_dir.expanduser().resolve()
    python = venv / "bin/python"
    if args.setup:
        subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
        subprocess.run([str(python), "-m", "pip", "install", "-r", str(root / "mcp_server/requirements.txt")], check=True)
        print("Launchpad dependencies installed. Run this command again without --setup.")
        return 0
    override = os.environ.get("STAR_SUITE_LAUNCHPAD_PYTHON")
    interpreter = override or (str(python) if python.is_file() else sys.executable)
    ready = subprocess.run([interpreter, "-c", "import fastmcp, yaml, uvicorn, multipart"],
                           capture_output=True, text=True)
    if ready.returncode:
        print("Launchpad dependencies are missing. Run star-suite-launchpad --setup (or this source script with --setup).", file=sys.stderr)
        return 2
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    if not source:
        env["PATH"] = str(prefix / "bin") + os.pathsep + env.get("PATH", "")
    config = (args.config.expanduser().resolve() if args.config else root / "mcp_server/config.yaml")
    print(f"Launchpad: http://{args.host}:{args.port}/launchpad/", flush=True)
    os.chdir(root)
    os.execvpe(interpreter, [interpreter, "-m", "mcp_server.app", "--host", args.host,
                           "--port", str(args.port), "--config", str(config)], env)


if __name__ == "__main__":
    raise SystemExit(main())
