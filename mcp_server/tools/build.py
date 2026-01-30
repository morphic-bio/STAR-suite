"""
Build tools for STAR-suite MCP server.

Provides safe, clean builds with change detection to prevent stale binary issues.
"""

import os
import subprocess
import hashlib
import json
import time
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

# Build state file location
BUILD_STATE_FILE = ".star_build_state.json"


def get_source_hash(source_dir: Path, extensions: List[str] = None) -> str:
    """
    Compute a hash of all source files in a directory.
    Used to detect when sources have changed since last build.
    """
    if extensions is None:
        extensions = ['.cpp', '.c', '.h', '.hpp']
    
    hasher = hashlib.sha256()
    
    # Sort files for deterministic hashing
    source_files = []
    for ext in extensions:
        source_files.extend(sorted(source_dir.rglob(f'*{ext}')))
    
    for src_file in sorted(source_files):
        try:
            # Hash file path and mtime for speed (avoid reading huge files)
            stat = src_file.stat()
            file_info = f"{src_file}:{stat.st_mtime}:{stat.st_size}"
            hasher.update(file_info.encode())
        except (OSError, IOError):
            continue
    
    return hasher.hexdigest()[:16]


def load_build_state(repo_root: Path) -> Dict[str, Any]:
    """Load the build state from disk."""
    state_file = repo_root / BUILD_STATE_FILE
    if state_file.exists():
        try:
            with open(state_file) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {}


def save_build_state(repo_root: Path, state: Dict[str, Any]):
    """Save the build state to disk."""
    state_file = repo_root / BUILD_STATE_FILE
    with open(state_file, 'w') as f:
        json.dump(state, f, indent=2)


def needs_rebuild(repo_root: Path, target: str) -> tuple[bool, str]:
    """
    Check if a target needs rebuilding based on source changes.
    
    Returns:
        (needs_rebuild: bool, reason: str)
    """
    state = load_build_state(repo_root)
    
    # Define source directories for each target
    source_dirs = {
        'core': [
            repo_root / 'core' / 'legacy' / 'source',
            repo_root / 'flex' / 'source',
            repo_root / 'slam' / 'source',
            repo_root / 'core' / 'features',
        ],
        'flex': [
            repo_root / 'flex' / 'source',
            repo_root / 'core' / 'features' / 'libscrna',
        ],
        'slam': [
            repo_root / 'slam' / 'source',
        ],
        'feature-tools': [
            repo_root / 'core' / 'features' / 'feature_barcodes',
            repo_root / 'core' / 'features' / 'process_features',
        ],
    }
    
    if target not in source_dirs:
        return True, f"Unknown target '{target}'"
    
    # Compute current source hash
    current_hash = ""
    for src_dir in source_dirs[target]:
        if src_dir.exists():
            current_hash += get_source_hash(src_dir)
    
    # Check if binary exists
    binary_paths = {
        'core': repo_root / 'core' / 'legacy' / 'source' / 'STAR',
        'flex': repo_root / 'flex' / 'source' / 'libflex' / 'libflex.a',
        'slam': repo_root / 'slam' / 'source' / 'slam_tools',
        'feature-tools': repo_root / 'core' / 'features' / 'feature_barcodes' / 'assignBarcodes',
    }
    
    binary = binary_paths.get(target)
    if binary and not binary.exists():
        return True, f"Binary not found: {binary}"
    
    # Check against saved state
    target_state = state.get(target, {})
    saved_hash = target_state.get('source_hash', '')
    
    if saved_hash != current_hash:
        return True, "Source files changed since last build"
    
    return False, "Build is up to date"


def build_target(
    repo_root: Path,
    target: str,
    clean: bool = False,
    force: bool = False,
    parallel: int = None,
) -> Dict[str, Any]:
    """
    Build a STAR-suite target with proper state tracking.
    
    Args:
        repo_root: Path to repository root
        target: Build target ('core', 'flex', 'slam', 'feature-tools')
        clean: If True, run make clean first
        force: If True, rebuild even if up to date
        parallel: Number of parallel jobs (default: nproc)
    
    Returns:
        Dict with build result information
    """
    result = {
        'target': target,
        'success': False,
        'clean': clean,
        'forced': force,
        'start_time': datetime.now().isoformat(),
        'duration_seconds': 0,
        'output': '',
        'error': '',
    }
    
    start_time = time.time()
    
    # Check if rebuild needed
    if not force and not clean:
        needs_build, reason = needs_rebuild(repo_root, target)
        if not needs_build:
            result['success'] = True
            result['skipped'] = True
            result['reason'] = reason
            result['duration_seconds'] = time.time() - start_time
            return result
    
    # Determine make target and working directory
    make_configs = {
        'core': {
            'target': 'STAR',
            'cwd': repo_root / 'core' / 'legacy' / 'source',
        },
        'flex': {
            'target': 'flex',
            'cwd': repo_root,
        },
        'slam': {
            'target': 'slam',
            'cwd': repo_root,
        },
        'feature-tools': {
            'target': 'feature-barcodes-tools',
            'cwd': repo_root,
        },
    }
    
    if target not in make_configs:
        result['error'] = f"Unknown target: {target}"
        return result
    
    config = make_configs[target]
    cwd = config['cwd']
    make_target = config['target']
    
    # Determine parallel jobs
    if parallel is None:
        try:
            parallel = os.cpu_count() or 4
        except:
            parallel = 4
    
    try:
        # Run make clean if requested
        if clean:
            clean_cmd = ['make', 'clean']
            result['clean_output'] = ''
            proc = subprocess.run(
                clean_cmd,
                cwd=cwd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            result['clean_output'] = proc.stdout + proc.stderr
            if proc.returncode != 0:
                result['error'] = f"make clean failed: {proc.stderr}"
                result['duration_seconds'] = time.time() - start_time
                return result
        
        # Run make build
        build_cmd = ['make', f'-j{parallel}', make_target]
        proc = subprocess.run(
            build_cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 min timeout
        )
        
        result['output'] = proc.stdout
        result['error'] = proc.stderr
        result['returncode'] = proc.returncode
        
        if proc.returncode == 0:
            result['success'] = True
            
            # Update build state
            state = load_build_state(repo_root)
            
            # Compute new source hash
            source_dirs_map = {
                'core': [
                    repo_root / 'core' / 'legacy' / 'source',
                    repo_root / 'flex' / 'source',
                    repo_root / 'slam' / 'source',
                    repo_root / 'core' / 'features',
                ],
                'flex': [
                    repo_root / 'flex' / 'source',
                ],
                'slam': [
                    repo_root / 'slam' / 'source',
                ],
                'feature-tools': [
                    repo_root / 'core' / 'features' / 'feature_barcodes',
                    repo_root / 'core' / 'features' / 'process_features',
                ],
            }
            
            current_hash = ""
            for src_dir in source_dirs_map.get(target, []):
                if src_dir.exists():
                    current_hash += get_source_hash(src_dir)
            
            state[target] = {
                'source_hash': current_hash,
                'build_time': datetime.now().isoformat(),
                'clean_build': clean,
            }
            save_build_state(repo_root, state)
        
    except subprocess.TimeoutExpired:
        result['error'] = "Build timed out after 30 minutes"
    except Exception as e:
        result['error'] = str(e)
    
    result['duration_seconds'] = time.time() - start_time
    return result


def ensure_clean_build(repo_root: Path, target: str = 'core') -> Dict[str, Any]:
    """
    Ensure a clean build exists. Always does make clean first.
    
    This is the recommended function to call before running tests
    to prevent stale binary issues.
    """
    return build_target(repo_root, target, clean=True, force=True)


def smart_build(repo_root: Path, target: str = 'core') -> Dict[str, Any]:
    """
    Smart build: only rebuild if sources changed.
    
    Use this for incremental development. Use ensure_clean_build()
    before running official test suites.
    """
    return build_target(repo_root, target, clean=False, force=False)


# MCP tool registration
def register_build_tools(app):
    """Register build tools with the MCP server."""
    
    @app.tool()
    async def build_star(
        target: str = "core",
        clean: bool = False,
        force: bool = False,
    ) -> Dict[str, Any]:
        """
        Build STAR-suite binaries.
        
        Args:
            target: Build target - 'core' (STAR), 'flex', 'slam', or 'feature-tools'
            clean: Run 'make clean' before building (recommended for test runs)
            force: Force rebuild even if sources haven't changed
        
        Returns:
            Build result with success status, output, and timing
        
        IMPORTANT: Always use clean=True before running test suites to
        prevent stale binary issues that can cause segfaults.
        """
        from ..config import get_config
        config = get_config()
        repo_root = Path(config['paths']['repo_root'])
        
        return build_target(repo_root, target, clean=clean, force=force)
    
    @app.tool()
    async def check_build_status(target: str = "core") -> Dict[str, Any]:
        """
        Check if a build target needs rebuilding.
        
        Args:
            target: Build target to check
        
        Returns:
            Status including whether rebuild is needed and why
        """
        from ..config import get_config
        config = get_config()
        repo_root = Path(config['paths']['repo_root'])
        
        needs_build, reason = needs_rebuild(repo_root, target)
        state = load_build_state(repo_root)
        target_state = state.get(target, {})
        
        return {
            'target': target,
            'needs_rebuild': needs_build,
            'reason': reason,
            'last_build': target_state.get('build_time'),
            'was_clean_build': target_state.get('clean_build'),
        }
    
    @app.tool()
    async def ensure_fresh_build(target: str = "core") -> Dict[str, Any]:
        """
        Ensure a fresh, clean build exists (recommended before test suites).
        
        This always runs 'make clean' first to prevent stale binary issues.
        Use this before running any test suite.
        
        Args:
            target: Build target - 'core', 'flex', 'slam', or 'feature-tools'
        
        Returns:
            Build result
        """
        from ..config import get_config
        config = get_config()
        repo_root = Path(config['paths']['repo_root'])
        
        return ensure_clean_build(repo_root, target)
