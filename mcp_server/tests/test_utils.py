"""Tests for utility functions."""

import os
import re
from pathlib import Path

import pytest

from mcp_server.config import load_config
from mcp_server.tools.utils import (
    generate_run_id,
    get_directory_size,
    get_disk_space_gb,
    is_path_allowed,
    validate_path,
)
import mcp_server.config as config_module


@pytest.fixture
def setup_config(sample_config_file: Path):
    """Load config for tests."""
    config = load_config(sample_config_file)
    yield config
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


class TestGenerateRunId:
    """Tests for run ID generation."""

    def test_format(self):
        """Test that run ID has correct format."""
        run_id = generate_run_id()
        # Format: YYYYMMDD_HHMMSS_<rand4>
        pattern = r"^\d{8}_\d{6}_[a-z0-9]{4}$"
        assert re.match(pattern, run_id), f"Run ID {run_id} doesn't match pattern"

    def test_uniqueness(self):
        """Test that run IDs are unique."""
        ids = [generate_run_id() for _ in range(100)]
        assert len(set(ids)) == len(ids), "Run IDs should be unique"


class TestIsPathAllowed:
    """Tests for path validation."""

    def test_allowed_path(self, setup_config, temp_dir: Path):
        """Test that paths within trusted roots are allowed."""
        assert is_path_allowed(temp_dir / "subdir")

    def test_disallowed_path(self, setup_config):
        """Test that paths outside trusted roots are rejected."""
        # /etc is not in trusted roots
        assert not is_path_allowed("/etc/passwd")

    def test_explicit_trusted_roots(self, temp_dir: Path):
        """Test with explicit trusted roots parameter."""
        trusted = [str(temp_dir)]
        assert is_path_allowed(temp_dir / "foo", trusted_roots=trusted)
        assert not is_path_allowed("/other/path", trusted_roots=trusted)

    def test_symlink_resolution(self, setup_config, temp_dir: Path):
        """Test that symlinks are resolved."""
        # Create a symlink within temp_dir pointing to temp_dir/target
        target = temp_dir / "target"
        target.mkdir()
        link = temp_dir / "link"
        link.symlink_to(target)

        # The symlink should resolve to a path within trusted roots
        assert is_path_allowed(link)


class TestValidatePath:
    """Tests for path validation with criteria."""

    def test_must_exist_pass(self, setup_config, temp_dir: Path):
        """Test must_exist with existing path."""
        existing = temp_dir / "existing"
        existing.touch()
        valid, error = validate_path(existing, must_exist=True)
        assert valid
        assert error is None

    def test_must_exist_fail(self, setup_config, temp_dir: Path):
        """Test must_exist with non-existing path."""
        valid, error = validate_path(temp_dir / "nonexistent", must_exist=True)
        assert not valid
        assert "does not exist" in error

    def test_must_be_file(self, setup_config, temp_dir: Path):
        """Test must_be_file validation."""
        file_path = temp_dir / "file.txt"
        file_path.touch()
        valid, error = validate_path(file_path, must_be_file=True)
        assert valid

        valid, error = validate_path(temp_dir, must_be_file=True)
        assert not valid

    def test_must_be_dir(self, setup_config, temp_dir: Path):
        """Test must_be_dir validation."""
        valid, error = validate_path(temp_dir, must_be_dir=True)
        assert valid

        file_path = temp_dir / "file.txt"
        file_path.touch()
        valid, error = validate_path(file_path, must_be_dir=True)
        assert not valid

    def test_outside_trusted_roots(self, setup_config):
        """Test that paths outside trusted roots fail."""
        valid, error = validate_path("/etc/passwd")
        assert not valid
        assert "outside trusted roots" in error

    def test_writable_with_nonexistent_ancestors(self, setup_config, temp_dir: Path):
        """Test that writable check finds first existing ancestor."""
        # Path with multiple non-existent ancestors
        deep_path = temp_dir / "a" / "b" / "c" / "d" / "output"
        valid, error = validate_path(deep_path, must_be_writable=True)
        # Should pass because temp_dir exists and is writable
        assert valid

    def test_writable_fails_if_ancestor_not_writable(self, setup_config, temp_dir: Path):
        """Test that writable check fails if ancestor is not writable."""
        # Create a read-only directory
        readonly_dir = temp_dir / "readonly"
        readonly_dir.mkdir()
        readonly_dir.chmod(0o555)

        try:
            deep_path = readonly_dir / "a" / "b" / "output"
            valid, error = validate_path(deep_path, must_be_writable=True)
            assert not valid
            assert "not writable" in error
        finally:
            # Restore permissions for cleanup
            readonly_dir.chmod(0o755)


class TestGetDiskSpaceGb:
    """Tests for disk space checking."""

    def test_returns_positive_value(self, temp_dir: Path):
        """Test that disk space check returns positive value."""
        space = get_disk_space_gb(temp_dir)
        assert space > 0

    def test_handles_nonexistent_path(self, temp_dir: Path):
        """Test that non-existent paths check parent."""
        space = get_disk_space_gb(temp_dir / "nonexistent" / "path")
        assert space > 0


class TestGetDirectorySize:
    """Tests for directory size calculation."""

    def test_empty_directory(self, temp_dir: Path):
        """Test size of empty directory."""
        empty_dir = temp_dir / "empty"
        empty_dir.mkdir()
        size = get_directory_size(empty_dir)
        assert size == 0

    def test_directory_with_files(self, temp_dir: Path):
        """Test size of directory with files."""
        dir_path = temp_dir / "with_files"
        dir_path.mkdir()
        (dir_path / "file1.txt").write_text("hello")
        (dir_path / "file2.txt").write_text("world")

        size = get_directory_size(dir_path)
        assert size == 10  # 5 + 5 bytes

    def test_nonexistent_directory(self, temp_dir: Path):
        """Test that non-existent directory returns 0."""
        size = get_directory_size(temp_dir / "nonexistent")
        assert size == 0
