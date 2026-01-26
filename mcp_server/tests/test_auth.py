"""Tests for authentication."""

from pathlib import Path

import pytest
import yaml

from mcp_server.config import load_config
from mcp_server.app import check_auth
# Import the underlying functions, not the decorated tools
from mcp_server.app import list_datasets as list_datasets_tool
from mcp_server.app import find_docs as find_docs_tool
import mcp_server.config as config_module


# Access the underlying functions from the FunctionTool wrappers
def list_datasets(**kwargs):
    """Call the list_datasets tool function."""
    return list_datasets_tool.fn(**kwargs)


def find_docs(**kwargs):
    """Call the find_docs tool function."""
    return find_docs_tool.fn(**kwargs)


@pytest.fixture
def config_with_auth(temp_dir: Path) -> Path:
    """Create a config file with auth enabled."""
    config_dict = {
        "server": {
            "auth_token": "test-secret-token",
        },
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
        },
        "trusted_roots": [str(temp_dir)],
        "datasets": [],
        "scripts": [],
        "test_suites": [],
    }
    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)
    return config_path


@pytest.fixture
def config_without_auth(temp_dir: Path) -> Path:
    """Create a config file without auth."""
    config_dict = {
        "server": {
            "auth_token": "",  # Empty = no auth
        },
        "paths": {
            "repo_root": str(temp_dir),
            "artifact_log_root": str(temp_dir / "artifacts"),
        },
        "trusted_roots": [str(temp_dir)],
        "datasets": [],
        "scripts": [],
        "test_suites": [],
    }
    config_path = temp_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)
    return config_path


@pytest.fixture
def setup_auth_config(config_with_auth: Path):
    """Load config with auth enabled."""
    config = load_config(config_with_auth)
    yield config
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


@pytest.fixture
def setup_no_auth_config(config_without_auth: Path):
    """Load config without auth."""
    config = load_config(config_without_auth)
    yield config
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


class TestCheckAuth:
    """Tests for check_auth function."""

    def test_auth_passes_with_correct_token(self, setup_auth_config):
        """Test that correct token passes auth."""
        result = check_auth("test-secret-token")
        assert result is None

    def test_auth_fails_with_wrong_token(self, setup_auth_config):
        """Test that wrong token fails auth."""
        result = check_auth("wrong-token")
        assert result is not None
        assert result.code == "AUTH_FAILED"

    def test_auth_fails_with_no_token(self, setup_auth_config):
        """Test that missing token fails auth."""
        result = check_auth(None)
        assert result is not None
        assert result.code == "AUTH_FAILED"
        assert "required" in result.message.lower()

    def test_auth_skipped_when_not_configured(self, setup_no_auth_config):
        """Test that auth is skipped when no token configured."""
        result = check_auth(None)
        assert result is None

        result = check_auth("any-token")
        assert result is None


class TestToolsWithAuth:
    """Tests for tools enforcing auth."""

    def test_list_datasets_requires_auth_when_not_public(self, temp_dir: Path):
        """Test that list_datasets requires auth when public_discovery is False."""
        # Create config with auth and public_discovery=False
        config_dict = {
            "server": {
                "auth_token": "test-secret-token",
                "public_discovery": False,
            },
            "paths": {
                "repo_root": str(temp_dir),
                "artifact_log_root": str(temp_dir / "artifacts"),
            },
            "trusted_roots": [str(temp_dir)],
            "datasets": [],
            "scripts": [],
            "test_suites": [],
        }
        config_path = temp_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        load_config(config_path)

        try:
            # Without token - should fail
            result = list_datasets()
            assert result.get("error") is True
            assert result.get("code") == "AUTH_FAILED"

            # With correct token - should pass
            result = list_datasets(auth_token="test-secret-token")
            assert result.get("error") is not True
            assert "datasets" in result
        finally:
            config_module._config = None
            config_module._config_path = None
            config_module._config_loaded_at = None

    def test_discovery_public_by_default(self, setup_auth_config):
        """Test that discovery tools work without auth when public_discovery is True (default)."""
        # Modify config to have public_discovery=True
        config_module._config.server.public_discovery = True

        # Without token - should work because discovery is public
        result = list_datasets()
        assert result.get("error") is not True
        assert "datasets" in result

    def test_find_docs_requires_auth_when_not_public(self, temp_dir: Path):
        """Test that find_docs requires auth when public_discovery is False."""
        # Create docs dir
        (temp_dir / "docs").mkdir(exist_ok=True)

        config_dict = {
            "server": {
                "auth_token": "test-secret-token",
                "public_discovery": False,
            },
            "paths": {
                "repo_root": str(temp_dir),
                "artifact_log_root": str(temp_dir / "artifacts"),
            },
            "trusted_roots": [str(temp_dir)],
            "datasets": [],
            "scripts": [],
            "test_suites": [],
        }
        config_path = temp_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        load_config(config_path)

        try:
            # Without token - should fail
            result = find_docs(topic="test")
            assert result.get("error") is True
            assert result.get("code") == "AUTH_FAILED"

            # With correct token - should pass
            result = find_docs(topic="test", auth_token="test-secret-token")
            assert result.get("error") is not True
            assert "matches" in result
        finally:
            config_module._config = None
            config_module._config_path = None
            config_module._config_loaded_at = None

    def test_tools_work_without_auth_when_not_configured(self, setup_no_auth_config):
        """Test that tools work without token when auth not configured."""
        result = list_datasets()
        assert result.get("error") is not True
        assert "datasets" in result
