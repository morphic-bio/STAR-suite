"""Tests for configuration loading."""

import os
from pathlib import Path

import pytest
import yaml

from mcp_server.config import (
    expand_env_vars,
    load_config,
    get_config,
    reload_config,
)
import mcp_server.config as config_module


class TestExpandEnvVars:
    """Tests for environment variable expansion."""

    def test_expand_simple_var(self, monkeypatch):
        """Test expanding a simple environment variable."""
        monkeypatch.setenv("TEST_VAR", "test_value")
        result = expand_env_vars("${TEST_VAR}")
        assert result == "test_value"

    def test_expand_missing_var(self):
        """Test that missing vars expand to empty string."""
        result = expand_env_vars("${NONEXISTENT_VAR_12345}")
        assert result == ""

    def test_expand_in_string(self, monkeypatch):
        """Test expanding variable within a string."""
        monkeypatch.setenv("NAME", "world")
        result = expand_env_vars("hello ${NAME}!")
        assert result == "hello world!"

    def test_expand_in_dict(self, monkeypatch):
        """Test expanding variables in a dictionary."""
        monkeypatch.setenv("HOST", "localhost")
        monkeypatch.setenv("PORT", "8080")
        data = {"server": {"host": "${HOST}", "port": "${PORT}"}}
        result = expand_env_vars(data)
        assert result == {"server": {"host": "localhost", "port": "8080"}}

    def test_expand_in_list(self, monkeypatch):
        """Test expanding variables in a list."""
        monkeypatch.setenv("PATH1", "/path/one")
        monkeypatch.setenv("PATH2", "/path/two")
        data = ["${PATH1}", "${PATH2}"]
        result = expand_env_vars(data)
        assert result == ["/path/one", "/path/two"]

    def test_no_expansion_needed(self):
        """Test that strings without vars are unchanged."""
        result = expand_env_vars("no vars here")
        assert result == "no vars here"


class TestLoadConfig:
    """Tests for config loading."""

    def test_load_valid_config(self, sample_config_file: Path):
        """Test loading a valid config file."""
        config = load_config(sample_config_file)
        assert config.server.host == "127.0.0.1"
        assert config.server.port == 9999
        assert len(config.datasets) == 1
        assert config.datasets[0].id == "test-dataset"

    def test_load_missing_file(self, temp_dir: Path):
        """Test that loading a missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_config(temp_dir / "nonexistent.yaml")

    def test_load_invalid_yaml(self, temp_dir: Path):
        """Test that invalid YAML raises an error."""
        bad_file = temp_dir / "bad.yaml"
        bad_file.write_text("this: is: not: valid: yaml: [")
        with pytest.raises(Exception):
            load_config(bad_file)

    def test_env_var_expansion_in_config(self, temp_dir: Path, monkeypatch):
        """Test that env vars are expanded when loading config."""
        monkeypatch.setenv("TEST_AUTH_TOKEN", "secret-token-123")

        config_dict = {
            "server": {"auth_token": "${TEST_AUTH_TOKEN}"},
            "paths": {
                "repo_root": str(temp_dir),
                "artifact_log_root": str(temp_dir / "artifacts"),
            },
            "trusted_roots": [str(temp_dir)],
        }
        config_path = temp_dir / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config_dict, f)

        config = load_config(config_path)
        assert config.server.auth_token == "secret-token-123"


class TestGetConfig:
    """Tests for getting the current config."""

    def test_get_config_after_load(self, loaded_config):
        """Test that get_config returns the loaded config."""
        config = get_config()
        assert config is loaded_config

    def test_get_config_before_load(self):
        """Test that get_config raises RuntimeError if not loaded."""
        # Reset global state
        config_module._config = None
        with pytest.raises(RuntimeError, match="Configuration not loaded"):
            get_config()


class TestReloadConfig:
    """Tests for config reloading."""

    def test_reload_config(self, sample_config_file: Path, temp_dir: Path):
        """Test reloading configuration."""
        # Initial load
        config = load_config(sample_config_file)
        assert config.server.port == 9999

        # Modify config file
        with open(sample_config_file) as f:
            config_dict = yaml.safe_load(f)
        config_dict["server"]["port"] = 8888
        with open(sample_config_file, "w") as f:
            yaml.dump(config_dict, f)

        # Reload
        new_config, loaded_at = reload_config()
        assert new_config.server.port == 8888
        assert loaded_at is not None

    def test_reload_preserves_old_on_error(self, sample_config_file: Path):
        """Test that invalid new config preserves old config."""
        config = load_config(sample_config_file)
        original_port = config.server.port

        # Write invalid config
        sample_config_file.write_text("this: is: invalid: [")

        # Reload should fail but preserve old config
        with pytest.raises(Exception):
            reload_config()

        # Old config should still be accessible
        current = get_config()
        assert current.server.port == original_port

    def test_reload_without_initial_load(self):
        """Test that reload fails if no config was loaded."""
        config_module._config = None
        config_module._config_path = None

        with pytest.raises(RuntimeError, match="No config previously loaded"):
            reload_config()
