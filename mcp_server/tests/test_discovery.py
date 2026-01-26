"""Tests for discovery tools."""

from pathlib import Path

import pytest

from mcp_server.config import load_config
from mcp_server.tools.discovery import (
    find_docs,
    find_tests,
    list_datasets,
    list_test_suites,
)
import mcp_server.config as config_module


@pytest.fixture
def setup_config(sample_config_file: Path, mock_repo_structure: Path):
    """Load config with mock repo structure."""
    # Update config to use mock repo
    import yaml

    with open(sample_config_file) as f:
        config_dict = yaml.safe_load(f)

    config_dict["paths"]["repo_root"] = str(mock_repo_structure)
    config_dict["datasets"][0]["path"] = str(mock_repo_structure / "datasets" / "test")

    with open(sample_config_file, "w") as f:
        yaml.dump(config_dict, f)

    config = load_config(sample_config_file)
    yield config

    # Cleanup
    config_module._config = None
    config_module._config_path = None
    config_module._config_loaded_at = None


class TestListDatasets:
    """Tests for list_datasets."""

    def test_list_datasets_returns_configured(self, setup_config, mock_repo_structure):
        """Test that list_datasets returns configured datasets."""
        result = list_datasets()
        assert len(result.datasets) == 1
        assert result.datasets[0].id == "test-dataset"

    def test_list_datasets_includes_path(self, setup_config, mock_repo_structure):
        """Test that dataset info includes path."""
        result = list_datasets()
        assert "datasets" in result.datasets[0].path

    def test_list_datasets_existing_path_has_mtime(self, setup_config, mock_repo_structure):
        """Test that existing datasets have mtime."""
        # Create the dataset directory
        ds_path = mock_repo_structure / "datasets" / "test"
        ds_path.mkdir(parents=True, exist_ok=True)

        result = list_datasets()
        # mtime should be set for existing paths
        assert result.datasets[0].mtime is not None


class TestListTestSuites:
    """Tests for list_test_suites."""

    def test_list_test_suites_returns_all(self, setup_config):
        """Test that list_test_suites returns all suites."""
        result = list_test_suites()
        assert len(result.suites) == 1
        assert result.suites[0].module == "test"

    def test_list_test_suites_filter_by_module(self, setup_config):
        """Test filtering by module."""
        result = list_test_suites(module="test")
        assert len(result.suites) == 1

        result = list_test_suites(module="nonexistent")
        assert len(result.suites) == 0

    def test_list_test_suites_shows_scripts(self, setup_config):
        """Test that suites include script info."""
        result = list_test_suites()
        assert len(result.suites[0].scripts) == 2
        script_names = [s.name for s in result.suites[0].scripts]
        assert "test_script" in script_names

    def test_list_test_suites_missing_fixtures_marked(self, setup_config):
        """Test that scripts with missing fixtures are marked."""
        result = list_test_suites()
        fixture_script = next(
            s for s in result.suites[0].scripts if s.name == "script_with_fixtures"
        )
        assert not fixture_script.runnable
        assert len(fixture_script.missing_fixtures) > 0


class TestFindDocs:
    """Tests for find_docs."""

    def test_find_docs_by_topic(self, setup_config, mock_repo_structure):
        """Test finding docs by topic."""
        result = find_docs("test")
        assert len(result.matches) >= 1

    def test_find_docs_case_insensitive(self, setup_config, mock_repo_structure):
        """Test that search is case insensitive."""
        result = find_docs("CBUB")
        cbub_matches = [m for m in result.matches if "CBUB" in m.path.upper()]
        assert len(cbub_matches) >= 1

    def test_find_docs_includes_snippet(self, setup_config, mock_repo_structure):
        """Test that matches include snippets."""
        result = find_docs("documentation")
        if result.matches:
            assert result.matches[0].snippet is not None

    def test_find_docs_no_matches(self, setup_config, mock_repo_structure):
        """Test that no matches returns empty list."""
        result = find_docs("xyznonexistent123")
        assert len(result.matches) == 0


class TestFindTests:
    """Tests for find_tests."""

    def test_find_tests_by_name(self, setup_config):
        """Test finding tests by name."""
        result = find_tests("test_script")
        assert len(result.matches) >= 1
        assert result.matches[0].name == "test_script"

    def test_find_tests_by_module(self, setup_config):
        """Test finding tests by module."""
        result = find_tests("test")
        assert len(result.matches) >= 1

    def test_find_tests_case_insensitive(self, setup_config):
        """Test that search is case insensitive."""
        result = find_tests("TEST")
        assert len(result.matches) >= 1

    def test_find_tests_no_matches(self, setup_config):
        """Test that no matches returns empty list."""
        result = find_tests("xyznonexistent123")
        assert len(result.matches) == 0

    def test_find_tests_includes_description(self, setup_config):
        """Test that matches include description."""
        result = find_tests("test_script")
        assert result.matches[0].description == "A test script"
