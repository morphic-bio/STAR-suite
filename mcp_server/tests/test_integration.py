"""Integration tests for MCP server end-to-end workflows.

Tests the full discover → preflight → run → collect flow.
"""

import time
from pathlib import Path
from unittest.mock import patch

import pytest

from mcp_server.config import get_config, load_config
from mcp_server.schemas.run_config import RunConfig
from mcp_server.tools.discovery import find_docs, find_tests, list_datasets, list_test_suites
from mcp_server.tools.executor import collect_outputs, get_executor, run_script
from mcp_server.tools.preflight import run_preflight


class TestDiscoverPreflightFlow:
    """Test discovery → preflight flow."""

    @pytest.fixture(autouse=True)
    def setup_config(self, loaded_config, mock_repo_structure):
        """Use mock config for all tests."""
        pass

    def test_list_test_suites_then_preflight(self):
        """Test that discovered test suites can be preflighted."""
        # Discovery
        suites_response = list_test_suites()
        assert suites_response.suites

        # Get first runnable script
        runnable_script = None
        for suite in suites_response.suites:
            for script in suite.scripts:
                if script.runnable:
                    runnable_script = script.name
                    break
            if runnable_script:
                break

        assert runnable_script, "No runnable scripts found in test suites"

        # Preflight
        run_config = RunConfig(script=runnable_script)
        preflight_result = run_preflight(run_config)

        # Should pass (may have warnings about missing binaries in test env)
        assert preflight_result.checks
        script_check = next(
            (c for c in preflight_result.checks if c.name == "script_allowed"), None
        )
        assert script_check is not None
        assert script_check.passed

    def test_find_tests_then_preflight(self):
        """Test that found tests can be preflighted."""
        # Search for a test
        search_response = find_tests(tag="test")
        
        if not search_response.matches:
            pytest.skip("No tests found matching 'test'")

        # Get first runnable match
        runnable_script = None
        config = get_config()
        for match in search_response.matches:
            script = config.get_script(match.name)
            if script and script.runnable:
                runnable_script = match.name
                break

        if not runnable_script:
            pytest.skip("No runnable tests found")

        # Preflight
        run_config = RunConfig(script=runnable_script)
        preflight_result = run_preflight(run_config)

        assert preflight_result.checks
        script_check = next(
            (c for c in preflight_result.checks if c.name == "script_allowed"), None
        )
        assert script_check is not None
        assert script_check.passed


class TestPreflightRunCollectFlow:
    """Test preflight → run → collect flow."""

    @pytest.fixture(autouse=True)
    def setup_config(self, loaded_config, mock_repo_structure, temp_dir):
        """Use mock config and temp directory for all tests."""
        self.temp_dir = temp_dir

    @pytest.fixture
    def echo_script(self, temp_dir: Path) -> str:
        """Create a simple echo script for testing."""
        scripts_dir = temp_dir / "tests"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        
        script_path = scripts_dir / "echo_test.sh"
        script_path.write_text(
            "#!/bin/bash\necho 'Hello from integration test'\nexit 0\n"
        )
        script_path.chmod(0o755)

        # Add script to config
        config = get_config()
        from mcp_server.schemas.config import ScriptConfig
        
        new_script = ScriptConfig(
            name="echo_test",
            path=str(script_path),
            module="test",
            description="Integration test echo script",
            timeout_seconds=10,
        )
        config.scripts.append(new_script)

        return "echo_test"

    def test_preflight_run_collect_success(self, echo_script):
        """Test successful preflight → run → collect flow."""
        run_config = RunConfig(
            script=echo_script,
            out_dir=str(self.temp_dir / "output"),
        )

        # Preflight should pass
        preflight_result = run_preflight(run_config)
        assert preflight_result.valid, f"Preflight failed: {preflight_result.errors}"

        # Run script
        run_response = run_script(run_config)
        
        # Check for error response
        if hasattr(run_response, "code"):
            pytest.fail(f"run_script returned error: {run_response.message}")

        assert run_response.run_id
        assert run_response.status in ("queued", "running")

        # Wait for completion (with timeout)
        executor = get_executor()
        max_wait = 15  # seconds
        waited = 0
        while waited < max_wait:
            state = executor.get_state(run_response.run_id)
            if state and state.status.value in ("completed", "failed", "timeout"):
                break
            time.sleep(0.5)
            waited += 0.5

        # Collect outputs
        collect_response = collect_outputs(run_response.run_id)

        # Check for error response
        if hasattr(collect_response, "code"):
            pytest.fail(f"collect_outputs returned error: {collect_response.message}")

        assert collect_response.run_id == run_response.run_id
        assert collect_response.status == "completed"
        assert collect_response.exit_code == 0
        assert collect_response.log_files
        assert "combined" in collect_response.log_files

    def test_preflight_failure_prevents_run(self):
        """Test that preflight failure prevents script execution."""
        run_config = RunConfig(
            script="nonexistent_script",
        )

        # Preflight should fail
        preflight_result = run_preflight(run_config)
        assert not preflight_result.valid

        # run_script should also fail (it runs preflight internally)
        run_response = run_script(run_config)
        
        assert hasattr(run_response, "code")
        assert run_response.code == "PREFLIGHT_FAILED"


class TestSecurityFlow:
    """Test security validation across the flow."""

    @pytest.fixture(autouse=True)
    def setup_config(self, loaded_config, mock_repo_structure, temp_dir):
        """Use mock config and temp directory for all tests."""
        self.temp_dir = temp_dir

    def test_untrusted_output_dir_rejected(self):
        """Test that output directories outside trusted roots are rejected."""
        config = get_config()
        
        # Get a runnable script
        runnable_script = None
        for script in config.scripts:
            if script.runnable:
                runnable_script = script.name
                break

        if not runnable_script:
            pytest.skip("No runnable scripts in config")

        run_config = RunConfig(
            script=runnable_script,
            out_dir="/etc/passwd",  # Definitely outside trusted roots
        )

        preflight_result = run_preflight(run_config)
        
        # Should fail due to untrusted output directory
        output_check = next(
            (c for c in preflight_result.checks if c.name == "output_directory"), None
        )
        assert output_check is not None
        assert not output_check.passed
        assert "trusted roots" in output_check.message.lower()

    def test_working_dir_validation(self, temp_dir: Path):
        """Test that working_dir is validated against trusted roots."""
        scripts_dir = temp_dir / "tests"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        
        script_path = scripts_dir / "wdir_test.sh"
        script_path.write_text("#!/bin/bash\necho 'test'\n")
        script_path.chmod(0o755)

        # Add script with untrusted working_dir
        config = get_config()
        from mcp_server.schemas.config import ScriptConfig
        
        new_script = ScriptConfig(
            name="wdir_test",
            path=str(script_path),
            module="test",
            description="Working dir test",
            working_dir="/etc",  # Outside trusted roots
        )
        config.scripts.append(new_script)

        run_config = RunConfig(script="wdir_test")
        preflight_result = run_preflight(run_config)

        # Should fail on working_dir check
        wdir_check = next(
            (c for c in preflight_result.checks if c.name == "working_dir"), None
        )
        assert wdir_check is not None
        assert not wdir_check.passed
        assert "trusted roots" in wdir_check.message.lower()


class TestFullEndToEnd:
    """Full end-to-end test simulating agent workflow."""

    @pytest.fixture(autouse=True)
    def setup_config(self, loaded_config, mock_repo_structure, temp_dir):
        """Use mock config and temp directory for all tests."""
        self.temp_dir = temp_dir

    def test_agent_workflow(self, temp_dir: Path):
        """Simulate a typical agent workflow:
        1. List available test suites
        2. Find tests matching a keyword
        3. Preflight the selected test
        4. Run the test
        5. Collect outputs
        """
        # Create a test script
        scripts_dir = temp_dir / "tests"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        
        script_path = scripts_dir / "agent_workflow_test.sh"
        script_path.write_text(
            "#!/bin/bash\n"
            "echo 'Starting agent workflow test'\n"
            "sleep 1\n"
            "echo 'Test completed successfully'\n"
            "exit 0\n"
        )
        script_path.chmod(0o755)

        # Add script to config
        config = get_config()
        from mcp_server.schemas.config import ScriptConfig, TestSuiteConfig
        
        new_script = ScriptConfig(
            name="agent_workflow",
            path=str(script_path),
            module="integration",
            description="Agent workflow integration test",
            timeout_seconds=30,
        )
        config.scripts.append(new_script)

        # Add test suite
        new_suite = TestSuiteConfig(
            module="integration",
            description="Integration tests",
            scripts=["agent_workflow"],
        )
        config.test_suites.append(new_suite)

        # Step 1: List test suites
        suites = list_test_suites(module="integration")
        assert suites.suites
        assert any(s.module == "integration" for s in suites.suites)

        # Step 2: Find tests (optional, but verifies search works)
        tests = find_tests(tag="workflow")
        # May or may not find depending on mock config setup

        # Step 3: Preflight
        run_config = RunConfig(
            script="agent_workflow",
            out_dir=str(temp_dir / "workflow_output"),
        )
        preflight = run_preflight(run_config)
        assert preflight.valid, f"Preflight failed: {preflight.errors}"

        # Step 4: Run
        result = run_script(run_config)
        if hasattr(result, "code"):
            pytest.fail(f"run_script failed: {result.message}")

        assert result.run_id

        # Step 5: Wait and collect
        executor = get_executor()
        max_wait = 35  # seconds (script has 1s sleep + overhead)
        waited = 0
        while waited < max_wait:
            state = executor.get_state(result.run_id)
            if state and state.status.value in ("completed", "failed", "timeout"):
                break
            time.sleep(0.5)
            waited += 0.5

        outputs = collect_outputs(result.run_id)
        if hasattr(outputs, "code"):
            pytest.fail(f"collect_outputs failed: {outputs.message}")

        assert outputs.status == "completed"
        assert outputs.exit_code == 0
        assert outputs.duration_seconds is not None
        assert outputs.duration_seconds >= 1  # Script has 1s sleep
