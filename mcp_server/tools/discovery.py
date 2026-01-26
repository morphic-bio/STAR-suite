"""Discovery tools for MCP server."""

import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import yaml

from ..config import get_config
from ..schemas.responses import (
    DatasetInfo,
    DocMatch,
    FindDocsResponse,
    FindTestsResponse,
    ListDatasetsResponse,
    ListTestSuitesResponse,
    ScriptInfo,
    TestMatch,
    TestSuiteInfo,
)
from .utils import get_directory_size


def list_datasets() -> ListDatasetsResponse:
    """List all configured datasets with metadata.

    Returns:
        ListDatasetsResponse with dataset information.
    """
    config = get_config()
    datasets: list[DatasetInfo] = []

    for ds_config in config.datasets:
        ds_path = Path(ds_config.path)

        # Get size and mtime if path exists
        size_bytes: Optional[int] = None
        mtime: Optional[datetime] = None
        metadata: Optional[dict[str, Any]] = None

        if ds_path.exists():
            try:
                size_bytes = get_directory_size(ds_path)
                stat = ds_path.stat()
                mtime = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc)
            except (OSError, PermissionError):
                pass

            # Load metadata file if specified
            if ds_config.metadata_file:
                metadata_path = ds_path / ds_config.metadata_file
                if metadata_path.exists():
                    try:
                        with open(metadata_path) as f:
                            metadata = yaml.safe_load(f)
                    except Exception:
                        pass

        datasets.append(
            DatasetInfo(
                id=ds_config.id,
                path=str(ds_config.path),
                description=ds_config.description,
                size_bytes=size_bytes,
                mtime=mtime,
                metadata=metadata,
            )
        )

    return ListDatasetsResponse(datasets=datasets)


def list_test_suites(module: Optional[str] = None) -> ListTestSuitesResponse:
    """List test suites with fixture availability.

    Args:
        module: Optional module filter.

    Returns:
        ListTestSuitesResponse with suite information.
    """
    config = get_config()
    suites: list[TestSuiteInfo] = []

    for suite_config in config.test_suites:
        # Filter by module if specified
        if module and suite_config.module != module:
            continue

        # Build script info for each script in the suite
        scripts: list[ScriptInfo] = []
        suite_runnable = suite_config.runnable

        for script_name in suite_config.scripts:
            script_config = config.get_script(script_name)
            if not script_config:
                continue

            # Check for missing fixtures
            missing_fixtures: list[str] = []
            script_runnable = script_config.runnable

            if script_config.fixtures:
                for fixture_path in script_config.fixtures:
                    if not Path(fixture_path).exists():
                        missing_fixtures.append(fixture_path)
                        script_runnable = False

            scripts.append(
                ScriptInfo(
                    name=script_config.name,
                    path=script_config.path,
                    description=script_config.description,
                    runnable=script_runnable and suite_runnable,
                    missing_fixtures=missing_fixtures,
                )
            )

        # Determine disk space requirement
        disk_space_min_gb = (
            suite_config.disk_space_min_gb
            if suite_config.disk_space_min_gb is not None
            else config.disk_space.default_min_gb
        )

        suites.append(
            TestSuiteInfo(
                module=suite_config.module,
                description=suite_config.description,
                disk_space_min_gb=disk_space_min_gb,
                runnable=suite_runnable and any(s.runnable for s in scripts),
                scripts=scripts,
            )
        )

    return ListTestSuitesResponse(suites=suites)


def find_docs(topic: str) -> FindDocsResponse:
    """Search documentation by topic/keyword.

    Args:
        topic: Search term.

    Returns:
        FindDocsResponse with matching documents.
    """
    config = get_config()
    docs_dir = config.paths.repo_root / "docs"
    matches: list[DocMatch] = []

    if not docs_dir.exists():
        return FindDocsResponse(matches=matches)

    # Case-insensitive search pattern
    pattern = re.compile(re.escape(topic), re.IGNORECASE)

    for doc_path in docs_dir.rglob("*.md"):
        try:
            content = doc_path.read_text()
        except (OSError, PermissionError):
            continue

        # Check if topic matches in filename or content
        filename_match = pattern.search(doc_path.name)
        content_match = pattern.search(content)

        if filename_match or content_match:
            # Extract title from first heading
            title: Optional[str] = None
            for line in content.split("\n"):
                if line.startswith("# "):
                    title = line[2:].strip()
                    break

            # Extract snippet around match
            snippet: Optional[str] = None
            if content_match:
                start = max(0, content_match.start() - 50)
                end = min(len(content), content_match.end() + 100)
                snippet = content[start:end].replace("\n", " ").strip()
                if start > 0:
                    snippet = "..." + snippet
                if end < len(content):
                    snippet = snippet + "..."

            matches.append(
                DocMatch(
                    path=str(doc_path.relative_to(config.paths.repo_root)),
                    title=title,
                    snippet=snippet,
                )
            )

    return FindDocsResponse(matches=matches)


def find_tests(tag: str) -> FindTestsResponse:
    """Search test scripts by tag/keyword.

    Args:
        tag: Search term (e.g., "cbub", "flex", "crispr").

    Returns:
        FindTestsResponse with matching scripts.
    """
    config = get_config()
    matches: list[TestMatch] = []

    # Case-insensitive search pattern
    pattern = re.compile(re.escape(tag), re.IGNORECASE)

    for script in config.scripts:
        # Search in name, path, description, and module
        if (
            pattern.search(script.name)
            or pattern.search(script.path)
            or pattern.search(script.description)
            or pattern.search(script.module)
        ):
            matches.append(
                TestMatch(
                    name=script.name,
                    path=script.path,
                    module=script.module,
                    description=script.description,
                )
            )

    return FindTestsResponse(matches=matches)
