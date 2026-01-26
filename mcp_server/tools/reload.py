"""Config reload tool for MCP server."""

from datetime import datetime, timezone

from ..config import get_config_path, reload_config as _reload_config
from ..schemas.responses import ReloadConfigResponse
from .executor import _executor


def reload_config() -> ReloadConfigResponse:
    """Reload configuration without restarting the server.

    This also updates the executor's concurrency and queue settings.

    Returns:
        ReloadConfigResponse with reload status.
    """
    config_path = get_config_path()

    try:
        new_config, loaded_at = _reload_config()
        
        # Update executor settings if it exists
        if _executor is not None:
            _executor.update_settings(
                max_concurrent=new_config.execution.max_concurrent_jobs,
                max_queue_size=new_config.execution.queue_max_size,
            )
        
        return ReloadConfigResponse(
            reloaded=True,
            config_path=str(config_path) if config_path else "",
            loaded_at=loaded_at or datetime.now(timezone.utc),
        )
    except Exception as e:
        return ReloadConfigResponse(
            reloaded=False,
            config_path=str(config_path) if config_path else "",
            loaded_at=datetime.now(timezone.utc),
            error=str(e),
        )
