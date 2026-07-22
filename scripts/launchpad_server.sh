#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/launchpad_server.sh up [--host HOST] [--port PORT] [--python PY] [--config PATH] [--foreground]
  bash scripts/launchpad_server.sh down [--force]
  bash scripts/launchpad_server.sh status
  bash scripts/launchpad_server.sh logs

Notes:
  - "up" starts the STAR Server (MCP + Launchpad) on one port.
  - "down" stops only the PID tracked by the pidfile (safe by default).
  - The "Run in shell" button is enabled only when you open Launchpad from
    the server host via localhost/127.0.0.1.

Examples:
  bash scripts/launchpad_server.sh up
  bash scripts/launchpad_server.sh up --host 0.0.0.0 --port 8765
  bash scripts/launchpad_server.sh status
  bash scripts/launchpad_server.sh logs
  bash scripts/launchpad_server.sh down
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pidfile="${repo_root}/plans/artifacts/launchpad_server.pid"
logfile="${repo_root}/plans/artifacts/launchpad_server.log"

cmd="${1:-}"
shift || true

python_bin="${PYTHON:-python3}"
host="${HOST:-127.0.0.1}"
port="${PORT:-8765}"
config_path="${CONFIG:-${repo_root}/mcp_server/config.yaml}"
foreground="0"
force="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      python_bin="${2:-}"; shift 2 ;;
    --host)
      host="${2:-}"; shift 2 ;;
    --port)
      port="${2:-}"; shift 2 ;;
    --config)
      config_path="${2:-}"; shift 2 ;;
    --foreground)
      foreground="1"; shift ;;
    --force)
      force="1"; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

mkdir -p "$(dirname "${pidfile}")"

pidfile_read() {
  local f="$1"
  PID=""
  PF_HOST=""
  PF_PORT=""
  PF_LOG=""
  [[ -f "$f" ]] || return 1
  while IFS='=' read -r k v; do
    case "$k" in
      pid) PID="$v" ;;
      host) PF_HOST="$v" ;;
      port) PF_PORT="$v" ;;
      logfile) PF_LOG="$v" ;;
    esac
  done < "$f"
  [[ "${PID}" =~ ^[0-9]+$ ]] || return 1
  return 0
}

pid_alive() {
  local pid="$1"
  kill -0 "$pid" 2>/dev/null
}

pid_cmdline() {
  local pid="$1"
  ps -p "$pid" -o args= 2>/dev/null || true
}

print_urls() {
  local h="$1"
  local p="$2"
  local h_show="$h"
  if [[ "$h" == "0.0.0.0" ]]; then
    h_show="127.0.0.1"
  fi
  echo "Launchpad: http://${h_show}:${p}/launchpad/"
  echo "Launchpad API: http://${h_show}:${p}/launchpad/api/workflows"
}

case "$cmd" in
  up)
    if pidfile_read "$pidfile" && pid_alive "$PID"; then
      echo "Already running (pid ${PID})."
      print_urls "${PF_HOST:-$host}" "${PF_PORT:-$port}"
      exit 0
    fi

    rm -f "$pidfile" 2>/dev/null || true

    if [[ "$foreground" == "1" ]]; then
      echo "Starting STAR Server in foreground (Ctrl+C to stop)."
      print_urls "$host" "$port"
      exec "$python_bin" -m mcp_server.app \
        --transport http \
        --host "$host" \
        --port "$port" \
        --config "$config_path"
    fi

    echo "Starting STAR Server in background..."
    echo "Log: ${logfile}"
    # Use --config explicitly so the pidfile uniquely refers to this repo.
    nohup "$python_bin" -m mcp_server.app \
      --transport http \
      --host "$host" \
      --port "$port" \
      --config "$config_path" \
      >"$logfile" 2>&1 &
    pid="$!"

    {
      echo "pid=${pid}"
      echo "host=${host}"
      echo "port=${port}"
      echo "logfile=${logfile}"
      echo "started_at=$(date -Is)"
    } > "$pidfile"

    # Basic health check: process still alive shortly after start.
    sleep 0.8
    if ! pid_alive "$pid"; then
      echo "Server failed to start (pid ${pid} exited). See log: ${logfile}" >&2
      rm -f "$pidfile" 2>/dev/null || true
      exit 1
    fi

    echo "Started (pid ${pid})."
    print_urls "$host" "$port"
    ;;

  down)
    if ! pidfile_read "$pidfile"; then
      echo "No pidfile found: ${pidfile}" >&2
      echo "If you started in foreground, stop it with Ctrl+C." >&2
      exit 1
    fi

    if ! pid_alive "$PID"; then
      echo "Not running (stale pidfile pid ${PID}). Removing pidfile."
      rm -f "$pidfile" 2>/dev/null || true
      exit 0
    fi

    args="$(pid_cmdline "$PID")"
    if [[ "$force" != "1" ]] && [[ "$args" != *"mcp_server.app"* ]]; then
      echo "Refusing to stop pid ${PID} (does not look like mcp_server.app)." >&2
      echo "Command line: ${args}" >&2
      echo "Re-run with --force if you're sure." >&2
      exit 2
    fi

    echo "Stopping STAR Server (pid ${PID})..."
    kill "$PID" 2>/dev/null || true
    for _ in $(seq 1 40); do
      if ! pid_alive "$PID"; then
        break
      fi
      sleep 0.25
    done
    if pid_alive "$PID"; then
      echo "Still running; sending SIGKILL." >&2
      kill -9 "$PID" 2>/dev/null || true
    fi
    rm -f "$pidfile" 2>/dev/null || true
    echo "Stopped."
    ;;

  status)
    if pidfile_read "$pidfile" && pid_alive "$PID"; then
      echo "Running (pid ${PID})."
      print_urls "${PF_HOST:-$host}" "${PF_PORT:-$port}"
      echo "Log: ${PF_LOG:-$logfile}"
      exit 0
    fi
    echo "Not running."
    exit 1
    ;;

  logs)
    if pidfile_read "$pidfile" && [[ -n "${PF_LOG}" ]] && [[ -f "${PF_LOG}" ]]; then
      tail -n 80 "${PF_LOG}"
    elif [[ -f "${logfile}" ]]; then
      tail -n 80 "${logfile}"
    else
      echo "No log found." >&2
      exit 1
    fi
    ;;

  ""|-h|--help)
    usage
    exit 0
    ;;

  *)
    echo "Unknown command: ${cmd}" >&2
    usage
    exit 2
    ;;
esac

