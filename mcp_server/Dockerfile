# MCP Server Dockerfile
# Build: docker build -t star-mcp-server .
# Run:   docker run -p 8765:8765 -v /storage:/storage:ro -v /mnt/pikachu/STAR-suite:/repo:ro star-mcp-server

FROM python:3.10-slim

LABEL maintainer="STAR-suite team"
LABEL description="MCP Server for STAR-suite agent workflows"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Use Docker-specific config
RUN if [ -f config.docker.yaml ]; then cp config.docker.yaml config.yaml; fi

# Create directories for logs and temp files
RUN mkdir -p /app/artifacts /tmp/mcp_runs

# Environment variables (can be overridden)
ENV MCP_AUTH_TOKEN=""
ENV MCP_HOST="0.0.0.0"
ENV MCP_PORT="8765"
ENV MCP_TRANSPORT="http"
ENV MCP_CONFIG_PATH="/app/config.yaml"

# Expose port
EXPOSE 8765

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8765/health')" || exit 1

# Run the server
CMD ["python", "-m", "mcp_server.app"]
