FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    bc \
    build-essential \
    ca-certificates \
    curl \
    g++ \
    git-lfs \
    jq \
    make \
    nodejs \
    npm \
    pigz \
    python3 \
    python3-pandas \
    python3-pip \
    python3-venv \
    samtools \
    sra-toolkit \
    time \
    unzip \
    zlib1g-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
 && rm -rf /var/lib/apt/lists/*
