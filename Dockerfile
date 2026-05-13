FROM python:3.11-slim

ARG DEBIAN_FRONTEND=noninteractive
ARG BINDASH_TAG=v2.6

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang \
    cmake \
    curl \
    g++ \
    gcc \
    git \
    make \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir uv

WORKDIR /workspace
COPY . /workspace

RUN uv sync
RUN make -C src
RUN BINDASH_TAG=${BINDASH_TAG} bash scripts/bootstrap.sh
RUN ln -sf /workspace/src/oddsketch /usr/local/bin/oddsketch

ENV PATH="/workspace/.venv/bin:/workspace/experiments/tools/bin:${PATH}"
ENV ODDSKETCH_BIN="/workspace/src/oddsketch"

CMD ["oddsketch", "--help"]
