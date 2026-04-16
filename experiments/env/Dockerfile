FROM python:3.11-slim

ARG DEBIAN_FRONTEND=noninteractive
ARG BINDASH_REF=v2.2.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    zlib1g-dev \
    ca-certificates \
    curl \
    clang \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir uv

WORKDIR /workspace
COPY . /workspace

RUN uv sync
RUN make -C src
RUN bash experiments/tools/scripts/build_tools.sh --method source --ref ${BINDASH_REF}

ENV PATH="/workspace/experiments/tools/bin:${PATH}"
ENV ODDSKETCH_BIN="/workspace/src/oddsketch"

CMD ["uv", "run", "python", "experiments/pair_task/scripts/project_runner.py", "--config", "experiments/pair_task/config.json"]
