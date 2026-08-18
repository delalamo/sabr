# syntax=docker/dockerfile:1

FROM python:3.11-slim AS build

ARG SOURCE_DIR=.
ARG SABR_VERSION=0.0.0

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /build
COPY ${SOURCE_DIR}/pyproject.toml ${SOURCE_DIR}/README.md \
    ${SOURCE_DIR}/LICENSE ${SOURCE_DIR}/constraints.txt ./
COPY ${SOURCE_DIR}/src ./src
RUN python -m pip install --upgrade pip && \
    SETUPTOOLS_SCM_PRETEND_VERSION="${SABR_VERSION}" \
    HATCH_VCS_PRETEND_VERSION="${SABR_VERSION}" \
    python -m pip install --constraint constraints.txt .

# -----------------------
# Stage 2: Runtime image
# -----------------------
FROM python:3.11-slim
COPY --from=build /usr/local /usr/local
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
ENTRYPOINT ["sabr"]
CMD ["--help"]
