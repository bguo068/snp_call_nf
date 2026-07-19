# syntax=docker/dockerfile:1

# =============================================================================
# Stage 1 — Build: install all conda packages via pixi, then strip cruft
# =============================================================================
FROM ghcr.io/prefix-dev/pixi:latest AS builder

WORKDIR /build
COPY pixi.toml ./

# Install packages (generates pixi.lock on first run if missing)
RUN pixi install

# Strip everything not needed at runtime **in the build stage** so the
# COPY in the runtime stage only pulls across what we actually need.
RUN rm -rf \
    /build/.pixi/envs/default/conda-meta \
    /build/.pixi/envs/default/include \
    /build/.pixi/envs/default/share/man \
    /build/.pixi/envs/default/share/doc \
    /build/.pixi/envs/default/share/gtk-doc \
    /build/.pixi/envs/default/share/info \
    /build/.pixi/envs/default/share/icons \
    /build/.pixi/envs/default/share/locale \
    /build/.pixi/envs/default/lib/pkgconfig \
    /build/.pixi/envs/default/lib/cmake \
    /build/.pixi/envs/default/lib/*.a \
    /build/.pixi/envs/default/lib/*.la \
    /build/.pixi/envs/default/man \
    /build/.pixi/envs/default/src.zip \
    && find /build/.pixi/envs/default -name "*.pyc" -delete \
    && find /build/.pixi/envs/default -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true \
    && find /build/.pixi/envs/default -name "*.h" -delete 2>/dev/null || true \
    && find /build/.pixi/envs/default -name "*-sources.jar" -delete 2>/dev/null || true \
    && find /build/.pixi/envs/default -name "*-javadoc.jar" -delete 2>/dev/null || true

# =============================================================================
# Stage 2 — Runtime: minimal Debian base with only the stripped environment
# =============================================================================
FROM debian:bookworm-slim

# Install only what's strictly needed at runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libgomp1 \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Copy the already-stripped conda environment
COPY --from=builder /build/.pixi/envs/default /opt/env

ENV PATH="/opt/env/bin:$PATH"
ENV CONDA_PREFIX="/opt/env"

WORKDIR /data
CMD ["/bin/bash"]
