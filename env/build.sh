#!/usr/bin/env bash
#
# build.sh — Build the snp_call_nf Docker image and Apptainer SIF
#
# Usage:
#   ./env/build.sh              # build both
#   ./env/build.sh docker       # Docker only
#   ./env/build.sh apptainer    # Apptainer SIF only (Docker image must exist)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

IMAGE_NAME="${IMAGE_NAME:-snp_call_nf}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
SIF_FILE="${PROJECT_DIR}/snp_call_nf.sif"

# ---------------------------------------------------------------------------
# build_docker
# ---------------------------------------------------------------------------
build_docker() {
    echo "=============================================="
    echo " Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"
    echo "=============================================="
    docker build --no-cache \
        -t "${IMAGE_NAME}:${IMAGE_TAG}" \
        -f "${SCRIPT_DIR}/Dockerfile" \
        "${SCRIPT_DIR}"

    echo ""
    echo "Docker image built:"
    docker images "${IMAGE_NAME}:${IMAGE_TAG}"
}

# ---------------------------------------------------------------------------
# build_apptainer
# ---------------------------------------------------------------------------
build_apptainer() {
    # Ensure the Docker image exists
    if ! docker image inspect "${IMAGE_NAME}:${IMAGE_TAG}" &>/dev/null; then
        echo "ERROR: Docker image '${IMAGE_NAME}:${IMAGE_TAG}' not found."
        echo "       Run './env/build.sh docker' first."
        exit 1
    fi

    echo "=============================================="
    echo " Building Apptainer SIF: ${SIF_FILE}"
    echo "=============================================="
    apptainer build --force "${SIF_FILE}" "docker-daemon://${IMAGE_NAME}:${IMAGE_TAG}"

    echo ""
    echo "SIF file built:"
    ls -lh "${SIF_FILE}"
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
case "${1:-all}" in
    docker)
        build_docker
        ;;
    apptainer)
        build_apptainer
        ;;
    all)
        build_docker
        echo ""
        build_apptainer
        echo ""
        echo "===== Done ====="
        ;;
    *)
        echo "Usage: $0 {docker|apptainer|all}"
        exit 1
        ;;
esac
