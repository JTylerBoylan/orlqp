#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

PROJECT_NAME="10-lib-cpp"

docker build -t ${PROJECT_NAME} -f "${SCRIPT_DIR}/docker/Dockerfile" "${SCRIPT_DIR}"

# Start CUDA docker
docker run -it --rm \
    -v $SCRIPT_DIR:/app \
    -w /app \
    ${PROJECT_NAME} \
    bash
