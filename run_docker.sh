#!/bin/bash

# Set image name
IMAGE_NAME="bemcpp"

# Rebuild image (optional: add --no-cache if needed)
echo "Building Docker image: $IMAGE_NAME"
docker build -t $IMAGE_NAME .

# Run container interactively with current project directory mounted
echo "Starting container with mounted project..."
docker run -it --rm \
    -v "$(pwd)":/app \
    -w /app \
    $IMAGE_NAME \
    /bin/bash
