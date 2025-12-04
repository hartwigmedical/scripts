#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: ./build_and_push.sh <path/to/Dockerfile> <image_name> <version>"
  exit 1
fi

DOCKERFILE_PATH="$1"
IMAGE_NAME="$2"
VERSION="$3"

REGISTRY="europe-west4-docker.pkg.dev"
PROJECT="hmf-build"
REPO="hmf-docker-crunch"

FULL_IMAGE="$REGISTRY/$PROJECT/$REPO/$IMAGE_NAME:$VERSION"

BUILD_CONTEXT="$(dirname "$DOCKERFILE_PATH")"
DOCKERFILE_NAME="$(basename "$DOCKERFILE_PATH")"

echo "Building image..."
docker buildx build \
  --platform linux/amd64 \
  -f "$DOCKERFILE_PATH" \
  -t "$IMAGE_NAME:$VERSION" \
  --load \
  "$BUILD_CONTEXT"

echo "Tagging image..."
docker tag "$IMAGE_NAME:$VERSION" "$FULL_IMAGE"

echo "Pushing image..."
docker push "$FULL_IMAGE"

echo "Done."
echo "Pushed: $FULL_IMAGE"