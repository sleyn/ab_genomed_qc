# Python environment for scripts

Environment is build based on official ubuntu:20.04 docker image.
Contains:

- Biopython
- Pandas

## Building image

The image `` was built with the following commands:

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/arm64,linux/amd64 -t semenleyn/ab-gen-qual-py-env:latest --push .
```