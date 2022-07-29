# MAFFT docker container

Docker container containing MAFFT v.7.505 alignment software.

## Building image

The image `` was built with the following commands:

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/arm64,linux/amd64 -t semenleyn/mafft:latest --push .
```

