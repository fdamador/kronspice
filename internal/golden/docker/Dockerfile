# Toolchain for the checks that need more than the Go distribution:
# the race detector (cgo, gcc) and ngspice for golden references.
#
#   docker build -t kronspice-tools internal/golden/docker
#   docker run --rm -v "$PWD":/src kronspice-tools go test -race ./...
#   docker run --rm -v "$PWD":/src kronspice-tools ngspice -b file.cir
#
# The source is mounted, not copied, so the image holds only tools.
FROM golang:1.24-bookworm
RUN apt-get update \
 && apt-get install -y --no-install-recommends ngspice \
 && rm -rf /var/lib/apt/lists/*
ENV CGO_ENABLED=1
WORKDIR /src
