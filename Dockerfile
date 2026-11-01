# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
#
# Runtime image: the kronspice and evolve commands only, as static binaries
# (pure Go, no cgo) on a minimal base, about 10 MB. For running GA
# experiments on servers and clusters, one container per seed or profile.
# The tools image for development (race detector, ngspice) is separate:
# internal/golden/docker/Dockerfile.
#
#   docker build -t kronspice .
#   docker run --rm -v "$PWD":/work -w /work kronspice testdata/netlists/rc_step.cir -o out.csv
#   docker run --rm -v "$PWD":/work -w /work --entrypoint /usr/local/bin/evolve kronspice -config examples/rc_lowpass_values/config.json
#
# Results are bit-identical for the same input on every machine only with
# the same Go toolchain, so the builder version is pinned.
FROM golang:1.24-bookworm AS build
WORKDIR /src
COPY go.mod go.sum* ./
RUN go mod download
COPY . .
ARG VERSION=dev
ENV CGO_ENABLED=0
RUN go build -trimpath -ldflags="-s -w -X main.version=${VERSION}" -o /out/kronspice ./cmd/kronspice \
 && go build -trimpath -ldflags="-s -w -X main.version=${VERSION}" -o /out/evolve ./cmd/evolve

FROM gcr.io/distroless/static-debian12:nonroot
ARG VERSION=dev
LABEL org.opencontainers.image.title="kronspice" \
      org.opencontainers.image.description="SPICE-compatible analog circuit simulator in pure Go, with a genetic-algorithm engine" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.licenses="GPL-3.0-only OR LicenseRef-kronspice-Commercial" \
      org.opencontainers.image.source="https://github.com/fdamador/kronspice"
# GPL-3.0 distribution: the license texts travel with the binaries, and the
# source is the repository named above at the image's version tag.
COPY LICENSE LICENSING.md COMMERCIAL-LICENSE.md /usr/share/doc/kronspice/
COPY --from=build /out/kronspice /out/evolve /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/kronspice"]
