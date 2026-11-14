#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
#
# Runs every ngspice twin in examples/bearden_meg/ngspice/ (the same
# circuits as netlists/, in ngspice syntax: each gyrator as two VCCS, each
# core as a B source carrying d/dt of its flux curve) and leaves NAME.dat
# beside it for TestNgspiceTwin. Run from the repository root inside the
# tools container:
#   docker build -t kronspice-tools internal/golden/docker
#   docker run --rm -v "$PWD":/src kronspice-tools sh examples/bearden_meg/ngspice_run.sh
#   go test -v -run NgspiceTwin ./examples/bearden_meg
set -eu
ver=$(ngspice -v 2>/dev/null | grep -o 'ngspice-[0-9.]*' | head -n 1)
root=$(pwd)
dir=examples/bearden_meg/ngspice
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
echo "$ver" >"$dir/VERSION"
for ref in "$dir"/*.cir; do
	name=$(basename "$ref" .cir)
	rm -f "$tmp/$name.dat"
	# ngspice -b can exit nonzero after wrdata/quit; judge by the .dat file.
	(cd "$tmp" && ngspice -b "$root/$ref" >log.txt 2>&1) || true
	if [ ! -s "$tmp/$name.dat" ]; then
		echo "$name: ngspice wrote no data" >&2
		tail -n 30 "$tmp/log.txt" >&2
		exit 1
	fi
	mv "$tmp/$name.dat" "$dir/$name.dat"
	echo "$dir/$name.dat: $(wc -l <"$dir/$name.dat") points ($ver)"
done
