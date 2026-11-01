#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
#
# Regenerates the transient golden tables whose reference is ngspice: each
# internal/golden/ref/ngspice/NAME.cir is an ngspice twin of
# testdata/netlists/NAME.cir (same circuit, ngspice syntax). Its .control
# block writes out.txt with wrdata after linearize; the "* columns:" line
# names the columns as kronspice reports them.
#
# Run from the repository root inside the tools container:
#   docker build -t kronspice-tools internal/golden/docker
#   docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_ref.sh
set -eu
ver=$(ngspice -v 2>/dev/null | grep -o 'ngspice-[0-9.]*' | head -n 1)
root=$(pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for ref in internal/golden/ref/ngspice/*.cir; do
	name=$(basename "$ref" .cir)
	cols=$(sed -n 's/^\* columns: *//p' "$ref" | tr ' ' ',')
	# ngspice -b exits nonzero for a deck without .print lines; judge by out.txt.
	rm -f "$tmp/out.txt"
	(cd "$tmp" && ngspice -b "$root/$ref" >log.txt 2>&1) || true
	[ -s "$tmp/out.txt" ] || { cat "$tmp/log.txt"; exit 1; }
	out=testdata/golden/$name.csv
	{
		echo "# reference: $ver transient of internal/golden/ref/ngspice/$name.cir (tight tolerances, small max step)"
		echo "time,$cols"
		# wrdata writes "time value time value ..."; keep the first time column.
		awk 'NF { printf "%s", $1; for (i = 2; i <= NF; i += 2) printf ",%s", $i; print "" }' "$tmp/out.txt"
	} >"$out"
	echo "$out: $(($(wc -l <"$out") - 2)) points"
done
