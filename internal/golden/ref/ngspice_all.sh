#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
#
# Second, ngspice-computed reference for every test netlist (CLAUDE.md §1
# goal 4). testdata/golden/NAME.csv holds an independent reference (exact
# solutions, Python re-implementations); this script runs the same
# testdata/netlists/NAME.cir through ngspice itself and writes
# testdata/golden/ngspice/NAME.csv with the same columns, for
# analysis.TestNgspiceTwin.
#
# The netlist is used unchanged: a .control block is appended that runs
# its analysis and writes the golden file's columns. Column names map to
# ngspice vectors as
#   v(node)                     -> v(node)
#   i(Vx), i(Lx), i(Ex), i(Hx)  -> i(name)        (branch current)
#   i(Rx), i(Cx)                -> @name[i]       (.options savecurrents)
#   i(Dx)                       -> @name[id]
# Currents of transistors (Q, M) have no single ngspice equivalent of the
# kronspice "+ terminal" current and are left out (logged).
#
# Netlists that already have an ngspice-syntax twin in
# internal/golden/ref/ngspice/ (ngspice_ref.sh) are skipped, and a netlist
# ngspice rejects is logged and skipped, so one failure never blocks the
# rest. Run from the repository root in the tools container:
#   docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_all.sh
# (VS Code task "kronspice: ngspice twin references").
set -u
ver=$(ngspice -v 2>/dev/null | grep -o 'ngspice-[0-9.]*' | head -n 1)
root=$(pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p testdata/golden/ngspice
ok=0
bad=0
for net in testdata/netlists/*.cir; do
	name=$(basename "$net" .cir)
	gold=testdata/golden/$name.csv
	[ -f "$gold" ] || { echo "$name: no golden file, skipped"; continue; }
	[ -f "internal/golden/ref/ngspice/$name.cir" ] && { echo "$name: has an ngspice twin (ngspice_ref.sh), skipped"; continue; }
	header=$(grep -v '^#' "$gold" | head -n 1 | tr -d '\r')
	kind=table
	if [ "$header" = "name,value" ]; then
		kind=op
		cols=$(grep -v '^#' "$gold" | tail -n +2 | cut -d, -f1 | tr -d '\r' | tr '\n' ' ')
		scale=""
	else
		scale=${header%%,*}
		cols=$(echo "$header" | cut -d, -f2- | tr ',' ' ')
	fi
	# Build the .control block: one "let gK = expr" per mappable column.
	lets=""
	names=""
	k=0
	for c in $cols; do
		lc=$(echo "$c" | tr 'A-Z' 'a-z')
		expr=""
		case "$lc" in
		v\(*\)) expr=$lc ;;
		vm\(*\)) expr="vm(${lc#vm(}" ;;
		vp\(*\)) expr="(180/pi)*vp(${lc#vp(}" ;;
		i\(*\))
			dev=${lc#i(}
			dev=${dev%)}
			case "$dev" in
			v* | l* | e* | h*) expr="i($dev)" ;;
			r* | c*) expr="@$dev[i]" ;;
			d*) expr="@$dev[id]" ;;
			*) echo "$name: column $c not mappable to ngspice, left out" ;;
			esac
			;;
		*) echo "$name: column $c not mappable to ngspice, left out" ;;
		esac
		[ -n "$expr" ] || continue
		k=$((k + 1))
		if [ "$kind" = table ] && [ "$scale" = freq ]; then
			lets="${lets}let g$k = real($expr)
"
		else
			lets="${lets}let g$k = $expr
"
		fi
		names="$names $lc"
	done
	[ "$k" -gt 0 ] || { echo "$name: no mappable columns, skipped"; continue; }
	gl=""
	i=1
	while [ $i -le $k ]; do
		gl="$gl g$i"
		i=$((i + 1))
	done
	deck=$tmp/$name.cir
	# Netlist without its .end, then savecurrents, tight tolerances (the
	# reference should be the accurate side), the control
	# block and .end.
	grep -iv '^[[:space:]]*\.end[[:space:]]*$' "$net" | tr -d '\r' >"$deck"
	{
		echo ".options savecurrents"
		# ngspice stops Newton at reltol=1e-3 by default, which leaves its
		# own solution off KCL by ~0.1 % on exponential devices (diode
		# clipper, zener). References converge tightly; transients keep a
		# looser reltol, since it also sets the time-step error control.
		if grep -qi '^[[:space:]]*\.tran' "$net"; then
			echo ".options reltol=1e-6 abstol=1e-15 vntol=1e-9 chgtol=1e-18"
		else
			echo ".options reltol=1e-9 abstol=1e-15 vntol=1e-12"
		fi
		echo ".control"
		echo "set numdgt=15"
		echo "run"
		# No "linearize": ngspice-39 zeroed the last 31 points of
		# mos1_inv_tran.cir after it (exact 0 V, then 5 V at the end). The
		# comparison interpolates the simulation onto the reference's own
		# time points, so the raw points serve.
		printf '%s' "$lets"
		if [ "$kind" = op ]; then
			echo "print$gl"
		else
			echo "wrdata $tmp/out.txt$gl"
		fi
		echo "quit"
		echo ".endc"
		echo ".end"
	} >>"$deck"
	rm -f "$tmp/out.txt"
	(cd "$tmp" && ngspice -b "$deck" >"$tmp/log.txt" 2>&1) || true
	out=testdata/golden/ngspice/$name.csv
	if [ "$kind" = op ]; then
		# "gK = value" lines from print.
		vals=$(grep -E '^g[0-9]+ = ' "$tmp/log.txt")
		if [ -z "$vals" ]; then
			echo "$name: ngspice gave no values:"
			tail -n 15 "$tmp/log.txt"
			bad=$((bad + 1))
			continue
		fi
		{
			echo "# reference: $ver operating point of testdata/netlists/$name.cir (internal/golden/ref/ngspice_all.sh)"
			echo "name,value"
			echo "$vals" | awk -v names="$names" 'BEGIN { n = split(names, nm, " ") }
				{ sub(/^g/, "", $1); printf "%s,%s\n", nm[$1 + 0], $3 }'
		} >"$out"
	else
		if [ ! -s "$tmp/out.txt" ]; then
			echo "$name: ngspice wrote no table:"
			tail -n 15 "$tmp/log.txt"
			bad=$((bad + 1))
			continue
		fi
		# wrdata writes, per vector, "scale value" (real) or a wider block
		# (complex); take the scale from the first field and each vector's
		# value from the first value field of its block.
		echo "$name: first wrdata row: $(head -n 1 "$tmp/out.txt")"
		{
			echo "# reference: $ver of testdata/netlists/$name.cir (internal/golden/ref/ngspice_all.sh)"
			echo "$scale,$(echo $names | tr ' ' ',')"
			awk -v k="$k" 'NF {
				b = NF / k
				printf "%s", $1
				for (j = 0; j < k; j++) {
					if (b == 2) f = 2 + 2 * j
					else if (b == 3) f = 2 + 3 * j
					else f = 3 + 4 * j
					printf ",%s", $f
				}
				print ""
			}' "$tmp/out.txt"
		} >"$out"
	fi
	echo "$out: $(($(wc -l <"$out") - 2)) rows"
	ok=$((ok + 1))
done
echo "ngspice twins: $ok written, $bad failed ($ver)"
