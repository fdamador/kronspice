// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strings"
)

// ErrBadModel is returned for missing, mistyped or out-of-range models.
var ErrBadModel = errors.New("circuit: invalid model")

// Model is a .model card: a named parameter set for a device type.
type Model struct {
	Name   string             // lower-cased
	Type   string             // lower-cased device type: "d"
	Params map[string]float64 // lower-cased parameter names, SI values
}

// paramRule describes one model parameter: whether it may be given and its
// valid range (inclusive unless open is set).
type paramRule struct {
	min, max         float64
	openMin, openMax bool
}

// Known diode parameters (SPICE3 level-1 diode, Quarles et al. SPICE3f3
// User's Manual §4.2 / ngspice manual §7.2). Ranges reject values the model
// equations cannot use.
var diodeParams = map[string]paramRule{
	"is":   {min: 0, max: math.Inf(1), openMin: true},
	"n":    {min: 0, max: math.Inf(1), openMin: true},
	"rs":   {min: 0, max: math.Inf(1)},
	"cjo":  {min: 0, max: math.Inf(1)},
	"cj0":  {min: 0, max: math.Inf(1)},
	"vj":   {min: 0, max: math.Inf(1), openMin: true},
	"m":    {min: 0, max: 1, openMax: true},
	"tt":   {min: 0, max: math.Inf(1)},
	"bv":   {min: 0, max: math.Inf(1), openMin: true},
	"ibv":  {min: 0, max: math.Inf(1), openMin: true},
	"fc":   {min: 0, max: 1, openMax: true},
	"eg":   {min: 0, max: math.Inf(1), openMin: true},
	"xti":  {min: math.Inf(-1), max: math.Inf(1)},
	"kf":   {min: 0, max: math.Inf(1)},
	"af":   {min: 0, max: math.Inf(1), openMin: true},
	"tnom": {min: -273.15, max: math.Inf(1), openMin: true},
}

// Known BJT parameters (SPICE3 Gummel–Poon; SPICE3f3 User's Manual §4.3,
// ngspice manual "BJT models"). All are accepted at the nominal temperature
// (27 °C); temperature-coefficient parameters are accepted and unused until
// circuit temperature is modeled (§3d).
var bjtParams = map[string]paramRule{
	"is": pos, "bf": pos, "nf": pos, "vaf": nonneg, "va": nonneg, "ikf": nonneg, "ik": nonneg,
	"ise": nonneg, "ne": pos, "br": pos, "nr": pos, "var": nonneg, "vb": nonneg,
	"ikr": nonneg, "isc": nonneg, "nc": pos, "rb": nonneg, "irb": nonneg, "rbm": nonneg,
	"re": nonneg, "rc": nonneg, "cje": nonneg, "vje": pos, "pe": pos, "mje": unitOpen, "me": unitOpen,
	"tf": nonneg, "xtf": nonneg, "vtf": nonneg, "itf": nonneg, "ptf": {min: 0, max: 0},
	"cjc": nonneg, "vjc": pos, "pc": pos, "mjc": unitOpen, "mc": unitOpen, "xcjc": {min: 0, max: 1},
	"tr": nonneg, "cjs": nonneg, "ccs": nonneg, "vjs": pos, "ps": pos, "mjs": unitOpen, "ms": unitOpen,
	"xtb": any, "eg": pos, "xti": any, "kf": nonneg, "af": pos, "fc": {min: 0, max: 1, openMax: true},
	"tnom": {min: -273.15, max: math.Inf(1), openMin: true},
}

var (
	pos      = paramRule{min: 0, max: math.Inf(1), openMin: true}
	nonneg   = paramRule{min: 0, max: math.Inf(1)}
	unitOpen = paramRule{min: 0, max: 1, openMax: true}
	any      = paramRule{min: math.Inf(-1), max: math.Inf(1)}
)

// Known MOSFET parameters, levels 1–3 (SPICE3f3 User's Manual §4.4,
// ngspice manual "MOSFET models (levels 1, 2, 3 and 6)"). LEVEL selects the
// drain-current equations; parameters a level does not use are accepted and
// ignored, as in SPICE.
var mosParams = map[string]paramRule{
	"level": {min: 1, max: 3}, "vto": any, "vt0": any, "kp": nonneg, "gamma": nonneg, "phi": pos,
	"lambda": nonneg, "rd": nonneg, "rs": nonneg, "cbd": nonneg, "cbs": nonneg, "is": nonneg, "pb": pos,
	"cgso": nonneg, "cgdo": nonneg, "cgbo": nonneg, "rsh": nonneg, "cj": nonneg, "mj": unitOpen,
	"cjsw": nonneg, "mjsw": unitOpen, "js": nonneg, "tox": pos, "nsub": pos, "nss": any, "nfs": nonneg,
	"tpg": {min: -1, max: 1}, "xj": nonneg, "ld": nonneg, "uo": pos, "u0": pos, "ucrit": pos, "uexp": nonneg,
	"utra": {min: 0, max: 0.5}, "vmax": nonneg, "neff": pos, "kf": nonneg, "af": pos, "fc": {min: 0, max: 1, openMax: true},
	"delta": nonneg, "theta": nonneg, "eta": nonneg, "kappa": nonneg, "xqc": {min: 0, max: 1},
	"wd": nonneg, "tnom": {min: -273.15, max: math.Inf(1), openMin: true},
}

var modelParams = map[string]map[string]paramRule{"d": diodeParams, "npn": bjtParams, "pnp": bjtParams,
	"nmos": mosParams, "pmos": mosParams}

func (m Model) validate() error {
	if m.Type == "machine" {
		_, err := ParseMachine(m)
		return err
	}
	rules, ok := modelParams[m.Type]
	if !ok {
		return fmt.Errorf("%w: %s has unsupported type %q", ErrBadModel, m.Name, m.Type)
	}
	keys := make([]string, 0, len(m.Params))
	for k := range m.Params {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		v := m.Params[k]
		r, ok := rules[k]
		if !ok {
			return fmt.Errorf("%w: %s: unknown %s parameter %q", ErrBadModel, m.Name, strings.ToUpper(m.Type), k)
		}
		if math.IsNaN(v) || math.IsInf(v, 0) ||
			v < r.min || (r.openMin && v == r.min) || v > r.max || (r.openMax && v == r.max) {
			return fmt.Errorf("%w: %s: %s = %g is out of range", ErrBadModel, m.Name, k, v)
		}
	}
	return nil
}

// kindModelType maps element kinds to the model type they need.
var kindModelType = map[Kind]string{Diode: "d"}
