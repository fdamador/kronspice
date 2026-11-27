// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strings"
)

// Switch model parameters (ngspice manual, "Switch Model (SW/CSW)"):
// threshold VT [V] or IT [A], hysteresis VH or IH, on and off resistance
// RON, ROFF [Ω]. A negative hysteresis selects ngspice's smooth switch,
// whose resistance moves continuously from ROFF to RON as the control
// passes from VT − |VH| to VT + |VH|.
var swParams = map[string]paramRule{"vt": any, "vh": any, "ron": pos, "roff": pos}
var cswParams = map[string]paramRule{"it": any, "ih": any, "ron": pos, "roff": pos}

func init() {
	modelParams["sw"] = swParams
	modelParams["csw"] = cswParams
	kindModelType[Switch] = "sw"
	kindModelType[CurrentSwitch] = "csw"
}

// SwitchDefaults are the ngspice defaults: VT = IT = 0, VH = IH = 0,
// RON = 1 Ω, ROFF = 1/GMIN = 1e12 Ω.
const (
	SwitchRonDefault  = 1.0
	SwitchRoffDefault = 1e12
)

// SwitchParams returns a switch model's threshold, hysteresis, on and off
// resistance with the defaults filled in.
func SwitchParams(m Model) (vt, vh, ron, roff float64) {
	p := m.Params
	get := func(k string, d float64) float64 {
		if v, ok := p[k]; ok {
			return v
		}
		return d
	}
	if m.Type == "csw" {
		vt, vh = get("it", 0), get("ih", 0)
	} else {
		vt, vh = get("vt", 0), get("vh", 0)
	}
	return vt, vh, get("ron", SwitchRonDefault), get("roff", SwitchRoffDefault)
}

func (c *Circuit) validateSwitch(e Element) error {
	m, ok := c.ModelFor(e)
	if !ok {
		return fmt.Errorf("%w: %s uses undefined model %q", ErrBadModel, e.Name, e.Model)
	}
	if want := kindModelType[e.Kind]; m.Type != want {
		return fmt.Errorf("%w: %s needs a %s model, %q is type %s", ErrBadModel, e.Name,
			strings.ToUpper(want), e.Model, strings.ToUpper(m.Type))
	}
	if _, vh, ron, roff := SwitchParams(m); vh < 0 && ron == roff {
		return fmt.Errorf("%w: %s: a smooth switch (negative hysteresis) needs RON ≠ ROFF", ErrBadModel, e.Model)
	}
	if e.HasIC && e.IC != 0 && e.IC != 1 {
		return fmt.Errorf("%w: %s initial state must be ON (1) or OFF (0)", ErrBadValue, e.Name)
	}
	if e.Wave.Kind != WaveNone || len(e.Extra) != 0 || math.IsNaN(e.Value) {
		return fmt.Errorf("%w: %s", ErrBadValue, e.Name)
	}
	return c.validateControl(e, e.Kind == Switch)
}

// switchCtrl returns the control nodes of a voltage-controlled switch.
func (c *Circuit) switchCtrl(e Element) []int {
	if e.Kind == Switch {
		return e.CtrlNodes[:]
	}
	return nil
}

// VSwitch adds a voltage-controlled switch (SPICE S) between a and b,
// controlled by V(cp) − V(cn), with a model of type "sw".
func (b *Builder) VSwitch(name, a, c, cp, cn, model string) *Builder {
	return b.Add(Element{Name: name, Kind: Switch, Model: strings.ToLower(model),
		CtrlNodes: [2]int{b.Node(cp), b.Node(cn)}}, a, c)
}

// ISwitch adds a current-controlled switch (SPICE W) between a and b,
// controlled by the current of voltage source vsrc, with a model of type
// "csw".
func (b *Builder) ISwitch(name, a, c, vsrc, model string) *Builder {
	return b.Add(Element{Name: name, Kind: CurrentSwitch, Model: strings.ToLower(model),
		CtrlSource: strings.ToLower(vsrc)}, a, c)
}
