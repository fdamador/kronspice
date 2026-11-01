// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

type denseStub struct {
	n int
	a []float64
}

func (d *denseStub) Add(i, j int, v float64) {
	if i < 0 || j < 0 {
		return
	}
	d.a[i*d.n+j] += v
}

func allDevices() []Device {
	return []Device{
		NewResistor(circuit.Element{Name: "R1", Kind: circuit.Resistor, Nodes: [2]int{1, 0}, Value: 100}),
		NewResistor(circuit.Element{Name: "RN", Kind: circuit.Resistor, Nodes: [2]int{1, 0}, Value: -100}),
		NewVoltageSource(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, Nodes: [2]int{1, 0}, Value: 5}),
		NewCurrentSource(circuit.Element{Name: "I1", Kind: circuit.CurrentSource, Nodes: [2]int{1, 0}, Value: 1e-3}),
		NewCapacitor(circuit.Element{Name: "C1", Kind: circuit.Capacitor, Nodes: [2]int{1, 0}, Value: 1e-6}),
		NewInductor(circuit.Element{Name: "L1", Kind: circuit.Inductor, Nodes: [2]int{1, 0}, Value: 1e-3}),
		NewBehavioral(circuit.Element{Name: "B1", Kind: circuit.Behavioral, Nodes: [2]int{1, 0}, BOut: 'I', BCtrl: 'V',
			CtrlNodes: [2]int{1, 0}, Curve: []float64{0, 0, 1, 1e-3, 2, 0}}, true),
		NewBehavioral(circuit.Element{Name: "B2", Kind: circuit.Behavioral, Nodes: [2]int{1, 0}, BOut: 'V', BCtrl: 'V',
			CtrlNodes: [2]int{1, 0}, Curve: []float64{0, 0, 1, 1}}, false),
		NewDiode(circuit.Element{Name: "D1", Kind: circuit.Diode, Nodes: [2]int{1, 0}, Model: "d"},
			circuit.Model{Name: "d", Type: "d", Params: map[string]float64{"rs": 5, "cjo": 1e-12, "tt": 1e-9}}),
		NewControlled(circuit.Element{Name: "E1", Kind: circuit.VCVS, Nodes: [2]int{1, 0}, Value: 2, CtrlNodes: [2]int{1, 0}}),
		NewControlled(circuit.Element{Name: "G1", Kind: circuit.VCCS, Nodes: [2]int{1, 0}, Value: 1e-3, CtrlNodes: [2]int{1, 0}}),
		NewControlled(circuit.Element{Name: "F1", Kind: circuit.CCCS, Nodes: [2]int{1, 0}, Value: 3, CtrlSource: "v1"}),
		NewControlled(circuit.Element{Name: "H1", Kind: circuit.CCVS, Nodes: [2]int{1, 0}, Value: 50, CtrlSource: "v1"}),
		NewBJT(circuit.Element{Name: "Q1", Kind: circuit.BJT, Nodes: [2]int{1, 0}, Extra: []int{1}, Model: "q"},
			circuit.Model{Name: "q", Type: "npn", Params: map[string]float64{"rb": 50, "rc": 5, "re": 1, "cje": 1e-12,
				"cjc": 1e-12, "tf": 1e-10, "cjs": 1e-12, "xcjc": 0.5}}),
		mustMOS(),
		NewMutual(circuit.Element{Name: "K1", Kind: circuit.Mutual, Value: 0.5, Couples: [2]string{"l1", "l2"}}, 1e-3, 4e-3),
		NewGyrator(circuit.Element{Name: "Y1", Kind: circuit.Gyrator, Nodes: [2]int{1, 0}, Extra: []int{2, 0}, Value: 0.5}),
		NewMachine(circuit.Element{Name: "N1", Kind: circuit.Machine, Nodes: [2]int{1, 0}, Extra: []int{2, 0, 3, 0}}, testMachineSpec()),
		NewNonlinearCapacitor(circuit.Element{Name: "CN", Kind: circuit.NonlinearCapacitor, Nodes: [2]int{1, 0},
			Curve: []float64{-1, -2e-6, 0, 0, 1, 1e-6}, Scale: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, 1e3}}}),
		NewNonlinearCapacitor(circuit.Element{Name: "CV", Kind: circuit.NonlinearCapacitor, Nodes: [2]int{1, 0}, Value: 1e-6,
			Scale: circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{1, 2}}}),
		NewNonlinearInductor(circuit.Element{Name: "LN", Kind: circuit.NonlinearInductor, Nodes: [2]int{1, 0},
			Curve: []float64{-1, -1e-3, 0, 0, 1, 1e-3, 2, 1.1e-3}, Scale: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, 1e3}}}),
		NewNonlinearInductor(circuit.Element{Name: "LV", Kind: circuit.NonlinearInductor, Nodes: [2]int{1, 0}, Value: 1e-3,
			Scale: circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{1, 2}}}),
		NewBearing(circuit.Element{Name: "P1", Kind: circuit.BearingMotor, Nodes: [2]int{1, 0}, Extra: []int{2, 0}, Value: 1e-3,
			Params: map[string]float64{"kappa": 1e-6, "w0": 10}}),
	}
}

// TestAllChannelsPresent enforces CLAUDE.md §3.1 rule 1 for every device.
func TestAllChannelsPresent(t *testing.T) {
	for _, d := range allDevices() {
		ports := d.Ports()
		if len(ports) != exchange.NumChannels {
			t.Errorf("%s: %d ports, want %d", d.Name(), len(ports), exchange.NumChannels)
		}
		for i, p := range ports {
			if p.Channel != exchange.ChannelID(i) {
				t.Errorf("%s: port %d is channel %v", d.Name(), i, p.Channel)
			}
		}
	}
}

// TestUnmodeledPortsAreZero enforces §3.1 rule 3.
func TestUnmodeledPortsAreZero(t *testing.T) {
	op := OperatingPoint{Voltage: 3, Current: 2, PowerAbsorbed: 6}
	for _, d := range allDevices() {
		for _, p := range d.Ports() {
			if p.Status == exchange.NotModeled {
				if got := d.ExteriorPower(p.Channel, op); got != 0 {
					t.Errorf("%s: unmodeled %v reports %g W", d.Name(), p.Channel, got)
				}
			}
		}
	}
}

func TestResistorStamp(t *testing.T) {
	r := NewResistor(circuit.Element{Name: "R1", Kind: circuit.Resistor, Nodes: [2]int{1, 2}, Value: 4})
	ctx := NewSetupContext(3)
	_ = r.Setup(ctx)
	m := &denseStub{n: 2, a: make([]float64, 4)}
	_ = r.Load(&LoadContext{M: m, RHS: make([]float64, 2)})
	want := []float64{0.25, -0.25, -0.25, 0.25}
	for i := range want {
		if m.a[i] != want[i] {
			t.Fatalf("stamp = %v, want %v", m.a, want)
		}
	}
	op := r.Operating([]float64{2, 0})
	if op.Current != 0.5 || op.PowerAbsorbed != 1 {
		t.Fatalf("op = %+v", op)
	}
	if got := r.ExteriorPower(exchange.Thermal, op); got != -1 {
		t.Fatalf("thermal out = %g, want -1", got)
	}
}

func TestVoltageSourceStampToGround(t *testing.T) {
	s := NewVoltageSource(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, Nodes: [2]int{1, 0}, Value: 5})
	ctx := NewSetupContext(2)
	_ = s.Setup(ctx)
	if ctx.Size() != 2 {
		t.Fatalf("size = %d", ctx.Size())
	}
	m := &denseStub{n: 2, a: make([]float64, 4)}
	rhs := make([]float64, 2)
	_ = s.Load(&LoadContext{M: m, RHS: rhs})
	if m.a[0*2+1] != 1 || m.a[1*2+0] != 1 || rhs[1] != 5 {
		t.Fatalf("stamp = %v rhs = %v", m.a, rhs)
	}
}

func TestLoadDoesNotAllocate(t *testing.T) {
	devs := allDevices()
	ctx := NewSetupContext(2)
	for _, d := range devs {
		_ = d.Setup(ctx)
	}
	n := ctx.Size()
	lc := &LoadContext{M: &denseStub{n: n, a: make([]float64, n*n)}, RHS: make([]float64, n)}
	allocs := testing.AllocsPerRun(100, func() {
		for _, d := range devs {
			_ = d.Load(lc)
		}
	})
	if allocs != 0 {
		t.Fatalf("Load allocates %v times", allocs)
	}
}

func mustMOS() Device {
	q, err := NewMOSFET(circuit.Element{Name: "MX", Kind: circuit.MOSFET, Nodes: [2]int{1, 0}, Extra: []int{1, 0},
		Model: "m", Params: map[string]float64{"ad": 1e-10, "as": 1e-10}},
		circuit.Model{Name: "m", Type: "nmos", Params: map[string]float64{"tox": 2e-8, "rd": 5, "rs": 5, "cj": 1e-4}})
	if err != nil {
		panic(err)
	}
	return q
}
