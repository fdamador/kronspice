// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
)

// ErrBadCheckpoint is returned for a checkpoint that cannot be resumed.
var ErrBadCheckpoint = errors.New("evo: invalid checkpoint")

// checkpointVersion is bumped when the State layout changes.
const checkpointVersion = 1

type checkpointFile struct {
	Version int    `json:"version"`
	State   *State `json:"state"`
}

// jsonFloat stores non-finite values (none are expected, but penalties and
// infinities must never make a checkpoint unwritable).
func sanitize(st *State) {
	fix := func(v *float64) {
		if math.IsNaN(*v) || math.IsInf(*v, 0) {
			*v = PenaltyOther
		}
	}
	for i := range st.Population {
		fix(&st.Population[i].Fitness)
	}
	fix(&st.Best.Fitness)
	for i := range st.Front {
		fix(&st.Front[i].Fitness)
	}
	for i := range st.History {
		h := &st.History[i]
		fix(&h.Best)
		fix(&h.Mean)
		fix(&h.Worst)
	}
}

// WriteCheckpoint writes st as JSON. float64 values round-trip exactly
// through encoding/json, so a resumed run is bit-identical.
func WriteCheckpoint(w io.Writer, st *State) error {
	cp := *st
	cp.Population = append([]Individual(nil), st.Population...)
	cp.History = append([]GenStats(nil), st.History...)
	cp.Front = append([]Individual(nil), st.Front...)
	sanitize(&cp)
	enc := json.NewEncoder(w)
	enc.SetIndent("", " ")
	return enc.Encode(checkpointFile{Version: checkpointVersion, State: &cp})
}

// SaveCheckpoint writes st to path atomically (temporary file + rename), so
// an interrupted run never leaves a truncated checkpoint.
func SaveCheckpoint(path string, st *State) error {
	tmp, err := os.CreateTemp(filepath.Dir(path), ".checkpoint-*")
	if err != nil {
		return err
	}
	if err := WriteCheckpoint(tmp, st); err != nil {
		tmp.Close()
		os.Remove(tmp.Name())
		return err
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmp.Name())
		return err
	}
	return os.Rename(tmp.Name(), path)
}

// ReadCheckpoint reads a checkpoint written by WriteCheckpoint.
func ReadCheckpoint(r io.Reader) (*State, error) {
	var cf checkpointFile
	dec := json.NewDecoder(r)
	if err := dec.Decode(&cf); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrBadCheckpoint, err)
	}
	if cf.Version != checkpointVersion || cf.State == nil || len(cf.State.Population) == 0 {
		return nil, fmt.Errorf("%w: version %d", ErrBadCheckpoint, cf.Version)
	}
	return cf.State, nil
}

// LoadCheckpoint reads a checkpoint file.
func LoadCheckpoint(path string) (*State, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ReadCheckpoint(f)
}
