# 2. Evolution

*Planned.* How kronspice serves as the fitness engine of a genetic
algorithm (CLAUDE.md §4): value genomes with E-series snapping, topology
genomes grown from an empty embryo, NSGA-II fronts, penalty classes for
every way a circuit can fail, determinism from a seed, and diakoptics
caching the pieces of a circuit that a mutation did not touch.

Examples: [amplifier_gain](../examples/amplifier_gain),
[filter_topology](../examples/filter_topology),
[mos_amplifier_pareto](../examples/mos_amplifier_pareto),
[torn_ga](../examples/torn_ga).

Planned sections: the value genome and the 2N3904 stage; why a failed
candidate gets a penalty and not a crash; the topology genome and the
elliptic low-pass grown from nothing; constrained dominance and the
gain / bandwidth / power front; tearing a circuit so only the mutated
piece is re-solved.
