# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H21 report: a ring of charge as one of Kron's machines, as a PDF.

    python examples/kron_electron/viz/make_report.py examples/kron_electron/kron_electron_report.pdf

All numbers come from examples/kron_electron/README.md (the example has no
run_output.txt); figures are viz/ring_field.png and the first frame of
viz/ring_motion.gif. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kron_electron", "A ring of charge as a Kron machine",
               "Do Kron's machine equations model the ring electron? — register entry H21",
               "register entry H21")

    r.h1("Summary")
    r.p("A closed ring of circulating charge is written as one of Gabriel Kron's machines: a winding whose "
        "current is i = q·dθ/dt/2π, with the ring's mass entering Kron's metric tensor as an inductance. The example "
        "first checks this machine against Kelvin's 1899 charged ring, solved independently with the Lorentz force, "
        "then fixes a ring by the electron's charge and spin and compares its moment and g-factor with the measured "
        "electron. The hypothesis (H21) is that the electron <i>is</i> such a ring; the standard prediction beside "
        "every row is the Dirac electron with QED (g = 2.002 319 304 36, μ = 1.001 16 μ<sub>B</sub>).")
    r.bullets([
        "<b>Twin test passes.</b> Kelvin's ring as a three-winding machine gives moments of momentum "
        "−/+1.11207×10<super>−10</super> kg m²/s (= eΦ/2π), the fluid at Larmor's rate 4.44827×10<super>−6</super> rad/s "
        "and kinetic energy 4.94677×10<super>−16</super> J, against 4.94683×10<super>−16</super> J from the "
        "Lorentz-force bench.",
        "<b>Torque from the metric is −μB sin θ.</b> The tilting ring's period is 19.873 s against the pendulum value "
        "4K(k)/ω<sub>0</sub> = 19.872 s (2×10<super>−5</super>); ledger mismatch 6×10<super>−10</super> J.",
        "<b>A ring whose mass circulates with its charge has g ≤ 1:</b> g = 1/(1 + L·e²/4π²·m·a²), so it precesses at "
        "Larmor's eB/2m, half the electron's spin precession.",
        "<b>Only a massless ring reaches the measured moment</b> at S = ħ/2, at radius 8.86×10<super>−15</super> m, "
        "with v = 43.6 c and a magnetic energy of 475 m<sub>e</sub>c².",
        "<b>Every ring is far larger than the 10<super>−18</super> m scattering bound</b> on the electron's size.",
        "<b>Status: bounded.</b> The machine reproduces the classical ring exactly and inherits its obstacles "
        "(g ≤ 1 or v &gt; c).",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis H21</b> (user hypothesis, 2026-10-31): a closed ring of circulating charge is a one-winding "
        "Kron machine, with self-inductance L, the circulating mass as rotor inertia and an external field entering "
        "as speed voltages; voltage and current are the pressure and flow of energy at that scale, and the electron "
        "is such a ring. The established part is Kron's own: his eq. 11 (<i>Equivalent Circuits of Electric "
        "Machinery</i>, 1951, Epilogue p. 259) states that with e = R = 0 the machine equation is the motion of a "
        "charged particle in an electromagnetic field.")
    r.p("<b>Standard prediction:</b> the Dirac electron with QED, g = 2.002 319 304 36 and μ = 1.001 16 μ<sub>B</sub>; "
        "for a classical ring, g = 1 when mass circulates with charge, a Compton-scale radius for μ<sub>B</sub> at "
        "v ≤ c, and no stable state. Scattering bounds the electron below 10<super>−18</super> m.")

    r.h1("The model")
    r.p("For a ring of charge q and mass m turning at dθ/dt, the current is i = q·dθ/dt/2π, so the electrical and the "
        "mechanical coordinate are the same axis. The kinetic energy ½·m·a²·dθ/dt² = ½·(m·a²·(2π/q)²)·i² is then a "
        "magnetic energy on that axis, and the inertia enters Kron's metric as an inductance m·a²·(2π/q)², beside "
        "the ring's own self-inductance L. Kron defined the metric as \"all the self- and mutual inductances and the "
        "rotor inertia\" (1951, Epilogue p. 254, eqs. 3–4). Everything uses the N machine element "
        "(device/machine.go) unchanged.")
    r.table([
        ["quantity", "formula / value", "source"],
        ["ring self-inductance", "L = μ<sub>0</sub>·a·(ln(8a/r) − 7/4)", "Maxwell, <i>Treatise</i> Art. 706"],
        ["ring capacitance", "thin-ring formula", "Smythe §5.17"],
        ["inertial inductance", "m·a²·(2π/q)²", "Kron's metric"],
        ["coil–ring mutual", "em.MutualCircles (elliptic form)", "Maxwell"],
        ["twin test: coil", "6.201×10<super>−7</super> H", "printed metric"],
        ["twin test: coil–ring mutual", "3.494×10<super>−8</super> H", "printed metric"],
        ["twin test: tube, fluid inertial inductance", "9.87×10<super>8</super> H each", "printed metric"],
        ["twin test: ring's own L", "3.1×10<super>−7</super> H (dropped as in Kelvin, or kept: no change)", ""],
    ], [0.32, 0.43, 0.25])

    r.h1("1. Twin test: Kelvin's ring in Kron's frame (standard)")
    r.p("Kelvin's ring (tube +e, fluid −e, each 10 g, radius 5 cm; coil of 20 kA-turns, radius 10 cm, switched on "
        "over 1 s, 5 cm below) becomes a three-winding machine: coil, tube and fluid. Tube and fluid windings are "
        "closed through 0 V sources (frictionless), so their flux linkages are conserved.")
    r.table([
        ["quantity", "Kron machine", "Lorentz-force bench"],
        ["flux through the ring", "6.98732×10<super>−4</super> Wb (Maxwell elliptic)",
         "6.98736×10<super>−4</super> Wb (Biot–Savart)"],
        ["moment of momentum, tube", "−1.11207×10<super>−10</super> kg m²/s", "−1.11207×10<super>−10</super>"],
        ["moment of momentum, fluid", "+1.11207×10<super>−10</super>", "+1.11207×10<super>−10</super>"],
        ["Kelvin's eΦ/2π", "1.11207×10<super>−10</super>", ""],
        ["fluid angular velocity", "4.44827×10<super>−6</super> rad/s", "4.44829×10<super>−6</super>"],
        ["Larmor e·B<sub>avg</sub>/2w", "4.44827×10<super>−6</super>", ""],
        ["kinetic energy", "4.94677×10<super>−16</super> J", "4.94683×10<super>−16</super> J"],
        ["source work against the ring",
         "= kinetic + mutual-field energy (−9.89×10<super>−16</super> J), residual 1.5×10<super>−13</super>",
         "residual 8×10<super>−16</super>"],
    ], [0.28, 0.42, 0.30])
    r.p("The 1.1×10<super>−5</super> difference in kinetic energy is the bench's Biot–Savart flux against Maxwell's "
        "elliptic value used here. The engine's own ledger also closes (source 124.028 J, machine's ½iᵀLi 124.02 J, "
        "mismatch 7.8×10<super>−3</super> J over the ramp's breakpoints), but it is dominated by the coil's own 124 J "
        "of field energy, which is why the bench's bookkeeping is reported.")

    r.h1("2. Torque from the metric (standard)")
    r.p("The same ring, spinning at 1000 rad/s (current 0.159 mA, moment 1.25×10<super>−6</super> A·m²), tilts about a "
        "diameter in 1 T. In Kron's frame I, M(θ) = B·π·a²·cos θ couples a field winding and the ring, and the tilt "
        "inertia ½·m·a² is a capacitor on the mechanical port. Kron's torque ½·iᵀ·(∂L/∂θ)·i is −μ·B·sin θ. Released "
        "at 0.05 rad, the oscillation has period 19.873 s against the pendulum value 4K(k)/ω<sub>0</sub> = 19.872 s "
        "(2×10<super>−5</super>). The ring current drifts by 6×10<super>−11</super> (flux conservation) and the ledger "
        "mismatch is 6×10<super>−10</super> J.")
    r.p("The ring current has no harmonics (a uniform ring's field is static). The tilt is a nonlinear pendulum, "
        "and the solver's ω(t) carries exactly the pendulum's odd series (Jacobi nome q = e<super>−πK′/K</super>, "
        "k = sin θ<sub>0</sub>/2), measured by leak-free Fourier sums over two periods:")
    r.table([
        ["θ<sub>0</sub>", "fundamental", "ω<sub>3</sub>/ω<sub>1</sub> run", "exact", "ω<sub>5</sub>/ω<sub>1</sub> run",
         "exact", "ω<sub>2</sub>/ω<sub>1</sub>"],
        ["0.05 rad", "0.05032 Hz", "3.906×10<super>−5</super>", "3.907×10<super>−5</super>", "2×10<super>−8</super>",
         "1.5×10<super>−9</super>", "2×10<super>−8</super>"],
        ["1 rad", "0.0472 Hz", "0.01658", "0.01658", "2.70×10<super>−4</super>", "2.71×10<super>−4</super>",
         "2×10<super>−8</super>"],
        ["2.5 rad", "0.03063 Hz", "0.1597", "0.1598", "0.0225", "0.0225", "8×10<super>−7</super>"],
    ], [0.11, 0.14, 0.15, 0.15, 0.15, 0.15, 0.15])
    r.p("Even harmonics are absent by the symmetry θ → −θ.")

    r.h1("3. The electron as the ring (hypothesis)")
    r.p("Charge e; the total angular momentum of the metric is set to ħ/2 and the rest follows (core radius a/100). "
        "KE is the kinetic energy, ½LI² the magnetic energy and W<sub>es</sub> the electrostatic energy, each in units "
        "of m<sub>e</sub>c².")
    r.table([
        ["ring", "a [m]", "v/c", "μ/μ<sub>B</sub>", "g", "KE", "½LI²", "W<sub>es</sub>"],
        ["Compton ring, electron mass circulating", "3.86×10<super>−13</super>", "0.494", "0.494", "0.989", "0.122",
         "0.0014", "0.0078"],
        ["Half-Compton ring, mass circulating", "1.93×10<super>−13</super>", "0.978", "0.489", "0.978", "0.478",
         "0.011", "0.0155"],
        ["Compton ring, massless (field inertia)", "3.86×10<super>−13</super>", "43.6", "43.6", "87.2", "0", "10.9",
         "0.0078"],
        ["Classical-radius ring, massless", "2.82×10<super>−15</super>", "43.6", "0.318", "0.637", "0", "1490", "1.06"],
        ["Massless ring fitted to the measured moment", "8.86×10<super>−15</super>", "43.6", "1.001", "2.002", "0",
         "475", "0.338"],
        ["<b>Measured (Dirac + QED)</b>", "&lt; 10<super>−18</super> (scattering)", "", "1.00116", "2.00232", "", "", ""],
    ], [0.28, 0.15, 0.07, 0.09, 0.09, 0.08, 0.08, 0.08], bold_rows=(6,))
    r.bullets([
        "A ring whose mass circulates with its charge has g = 1/(1 + L·e²/4π²·m·a²) ≤ 1: its precession in a field is "
        "Larmor's eB/2m, half the electron's spin precession.",
        "Only a massless ring (inertia from its own field) reaches the measured moment at S = ħ/2, at 8.9 fm with "
        "v = 43.6 c and a magnetic energy of 475 m<sub>e</sub>c². Every massless thin ring with S = ħ/2 moves at "
        "v/c = π/(2α(ln(8a/r) − 7/4)), fixed by the fine-structure constant α; at v ≤ c (a ≥ ħ/2m<sub>e</sub>c) its "
        "moment is far below μ<sub>B</sub>.",
        "Every ring radius lies far above the 10<super>−18</super> m scattering bound.",
    ])

    r.h2("Natural frequencies of the rings (S = ħ/2, B = 1 T)")
    r.table([
        ["ring", "rotation of the charge", "precession (ESR)", "tilt libration", "LC self-resonance",
         "wave once round at c"],
        ["Compton ring, mass circulating", "6.1×10<super>19</super> Hz", "13.8 GHz", "1.3×10<super>15</super> Hz",
         "2.3×10<super>19</super> Hz", "1.24×10<super>20</super> Hz"],
        ["Half-Compton ring, mass circulating", "2.4×10<super>20</super>", "13.7 GHz", "2.6×10<super>15</super>",
         "4.6×10<super>19</super>", "2.47×10<super>20</super>"],
        ["Compton ring, massless", "5.4×10<super>21</super>", "1.22 THz", "—", "2.3×10<super>19</super>",
         "1.24×10<super>20</super>"],
        ["Classical-radius ring, massless", "7.4×10<super>23</super>", "8.9 GHz", "—", "3.1×10<super>21</super>",
         "1.7×10<super>22</super>"],
        ["Massless ring fitted to the moment", "2.4×10<super>23</super>", "28.0 GHz", "—", "1.0×10<super>21</super>",
         "5.4×10<super>21</super>"],
        ["<b>Measured / reference</b>", "", "<b>28.025 GHz</b> (ESR)", "", "",
         "<b>1.236×10<super>20</super> Hz</b> (m<sub>e</sub>c²/h)"],
    ], [0.27, 0.14, 0.15, 0.13, 0.15, 0.16], bold_rows=(6,))
    r.bullets([
        "<b>Precession</b> is the spin's magnetic-resonance frequency g·eB/4πm, measured in examples/kron_nmr; only "
        "the fitted massless ring reaches the measured value, by construction.",
        "<b>Rotation</b> is the charge's turning rate, the current's fundamental: about half the Compton frequency for "
        "the Compton ring (v ≈ c/2).",
        "<b>Tilt libration</b> √(μB/J)/2π is optical (10<super>15</super> Hz) for an electron-mass ring in 1 T.",
        "<b>Wave once round at c</b> equals the Compton frequency m<sub>e</sub>c²/h for the Compton ring and "
        "2m<sub>e</sub>c²/h (Zitterbewegung) for the half-Compton ring: identities of the chosen radius, not predictions.",
        "There is no radial breathing mode: the model has no restoring force for the ring's radius.",
        "Kelvin's 5 cm ring of part 2 (1 µC, 10 g, 1000 rad/s, 1 T): rotation 159 Hz, precession eB/4πm = "
        "8×10<super>−6</super> Hz, libration 0.0503 Hz (the 19.9 s period), LC self-resonance 177 MHz.",
    ])

    r.h2("Against existing measurements")
    r.table([
        ["quantity in the model", "measured value", "match?"],
        ["Precession of a ring with circulating mass, eB/4πm = 14.0 GHz/T (g = 1)",
         "Normal Zeeman effect: orbital moments precess at Larmor's eB/4πm, g<sub>L</sub> = 1 (Zeeman 1897; Larmor 1897)",
         "<b>Yes</b> for orbital moments; not the electron's spin (g = 2.0023)."],
        ["Precession of the fitted massless ring, 28.02 GHz/T",
         "Electron spin resonance; g/2 = 1.001 159 652 180 59(13) (Fan et al. 2023; Hanneke et al. 2008)",
         "By construction only; needs v = 43.6 c."],
        ["Pendulum period 4K(k)/ω<sub>0</sub> and odd harmonics of the tilt",
         "Large-amplitude pendulum measurements agree with the elliptic-integral period to about "
         "10<super>−4</super> (Nelson &amp; Olsson 1986)",
         "<b>Yes</b>, for the classical pendulum."],
        ["Wave once round the Compton ring, 1.236×10<super>20</super> Hz; 2m<sub>e</sub>c²/h for the half-Compton ring",
         "Identity of the radius. A channeling dip at m<sub>e</sub>c²/h was reported (Gouanère et al. 2008), not "
         "reproduced; Zitterbewegung simulated with a trapped ion (Gerritsma et al. 2010), not seen in an electron",
         "No measurement on the electron."],
        ["Charge rotation, tilt libration, LC self-resonance, modal series n·c/2πa",
         "None: a uniform ring current radiates nothing at these frequencies", "No measurement."],
    ], [0.32, 0.45, 0.23])
    r.p("The two genuine matches are classical: the ring's precession is the orbital Larmor frequency of the normal "
        "Zeeman effect, and its tilt is the measured large-amplitude pendulum. Neither bears on the electron's spin, "
        "which is what the g ≤ 1 bound says.")

    r.h1("4. Figures")
    r.figure(os.path.join(HERE, "ring_field.png"),
             "Figure 1. Field lines of the ring current from Biot–Savart. Left: the meridional plane, the wire as two "
             "red dots and the axis vertical. Right: oblique view at eight azimuths, the ring in red. Every line lies in "
             "a meridional plane and threads the ring once: dipole far field, near field circling the wire. Lines "
             "started at 0.3, 0.6 and 0.88 of the radius close to 10<super>−9</super> (TestFieldLinesClose).")
    r.p("A single ring current cannot produce field lines that wind around a sphere, as drawn in some toroidal "
        "electron pictures; that needs a current component along the ring's axis, a second winding in Kron's frame.")
    r.figure(os.path.join(HERE, "ring_motion.gif"),
             "Figure 2. First frame of viz/ring_motion.gif (90 frames). Left: Kelvin's switch-on seen from above, tube "
             "marker red and fluid marker blue (angles scaled; the real ones are microradians), coil current as the "
             "green bar. Middle: the tilting ring edge-on (black bar, wire ends red) in the vertical field (grey "
             "arrows), its magnetic moment as the red arrow; the 0.05 rad swing drawn ×8. Right: its dipole field "
             "lines, which rock with it at the pendulum period. The animation is drawn from the solver's own time "
             "points.", width_frac=0.95)

    r.h1("Validation")
    r.table([
        ["test", "check", "held to"],
        ["TestKelvinTwin", "moments −/+eΦ/2π, fluid at Larmor's rate, bench ledger closed, with and without the ring's "
         "self-inductance; against the Lorentz-force bench", "1e-6; bench 1e-5, kinetic energy 2e-5"],
        ["TestTiltPeriod", "period against 4K(k)/ω<sub>0</sub>; ring current conserved", "2e-4; drift 1e-9"],
        ["TestTiltHarmonics", "3rd and 5th harmonics against the Jacobi series at θ<sub>0</sub> = 1 rad; no even "
         "harmonics; small-angle 3rd = θ<sub>0</sub>²/64", "2e-4; 1e-3"],
        ["TestElectronRing", "g·(1 + Le²/4π²ma²) = 1; precession = g × Larmor; Compton ring below c; fitted ring "
         "carries the measured moment and is superluminal; every scenario above the size bound", "1e-12; 1e-9"],
        ["TestMasslessSpeed", "v/c = π/(2α(ln(8a/r) − 7/4)) at a = 1 fm, ħ/m<sub>e</sub>c and 0.1 nm", "1e-9"],
        ["TestFrequencies", "c/2πa = m<sub>e</sub>c²/h for the Compton ring; precession = g·eB/4πm; fitted ring at "
         "the measured ESR; Kelvin ring libration period 19.869 s", "1e-9; 1e-12; 1e-4; 0.01 s"],
        ["TestFieldLinesClose", "field lines of the ring current close", "1e-9"],
    ], [0.22, 0.56, 0.22])

    r.h1("Limits of the model")
    r.bullets([
        "The N element has one angle, so the ring cannot precess in 3-D; the Zeeman polarisation test stays with the "
        "3-D bench of examples/kelvin_ring.",
        "No restoring force holds the ring's radius (the electrostatic hoop stress is unbalanced); a breathing mode "
        "would need a declared stiffness.",
        "Thin-ring formulas for L and C, with a core radius of a/100.",
    ])

    r.h1("Status")
    r.p("<b>H21: bounded.</b> Step 1 (the twin test) passes; step 2 is computed: no ring with the electron's "
        "charge-to-mass ratio reaches g = 2 (g ≤ 1, Larmor precession eB/2m), and the massless ring that carries the "
        "measured moment at S = ħ/2 needs v = 43.6 c at a = 8.9 fm. The standard prediction (Dirac + QED) is given "
        "beside every row. Tewari's space-vortex electron (H13) on the same bench gives g = 15/8.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/kron_electron -png viz/ring_field.png\n"
           "go run ./examples/kron_electron -gif viz/ring_motion.gif\n"
           "go test ./examples/kron_electron\n"
           "KRONSPICE_EXPORT=1 go test ./examples/kron_electron -run TestExport\n"
           "python examples/kron_electron/viz/make_report.py examples/kron_electron/kron_electron_report.pdf")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kron_electron_report.pdf"))
