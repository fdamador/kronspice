# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""One rotor in three gravitational fields, as a PDF report.

Run from the repository root:
    python examples/machine_spacetime/viz/make_report.py examples/machine_spacetime/machine_spacetime_report.pdf

Numbers come from examples/machine_spacetime/run_output.txt (the program's
own output) and README.md. The example has no figures. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/machine_spacetime", "One rotor in three gravitational fields",
               "What the clocks, spin axis and a light loop of a fast rotor show in Earth's field, in orbit and near a "
               "neutron star, and how much the rotor itself curves spacetime. CLAUDE.md §3l",
               "CLAUDE.md §3l",
               "T1000G disk, 0.5 m radius, at half its burst speed")

    r.h1("Summary")
    r.p("A 0.5 m carbon-fibre (Toray T1000G) disk, 2 cm thick, spins at half its burst speed and carries a 1 m² fibre "
        "loop. The example combines rotation (special relativity) and curvature (general relativity) for this one "
        "machine at three sites, using the <i>spacetime/</i> package, and then asks the reverse question: how much does "
        "the machine itself curve and twist spacetime? It is a standard-physics example; the one hypothesis on this "
        "bench, the rotor energy-shift anomaly H17, is computed in <i>spacetime.RotorShift</i> with the standard "
        "k = 1/2 beside a declared k.")
    r.bullets([
        "<b>The rim loses 4.25×10<super>−7</super> s/day (0.425 µs/day) against the hub at every site</b>, because the "
        "lag depends only on the rim speed (941 m/s); the hub rate carries the gravity.",
        "<b>Gravity Probe B check:</b> averaged over the 642 km polar orbit the model gives 6621.0 mas/yr geodetic and "
        "40.85 mas/yr frame dragging, 39.1 mas/yr of it seen against GP-B's guide star; measured 6601.8 ± 18.3 and "
        "37.2 ± 7.2 mas/yr.",
        "<b>The material limits the rotor first:</b> burst margin 2.00 at every site; at the burst speed (v/c)² is "
        "4×10<super>−11</super>.",
        "<b>Near the illustrative neutron star</b> the hub clock runs at 0.997930583704828 (GM/rc² = 0.00207, a strong "
        "field) and frame dragging reaches 4.87×10<super>11</super> mas/yr; tidal stress 3.34×10<super>5</super> Pa still "
        "leaves a margin of 1.9×10<super>4</super>.",
        "<b>The machine's own spacetime is tiny:</b> 1 m away its clock shift is 2.1×10<super>−26</super> (3×10<super>−17</super> "
        "of Earth's) and its frame dragging 6.43×10<super>−8</super> mas/yr (3.7×10<super>−10</super> of Earth's at 45°); "
        "only its tidal pull (0.0012 of Earth's) is within reach of torsion balances.",
        "<b>No energy enters a circuit from spacetime curvature:</b> the gravitational channel carries clock rates and "
        "frames, not power.",
    ])

    r.h1("The construction")
    r.p("CLAUDE.md §3l keeps three \"curvatures\" apart; this example uses the first two.")
    r.bullets([
        "<b>Gravitational curvature.</b> The hub sits in the Schwarzschild field of the body at its site, moving with the "
        "site: clock rate √(1 − r<sub>s</sub>/r − v²/c²) (r<sub>s</sub> the Schwarzschild radius, r the distance, v the "
        "site's speed), exact form in strong fields, with tidal stretching and the horizon and photon-sphere limits as "
        "typed errors.",
        "<b>Rotation in flat spacetime.</b> In the hub's local inertial frame the rim is ordinary special relativity (the "
        "Langevin metric): rim clocks run at √(1 − ω²r²/c²), a loop on the rotor shows the Sagnac delay 4Aω/c² (A the "
        "loop area, ω the spin rate), and the rim speed is bounded by the material's burst speed √(σ/ρ).",
        "By the equivalence principle the two compose: the exact rim rate equals the first-order weak-field metric's rate "
        "to second order in the small quantities. The spin axis precesses by Thomas + geodetic + Lense–Thirring "
        "(Schiff 1960); the loop's Sagnac delay is taken relative to the local gyroscopes, with a \"local-frame "
        "correction\" for the local frame's precession.",
    ])
    r.h2("The rotor")
    r.p("Toray T1000G fibre, 6,370 MPa and 1.80 g/cm³, bursts as a thin hoop at √(σ/ρ) = 1,881 m/s (a solid disk "
        "bursts lower; the fibre value is an upper bound for a composite rotor). The example runs at half that, "
        "941 m/s at 0.5 m, so ω = 1,881 rad/s. The disk's mass is 1800·π·0.5²·0.02 = 28.3 kg and its angular momentum "
        "½MR²ω = 6.65×10<super>3</super> kg·m²/s.")
    r.h2("The three sites")
    r.table([
        ["site", "body", "position and motion"],
        ["Earth laboratory", "Earth (GM, J2-free, rotating)", "surface at 45° latitude, moving with the ground"],
        ["Gravity Probe B orbit", "Earth", "circular polar orbit at 642 km above the equatorial radius, inertially pointed"],
        ["neutron star (illustrative)", "1.4 solar masses, 12 km radius, 10 ms spin, I = 0.4·M·R² (assumed values, not "
         "a measured object)", "hovering 1,000 km from the centre, held up by a rocket"],
    ], [0.22, 0.38, 0.4])
    r.small("Solar GM is the IAU 2015 Resolution B3 nominal value.")

    r.h1("Results")
    r.h2("The rotor at each site")
    r.table([
        ["", "Earth lab, 45°", "GP-B orbit, 642 km", "neutron star, 1,000 km (illustrative)"],
        ["hub clock rate", "0.999999999304050", "0.999999999052363", "0.997930583704828"],
        ["rim clock rate", "0.999999999299128", "0.999999999047441", "0.997930583699917"],
        ["rim loses against hub", "4.25×10<super>−7</super> s/day", "4.25×10<super>−7</super> s/day",
         "4.25×10<super>−7</super> s/day"],
        ["GM/rc²", "6.95×10<super>−10</super>", "6.32×10<super>−10</super>", "0.00207"],
        ["v²/c² at the rim", "9.84×10<super>−12</super>", "9.84×10<super>−12</super>", "9.84×10<super>−12</super>"],
        ["burst margin", "2.00", "2.00", "2.00"],
        ["tidal stress on the disk (margin)", "0.00277 Pa (2.3×10<super>12</super>)", "0.00207 Pa (3.07×10<super>12</super>)",
         "3.34×10<super>5</super> Pa (1.9×10<super>4</super>)"],
        ["Thomas precession, here", "116 mas/yr", "0", "0"],
        ["geodetic precession, here", "350 mas/yr", "6.62×10<super>3</super> mas/yr", "0"],
        ["frame dragging, here", "172 mas/yr", "81.7 mas/yr", "4.87×10<super>11</super> mas/yr"],
        ["Sagnac on the 1 m² loop", "8.372×10<super>−14</super> s per pass", "8.372×10<super>−14</super> s",
         "8.372×10<super>−14</super> s"],
        ["local-frame correction to Sagnac", "−1.05×10<super>−30</super> s", "4.53×10<super>−29</super> s", "0 s"],
    ], [0.27, 0.22, 0.22, 0.29])
    r.bullets([
        "Averaged over the GP-B orbit the model gives 6621.0 mas/yr geodetic and 40.85 mas/yr frame dragging along "
        "Earth's axis, of which 39.1 mas/yr is seen against the mission's guide star IM Pegasi (declination +16°50′; "
        "GP-B measured the component perpendicular to the star line). Measured: 6601.8 ± 18.3 and 37.2 ± 7.2 mas/yr.",
        "At the neutron star the program warns that GM/rc² = 0.00207 is a strong field: the clock rates use the exact "
        "Schwarzschild form, and the precession and Sagnac terms are first-order and only indicative.",
    ])
    r.h2("The machine's own spacetime")
    r.p("1 m from the 28.3 kg rotor (J = 6.65×10<super>3</super> kg·m²/s), on its axis:")
    r.table([
        ["", "rotor", "Earth (surface, or at 45°)", "ratio"],
        ["clock shift GM/dc²", "2.1×10<super>−26</super>", "6.95×10<super>−10</super>", "3×10<super>−17</super>"],
        ["tidal curvature GM/c²d³", "2.1×10<super>−26</super> /m²", "1.71×10<super>−23</super> /m²", "0.0012"],
        ["frame dragging on the axis", "9.87×10<super>−24</super> rad/s = 6.43×10<super>−8</super> mas/yr",
         "172 mas/yr", "3.7×10<super>−10</super>"],
    ], [0.27, 0.33, 0.24, 0.16])
    r.p("A machine changes curvature through its mass and energy, and its spin drags frames, but at laboratory scale "
        "the tidal curvature is the only part within reach (torsion balances weigh that pull), and the frame dragging "
        "is ten orders below Earth's own. The rotor's rim speed is limited by its material long before any relativistic "
        "effect becomes large.")

    r.h1("Validation")
    r.p("The example's own test checks the self-field arithmetic (a 1 kg point mass at 1 m shifts clocks by "
        "G/c² = 7.4261×10<super>−28</super>, and a disk's J is ½MR²ω) and that the program runs. The physics is held in "
        "<i>spacetime/</i>:")
    r.table([
        ["test", "check", "held to"],
        ["TestGravityProbeB", "geodetic 6,606.1 and frame dragging 39.2 mas/yr (GR, 642 km polar orbit)",
         "0.5 % and 3 % (spherical Earth, no J2); within 2σ of the measured 6,601.8 ± 18.3 and 37.2 ± 7.2"],
        ["TestRimRateTwin", "exact rim rate = weak-field first-order rate", "second order in the small quantities"],
        ["TestSagnacEarthGyro", "1 m² loop on the ground at 45°: 4Aω sin φ/c²; local-frame correction "
         "10<super>−11</super> to 10<super>−8</super> of it", "1e-12"],
        ["TestGPSClock", "GPS clock net +38.6 µs/day against the rotating Earth (model 38.54 without J2); gravitational "
         "part 45.7", "Ashby 2003"],
        ["TestPoundRebka", "gh/c² over 22.5 m = 2.45×10<super>−15</super>", "Pound and Snider 1965"],
        ["TestShapiroDelay, TestSolarDeflection", "closed form against the Plebanski index path integral; solar-limb "
         "deflection 1.751″", "Shapiro 1964; Will 2014"],
        ["TestHorizonAndOrbitErrors, TestTidal", "r = 3r<sub>s</sub> static rate √(2/3); typed errors at the horizon "
         "and inside the photon sphere; Earth-surface tidal gradient 3.08×10<super>−6</super> s<super>−2</super>", "exact"],
        ["TestRotorShift, TestSagnac, TestBurstAndRotor", "rotor shift with k = 1/2 equals exact rim dilation (H17's k "
         "beside it); exact Sagnac equals 4Aω/c²; T1000G burst speed", "exact"],
    ], [0.27, 0.45, 0.28])

    r.h1("Limits")
    r.bullets([
        "Earth is spherical with no J2 term; the GP-B test tolerances (0.5 % geodetic, 3 % frame dragging) are stated with that model.",
        "The neutron star uses assumed values, not a measured object; in its strong field the precession and Sagnac "
        "terms are first-order and indicative only.",
        "The thin-hoop burst speed is an upper bound; a solid composite disk bursts lower.",
        "H17 (rotor Mössbauer anomaly, k ≈ 0.68 reported against the relativistic 1/2) is carried as a parameter in "
        "<i>spacetime.RotorShift</i>; this run does not print H17 values. Its register status is \"deviation reported, "
        "not independently reproduced\".",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/machine_spacetime\n"
           "go test ./examples/machine_spacetime ./spacetime")
    r.small("Sources: C. W. F. Everitt et al., PRL 106, 221101 (2011); N. Bartel et al., arXiv:1509.07529; L. I. Schiff, "
            "PRL 4, 215 (1960); N. Ashby, Living Rev. Relativ. 6, 1 (2003); Misner, Thorne and Wheeler, Gravitation "
            "(1973); J. Plebanski, Phys. Rev. 118, 1396 (1960); C. M. Will, Living Rev. Relativ. 17, 4 (2014); "
            "G. Sagnac (1913); E. J. Post, Rev. Mod. Phys. 39, 475 (1967); J. Bailey et al., Nature 268, 301 (1977); "
            "B. Mashhoon, Phys. Lett. A 145, 147 (1990); IAU 2015 Resolution B3; Toray T1000G data sheet (2024); "
            "S. Timoshenko, Strength of Materials II (1956).")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "machine_spacetime_report.pdf"))
