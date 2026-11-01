# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent references for the M9 MOSFET golden files.

As for the BJT files (bjt_ref.py, whose nodal framework this reuses), the
model is written out again from the published equations as terminal
currents and charges only; every Jacobian entry, conductance and junction
capacitance comes from central finite differences, and the solvers are
independent of the Go simulator (damped Newton; (G + jwC) x = b; fixed-step
trapezoidal rule).

Gate capacitances follow Meyer (RCA Review 1971) as SPICE applies them: the
gate charge advances by the step-averaged capacitance (plus overlap) times
the voltage step, and starts from the full capacitance times the voltage at
the operating point. Meyer's model is not charge-conserving, so this is
reproduced as a definition, not derived from a charge function.

Sources: Shichman & Hodges, IEEE JSSC 1968 (level 1); Vladimirescu & Liu,
UCB/ERL M80/7, 1980 (levels 1-3, process parameters); SPICE3f3 User's
Manual 4.4 (defaults).

Run from the repository root:  python3 internal/golden/ref/mos_ref.py [name ...]
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from bjt_ref import (Circuit, Element, R, Cap, Vsrc, V, addf, jac, newton, op,  # noqa: E402
                     write, write_op, frange, VT, GMIN)

EPS_OX = 3.9 * 8.854214871e-12
EPS_SI = 11.7 * 8.854214871e-12
Q = 1.602176634e-19
NI = 1.45e16
KB = 1.380649e-23
T = 300.15

DEF = dict(LEVEL=1, VTO=None, KP=None, GAMMA=None, PHI=None, LAMBDA=0.0, RD=0.0, RS=0.0, CBD=None, CBS=None,
           IS=1e-14, PB=0.8, CGSO=0.0, CGDO=0.0, CGBO=0.0, RSH=0.0, CJ=None, MJ=0.5, CJSW=0.0, MJSW=None,
           JS=0.0, TOX=None, NSUB=None, NSS=0.0, TPG=1.0, LD=0.0, UO=600.0, FC=0.5,
           NFS=0.0, XJ=0.0, VMAX=0.0, THETA=0.0, ETA=0.0, KAPPA=None, DELTA=0.0, WD=0.0,
           UCRIT=1e4, UEXP=0.0, UTRA=0.0, NEFF=1.0)


def depletion_q(v, cz, pb, m, fc):
    if cz == 0:
        return 0.0
    if v < fc * pb:
        return cz * pb * (1 - (1 - v / pb) ** (1 - m)) / (1 - m)
    v0 = fc * pb
    q0 = cz * pb * (1 - (1 - fc) ** (1 - m)) / (1 - m)
    k = cz / (1 - fc) ** (1 + m)
    return q0 + k * ((1 - fc * (1 + m)) * (v - v0) + m / (2 * pb) * (v * v - v0 * v0))


def meyer(vgs, vgd, von, vdsat, phi, cox):
    """Full Meyer capacitances (cgs, cgd, cgb), normal mode."""
    vgst = vgs - von
    vdsat = max(vdsat, 0.0)
    if vgst <= -phi:
        return 0.0, 0.0, cox
    if vgst <= -phi / 2:
        return 0.0, 0.0, -vgst * cox / phi
    if vgst <= 0:
        return 2 * (vgst * cox / (1.5 * phi) + cox / 3), 0.0, -vgst * cox / phi
    vds = vgs - vgd
    if vdsat <= vds:
        return 2 * cox / 3, 0.0, 0.0
    a = 2 * vdsat - vds
    return (2 * cox * (1 - ((vdsat - vds) / a) ** 2) / 3, 2 * cox * (1 - (vdsat / a) ** 2) / 3, 0.0)


class MOS(Element):
    def __init__(self, c, name, d, g, s, b, model, typ, L=100e-6, W=100e-6, AD=0.0, AS=0.0, PD=0.0, PS=0.0,
                 NRD=0.0, NRS=0.0):
        p = dict(DEF)
        p.update(model)
        self.name, self.t = name, (1 if typ == "nmos" else -1)
        self.level = int(p["LEVEL"])
        if p["MJSW"] is None:
            p["MJSW"] = 0.33 if self.level >= 2 else 0.5
        if p["KAPPA"] is None:
            p["KAPPA"] = 0.2
        if p["TOX"] is None and self.level >= 2:
            p["TOX"] = 1e-7
        self.d, self.g, self.s, self.b = c.node(d), c.node(g), c.node(s), c.node(b)
        rd = p["RD"] if p["RD"] > 0 else p["RSH"] * NRD
        rs = p["RS"] if p["RS"] > 0 else p["RSH"] * NRS
        self.rd, self.rs = rd, rs
        self.dp = c.internal(name + "d") if rd > 0 else self.d
        self.sp = c.internal(name + "s") if rs > 0 else self.s
        self.leff = L - 2 * p["LD"]
        self.w = W - 2 * p["WD"] if self.level == 3 else W
        cox = EPS_OX / p["TOX"] if p["TOX"] else 0.0
        self.cox = cox
        if p["KP"] is None:
            p["KP"] = p["UO"] * cox * 1e-4 if cox else 2e-5
        vt = VT
        if cox and p["NSUB"]:
            nsub = p["NSUB"] * 1e6
            ni = NI
            if self.level == 3:
                # ngspice mos3temp.c scales ni from 300 K to TNOM (level 3 only).
                tn = T
                eg = 1.16 - 7.02e-4 * tn * tn / (tn + 1108)
                ni *= (tn / 300) ** 1.5 * math.exp(0.5 * eg * (1 / 300 - 1 / tn) / (KB / Q))
            if nsub <= ni:
                raise ValueError("NSUB below the intrinsic density")
            if p["PHI"] is None:
                p["PHI"] = max(0.1, 2 * vt * math.log(nsub / ni))
            if p["GAMMA"] is None:
                p["GAMMA"] = math.sqrt(2 * EPS_SI * Q * nsub) / cox
            if p["VTO"] is None:
                eg = 1.16 - 7.02e-4 * T * T / (T + 1108)
                fermis = self.t * 0.5 * p["PHI"]
                wkfng = 3.2
                if p["TPG"] != 0:
                    wkfng = 3.25 + 0.5 * eg - self.t * p["TPG"] * 0.5 * eg
                vfb = wkfng - (3.25 + 0.5 * eg + fermis) - p["NSS"] * 1e4 * Q / cox
                p["VTO"] = vfb + self.t * (p["GAMMA"] * math.sqrt(p["PHI"]) + p["PHI"])
            if p["CJ"] is None:
                p["CJ"] = math.sqrt(EPS_SI * Q * nsub / (2 * p["PB"]))
        for k, v in (("PHI", 0.6), ("GAMMA", 0.0), ("VTO", 0.0), ("CJ", 0.0)):
            if p[k] is None:
                p[k] = v
        self.p = p
        self.beta = p["KP"] * self.w / self.leff
        self.ocap = cox * self.leff * self.w
        # JS x area only when JS, AD and AS are all non-zero (SPICE)
        both = p["JS"] > 0 and AD > 0 and AS > 0
        self.isd = p["JS"] * AD if both else p["IS"]
        self.iss = p["JS"] * AS if both else p["IS"]
        self.czbd = p["CBD"] if p["CBD"] is not None else p["CJ"] * AD
        self.czbs = p["CBS"] if p["CBS"] is not None else p["CJ"] * AS
        self.czbdsw, self.czbssw = p["CJSW"] * PD, p["CJSW"] * PS
        self.ov = (p["CGSO"] * self.w, p["CGDO"] * self.w, p["CGBO"] * self.leff)
        self.vto = self.t * p["VTO"]
        self.vbi = self.vto - p["GAMMA"] * math.sqrt(p["PHI"])
        self.nsub_given = bool(cox and p["NSUB"])
        # depletion-width coefficient sqrt(2 eps_Si / (q NSUB))
        self.xd = math.sqrt(2 * EPS_SI / (Q * p["NSUB"] * 1e6)) if self.nsub_given else 0.0

    # ---- drain current, normal mode (vds >= 0) ----
    def ids(self, vgs, vds, vbs):
        """Returns (id, von, vdsat)."""
        if self.level == 2:
            return self.ids2(vgs, vds, vbs)
        if self.level == 3:
            return self.ids3(vgs, vds, vbs)
        p = self.p
        if vbs <= 0:
            sarg = math.sqrt(p["PHI"] - vbs)
        else:
            sarg = math.sqrt(p["PHI"])
            sarg = max(0.0, sarg - vbs / (2 * sarg))
        von = self.vbi + p["GAMMA"] * sarg
        vgst = vgs - von
        vdsat = max(vgst, 0.0)
        if vgst <= 0:
            return 0.0, von, vdsat
        betap = self.beta * (1 + p["LAMBDA"] * vds)
        if vgst <= vds:
            return betap * vgst * vgst / 2, von, vdsat
        return betap * vds * (vgst - vds / 2), von, vdsat

    def ids2(self, vgs, vds, vbs):
        """Level 2, Grove-Frohman (Vladimirescu & Liu 1980, 2.3; SPICE mos2 regions)."""
        p, L, W, cox, vt = self.p, self.leff, self.w, self.cox, VT
        phi = p["PHI"]
        sphi = math.sqrt(phi)
        gamma = p["GAMMA"]
        pv = phi - vbs  # phi - vbs

        def root(u, vb):
            # sqrt(u) with the SPICE continuation sqrt(phi)/(1 + vb/(2 phi)) for vb > 0
            return math.sqrt(u) if vb <= 0 else sphi / (1 + 0.5 * vb / phi)

        def droot(u, vb):
            # d root / d vbs
            return -0.5 / math.sqrt(u) if vb <= 0 else -0.5 * root(u, vb) ** 2 / (phi * sphi)

        s_ = root(pv, vbs)
        b_ = root(pv + vds, vbs - vds)
        fac = p["DELTA"] * math.pi * EPS_SI / (4 * cox * W)
        eta = 1 + fac
        vbin = self.vbi + fac * pv
        gs, dgs = gamma, 0.0  # gamma with junction-depth correction, its vbs slope
        if (gamma > 0 or self.nsub_given) and p["XJ"] > 0:
            xj, xd = p["XJ"], self.xd

            def gam(vb):
                ss = root(phi - vb, vb)
                bb = root(phi - vb + vds, vb - vds)
                a_s = 0.5 * xj / L * (math.sqrt(1 + 2 * xd * ss / xj) - 1)
                a_d = 0.5 * xj / L * (math.sqrt(1 + 2 * xd * bb / xj) - 1)
                return gamma * (1 - a_s - a_d)
            gs = gam(vbs)
            # analytic slope of gamma_sd in vbs (enters the weak-inversion factor xn)
            ss, bb = s_, b_
            dgs = -gamma * (0.5 / L * xd * droot(pv, vbs) / math.sqrt(1 + 2 * xd * ss / xj)
                            + 0.5 / L * xd * droot(pv + vds, vbs - vds) / math.sqrt(1 + 2 * xd * bb / xj))
        von = vbin + gs * s_
        nfs = p["NFS"] != 0
        if nfs:
            xn = 1 + Q * p["NFS"] * 1e4 / cox - (gs * droot(pv, vbs) + dgs * s_) + fac
            von += vt * xn
        elif vgs <= vbin:
            return 0.0, von, 0.0
        vgst = vgs - von
        crit = p["UCRIT"] * 100 * EPS_SI / cox
        uf = (crit / vgst) ** p["UEXP"] if vgst > crit else 1.0
        ueff = p["UO"] * 1e-4 * uf
        vgx = max(vgs, von) if nfs else vgs
        gd = gs / eta
        u = (vgx - vbin) / eta
        vdsat = 0.0
        if gd > 0:
            if u + pv > 0:
                vdsat = max(u + gd * gd * (1 - math.sqrt(1 + 4 * (u + pv) / (gd * gd))) / 2, 0.0)
        else:
            vdsat = max(u, 0.0)
        s3 = s_ ** 3
        if p["VMAX"] > 0:
            xv = p["VMAX"] * L / ueff
            # x = sqrt(vdsat + phi - vbs) solves x^4 + a1 x^3 + b1 x^2 + c1 x + d1 = 0 (Baum-Beneking)
            v1, v2 = u + pv, pv
            co = [1.0, gd / 0.75, -2 * (v1 + xv), -2 * gd * xv, 2 * v1 * (v2 + xv) - v2 * v2 - 4 / 3 * gd * s3]
            good = [r.real for r in np.roots(co) if abs(r.imag) < 1e-9 * max(1, abs(r)) and r.real > 0]
            good = [r for r in good if abs(np.polyval(co, r)) <= 1e-6]
            if good:
                x = min(good)
                vdsat = x * x - pv
        xlam = p["LAMBDA"]
        if vds != 0:
            bs = root(vdsat + pv, vbs - vdsat)
            bodys = bs ** 3 - s3
            if self.nsub_given and not xlam > 0:
                if p["VMAX"] <= 0:
                    a = (vds - vdsat) / 4
                    xlam = self.xd / (L * vds) * math.sqrt(a + math.sqrt(1 + a * a))
                else:
                    xdv = self.xd / math.sqrt(p["NEFF"])
                    xlv = p["VMAX"] * xdv / (2 * ueff)
                    xlam = xdv / (L * vds) * (math.sqrt(xlv * xlv + max(vds - vdsat, 0)) - xlv)
        xwb = self.xd * math.sqrt(p["PB"])
        xld = L - xwb
        clf = 1 - xlam * vds
        if not self.nsub_given:
            xwb = 0.25e-6
        if L * clf < xwb:  # punch-through
            clf = xwb / (1 + (xlam * vds * L - xld) / xwb) / L
        b1 = self.beta * uf / clf
        if vds <= 1e-10:
            return 0.0, von, vdsat
        body = b_ ** 3 - s3
        if nfs and vgs <= von:
            if vdsat <= 0:
                return 0.0, von, vdsat
            vo = min(vdsat, vds)
            bd = bodys if vds > vdsat else body
            cd = b1 * ((von - vbin - eta * vo / 2) * vo - gs * bd / 1.5)
            return cd * math.exp((vgs - von) / (vt * xn)), von, vdsat
        if vds <= vdsat:
            return b1 * ((vgs - vbin - eta * vds / 2) * vds - gs * body / 1.5), von, vdsat
        return b1 * ((vgs - vbin - eta * vdsat / 2) * vdsat - gs * bodys / 1.5), von, vdsat

    def ids3(self, vgs, vds, vbs):
        """Level 3, semi-empirical (Vladimirescu & Liu 1980, 2.4; ngspice mos3 regions)."""
        p, L, W, cox, vt = self.p, self.leff, self.w, self.cox, VT
        phi = p["PHI"]
        if vbs <= 0:
            phibs = phi - vbs
            sq = math.sqrt(phibs)
        else:
            sq = math.sqrt(phi) / (1 + vbs / (2 * phi))
            phibs = sq * sq
        fsh = 1.0
        if p["XJ"] != 0 and self.xd != 0:
            xj = p["XJ"]
            wp = self.xd * sq / xj
            wc = 0.0631353 + 0.8013292 * wp - 0.01110777 * wp * wp
            ac = wp / (1 + wp)
            fsh = 1 - xj / L * ((wc + p["LD"] / xj) * math.sqrt(1 - ac * ac) - p["LD"] / xj)
        gms = p["GAMMA"] * fsh
        nar = p["DELTA"] * 0.5 * math.pi * EPS_SI / cox / W
        fb = gms / (4 * sq) + nar
        qb = gms * sq + nar * phibs
        eta = p["ETA"] * 8.15e-22 / (cox * L ** 3)
        vth = self.vbi - eta * vds + qb
        von = vth
        nfs = p["NFS"] != 0
        if nfs:
            xn = 1 + Q * p["NFS"] * 1e4 / cox + qb / (2 * phibs)
            von = vth + vt * xn
        elif vgs <= von:
            return 0.0, von, 0.0
        vgx = max(vgs, von)
        fg = 1 / (1 + p["THETA"] * (vgx - vth))
        us = p["UO"] * 1e-4 * fg
        vdsat = (vgx - vth) / (1 + fb)
        vmax = p["VMAX"]
        if vmax > 0:
            vc = L * vmax / us
            vdsat = vdsat + vc - math.sqrt(vdsat * vdsat + vc * vc)
        vx = min(vds, vdsat)
        if vx == 0:
            return 0.0, von, vdsat
        cd = self.beta * fg * (vgx - vth - 0.5 * (1 + fb) * vx) * vx
        fd = 1.0
        if vmax > 0:
            fd = 1 / (1 + vx / vc)
            cd *= fd
        alpha, kap = self.xd ** 2, p["KAPPA"]
        dl = None
        if vds <= vdsat:
            if vmax <= 0 and alpha != 0:
                dl = math.sqrt(kap * alpha * vdsat / 8) * (vds / vdsat) ** 4
        elif vmax <= 0:
            dl = math.sqrt(kap * alpha * (vds - vdsat + vdsat / 8))
        elif alpha != 0:
            gdsat = max(cd * (1 - fd) / vc, 1e-12)
            em = kap * cd / (L * gdsat)
            a = 0.5 * em * alpha
            dl = math.sqrt(a * a + kap * alpha * (vds - vdsat)) - a
        if dl is not None:
            if dl > 0.5 * L:
                dl = L - L * L / (4 * dl)
            cd /= 1 - dl / L
        if nfs and vgs < von:
            cd *= math.exp((vgs - von) / (vt * xn))
        return cd, von, vdsat

    def channel(self, vgs, vds, vbs):
        """Channel current drain' -> source' (normalized), von, vdsat, mode."""
        if vds >= 0:
            i, von, vdsat = self.ids(vgs, vds, vbs)
            return i, von, vdsat, 1
        i, von, vdsat = self.ids(vgs - vds, -vds, vbs - vds)
        return -i, von, vdsat, -1

    def volts(self, x):
        t = self.t
        vg, vb, vd, vs = V(x, self.g), V(x, self.b), V(x, self.dp), V(x, self.sp)
        return t * (vg - vs), t * (vd - vs), t * (vb - vs)

    def diode(self, v, i_s):
        if v <= -3 * VT:
            return GMIN * v - i_s
        return i_s * (math.exp(min(v / VT, 700)) - 1) + GMIN * v

    def static(self, x, t_, f, s):
        t = self.t
        vgs, vds, vbs = self.volts(x)
        vbd = vbs - vds
        ich, _, _, _ = self.channel(vgs, vds, vbs)
        ibs, ibd = self.diode(vbs, self.iss), self.diode(vbd, self.isd)
        addf(f, self.dp, t * (ich - ibd))
        addf(f, self.sp, t * (-ich - ibs))
        addf(f, self.b, t * (ibs + ibd))
        for a, b2, r in ((self.d, self.dp, self.rd), (self.s, self.sp, self.rs)):
            if r > 0:
                i = (V(x, a) - V(x, b2)) / r
                addf(f, a, i)
                addf(f, b2, -i)

    def charge(self, x, qv):
        """Junction charges only (conservative)."""
        t = self.t
        _, vds, vbs = self.volts(x)
        p = self.p
        for (a, v, cz, czsw) in ((self.sp, vbs, self.czbs, self.czbssw), (self.dp, vbs - vds, self.czbd, self.czbdsw)):
            qq = depletion_q(v, cz, p["PB"], p["MJ"], p["FC"]) + depletion_q(v, czsw, p["PB"], p["MJSW"], p["FC"])
            addf(qv, self.b, t * qq)
            addf(qv, a, -t * qq)

    def meyer_caps(self, x):
        """Full gate capacitances (gs, gd, gb) in the terminal frame, without overlap."""
        vgs, vds, vbs = self.volts(x)
        _, von, vdsat, mode = self.channel(vgs, vds, vbs)
        if self.ocap == 0:
            return (0.0, 0.0, 0.0)
        vgd = vgs - vds
        if mode > 0:
            return meyer(vgs, vgd, von, vdsat, self.p["PHI"], self.ocap)
        a, b, cgb = meyer(vgd, vgs, von, vdsat, self.p["PHI"], self.ocap)
        return (b, a, cgb)

    def gate_v(self, x):
        vgs, vds, vbs = self.volts(x)
        return (vgs, vgs - vds, vgs - vbs)

    def gate_ends(self):
        return ((self.g, self.sp), (self.g, self.dp), (self.g, self.b))


def mos_ac(c, x, freqs):
    """(G + jwC) with junction C by finite differences and Meyer + overlap C at x."""
    G = jac(lambda y: c.residual(y, None), x, 1e-6)
    C = jac(lambda y: c.charge(y), x, 1e-6)
    for e in c.elems:
        if isinstance(e, MOS):
            mc = e.meyer_caps(x)
            for (a, b), cm, ov in zip(e.gate_ends(), mc, e.ov):
                cc = cm + ov
                for i, j, sgn in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
                    if i >= 0 and j >= 0:
                        C[i, j] += sgn * cc
    b = c.ac_rhs()
    return [np.linalg.solve(G + 2j * math.pi * f * C, b) for f in freqs]


def mos_tran(c, stop, h, record_every):
    """Trapezoidal rule at a fixed step; Meyer gate charges advance by the
    step-averaged capacitance times the voltage step (SPICE convention)."""
    x = op(c, t=0.0)
    mos = [e for e in c.elems if isinstance(e, MOS)]
    n = c.size()

    def gate_q_vector(qg):
        qv = np.zeros(n)
        for e, qs in zip(mos, qg):
            for (a, b), qq in zip(e.gate_ends(), qs):
                addf(qv, a, e.t * qq)
                addf(qv, b, -e.t * qq)
        return qv

    # Initial gate charges: (C + overlap) * v at the operating point.
    caps0 = [e.meyer_caps(x) for e in mos]
    vg0 = [e.gate_v(x) for e in mos]
    qg = [tuple((cm + ov) * v for cm, ov, v in zip(cs, e.ov, vs)) for e, cs, vs in zip(mos, caps0, vg0)]
    qtot = c.charge(x) + gate_q_vector(qg)
    iq = np.zeros(n)
    rows = [(0.0, x.copy())]
    caps_prev, v_prev = caps0, vg0
    nsteps = int(round(stop / h))
    for k in range(1, nsteps + 1):
        t = k * h

        def new_gate_q(y):
            out = []
            for e, qs, cp, vp in zip(mos, qg, caps_prev, v_prev):
                cn, vn = e.meyer_caps(y), e.gate_v(y)
                out.append(tuple(qo + ((a + b2) / 2 + ov) * (v - v0)
                                 for qo, a, b2, ov, v, v0 in zip(qs, cn, cp, e.ov, vn, vp)))
            return out

        def fun(y):
            qn = c.charge(y) + gate_q_vector(new_gate_q(y))
            return c.residual(y, t) + 2 * (qn - qtot) / h - iq

        x = newton(fun, x.copy(), vmax=0.5)
        qg = new_gate_q(x)
        qnew = c.charge(x) + gate_q_vector(qg)
        iq = 2 * (qnew - qtot) / h - iq
        qtot = qnew
        caps_prev, v_prev = [e.meyer_caps(x) for e in mos], [e.gate_v(x) for e in mos]
        if k % record_every == 0:
            rows.append((t, x.copy()))
    return rows


def drain_current(e, x):
    """External drain current into the device (DC)."""
    if e.rd > 0:
        return (V(x, e.d) - V(x, e.dp)) / e.rd
    f = np.zeros(len(x))
    e.static(x, None, f, 1.0)
    return f[e.d] if e.d >= 0 else 0.0


# --------------------------------------------------------------- circuits --

NPROC = dict(TOX=20e-9, NSUB=3e16, UO=550, LD=0.1e-6, LAMBDA=0.03, RSH=20, CGSO=2e-10, CGDO=2e-10, CGBO=1e-10,
             CJ=4e-4, MJ=0.45, CJSW=3e-10, MJSW=0.3, PB=0.8, JS=1e-4)
PPROC = dict(TOX=20e-9, NSUB=1e16, UO=200, LD=0.1e-6, LAMBDA=0.05, RSH=40, CGSO=2e-10, CGDO=2e-10, CGBO=1e-10,
             CJ=5e-4, MJ=0.5, CJSW=3.5e-10, MJSW=0.33, PB=0.85, JS=1e-4)


SUFFIX = (("meg", 1e6), ("t", 1e12), ("g", 1e9), ("k", 1e3), ("m", 1e-3), ("u", 1e-6), ("n", 1e-9),
          ("p", 1e-12), ("f", 1e-15))


def num(tok):
    t = tok.lower()
    for sfx, mul in SUFFIX:
        if t.endswith(sfx) and not t.endswith("e" + sfx):
            return float(t[:-len(sfx)]) * mul
    return float(t)


def netlist_model(fname, mname):
    """The .model card of a test netlist as a parameter dict (input data, not simulator code)."""
    here = os.path.dirname(__file__)
    path = os.path.join(here, "..", "..", "..", "testdata", "netlists", fname)
    for line in open(path):
        w = line.split(None, 2)
        if len(w) == 3 and w[0].lower() == ".model" and w[1].lower() == mname:
            body = w[2][w[2].index("(") + 1:w[2].rindex(")")]
            return {k.upper(): num(v) for k, v in (kv.split("=") for kv in body.split())}
    raise KeyError(mname)


def add_pair(c, nm, pm, nin, nout, vdd, wn=10e-6, wp=25e-6, l=1e-6):
    # AD = AS = W x 5 um, PD = PS = 2 W + 10 um (as in the netlists)
    c.add(MOS(c, "MN", nout, nin, "0", "0", nm, "nmos", L=l, W=wn, AD=wn * 5e-6, AS=wn * 5e-6, PD=2 * wn + 10e-6,
              PS=2 * wn + 10e-6, NRD=0.5, NRS=0.5))
    c.add(MOS(c, "MP", nout, nin, vdd, vdd, pm, "pmos", L=l, W=wp, AD=wp * 5e-6, AS=wp * 5e-6, PD=2 * wp + 10e-6,
              PS=2 * wp + 10e-6, NRD=0.5, NRS=0.5))


def mos1_op():
    """Body-biased NMOS source follower and a PMOS current source, from process parameters."""
    c = Circuit()
    vdd = c.add(Vsrc(c, "VDD", "vdd", "0", 5))
    c.add(Vsrc(c, "VG", "g", "0", 3.2))
    c.add(Vsrc(c, "VB", "vb", "0", -1))
    m1 = c.add(MOS(c, "M1", "vdd", "g", "s", "vb", NPROC, "nmos", L=2e-6, W=20e-6, AD=1e-10, AS=1e-10, PD=30e-6,
                   PS=30e-6, NRD=1, NRS=1))
    c.add(R(c, "RS1", "s", "0", 10e3))
    c.add(Vsrc(c, "VGP", "gp", "0", 3.0))
    m2 = c.add(MOS(c, "M2", "dp", "gp", "vdd", "vdd", PPROC, "pmos", L=3e-6, W=30e-6, NRD=1, NRS=1))
    c.add(R(c, "RL", "dp", "0", 5e3))
    c.finish()
    x = op(c)
    vals = [("v(%s)" % n, V(x, c.nodes[n])) for n in ("vdd", "g", "s", "vb", "gp", "dp")]
    vals += [("i(vdd)", vdd.current(x)), ("i(m1)", drain_current(m1, x)), ("i(m2)", drain_current(m2, x))]
    write_op("mos1_op.csv", vals, "damped Newton on the level-1 equations, finite-difference Jacobian")


def mos1_idvd():
    """Output characteristic through source/drain interchange (negative VDS)."""
    m = dict(VTO=0.7, KP=60e-6, GAMMA=0.5, PHI=0.7, LAMBDA=0.04, RD=20, RS=20, IS=1e-14)
    rows, x = [], None
    for vds in frange(-2, 5, 0.1):
        c = Circuit()
        vd = c.add(Vsrc(c, "VDS", "d", "0", vds))
        c.add(Vsrc(c, "VGS", "g", "0", 2.5))
        c.add(Vsrc(c, "VBS", "b", "0", -0.5))
        c.add(MOS(c, "M1", "d", "g", "0", "b", m, "nmos", L=2e-6, W=10e-6))
        c.finish()
        x = op(c, x0=x) if x is not None else op(c)
        rows.append((vds, vd.current(x)))
    write("mos1_idvd_dc.csv", ["sweep", "i(vds)"], rows, "damped Newton on the level-1 equations at each point")


def idvd(name, fname, mname, L, W, note):
    """Output characteristic through source/drain interchange (VGS 2.5 V, VBS -0.5 V)."""
    m = netlist_model(fname, mname)
    rows, x = [], None
    for vds in frange(-2, 5, 0.1):
        c = Circuit()
        vd = c.add(Vsrc(c, "VDS", "d", "0", vds))
        c.add(Vsrc(c, "VGS", "g", "0", 2.5))
        c.add(Vsrc(c, "VBS", "b", "0", -0.5))
        c.add(MOS(c, "M1", "d", "g", "0", "b", m, "nmos", L=L, W=W, NRD=1, NRS=1))
        c.finish()
        x = op(c, x0=x) if x is not None else op(c)
        rows.append((vds, vd.current(x)))
    write(name, ["sweep", "i(vds)"], rows, note)


def subvt(name, fname, mname, L, W, vds, note):
    """Transfer characteristic VGS 0..2 V at fixed VDS, VBS -1 V."""
    m = netlist_model(fname, mname)
    rows, x = [], None
    for vgs in frange(0, 2, 0.025):
        c = Circuit()
        vd = c.add(Vsrc(c, "VDS", "d", "0", vds))
        c.add(Vsrc(c, "VGS", "g", "0", vgs))
        c.add(Vsrc(c, "VBS", "b", "0", -1))
        c.add(MOS(c, "M1", "d", "g", "0", "b", m, "nmos", L=L, W=W))
        c.finish()
        x = op(c, x0=x) if x is not None else op(c)
        rows.append((vgs, vd.current(x)))
    write(name, ["sweep", "i(vds)"], rows, note)


def cmos_vtc(name, level_n, level_p, note, l=1e-6):
    rows, x = [], None
    for vin in frange(0, 5, 0.05):
        c = Circuit()
        c.add(Vsrc(c, "VDD", "vdd", "0", 5))
        c.add(Vsrc(c, "VIN", "in", "0", vin))
        add_pair(c, level_n, level_p, "in", "out", "vdd", l=l)
        c.add(R(c, "RL", "out", "0", 1e6))
        c.finish()
        x = op(c, x0=x) if x is not None else op(c)
        rows.append((vin, V(x, c.nodes["out"])))
    write(name, ["sweep", "v(out)"], rows, note)


def mos1_cs_ac(name="mos1_cs_ac.csv", model=NPROC, vin=1.4):
    """Common-source stage with Meyer, overlap and junction capacitances."""
    c = Circuit()
    c.add(Vsrc(c, "VDD", "vdd", "0", 5))
    c.add(Vsrc(c, "VIN", "in", "0", vin, acmag=1))
    c.add(R(c, "RG", "in", "g", 1e3))
    c.add(R(c, "RD", "vdd", "d", 10e3))
    c.add(Cap(c, "CL", "d", "0", 0.2e-12))
    c.add(MOS(c, "M1", "d", "g", "0", "0", model, "nmos", L=1e-6, W=20e-6, AD=1e-10, AS=1e-10, PD=30e-6,
              PS=30e-6, NRD=0.5, NRS=0.5))
    c.finish()
    x = op(c)
    freqs = [10 ** (3 + k / 10) for k in range(0, 71)]  # .ac dec 10 1k 10g
    xs = mos_ac(c, x, freqs)
    n = c.nodes["d"]
    rows = [(f, abs(v[n]), math.degrees(math.atan2(v[n].imag, v[n].real))) for f, v in zip(freqs, xs)]
    write(name, ["freq", "vm(d)", "vp(d)"], rows,
          "(G + jwC) at the operating point: G and junction C by finite differences, Meyer C from its formula")


def mos1_inv_tran(name="mos1_inv_tran.csv", nm=NPROC, pm=PPROC):
    c = Circuit()
    c.add(Vsrc(c, "VDD", "vdd", "0", 5))
    c.add(Vsrc(c, "VIN", "in", "0", 0, pulse_args=(0, 5, 2e-9, 1e-9, 1e-9, 10e-9, 25e-9)))
    add_pair(c, nm, pm, "in", "out", "vdd")
    c.add(Cap(c, "CL", "out", "0", 50e-15))
    c.finish()
    rows = mos_tran(c, 30e-9, 1e-12, 100)  # record every 0.1 ns
    n = c.nodes["out"]
    write(name, ["time", "v(out)"], [(t, x[n]) for t, x in rows],
          "trapezoidal rule at a fixed 1 ps step; Meyer charge stepped as in SPICE")


if __name__ == "__main__":
    only = set(sys.argv[1:])
    jobs = {
        "mos1_op": mos1_op,
        "mos1_idvd": mos1_idvd,
        "mos1_vtc": lambda: cmos_vtc("mos1_cmos_vtc.csv", NPROC, PPROC, "damped Newton on the level-1 equations"),
        "mos1_cs_ac": mos1_cs_ac,
        "mos1_inv_tran": mos1_inv_tran,
        "mos2_idvd": lambda: idvd("mos2_idvd_dc.csv", "mos2_idvd_dc.cir", "n2", 2e-6, 10e-6,
                                  "damped Newton on the level-2 equations at each point (quartic vdsat by numpy.roots)"),
        "mos2_subvt": lambda: subvt("mos2_subvt_dc.csv", "mos2_subvt_dc.cir", "n2", 2e-6, 10e-6, 0.1,
                                    "damped Newton on the level-2 equations at each point"),
        "mos2_vtc": lambda: cmos_vtc("mos2_cmos_vtc.csv", netlist_model("mos2_cmos_vtc.cir", "n2"),
                                     netlist_model("mos2_cmos_vtc.cir", "p2"),
                                     "damped Newton on the level-2 equations", l=2e-6),
        "mos3_idvd": lambda: idvd("mos3_idvd_dc.csv", "mos3_idvd_dc.cir", "n3", 1.2e-6, 10e-6,
                                  "damped Newton on the level-3 equations at each point"),
        "mos3_subvt": lambda: subvt("mos3_subvt_dc.csv", "mos3_subvt_dc.cir", "n3", 1.2e-6, 10e-6, 3.0,
                                    "damped Newton on the level-3 equations at each point"),
        "mos3_vtc": lambda: cmos_vtc("mos3_cmos_vtc.csv", netlist_model("mos3_cmos_vtc.cir", "n3"),
                                     netlist_model("mos3_cmos_vtc.cir", "p3"), "damped Newton on the level-3 equations"),
        "mos3_cs_ac": lambda: mos1_cs_ac("mos3_cs_ac.csv", netlist_model("mos3_cs_ac.cir", "n3"), 1.2),
        "mos3_inv_tran": lambda: mos1_inv_tran("mos3_inv_tran.csv", netlist_model("mos3_inv_tran.cir", "n3"),
                                               netlist_model("mos3_inv_tran.cir", "p3")),
    }
    for name, fn in jobs.items():
        if not only or name in only:
            fn()
