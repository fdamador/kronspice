"""Graphs, diagrams and flux animations for examples/bearden_meg (KronSpice runs).

Run from the repository root after `go run ./examples/bearden_meg` (which writes
results.csv):

    go build -o kronspice ./cmd/kronspice
    KRON=./kronspice python3 examples/bearden_meg/viz/make_viz.py

Needs numpy, matplotlib and ffmpeg (for the .mp4 files).
"""
import csv, os, subprocess, tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib import animation, cm, colors

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.dirname(HERE)
OUT = HERE
RES = os.path.join(EX, 'results.csv')
KRON = os.environ.get('KRON', 'kronspice')
TMP = tempfile.mkdtemp()

mu0 = 4e-7 * np.pi
H = np.array([-1e5, -2000, -1000, -400, -200, -100, 0, 100, 200, 400, 1000, 2000, 1e5])
B = np.array([-1.6 - mu0*98000, -1.6, -1.5, -1.3, -1.0, -0.6, 0, 0.6, 1.0, 1.3, 1.5, 1.6, 1.6 + mu0*98000])
LEG, SEG = 0.05, 0.15 / 4          # side legs / return yokes; knife bar segments (gap 15 cm)


def bf(mmf, length=LEG):
    return np.interp(mmf / length, H, B)


def load(name):
    """Run netlists/NAME.cir through kronspice and read its transient CSV."""
    csvp = os.path.join(TMP, name + '.csv')
    subprocess.run([KRON, '-o', csvp, os.path.join(EX, 'netlists', name + '.cir')],
                   check=True, stdout=subprocess.DEVNULL)
    return np.genfromtxt(csvp, delimiter=',', names=True, deletechars='')


plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': 0.25, 'figure.dpi': 130})
BLUE, ORANGE, GREY, GREEN, RED, PURPLE = '#2c6fbb', '#d9822b', '#888888', '#3a9a5b', '#c0392b', '#6c3483'

rows = list(csv.DictReader(open(RES)))
for r in rows:
    for k in ('V', 'f_Hz', 'RL_ohm', 'mismatch', 'Cres_F', 'P_drive_W', 'P_magnets_W', 'P_load_W', 'COP'):
        r[k] = float(r[k])


def find(sweep, **kw):
    for r in rows:
        if r['sweep'] == sweep and all(abs(r[k] - v) < 1e-12 * max(1, abs(v)) for k, v in kw.items()):
            return r
    raise KeyError(sweep, kw)


kb = find('knife: drive 0.5-100 V', V=2)
mb = find('meg: drive 1-20 V', V=5)


def window(d, t0, t1):
    m = (d['time'] >= t0) & (d['time'] <= t1)
    return m, (d['time'][m] - t0) * 1e6


# ---------------------------------------------------------------- MEG
meg = load('meg_baseline')
m, tus = window(meg, 2.0e-3, 2.3e-3)
megBL, megBR = bf(meg['v(l2)']), bf(meg['v(r2)'])
fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
ax[0].plot(tus, meg['v(dl)'][m], color=BLUE, label='left control pulse')
ax[0].plot(tus, meg['v(dr)'][m], color=ORANGE, label='right control pulse')
ax[0].set_ylabel('drive [V]'); ax[0].legend(loc='upper right', fontsize=8)
ax[1].plot(tus, meg['i(rsl)'][m], color=BLUE, label='left control')
ax[1].plot(tus, meg['i(rsr)'][m], color=ORANGE, label='right control')
ax[1].set_ylabel('control current [A]'); ax[1].legend(loc='upper right', fontsize=8)
ax[2].plot(tus, megBL[m], color=BLUE, label='left leg')
ax[2].plot(tus, megBR[m], color=ORANGE, label='right leg')
ax[2].plot(tus, megBL[m] + megBR[m], color=GREY, ls='--', label='magnet leg (sum)')
ax[2].set_ylabel('flux density [T]'); ax[2].legend(loc='center right', fontsize=8)
ax[3].plot(tus, meg['v(ol)'][m], color=BLUE, label='left output')
ax[3].plot(tus, meg['v(or)'][m], color=ORANGE, label='right output')
ax[3].set_ylabel('output [V]'); ax[3].set_xlabel('time [µs] (three switching cycles in steady state)')
ax[3].legend(loc='upper right', fontsize=8)
fig.suptitle(f"MEG-style core, 5 V pulses at 10 kHz, 100 Ω loads: {mb['P_drive_W']:.3g} W in, "
             f"{mb['P_load_W']:.3g} W out (COP {mb['COP']:.2f})", fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/meg_waveforms.png'); plt.close(fig)

# ---------------------------------------------------------------- Knife edge (bucking magnets)
kn = load('knife_baseline')
m2, tus2 = window(kn, 11.0e-3, 11.3e-3)
kA = bf(kn['v(al)'] - kn['v(jl)'], SEG)        # + = toward the null (rightward)
kB = bf(kn['v(br)'] - kn['v(jr)'], SEG)        # + = toward the null (leftward)
kBr = bf(kn['v(y1)'] - kn['v(n)'], SEG)        # bridge at the null, + = rightward
kRL = bf(kn['v(jl)'], LEG)                     # left return yoke
kRR = bf(kn['v(jr)'], LEG)                     # right return yoke
fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
ax[0].plot(tus2, kn['v(ds)'][m2], color=BLUE); ax[0].set_ylabel('drive [V]')
ax[1].plot(tus2, -kn['i(vd)'][m2] * 1e3, color=BLUE); ax[1].set_ylabel('drive current [mA]\n(both coils, series)')
ax[2].plot(tus2, kA[m2], color=BLUE, label='bar at magnet A')
ax[2].plot(tus2, kB[m2], color=ORANGE, ls='--', label='bar at magnet B')
ax[2].plot(tus2, kRL[m2], color=PURPLE, label='return yokes')
ax[2].plot(tus2, kBr[m2], color=GREEN, lw=2.5, label='null (between the drive coils)')
ax[2].set_ylabel('flux density [T]'); ax[2].legend(loc='center right', fontsize=8)
ax[3].plot(tus2, kn['v(o1)'][m2], color=BLUE); ax[3].set_ylabel('output [V]')
ax[3].set_xlabel('time [µs] (three drive cycles in steady state)')
fig.suptitle(f"Bucking magnets, split drive coil at the null, 2 V at 10 kHz, 100 Ω: "
             f"{kb['P_drive_W']:.3g} W in, {kb['P_load_W']:.3g} W out (COP {kb['COP']:.2f})", fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/knife_waveforms.png'); plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 3.4))
ax.plot(tus2, kBr[m2] * 1e3, color=GREEN, lw=2)
ax.set_ylabel('B at the null [mT]'); ax.set_xlabel('time [µs]')
sw = (kBr[m2].max() - kBr[m2].min()) / 2 * 1e3
ax.set_title(f'The null barely moves: swing ±{sw:.1f} mT, against 1.3 T at the knee '
             '(the loaded output coil holds it)', fontsize=10)
fig.tight_layout(); fig.savefig(f'{OUT}/knife_null_zoom.png'); plt.close(fig)

# B-H paths
fig, ax = plt.subplots(figsize=(8, 5))
hh = np.linspace(-2500, 2500, 400)
ax.plot(hh, np.interp(hh, H, B), color=GREY, lw=2, label='core B-H curve')
ax.plot(meg['v(l2)'][m] / LEG, megBL[m], color=BLUE, lw=4, label='MEG side legs')
ax.plot((kn['v(al)'] - kn['v(jl)'])[m2] / SEG, kA[m2], color=ORANGE, lw=4, label='knife: bar at the magnets')
ax.plot(kn['v(jl)'][m2] / LEG, kRL[m2], color=PURPLE, lw=4, label='knife: return yokes')
ax.plot((kn['v(y1)'] - kn['v(n)'])[m2] / SEG, kBr[m2], color=GREEN, lw=5, label='knife: null (drive coils)')
ax.axhline(1.3, color=RED, lw=0.8, ls=':'); ax.text(-2400, 1.33, 'knee ≈ 1.3 T', color=RED, fontsize=8)
ax.set_xlabel('H [A/m]'); ax.set_ylabel('B [T]'); ax.set_xlim(-2500, 2500)
ax.set_title('Where each part operates on the B-H curve (each path is a dot: the swing is tiny)', fontsize=10)
ax.legend(fontsize=8, loc='lower right'); fig.tight_layout(); fig.savefig(f'{OUT}/bh_operating_paths.png'); plt.close(fig)

# energy split
fig, ax = plt.subplots(figsize=(8.5, 3.6))
best = max(rows, key=lambda r: r['COP'])
base = [('MEG-style\n5 V, 10 kHz, 100 Ω', mb), ('Bucking magnets\n2 V, 10 kHz, 100 Ω', kb),
        (f"Bucking magnets\n2 V, 10 kHz, {best['RL_ohm']/1e3:g} kΩ (best)", best)]
for i, (lab, r) in enumerate(base):
    pd, pl = r['P_drive_W'], r['P_load_W']
    ax.barh(i, pl / pd * 100, color=GREEN, label='to the load' if i == 0 else None)
    ax.barh(i, (pd - pl) / pd * 100, left=pl / pd * 100, color=GREY, alpha=0.6, label='lost in windings and leaks' if i == 0 else None)
    ax.text(101.5, i, f'{pl:.3g} W of {pd:.3g} W', va='center', fontsize=9)
ax.axvline(100, color=RED, ls='--', lw=1.2); ax.text(100, 2.45, ' break-even', color=RED, fontsize=8)
ax.set_yticks(range(len(base))); ax.set_yticklabels([b[0] for b in base]); ax.set_xlim(0, 140)
ax.set_xlabel('% of drive energy'); ax.legend(loc='lower right', fontsize=8)
ax.set_title('Where the drive energy goes (magnets contribute ≈ 0 net)', fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/energy_split.png'); plt.close(fig)

# resonance
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for RL, c in ((1000, BLUE), (10000, ORANGE)):
    rr = sorted([r for r in rows if r['sweep'].startswith('knife: series') and r['RL_ohm'] == RL], key=lambda r: r['Cres_F'])
    C = [r['Cres_F'] * 1e9 for r in rr]
    ax[0].plot(C, [r['P_drive_W'] for r in rr], 'o-', color=c, label=f'drive, {RL/1e3:g} kΩ load')
    ax[0].plot(C, [r['P_load_W'] for r in rr], 's--', color=c, alpha=0.7, label=f'load, {RL/1e3:g} kΩ load')
    ax[1].plot(C, [r['COP'] for r in rr], 'o-', color=c, label=f'{RL/1e3:g} kΩ load')
ax[0].set_yscale('log'); ax[0].set_xlabel('series capacitor [nF]'); ax[0].set_ylabel('power [W]')
ax[0].set_title('Resonance multiplies the power drawn…', fontsize=10); ax[0].legend(fontsize=7)
ax[1].axhline(1, color=RED, ls='--'); ax[1].text(35, 1.02, 'break-even', color=RED, fontsize=8)
ax[1].set_ylim(0.9, 1.04); ax[1].set_xlabel('series capacitor [nF]'); ax[1].set_ylabel('COP')
ax[1].set_title('…but not the share reaching the load', fontsize=10); ax[1].legend(fontsize=8, loc='lower right')
fig.tight_layout(); fig.savefig(f'{OUT}/resonance.png'); plt.close(fig)

# null leakage
rr = [r for r in rows if r['sweep'] == 'knife: air leakage at the null']
fig, ax = plt.subplots(figsize=(7, 3.6))
L = [1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
ax.semilogx(L, [r['COP'] for r in rr], 'o-', color=BLUE)
ax.axhline(kb['COP'], color=GREY, ls=':'); ax.text(1.1e-8, kb['COP'] + 0.003, 'no leakage', color=GREY, fontsize=8)
ax.set_xlabel('air-leakage permeance at the null [Wb/A]'); ax.set_ylabel('COP')
ax.set_title('Flux leaking out at the null bypasses the output coil', fontsize=10)
fig.tight_layout(); fig.savefig(f'{OUT}/null_leakage.png'); plt.close(fig)

# every run
fig, ax = plt.subplots(figsize=(9, 5.6))
groups = []
for r in rows:
    if r['sweep'] not in groups:
        groups.append(r['sweep'])
for r in rows:
    y = groups.index(r['sweep'])
    ax.plot(r['COP'], y, 'o', ms=7, color=BLUE if r is best else GREY, alpha=1 if r is best else 0.6)
ax.annotate(f"best {best['COP']:.3f}\n({best['RL_ohm']/1e3:g} kΩ load)", (best['COP'], groups.index(best['sweep'])),
            xytext=(0.80, groups.index(best['sweep']) + 1.3), fontsize=8, color=BLUE, arrowprops=dict(arrowstyle='-', color=BLUE))
ax.axvline(1, color=RED, ls='--'); ax.text(1.01, len(groups) - 0.6, 'break-even', color=RED, fontsize=8)
ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups, fontsize=8); ax.invert_yaxis()
ax.set_xlim(0, 1.15); ax.set_xlabel('COP = load power ÷ drive power')
ax.set_title(f'All {len(rows)} KronSpice runs: none reached break-even', fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/cop_all_runs.png'); plt.close(fig)

# ---------------------------------------------------------------- drawings
NORM = colors.Normalize(0, 1.6)
CMAP = matplotlib.colormaps['YlOrRd']
W = 0.7


def arrow(ax, x, y, dx, dy, mag):
    ax.add_patch(FancyArrow(x, y, dx, dy, width=0.04 + 0.16 * mag, head_width=0.12 + 0.3 * mag,
                            head_length=0.2, length_includes_head=True, color='#1a5276'))


def coil_v(ax, x, y, label, col, side):
    for k in range(5):
        ax.add_patch(Rectangle((x - W/2 - 0.12, y + k * 0.13), W + 0.24, 0.07, color=col))
    ax.text(x + side * 0.95, y + 0.3, label, ha='center', va='center', fontsize=7, color=col)


def coil_h(ax, x, y, label, col, below=False):
    for k in range(4):
        ax.add_patch(Rectangle((x + k * 0.13, y - W/2 - 0.12), 0.07, W + 0.24, color=col))
    ax.text(x + 0.22, y + (-0.75 if below else 0.72), label, ha='center', va='center', fontsize=7, color=col)


def draw_meg(ax, b, title):
    ax.clear(); ax.set_xlim(-0.6, 6.6); ax.set_ylim(-0.4, 5.2); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=10)
    ax.add_patch(Rectangle((0.65, 4.0), 4.7, 0.6, color='#bbbbbb'))
    ax.add_patch(Rectangle((0.65, 0.4), 4.7, 0.6, color='#bbbbbb'))
    X = {'L': 1.0, 'C': 3.0, 'R': 5.0}
    for leg, x in X.items():
        ax.add_patch(Rectangle((x - W/2, 1.0), W, 3.0, color=CMAP(NORM(abs(b[leg]))), ec='k', lw=0.6))
    ax.add_patch(Rectangle((3 - W/2 - 0.05, 2.1), W + 0.1, 0.8, color='#555', ec='k'))
    ax.text(3, 2.72, 'N', ha='center', va='center', color='w', fontsize=9, weight='bold')
    ax.text(3, 2.28, 'S', ha='center', va='center', color='w', fontsize=9, weight='bold')
    coil_v(ax, 1.0, 3.1, 'control', BLUE, -1); coil_v(ax, 1.0, 1.3, 'output', GREEN, -1)
    coil_v(ax, 5.0, 3.1, 'control', ORANGE, 1); coil_v(ax, 5.0, 1.3, 'output', GREEN, 1)
    for leg, x in X.items():
        v = b[leg]; mg = min(abs(v) / 1.6, 1)
        arrow(ax, x + (0.55 if leg != 'R' else -0.55), 2.5 - 0.45 * np.sign(v), 0, 0.9 * np.sign(v), mg)
        ax.text(x, 4.75, f'{abs(v):.2f} T', ha='center', fontsize=8)


def draw_knife(ax, b, title, null_scale=1):
    """b: A, B (bar at magnets, + toward the null), RL, RR (returns), N (null, + rightward)."""
    ax.clear(); ax.set_xlim(-0.2, 6.2); ax.set_ylim(-0.3, 5.3); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=10)
    grey = '#bbbbbb'
    ax.add_patch(Rectangle((0.0, 4.2), 6.0, 0.5, color=grey))          # top yoke
    ax.add_patch(Rectangle((0.0, 0.3), 6.0, 0.5, color=grey))          # bottom yoke
    ax.add_patch(Rectangle((0.0, 0.3), 0.5, 4.4, color=grey))          # left side (A's S end)
    ax.add_patch(Rectangle((5.5, 0.3), 0.5, 4.4, color=grey))          # right side (B's S end)
    yb = 2.5
    # return posts from the junctions to the yokes
    for x, v in ((1.9, b['RL']), (4.1, b['RR'])):
        ax.add_patch(Rectangle((x - 0.25, 0.8), 0.5, 3.4, color=CMAP(NORM(abs(v))), ec='k', lw=0.5))
    # bar segments
    ax.add_patch(Rectangle((1.3, yb - W/2), 0.35, W, color=CMAP(NORM(abs(b['A']))), ec='k', lw=0.5))
    ax.add_patch(Rectangle((4.35, yb - W/2), 0.35, W, color=CMAP(NORM(abs(b['B']))), ec='k', lw=0.5))
    ax.add_patch(Rectangle((2.15, yb - W/2), 1.7, W, color=CMAP(NORM(min(abs(b['N']) * null_scale, 1.6))), ec='k', lw=0.6))
    # magnets
    for x0, lab, left_N in ((0.5, 'A', False), (4.7, 'B', True)):
        ax.add_patch(Rectangle((x0, yb - W/2 - 0.05), 0.8, W + 0.1, color='#555', ec='k'))
        ax.text(x0 + 0.2, yb, 'N' if left_N else 'S', ha='center', va='center', color='w', fontsize=9, weight='bold')
        ax.text(x0 + 0.6, yb, 'S' if left_N else 'N', ha='center', va='center', color='w', fontsize=9, weight='bold')
        ax.text(x0 + 0.4, yb + 0.62, lab, ha='center', fontsize=9)
    coil_h(ax, 2.25, yb, 'drive 1', BLUE); coil_h(ax, 3.25, yb, 'drive 2', BLUE)
    coil_h(ax, 2.78, yb, 'output', GREEN, below=True)
    ax.text(3.0, yb + 1.05, 'null', ha='center', fontsize=8, style='italic')
    # flux arrows
    arrow(ax, 0.9 - 0.3 * np.sign(b['A']), yb - 0.75, 0.6 * np.sign(b['A']), 0, min(abs(b['A']) / 1.6, 1))
    arrow(ax, 5.1 + 0.3 * np.sign(b['B']), yb - 0.75, -0.6 * np.sign(b['B']), 0, min(abs(b['B']) / 1.6, 1))
    for x, v in ((1.9, b['RL']), (4.1, b['RR'])):
        mg = min(abs(v) / 1.6, 1)
        arrow(ax, x + (-0.45 if x < 3 else 0.45), 3.0, 0, 0.8 * np.sign(v), mg)
        arrow(ax, x + (-0.45 if x < 3 else 0.45), 2.0, 0, -0.8 * np.sign(v), mg)
    vn = b['N'] * null_scale; mg = min(abs(vn) / 1.6, 1)
    if mg > 0.01:
        arrow(ax, 3.0 - 0.35 * np.sign(vn), yb - 1.0, 0.7 * np.sign(vn), 0, mg)
    lab = f"null {abs(b['N'])*1e3:.1f} mT" + (f' (colour ×{null_scale})' if null_scale != 1 else '')
    ax.text(3.0, 0.05, lab, ha='center', fontsize=8)
    ax.text(1.9, 4.85, f"{abs(b['RL']):.2f} T", ha='center', fontsize=8)
    ax.text(4.1, 4.85, f"{abs(b['RR']):.2f} T", ha='center', fontsize=8)


def colorbar(fig, ax):
    sm = cm.ScalarMappable(norm=NORM, cmap=CMAP); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02); cb.set_label('|B| [T]', fontsize=8)


k0 = np.argmin(np.abs(meg['time'] - 2.02e-3))
fig, ax = plt.subplots(figsize=(6.4, 5.4))
draw_meg(ax, {'L': -megBL[k0], 'C': megBL[k0] + megBR[k0], 'R': -megBR[k0]},
         'MEG-style core: magnet in the centre leg,\ncontrol and output windings on each side leg')
colorbar(fig, ax); fig.tight_layout(); fig.savefig(f'{OUT}/meg_diagram.png'); plt.close(fig)

k1 = np.argmin(np.abs(kn['time'] - 11.025e-3))
kd = lambda k: {'A': kA[k], 'B': kB[k], 'RL': kRL[k], 'RR': kRR[k], 'N': kBr[k]}
fig, ax = plt.subplots(figsize=(6.6, 5.6))
draw_knife(ax, kd(k1), 'Bucking magnets (N facing N): each returns through its own yoke;\n'
           'two drive coils either side of the null, output at the null')
colorbar(fig, ax); fig.tight_layout(); fig.savefig(f'{OUT}/knife_diagram.png'); plt.close(fig)


# ---------------------------------------------------------------- animations
def animate(draw, d, legs, traces, labels, ylabel, title, fname, t0, t1, frames=160, **kw):
    m = (d['time'] >= t0) & (d['time'] <= t1)
    t = d['time'][m]
    ts = np.linspace(t0, t1, frames)
    fig = plt.figure(figsize=(11.6, 5.2))
    axc = fig.add_axes([0.0, 0.02, 0.46, 0.88]); axw = fig.add_axes([0.54, 0.14, 0.43, 0.70])
    for tr, lab, c in zip(traces, labels, [BLUE, ORANGE, GREEN, GREY]):
        axw.plot((t - t0) * 1e6, tr[m], color=c, label=lab)
    axw.set_xlabel('time [µs]'); axw.set_ylabel(ylabel); axw.legend(fontsize=8, loc='upper right')
    cur = axw.axvline(0, color='k', lw=1)
    fig.suptitle(title, fontsize=11)

    def frame(i):
        k = np.searchsorted(d['time'], ts[i])
        draw(axc, legs(k), f't = {(ts[i]-t0)*1e6:5.1f} µs', **kw)
        cur.set_xdata([(ts[i] - t0) * 1e6] * 2)
        return []
    ani = animation.FuncAnimation(fig, frame, frames=frames)
    ani.save(f'{OUT}/{fname}', writer=animation.FFMpegWriter(fps=24, bitrate=2400, codec='libx264',
             extra_args=['-pix_fmt', 'yuv420p', '-preset', 'veryfast']))
    plt.close(fig)


animate(draw_meg, meg, lambda k: {'L': -megBL[k], 'C': megBL[k] + megBR[k], 'R': -megBR[k]},
        [megBL, megBR, meg['v(ol)'] / 10], ['left leg B [T]', 'right leg B [T]', 'left output ÷ 10 [V]'], 'B [T] / scaled V',
        'MEG-style core: control pulses steer the magnet flux between the side legs (KronSpice)', 'meg_flux.mp4', 2.0e-3, 2.2e-3)
animate(draw_knife, kn, kd, [kBr * 1e3, kn['v(o1)'], -kn['i(vd)'] * 10],
        ['null B [mT]', 'output [V]', 'drive current × 10 [A]'], 'mT / V / A',
        'Bucking magnets: the split drive coil rocks the null; the output coil reads it (KronSpice)', 'knife_flux.mp4',
        11.0e-3, 11.2e-3, null_scale=50)
print('done', sorted(os.listdir(OUT)))
