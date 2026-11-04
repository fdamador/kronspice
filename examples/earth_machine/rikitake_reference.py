import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
mu, A = 1.0, 1.0
def f(t, s):
    x, y, z = s
    return [-mu*x + y*z, -mu*y + (z - A)*x, 1 - x*y]
sol = solve_ivp(f, [0, 1000], [1, 0.5, 1], max_step=0.05, rtol=1e-8)
t = sol.t; x, y, z = sol.y; m = t > 50
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
for a, v, l in zip(ax, [x, y, z], ["I1", "I2", "omega"]):
    a.plot(t[m], v[m], lw=0.6); a.set_ylabel(l); a.grid(alpha=.3)
ax[-1].set_xlabel("Dimensionless time"); fig.tight_layout()
fig.savefig("fig1_timeseries.png", dpi=300)
fig = plt.figure(figsize=(8, 7)); a = fig.add_subplot(projection="3d")
a.plot(x[m], y[m], z[m], lw=0.3); fig.savefig("fig2_attractor.png", dpi=300)
s = np.sign(x[m]); idx = np.where(s[1:] != s[:-1])[0]; iv = np.diff(t[m][idx])
fig, a = plt.subplots(figsize=(8, 5)); a.hist(iv, bins=20, edgecolor="k")
a.set_xlabel("Interval between reversals"); fig.savefig("fig3_reversals.png", dpi=300)
print("reversals:", len(idx), "interval min/mean/max:", iv.min(), iv.mean(), iv.max())
