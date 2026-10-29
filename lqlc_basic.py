"""
Basic LQLC single-point positioning
-----------------------------------
Unknowns per epoch (4):  x, y, z  (ECEF, m)   and   cdt  (receiver clock x c, m)
Input CSV columns used:  epoch, tow, sat_const, sat_x, sat_y, sat_z, cor_pseudo
Constellation:           GPS only (sat_const == 1), so one clock is enough

Cost:     J(theta) = sum_i  ln cosh( r_i / (2 s) )
Solver:   IRLS = Gauss-Newton where every step is a weighted LS
          weight  w_i = tanh(x_i) / x_i ,  x_i = |r_i| / (2 s)

Usage:    python lqlc_basic.py sppStruct_user.csv
"""
import sys
import numpy as np
import pandas as pd

# ---------------- settings ----------------
CSV        = sys.argv[1] if len(sys.argv) > 1 else "data.csv"
S          = 3.0     # LQLC scale s (m). Large s -> behaves like LS
ELEV_MASK  = 10.0    # degrees; elevation is computed from OUR position, not the file
MAX_ITER   = 30
TOL        = 1e-3    # stop when the position step is < 1 mm


# ---------------- helpers ----------------
def ecef_to_latlon(p):
    """ECEF (m) -> geodetic lat, lon (rad), height (m), WGS-84."""
    a, e2 = 6378137.0, 6.69437999014e-3
    x, y, z = p
    lon = np.arctan2(y, x)
    r = np.hypot(x, y)
    lat = np.arctan2(z, r * (1 - e2))
    for _ in range(5):
        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = r / np.cos(lat) - N
        lat = np.arctan2(z, r * (1 - e2 * N / (N + h)))
    return lat, lon, h


def enu_matrix(lat, lon):
    """Rotation ECEF difference -> local East, North, Up."""
    sl, cl, so, co = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    return np.array([[-so,       co,      0],
                     [-sl * co, -sl * so, cl],
                     [ cl * co,  cl * so, sl]])


def elevation(sat_xyz, rx_xyz):
    """Elevation (deg) of each satellite seen from rx_xyz."""
    lat, lon, _ = ecef_to_latlon(rx_xyz)
    enu = (sat_xyz - rx_xyz) @ enu_matrix(lat, lon).T
    return np.degrees(np.arctan2(enu[:, 2], np.hypot(enu[:, 0], enu[:, 1])))


def lqlc_weight(res, s):
    """w = tanh(x)/x with x = |r|/(2s);  w -> 1 for small residuals."""
    x = np.abs(res) / (2 * s)
    return np.where(x < 1e-8, 1.0, np.tanh(x) / np.maximum(x, 1e-12))


def solve_epoch(sat_xyz, rho, theta0, s, robust=True):
    """
    Gauss-Newton / IRLS for 4 unknowns.
    Phase 1: plain LS until converged (gets a sane position from any start).
    Phase 2: LQLC re-weighting until converged.
    Returns theta = [x, y, z, cdt], residuals, weights, iterations used.
    """
    theta = theta0.copy()
    w = np.ones(len(rho))
    phase = "LS"
    it_total = 0
    for it in range(2 * MAX_ITER):
        it_total += 1
        # 1) predicted pseudorange at the current guess
        vec = sat_xyz - theta[:3]
        rng = np.linalg.norm(vec, axis=1)
        pred = rng + theta[3]
        res = rho - pred                          # residual r_i

        # 2) LQLC weights (only in phase 2)
        if phase == "LQLC":
            w = lqlc_weight(res, s)

        # 3) linearise: H row = [-unit LOS vector, 1]
        H = np.column_stack([-vec / rng[:, None], np.ones(len(rho))])

        # 4) weighted LS step:  dtheta = (H^T W H)^-1 H^T W r
        W = np.diag(w)
        dtheta = np.linalg.solve(H.T @ W @ H, H.T @ W @ res)
        theta += dtheta

        # 5) convergence check
        if np.linalg.norm(dtheta[:3]) < TOL:
            if phase == "LS" and robust:
                phase = "LQLC"                    # switch on robust weights
            else:
                break
    vec = sat_xyz - theta[:3]
    res = rho - (np.linalg.norm(vec, axis=1) + theta[3])
    return theta, res, w, it_total


# ---------------- main loop ----------------
data = pd.read_csv(CSV)
data = data[data.sat_const == 1]                  # GPS only -> one clock

results = []
theta = np.zeros(4)                               # first epoch starts at Earth's centre
for ep, g in data.groupby("epoch"):
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    rho = g["cor_pseudo"].values

    # elevation mask from the previous fix (skip on the very first epoch)
    if np.linalg.norm(theta[:3]) > 6.0e6:
        keep = elevation(sat_xyz, theta[:3]) >= ELEV_MASK
        sat_xyz, rho = sat_xyz[keep], rho[keep]
    if len(rho) < 5:                              # need > 4 to have any redundancy
        continue

    for method, robust in [("LS", False), ("LQLC", True)]:
        th, res, w, nit = solve_epoch(sat_xyz, rho, theta, S, robust)
        results.append(dict(epoch=ep, tow=g.tow.iloc[0], method=method,
                            x=th[0], y=th[1], z=th[2], cdt=th[3],
                            nsat=len(rho), iters=nit, min_w=w.min()))
        if method == "LQLC":
            theta = th                            # next epoch starts from this fix

out = pd.DataFrame(results)
out.to_csv("lqlc_basic_out.csv", index=False)

# ---------------- quick report ----------------
# no surveyed truth in the file: report spread around each method's own mean
for method, r in out.groupby("method"):
    P = r[["x", "y", "z"]].values
    mean = P.mean(axis=0)
    lat, lon, h = ecef_to_latlon(mean)
    enu = (P - mean) @ enu_matrix(lat, lon).T
    hz = np.hypot(enu[:, 0], enu[:, 1])
    print(f"{method:5s} mean lat {np.degrees(lat):.6f} lon {np.degrees(lon):.6f} h {h:6.1f} m | "
          f"horiz spread RMS {np.sqrt((hz**2).mean()):5.2f} m, 95% {np.percentile(hz,95):5.2f} m | "
          f"vert std {enu[:,2].std():5.2f} m | sats {r.nsat.mean():.1f} | iters {r.iters.mean():.1f}")
print("saved lqlc_basic_out.csv")
