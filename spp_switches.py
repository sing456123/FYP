"""
SPP with correction switches: LS / WLS / LQLC
=============================================
Flow per epoch (as discussed):
  1) once per epoch : rho0 = raw (+ fixed satellite terms: clock, relativity, TGD)
  2) every iteration: from current theta_n -> el, az -> tropo, iono, Sagnac, mask, weights
                      rho_c = rho0 - iono - tropo        (always fresh from rho0)
                      y = rho_c - pred(theta_n);  dtheta = (H^T W H)^-1 H^T W y
  3) next epoch starts from this epoch's fix

Switches are in CONFIGS below. Each config is one run; a comparison table is printed.

What can be computed from sppStruct_user.csv alone:
  friend_cor=False -> raw_pseudo, NO corrections at all (baseline)
  friend_cor=True  -> cor_pseudo, the friend's tool corrections (sat clock+iono+tropo, all-in-one)
  tropo, sagnac, mask, weights, estimator -> computed here
Needs the broadcast navigation file (BRDC, 2022 day 173) -> not available yet:
  sat_clock, relativity, tgd, iono   (stubs raise a clear error until NAV is loaded)

Usage:  python spp_switches.py sppStruct_user.csv
"""
import sys
import numpy as np
import pandas as pd

C = 299792458.0
OMEGA_E = 7.2921151467e-5
SYS_NAME = {1: "GPS", 4: "GLO", 8: "GAL", 32: "BDS"}

# ============================== SWITCHES ==============================
BASE = dict(
    friend_cor=True,     # friend's lumped sat clock+iono+tropo (uses cor_pseudo); False = raw_pseudo
    sat_clock=False,     # needs NAV  (turn friend_cor off first)
    relativity=False,    # needs NAV
    tgd=False,           # needs NAV
    iono=False,          # needs NAV  (Klobuchar alpha/beta)
    tropo=False,         # Saastamoinen, recomputed every iteration
    sagnac=False,        # Earth-rotation range term, recomputed every iteration
    systems=[1],         # [1]=GPS only (4 unknowns) | [1,4,8,32]=all (one clock per system)
    mask_deg=10.0,       # elevation mask, recomputed every iteration
    weight="none",       # 'none' | 'elev' (sin^2 el) | 'snr' (10^(snr/10))
    estimator="LS",      # 'LS' | 'LQLC'
    s=3.0,               # LQLC scale (m)
)
CONFIGS = {   # name: changes from BASE
    "0 no corrections, LS":       dict(friend_cor=False),
    "1 corrected, LS":            dict(),
    "2 corrected, WLS elev":      dict(weight="elev"),
    "3 corrected, LQLC":          dict(estimator="LQLC"),
    "4 corrected, LQLC + elev":   dict(estimator="LQLC", weight="elev"),
    "5 all GNSS, LS":             dict(systems=[1, 4, 8, 32]),
    "6 all GNSS, LQLC":           dict(systems=[1, 4, 8, 32], estimator="LQLC"),
}
REF_XYZ = None    # put surveyed antenna ECEF (m) here, e.g. np.array([x, y, z]); None = no truth
MAX_ITER, TOL = 30, 1e-3
NAV = None        # will hold broadcast ephemeris once a BRDC file is provided
# ======================================================================


# ---------------- geometry helpers ----------------
def ecef_to_lla(p):
    a, e2 = 6378137.0, 6.69437999014e-3
    x, y, z = p
    lon = np.arctan2(y, x); r = np.hypot(x, y); lat = np.arctan2(z, r * (1 - e2))
    for _ in range(5):
        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = r / np.cos(lat) - N
        lat = np.arctan2(z, r * (1 - e2 * N / (N + h)))
    return lat, lon, h


def enu_matrix(lat, lon):
    sl, cl, so, co = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    return np.array([[-so, co, 0], [-sl * co, -sl * so, cl], [cl * co, cl * so, sl]])


def el_az(sat_xyz, rx):
    lat, lon, _ = ecef_to_lla(rx)
    enu = (sat_xyz - rx) @ enu_matrix(lat, lon).T
    el = np.arctan2(enu[:, 2], np.hypot(enu[:, 0], enu[:, 1]))
    az = np.arctan2(enu[:, 0], enu[:, 1]) % (2 * np.pi)
    return el, az


# ---------------- correction models ----------------
def tropo_saastamoinen(el, rx):
    """Saastamoinen with standard atmosphere (RTKLIB form). el in rad. Returns metres."""
    lat, _, h = ecef_to_lla(rx)
    h = max(h, 0.0)
    pres = 1013.25 * (1 - 2.2557e-5 * h) ** 5.2568
    temp = 15.0 - 6.5e-3 * h + 273.16
    e = 6.108 * 0.7 * np.exp((17.15 * temp - 4684.0) / (temp - 38.45))
    sinel = np.sin(np.maximum(el, np.radians(1.0)))
    zhd = 0.0022768 * pres / (1 - 0.00266 * np.cos(2 * lat) - 0.00028 * h / 1e3)
    zwd = 0.002277 * (1255.0 / temp + 0.05) * e
    return (zhd + zwd) / sinel


def sagnac(sat_xyz, rx):
    """Earth rotation during signal travel, as a range term (m)."""
    return OMEGA_E / C * (sat_xyz[:, 0] * rx[1] - sat_xyz[:, 1] * rx[0])


def need_nav(name):
    raise RuntimeError(f"'{name}' needs the broadcast navigation file (BRDC 2022 day 173). "
                       "Upload it and this switch will be implemented and tested.")


def fixed_sat_terms(g, cfg):
    """Once per epoch: rho0 = pseudorange + fixed satellite terms."""
    rho0 = (g.cor_pseudo if cfg["friend_cor"] else g.raw_pseudo).values.astype(float)
    if cfg["friend_cor"] and any(cfg[k] for k in ("sat_clock", "iono", "tropo", "tgd", "relativity")):
        raise ValueError("friend_cor already contains sat clock/iono/tropo -> would double count")
    if cfg["sat_clock"]:  need_nav("sat_clock")    # rho0 += c*dts
    if cfg["relativity"]: need_nav("relativity")   # rho0 += c*dtr
    if cfg["tgd"]:        need_nav("tgd")          # rho0 -= c*TGD
    return rho0


def iter_terms(sat_xyz, theta, cfg, have_pos):
    """Every iteration: position-dependent terms from the CURRENT theta."""
    n = len(sat_xyz)
    if not have_pos:                      # theta still at Earth's centre: nothing computable
        return np.zeros(n), np.ones(n, bool), np.full(n, np.nan)
    rx = theta[:3]
    el, az = el_az(sat_xyz, rx)
    corr = np.zeros(n)
    if cfg["tropo"]:  corr += tropo_saastamoinen(el, rx)
    if cfg["iono"]:   need_nav("iono")     # corr += klobuchar(el, az, rx, alpha, beta)
    if cfg["sagnac"]: corr += sagnac(sat_xyz, rx)   # added to the prediction side
    keep = el >= np.radians(cfg["mask_deg"])
    return corr, keep, el


def prior_weight(el, snr, cfg):
    if cfg["weight"] == "elev":
        return np.where(np.isnan(el), 1.0, np.sin(np.nan_to_num(el, nan=np.pi / 2)) ** 2)
    if cfg["weight"] == "snr":
        return 10 ** (snr / 10) / 10 ** (45 / 10)
    return np.ones(len(snr))


def lqlc_weight(y, s):
    x = np.abs(y) / (2 * s)
    return np.where(x < 1e-8, 1.0, np.tanh(x) / np.maximum(x, 1e-12))


# ---------------- one epoch ----------------
def solve_epoch(g, theta0, cfg):
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    sysid = g.sat_const.values
    snr = g.sat_snr.values.astype(float)
    systems = [s for s in cfg["systems"] if (sysid == s).any()]
    col = {s: 3 + k for k, s in enumerate(systems)}

    rho0 = fixed_sat_terms(g, cfg)                                # (1) once per epoch
    theta = np.r_[theta0[:3], np.zeros(len(systems))]
    if len(theta0) == len(theta): theta = theta0.copy()
    robust_on = False
    for it in range(2 * MAX_ITER):
        have_pos = np.linalg.norm(theta[:3]) > 6.0e6
        corr, keep, el = iter_terms(sat_xyz, theta, cfg, have_pos)  # (2) every iteration
        vec = sat_xyz - theta[:3]
        rng = np.linalg.norm(vec, axis=1)
        clk = np.array([theta[col[s]] for s in sysid])
        y = rho0 - (rng + clk + corr)                               # tropo/sagnac move with theta
        H = np.zeros((len(y), len(theta)))
        H[:, :3] = -vec / rng[:, None]
        for i, s in enumerate(sysid): H[i, col[s]] = 1.0

        w = prior_weight(el, snr, cfg) * keep
        if robust_on:
            w = w * lqlc_weight(y, cfg["s"])                        # prior trust x data trust
        if (keep.sum() <= len(theta)) or any(w[sysid == s].sum() == 0 for s in systems):
            return None
        HW = H * w[:, None]
        dtheta = np.linalg.solve(H.T @ HW, HW.T @ y)                # (H^T W H)^-1 H^T W y
        theta = theta + dtheta
        if np.linalg.norm(dtheta[:3]) < TOL:
            if cfg["estimator"] == "LQLC" and not robust_on:
                robust_on = True                                    # LS first, then LQLC
            else:
                return dict(theta=theta, nsat=int(keep.sum()), iters=it + 1, converged=True)
    return dict(theta=theta, nsat=int(keep.sum()), iters=it + 1, converged=False)


# ---------------- run all configs ----------------
def run(cfg, data):
    d = data[data.sat_const.isin(cfg["systems"])]
    theta = np.zeros(3)
    out = []
    for ep, g in d.groupby("epoch"):
        r = solve_epoch(g, theta, cfg)
        if r is None: continue
        out.append(dict(epoch=ep, x=r["theta"][0], y=r["theta"][1], z=r["theta"][2],
                        nsat=r["nsat"], iters=r["iters"], conv=r["converged"]))
        if r["converged"] and np.linalg.norm(r["theta"][:3]) > 6.0e6:
            theta = r["theta"]
    return pd.DataFrame(out)


def score(df, ref):
    P = df[["x", "y", "z"]].values
    mean = P.mean(0)
    lat, lon, h = ecef_to_lla(mean)
    R = enu_matrix(lat, lon)
    sp = (P - mean) @ R.T
    hs = np.hypot(sp[:, 0], sp[:, 1])
    row = dict(epochs=len(df), conv_pct=100 * df.conv.mean(), sats=df.nsat.mean(), iters=df.iters.mean(),
               h_spread=np.sqrt((hs ** 2).mean()), h95_spread=np.percentile(hs, 95), v_spread=sp[:, 2].std(),
               mean_h=h)
    if ref is not None:
        e = (P - ref) @ R.T
        hz = np.hypot(e[:, 0], e[:, 1])
        row.update(hRMS=np.sqrt((hz ** 2).mean()), h95=np.percentile(hz, 95), vRMS=np.sqrt((e[:, 2] ** 2).mean()),
                   biasE=e[:, 0].mean(), biasN=e[:, 1].mean(), biasU=e[:, 2].mean())
    else:
        lat0, lon0, _ = ecef_to_lla(mean)
        row.update(lat=np.degrees(lat0), lon=np.degrees(lon0))
    return row


if __name__ == "__main__":
    data = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "data.csv")
    rows = []
    for name, change in CONFIGS.items():
        cfg = {**BASE, **change}
        try:
            df = run(cfg, data)
            rows.append(dict(config=name, **score(df, REF_XYZ)))
            df.to_csv(f"out_{name.split()[0]}.csv", index=False)
        except Exception as ex:
            rows.append(dict(config=name, error=str(ex)[:80]))
        print("done:", name)
    T = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(T.round(2).to_string(index=False))
    T.to_csv("comparison.csv", index=False)
