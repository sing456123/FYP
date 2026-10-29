"""
SPP + per-satellite noise (RTKLIB model) inside LQLC   (copy of spp_switches.py with sigma_i added)
=================================================================================
One shared s treats every satellite as equally noisy. They are not: in feas1 the robust noise level
was ~3.2 m at 60-90 deg vs ~7.6 m at 0-15 deg, and GLONASS ~1.5x the others. RTKLIB's pntpos uses
the same idea for its weights (varerr), with default settings and no fitting to this data:

  sigma_i = F_sys * 100 * sqrt(0.003^2 + 0.003^2 / sin(el_i))     F = 1 (GPS/GAL/BDS), 1.5 (GLO)
  scaled to sigma_i / sigma_ref, where sigma_ref = GPS at zenith

Two switches, usable together:
  sig_w : prior weight 1/sigma_i^2               (RTKLIB-style weighted LS)
  sig_s : LQLC judges y_i / sigma_i               (= a per-satellite scale s_i = s * sigma_i/sigma_ref)

Usage:  python spp_sigma.py data.csv          -> results go to the next free run<N>/ folder
        python spp_sigma.py data.csv run1     -> re-run run1/ from its own config.json
"""
import sys
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

C = 299792458.0
OMEGA_E = 7.2921151467e-5
SYS_NAME = {1: "GPS", 4: "GLO", 8: "GAL", 32: "BDS"}
FREQ = {1: 1575.42e6, 4: 1602.0e6, 8: 1575.42e6, 32: 1561.098e6}   # L1 / G1 (nominal) / E1 / B1I
EFACT = {1: 1.0, 4: 1.5, 8: 1.0, 32: 1.0}                          # RTKLIB EFACT_GPS/GLO/GAL/CMP
ION_DEFAULT = [0.1118e-07, -0.7451e-08, -0.5961e-07, 0.1192e-06,  # alpha  (RTKLIB ion_default)
               0.1167e+06, -0.2294e+06, -0.1311e+06, 0.1049e+07]  # beta

# ============================== SWITCHES ==============================
BASE = dict(
    sagnac=False,        # keep False: already in cor_pseudo (only for raw_pseudo)
    tropo=False,         # keep False: already in cor_pseudo (Saastamoinen, only for raw_pseudo)
    iono=False,          # keep False: already in cor_pseudo (Klobuchar, only for raw_pseudo)
    systems=[1,4,8,32],         # [1]=GPS only (4 unknowns) | [1,4,8,32]=all (one clock per system)
    mask_deg=5.0,       # elevation mask (deg)
    estimator="LQLC",    # 'LS' | 'LQLC'
    s=5.0,               # LQLC scale (m); with sig_s it is the scale for a GPS satellite at zenith
    sig_w=False,         # prior weight 1/sigma_i^2 (RTKLIB model)
    sig_s=False,         # LQLC on y_i/sigma_i, i.e. per-satellite s_i
)
CONFIGS = {             # name: changes from BASE
    "LS":            dict(estimator="LS"),
    "WLS_rtk":       dict(estimator="LS", sig_w=True),          # what RTKLIB SPP effectively does
    "LQLC":          dict(),                                     # baseline: one s for everyone
    "LQLC_sig":      dict(sig_s=True),
    "LQLC_sig_w":    dict(sig_s=True, sig_w=True),
    "LQLC_sig_s3":   dict(sig_s=True, s=3.0),                    # s_i grows with sigma_i, so try a smaller s too
    "LQLC_sig_w_s3": dict(sig_s=True, sig_w=True, s=3.0),
}
REF_XYZ = None    # put surveyed antenna ECEF (m) here, e.g. np.array([x, y, z]); None = no truth
MAX_ITER, TOL = 30, 1e-3
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
    if not -100.0 < h < 1e4:              # still converging from Earth's centre: model invalid
        return np.zeros(len(el))
    h = max(h, 0.0)
    pres = 1013.25 * (1 - 2.2557e-5 * h) ** 5.2568
    temp = 15.0 - 6.5e-3 * h + 273.16
    e = 6.108 * 0.7 * np.exp((17.15 * temp - 4684.0) / (temp - 38.45))
    sinel = np.sin(np.maximum(el, np.radians(1.0)))
    zhd = 0.0022768 * pres / (1 - 0.00266 * np.cos(2 * lat) - 0.00028 * h / 1e3)
    zwd = 0.002277 * (1255.0 / temp + 0.05) * e
    return (zhd + zwd) / sinel


def iono_klobuchar(el, az, rx, tow, ion=ION_DEFAULT):
    """Klobuchar L1 delay (m), RTKLIB ionmodel. el, az in rad, tow = GPS seconds of week."""
    lat, lon, _ = ecef_to_lla(rx)
    E = el / np.pi                                           # semicircles from here on
    psi = 0.0137 / (E + 0.11) - 0.022                        # earth-centred angle
    phi = np.clip(lat / np.pi + psi * np.cos(az), -0.416, 0.416)   # pierce-point latitude
    lam = lon / np.pi + psi * np.sin(az) / np.cos(phi * np.pi)     # pierce-point longitude
    phi = phi + 0.064 * np.cos((lam - 1.617) * np.pi)        # geomagnetic latitude
    t = (43200.0 * lam + tow) % 86400.0                      # local time (s)
    F = 1.0 + 16.0 * (0.53 - E) ** 3                         # slant factor
    amp = np.maximum(np.polyval(ion[3::-1], phi), 0.0)
    per = np.maximum(np.polyval(ion[:3:-1], phi), 72000.0)
    x = 2 * np.pi * (t - 50400.0) / per
    return C * F * np.where(np.abs(x) < 1.57, 5e-9 + amp * (1 - x ** 2 / 2 + x ** 4 / 24), 5e-9)


def sagnac(sat_xyz, rx):
    """Earth rotation during signal travel, as a range term (m)."""
    return OMEGA_E / C * (sat_xyz[:, 0] * rx[1] - sat_xyz[:, 1] * rx[0])


def iter_terms(sat_xyz, sysid, tow, theta, cfg, have_pos):
    """Every iteration: position-dependent terms from the CURRENT theta."""
    n = len(sat_xyz)
    if not have_pos:                      # theta still at Earth's centre: nothing computable
        return np.zeros(n), np.ones(n, bool)
    rx = theta[:3]
    el, az = el_az(sat_xyz, rx)
    corr = np.zeros(n)
    if cfg["sagnac"]: corr += sagnac(sat_xyz, rx)
    if cfg["tropo"]:  corr += tropo_saastamoinen(el, rx)
    if cfg["iono"]:   corr += iono_klobuchar(el, az, rx, tow) * (FREQ[1] / np.vectorize(FREQ.get)(sysid)) ** 2
    keep = el >= np.radians(cfg["mask_deg"])
    return corr, keep


def sigma_rel(el, sysid):
    """RTKLIB pntpos varerr() with default err = (100, 0.003, 0.003), divided by GPS-at-zenith."""
    sin_el = np.sin(np.maximum(el, np.radians(5.0)))
    sig = np.vectorize(EFACT.get)(sysid) * 100 * np.sqrt(0.003 ** 2 + 0.003 ** 2 / sin_el)
    return sig / (100 * np.sqrt(2 * 0.003 ** 2))


def lqlc_weight(y, s):
    x = np.abs(y) / (2 * s)
    return np.where(x < 1e-8, 1.0, np.tanh(x) / np.maximum(x, 1e-12))


# ---------------- one epoch ----------------
def solve_epoch(g, theta0, cfg):
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    sysid = g.sat_const.values
    tow = float(g.tow.iloc[0])
    rho = g.cor_pseudo.values.astype(float)
    systems = [s for s in cfg["systems"] if (sysid == s).any()]
    col = {s: 3 + k for k, s in enumerate(systems)}

    theta = np.r_[theta0[:3], np.zeros(len(systems))]
    if len(theta0) == len(theta): theta = theta0.copy()
    robust_on = False
    for it in range(2 * MAX_ITER):
        have_pos = np.linalg.norm(theta[:3]) > 6.0e6
        corr, keep = iter_terms(sat_xyz, sysid, tow, theta, cfg, have_pos)
        vec = sat_xyz - theta[:3]
        rng = np.linalg.norm(vec, axis=1)
        clk = np.array([theta[col[s]] for s in sysid])
        y = rho - (rng + clk + corr)                                # corrections move with theta
        H = np.zeros((len(y), len(theta)))
        H[:, :3] = -vec / rng[:, None]
        for i, s in enumerate(sysid): H[i, col[s]] = 1.0

        sig = sigma_rel(el_az(sat_xyz, theta[:3])[0], sysid) if have_pos else np.ones(len(y))
        w = keep / sig ** 2 if cfg["sig_w"] else keep.astype(float)    # prior trust per satellite
        if robust_on:
            w = w * lqlc_weight(y / sig if cfg["sig_s"] else y, cfg["s"])
        if (keep.sum() <= len(theta)) or any(w[sysid == s].sum() == 0 for s in systems):
            return None
        HW = H * w[:, None]
        dtheta = np.linalg.solve(H.T @ HW, HW.T @ y)                # (H^T W H)^-1 H^T W y
        theta = theta + dtheta
        if np.linalg.norm(dtheta[:3]) < TOL:
            if cfg["estimator"] == "LQLC" and not robust_on:
                robust_on = True                                    # LS first, then LQLC
            else:
                return result(theta, y, keep, it, True)
    return result(theta, y, keep, it, False)


def result(theta, y, keep, it, converged):
    r = y[keep]                                                     # post-fit residuals of used sats
    return dict(theta=theta, nsat=int(keep.sum()), iters=it + 1, converged=converged,
                res_rms=np.sqrt(np.mean(r ** 2)), res_max=np.abs(r).max())


# ---------------- run all configs ----------------
def run(cfg, data):
    d = data[data.sat_const.isin(cfg["systems"])]
    theta = np.zeros(3)
    out = []
    for ep, g in d.groupby("epoch"):
        r = solve_epoch(g, theta, cfg)
        if r is None: continue
        out.append(dict(epoch=ep, x=r["theta"][0], y=r["theta"][1], z=r["theta"][2],
                        nsat=r["nsat"], iters=r["iters"], conv=r["converged"],
                        res_rms=r["res_rms"], res_max=r["res_max"]))
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
               mean_h=h, res_rms=df.res_rms.mean(), res_max=df.res_max.mean())
    if ref is not None:
        e = (P - ref) @ R.T
        hz = np.hypot(e[:, 0], e[:, 1])
        row.update(hRMS=np.sqrt((hz ** 2).mean()), h95=np.percentile(hz, 95), vRMS=np.sqrt((e[:, 2] ** 2).mean()),
                   biasE=e[:, 0].mean(), biasN=e[:, 1].mean(), biasU=e[:, 2].mean())
    else:
        row.update(lat=np.degrees(lat), lon=np.degrees(lon))
    return row


if __name__ == "__main__":
    data = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "data.csv")
    if len(sys.argv) > 2:                 # re-run an existing run<N>/ from its own config.json
        out = Path(sys.argv[2])
        configs = json.load(open(out / "config.json"))
    else:
        n = 1 + max([int(p.name[3:]) for p in Path(".").glob("run*") if p.name[3:].isdigit()], default=0)
        out = Path(f"run{n}"); out.mkdir()
        configs = {name: {**BASE, **change} for name, change in CONFIGS.items()}
        json.dump(configs, open(out / "config.json", "w"), indent=1)
    rows = []
    for name, cfg in configs.items():
        df = run(cfg, data)
        rows.append(dict(config=name, **score(df, REF_XYZ)))
        df.to_csv(out / f"out_{name}.csv", index=False)
        print("done:", name)
    T = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(T.round(2).to_string(index=False))
    T.to_csv(out / "comparison.csv", index=False)
    print("saved to", out)
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_viewer.py")), str(out)], check=True)
