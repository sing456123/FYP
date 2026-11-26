"""
LQLC with an adaptive scale s   (s computed each epoch from the residuals, instead of fixed)
==============================
Uses spp_switches.py as a library (same corrections, mask, scoring). Per epoch:

  LS phase -> residuals e of the used satellites
           -> s = max(S_MIN, k * 1.4826 * MAD(e))      MAD = median(|e - median(e)|), ignores up to half bad sats
           -> LQLC phase with that s
  S_MODE = "once": s stays fixed for the rest of the epoch
  S_MODE = "iter": s is recomputed from the LQLC residuals at every LQLC iteration

1.4826 * MAD is a robust standard deviation: it measures how noisy the GOOD satellites are this epoch.
Only k is tuned (unit-free), so one k may work in every environment, unlike one fixed s.
Fixed-s LQLC runs are included as baselines.

Usage:  python spp_adaptive.py data.csv     -> next free run<N>/ with comparison.csv and viewer.html
"""
import sys
import json
import subprocess
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
import spp_switches as spp

MASK_DEG = 10.0                      # elevation mask (deg), same for every config
KS = [0.4, 0.5, 0.6, 0.75, 0.9]      # k values to try for the adaptive s
FIXED_S = [1.0, 5.0]                 # fixed-s LQLC baselines (best on mid-urban / urban)
S_MIN = 0.3                          # lower limit on the adaptive s (m)
S_MODE = "iter"                      # "once": s from the LS residuals only | "iter": recomputed every LQLC iteration
REFS = {                             # true position per data file (lat deg, lon deg, h m); missing = precision only
    "data.csv": (22.313095, 114.19142326, 3.0),
}


def configs():
    base = {**spp.BASE, "mask_deg": MASK_DEG, "k": None}
    out = {"LS": {**base, "estimator": "LS"}}
    for s in FIXED_S:
        out[f"LQLC_s{s:g}"] = {**base, "estimator": "LQLC", "s": s}
    for k in KS:
        out[f"LQLC_ad{S_MODE}_k{k:g}"] = {**base, "estimator": "LQLC", "k": k, "s_mode": S_MODE}
    return out


def robust_std(e):
    return 1.4826 * np.median(np.abs(e - np.median(e)))


def solve_epoch(g, theta0, cfg):
    """spp.solve_epoch, but when LQLC switches on, s is set from this epoch's residuals (if cfg['k'] is set)."""
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    sysid = g.sat_const.values
    tow = float(g.tow.iloc[0])
    rho = g.cor_pseudo.values.astype(float)
    systems = [s for s in cfg["systems"] if (sysid == s).any()]
    col = {s: 3 + k for k, s in enumerate(systems)}

    theta = np.r_[theta0[:3], np.zeros(len(systems))]
    if len(theta0) == len(theta): theta = theta0.copy()
    robust_on, s = False, cfg["s"]
    for it in range(2 * spp.MAX_ITER):
        have_pos = np.linalg.norm(theta[:3]) > 6.0e6
        corr, keep = spp.iter_terms(sat_xyz, sysid, tow, theta, cfg, have_pos)
        vec = sat_xyz - theta[:3]
        rng = np.linalg.norm(vec, axis=1)
        clk = np.array([theta[col[c]] for c in sysid])
        y = rho - (rng + clk + corr)
        H = np.zeros((len(y), len(theta)))
        H[:, :3] = -vec / rng[:, None]
        for i, c in enumerate(sysid): H[i, col[c]] = 1.0

        w = keep.astype(float)
        if robust_on:
            if cfg["k"] is not None and cfg["s_mode"] == "iter":   # s from the current (LQLC) residuals
                s = max(S_MIN, cfg["k"] * robust_std(y[keep]))
            w = w * spp.lqlc_weight(y, s)
        if (keep.sum() <= len(theta)) or any(w[sysid == c].sum() == 0 for c in systems):
            return None
        HW = H * w[:, None]
        dtheta = np.linalg.solve(H.T @ HW, HW.T @ y)                # (H^T W H)^-1 H^T W y
        theta = theta + dtheta
        if np.linalg.norm(dtheta[:3]) < spp.TOL:
            if cfg["estimator"] == "LQLC" and not robust_on:
                robust_on = True                                    # LS first, then LQLC
                if cfg["k"] is not None:                            # adaptive: s from the LS residuals
                    s = max(S_MIN, cfg["k"] * robust_std(y[keep]))
            else:
                return {**spp.result(theta, y, keep, it, True), "s": s}
    return {**spp.result(theta, y, keep, it, False), "s": s}


def run(cfg, data):
    """spp.run with this file's solve_epoch; also records the s used each epoch."""
    d = data[data.sat_const.isin(cfg["systems"])]
    theta, out = np.zeros(3), []
    for ep, g in d.groupby("epoch"):
        r = solve_epoch(g, theta, cfg)
        if r is None: continue
        out.append(dict(epoch=ep, x=r["theta"][0], y=r["theta"][1], z=r["theta"][2],
                        nsat=r["nsat"], iters=r["iters"], conv=r["converged"],
                        res_rms=r["res_rms"], res_max=r["res_max"], s=r["s"]))
        if r["converged"] and np.linalg.norm(r["theta"][:3]) > 6.0e6:
            theta = r["theta"]
    return pd.DataFrame(out)


def load(path):                      # runs once in each worker process
    global DATA, REF
    DATA = pd.read_csv(path)
    lla = REFS.get(Path(path).name)
    REF = spp.lla_to_ecef(*lla) if lla else None


def job(item):
    name, cfg = item
    df = run(cfg, DATA)
    sc = spp.score(df, REF)
    if cfg["estimator"] == "LQLC":
        sc.update(s_median=df.s.median(), s_p10=df.s.quantile(.1), s_p90=df.s.quantile(.9))
    return name, df, sc


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data.csv"
    n = 1 + max([int(p.name[3:]) for p in Path(".").glob("run*") if p.name[3:].isdigit()], default=0)
    out = Path(f"run{n}"); out.mkdir()
    cfgs = configs()
    json.dump(dict(data=path, ref=REFS.get(Path(path).name), S_MIN=S_MIN, configs=cfgs),
              open(out / "settings.json", "w"), indent=1)
    json.dump(cfgs, open(out / "config.json", "w"), indent=1)

    rows = {}
    with Pool(initializer=load, initargs=(path,)) as pool:
        for name, df, sc in pool.imap_unordered(job, cfgs.items()):
            df.to_csv(out / f"out_{name}.csv", index=False)
            rows[name] = dict(config=name, **sc)
            print("done:", name, flush=True)
    T = pd.DataFrame([rows[k] for k in cfgs])
    T.to_csv(out / "comparison.csv", index=False)
    cols = [c for c in ["config", "h_spread", "h95_spread", "v_spread", "hRMS", "h95", "vRMS", "biasU",
                        "res_rms", "conv_pct", "s_median", "s_p10", "s_p90"] if c in T]
    pd.set_option("display.width", 220)
    print(T[cols].round(2).to_string(index=False))
    print("saved to", out)
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_viewer.py")), str(out)], check=True)
