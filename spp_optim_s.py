"""
LQLC with s chosen per epoch by minimum predicted variance   (no true position needed)
===========================================================
Uses spp_switches.py as a library (same corrections, mask, scoring). Per epoch, for every s in S_GRID:
  solve LQLC with that s -> residuals e_i of the used satellites -> predicted variance of the position

  V(s) = mean(psi(e)^2) / mean(psi'(e))^2 * trace( (H^T H)^-1 )_xyz        (M-estimator asymptotic covariance)
  LQLC (cost 4s^2 ln cosh(e/2s)):  psi(e) = 2s tanh(e/2s),  psi'(e) = sech^2(e/2s)

and keep the solution of the s with the smallest V(s).
Limits: s -> large gives the LS variance mean(e^2); s -> small gives the L1 (median) variance.
Assumes bad satellites are random outliers; several reflected signals that agree can fool it.

Usage:  python spp_optim_s.py data.csv     -> next free run<N>/ with comparison.csv and viewer.html
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
S_GRID = [0.5, 1, 2, 3, 5, 10]       # candidate s values (m) tried every epoch
FIXED_S = [1.0, 5.0]                 # fixed-s LQLC baselines (best on mid-urban / urban)
REFS = {                             # true position per data file (lat deg, lon deg, h m); missing = precision only
    "data.csv": (22.313095, 114.19142326, 3.0),
}


def configs():
    base = {**spp.BASE, "mask_deg": MASK_DEG, "optim": False}
    out = {"LS": {**base, "estimator": "LS"}}
    for s in FIXED_S:
        out[f"LQLC_s{s:g}"] = {**base, "estimator": "LQLC", "s": s}
    out["LQLC_optS"] = {**base, "estimator": "LQLC", "optim": True}
    return out


def predicted_var(g, theta, cfg, s):
    """V(s) at the solution theta: robust variance factor x geometry (position part of (H^T H)^-1)."""
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    sysid = g.sat_const.values
    systems = [c for c in cfg["systems"] if (sysid == c).any()]
    col = {c: 3 + k for k, c in enumerate(systems)}
    corr, keep = spp.iter_terms(sat_xyz, sysid, float(g.tow.iloc[0]), theta, cfg, True)
    vec = sat_xyz - theta[:3]
    rng = np.linalg.norm(vec, axis=1)
    e = (g.cor_pseudo.values - (rng + theta[[col[c] for c in sysid]] + corr))[keep]
    H = np.zeros((len(sysid), len(theta)))
    H[:, :3] = -vec / rng[:, None]
    for i, c in enumerate(sysid): H[i, col[c]] = 1.0
    H = H[keep]
    u = e / (2 * s)
    psi, dpsi = 2 * s * np.tanh(u), 1 / np.cosh(u) ** 2
    factor = np.mean(psi ** 2) / np.mean(dpsi) ** 2
    return factor * np.trace(np.linalg.inv(H.T @ H)[:3, :3])


def solve_epoch(g, theta0, cfg):
    """Fixed s: spp.solve_epoch. Optimised: solve for every s in S_GRID and keep the smallest V(s)."""
    if not cfg["optim"]:
        r = spp.solve_epoch(g, theta0, cfg)
        return r and {**r, "s": cfg["s"]}
    best = None
    for s in S_GRID:
        r = spp.solve_epoch(g, theta0, {**cfg, "s": s})
        if r is None or not r["converged"] or np.linalg.norm(r["theta"][:3]) < 6.0e6:
            continue
        v = predicted_var(g, r["theta"], cfg, s)
        if best is None or v < best["V"]:
            best = {**r, "s": s, "V": v}
    return best


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
    if cfg["optim"]:
        share = df.s.value_counts(normalize=True)
        sc.update(s_median=df.s.median(), **{f"pct_s{s:g}": 100 * share.get(s, 0.0) for s in S_GRID})
    return name, df, sc


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data.csv"
    n = 1 + max([int(p.name[3:]) for p in Path(".").glob("run*") if p.name[3:].isdigit()], default=0)
    out = Path(f"run{n}"); out.mkdir()
    cfgs = configs()
    json.dump(cfgs, open(out / "config.json", "w"), indent=1)
    json.dump(dict(data=path, ref=REFS.get(Path(path).name), S_GRID=S_GRID), open(out / "settings.json", "w"), indent=1)

    rows = {}
    with Pool(initializer=load, initargs=(path,)) as pool:
        for name, df, sc in pool.imap_unordered(job, cfgs.items()):
            df.to_csv(out / f"out_{name}.csv", index=False)
            rows[name] = dict(config=name, **sc)
            print("done:", name, flush=True)
    T = pd.DataFrame([rows[k] for k in cfgs])
    T.to_csv(out / "comparison.csv", index=False)
    cols = [c for c in ["config", "h_spread", "h95_spread", "v_spread", "hRMS", "h95", "vRMS", "biasU",
                        "res_rms", "conv_pct", "s_median"] + [f"pct_s{s:g}" for s in S_GRID] if c in T]
    pd.set_option("display.width", 240)
    print(T[cols].round(2).to_string(index=False))
    print("saved to", out)
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_viewer.py")), str(out)], check=True)
