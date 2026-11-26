"""
RTKLIB elevation/system noise model with LS and LQLC   (ideas 8, 9 and 8+9 in README "What we tried")
=====================================================
Uses spp_switches.py as a library (same corrections, mask, scoring); only the weight changes.

  sigma_i = F_sys * 100 * sqrt(0.003^2 + 0.003^2 / sin el_i)    RTKLIB varerr() defaults, F = 1.5 GLONASS, 1 others
            divided by GPS at zenith, so sigma_i = 1 for the cleanest satellite

Two switches:
  sig_w : prior trust 1/sigma_i^2                    (idea 8, RTKLIB-style)
  sig_s : LQLC judges e_i/sigma_i, i.e. scale s*sigma_i (idea 9, per-satellite s)

  W_ii = m_i * [1/sigma_i^2 if sig_w] * [tanh(x_i)/x_i if LQLC],   x_i = |e_i| / (2 s [sigma_i if sig_s])

Usage:  python spp_noise.py data.csv     -> next free run<N>/ with comparison.csv and viewer.html
"""
import sys
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import spp_switches as spp

EFACT = {1: 1.0, 4: 1.5, 8: 1.0, 32: 1.0}                  # RTKLIB EFACT_GPS/GLO/GAL/CMP

MASK_DEG = 15.0                                             # elevation mask (deg)
S = 1                                                     # LQLC scale s (m); with sig_s, satellite i uses S * sigma_i

REF_LLA = (22.313095, 114.19142326, 3.0)   # true position (lat deg, lon deg, h m) for scoring; None = precision only
                                           # (this is the urban data.csv site: set None or the right point for other data)

BASE = {**spp.BASE, "mask_deg": MASK_DEG, "s": S, "sig_w": False, "sig_s": False}
CONFIGS = {                                                 # name: changes from BASE
    "LS":            dict(estimator="LS"),                  # baseline
    "LQLC":          dict(estimator="LQLC"),                # baseline (idea 1)
    "LS_sig_w":      dict(estimator="LS", sig_w=True),      # idea 8: RTKLIB weighted LS
    "LQLC_sig_s":    dict(estimator="LQLC", sig_s=True),    # idea 9: per-satellite s
    "LQLC_sig_w_s":  dict(estimator="LQLC", sig_w=True, sig_s=True),   # idea 8+9: full fusion
}


def sigma_rel(el, sysid):
    """RTKLIB noise model, relative to GPS at zenith. el in rad, clamped at 5 deg."""
    sin_el = np.sin(np.maximum(el, np.radians(5.0)))
    sig = np.vectorize(EFACT.get)(sysid) * 100 * np.sqrt(0.003 ** 2 + 0.003 ** 2 / sin_el)
    return sig / (100 * np.sqrt(2 * 0.003 ** 2))


def solve_epoch(g, theta0, cfg):
    """spp.solve_epoch with sigma_i in the weight; everything else identical."""
    sat_xyz = g[["sat_x", "sat_y", "sat_z"]].values
    sysid = g.sat_const.values
    tow = float(g.tow.iloc[0])
    rho = g.cor_pseudo.values.astype(float)
    systems = [s for s in cfg["systems"] if (sysid == s).any()]
    col = {s: 3 + k for k, s in enumerate(systems)}

    theta = np.r_[theta0[:3], np.zeros(len(systems))]
    if len(theta0) == len(theta): theta = theta0.copy()
    robust_on = False
    for it in range(2 * spp.MAX_ITER):
        have_pos = np.linalg.norm(theta[:3]) > 6.0e6
        corr, keep = spp.iter_terms(sat_xyz, sysid, tow, theta, cfg, have_pos)
        vec = sat_xyz - theta[:3]
        rng = np.linalg.norm(vec, axis=1)
        clk = np.array([theta[col[s]] for s in sysid])
        y = rho - (rng + clk + corr)
        H = np.zeros((len(y), len(theta)))
        H[:, :3] = -vec / rng[:, None]
        for i, s in enumerate(sysid): H[i, col[s]] = 1.0

        sig = sigma_rel(spp.el_az(sat_xyz, theta[:3])[0], sysid) if have_pos else np.ones(len(y))
        w = keep / sig ** 2 if cfg["sig_w"] else keep.astype(float)        # m_i [/ sigma_i^2]
        if robust_on:
            w = w * spp.lqlc_weight(y / sig if cfg["sig_s"] else y, cfg["s"])   # x_i = |e_i| / 2s[sigma_i]
        if (keep.sum() <= len(theta)) or any(w[sysid == s].sum() == 0 for s in systems):
            return None
        HW = H * w[:, None]
        dtheta = np.linalg.solve(H.T @ HW, HW.T @ y)                        # (H^T W H)^-1 H^T W y
        theta = theta + dtheta
        if np.linalg.norm(dtheta[:3]) < spp.TOL:
            if cfg["estimator"] == "LQLC" and not robust_on:
                robust_on = True                                            # LS first, then LQLC
            else:
                return spp.result(theta, y, keep, it, True)
    return spp.result(theta, y, keep, it, False)


def ref_xyz():
    return spp.lla_to_ecef(*REF_LLA) if REF_LLA is not None else None


def run(cfg, data):
    """spp.run with this file's solve_epoch."""
    d = data[data.sat_const.isin(cfg["systems"])]
    theta, out = np.zeros(3), []
    for ep, g in d.groupby("epoch"):
        r = solve_epoch(g, theta, cfg)
        if r is None: continue
        out.append(dict(epoch=ep, x=r["theta"][0], y=r["theta"][1], z=r["theta"][2],
                        nsat=r["nsat"], iters=r["iters"], conv=r["converged"],
                        res_rms=r["res_rms"], res_max=r["res_max"]))
        if r["converged"] and np.linalg.norm(r["theta"][:3]) > 6.0e6:
            theta = r["theta"]
    return pd.DataFrame(out)


if __name__ == "__main__":
    data = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "data.csv")
    n = 1 + max([int(p.name[3:]) for p in Path(".").glob("run*") if p.name[3:].isdigit()], default=0)
    out = Path(f"run{n}"); out.mkdir()
    configs = {name: {**BASE, **change} for name, change in CONFIGS.items()}
    json.dump(configs, open(out / "config.json", "w"), indent=1)
    rows = []
    for name, cfg in configs.items():
        df = run(cfg, data)
        rows.append(dict(config=name, **spp.score(df, ref_xyz())))
        df.to_csv(out / f"out_{name}.csv", index=False)
        print("done:", name, flush=True)
    T = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(T.round(2).to_string(index=False))
    T.to_csv(out / "comparison.csv", index=False)
    print("saved to", out)
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_viewer.py")), str(out)], check=True)
