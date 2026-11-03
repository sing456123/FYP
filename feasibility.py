"""
Feasibility check for the MOMENT proposal: does SNR behaviour predict bad pseudoranges?
Uses spp_switches.py as a library (not modified). Solves LQLC as usual, then records every used
satellite's post-fit residual and LQLC weight next to its SNR, and tests three SNR indicators:
  (a) low SNR
  (b) SNR below normal for its elevation (per system)
  (c) sudden SNR change vs the satellite's own recent history   <- what MOMENT looks for

Usage:  python feasibility.py data.csv     -> next free feas<N>/ with feasibility.txt, sats.csv.gz
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import spp_switches as spp

CFG = {**spp.BASE, "estimator": "LQLC", "mask_deg": 5.0, "s": 5.0}   # best horizontal spread in run4
BAD = 10.0        # |residual| above this (m) counts as a "bad" measurement
FLAG = 0.10       # each indicator flags its worst 10% of measurements
HIST = 30         # epochs of history for the "sudden change" indicator


def per_satellite(data, cfg):
    """Same epoch loop as spp.run(), but keeps every used satellite's residual and weight."""
    d = data[data.sat_const.isin(cfg["systems"])]
    theta, recs = np.zeros(3), []
    for ep, g in d.groupby("epoch"):
        r = spp.solve_epoch(g, theta, cfg)
        if r is None or not r["converged"]:
            continue
        th = theta = r["theta"]
        sysid = g.sat_const.values
        systems = [s for s in cfg["systems"] if (sysid == s).any()]
        sat = g[["sat_x", "sat_y", "sat_z"]].values
        corr, k = spp.iter_terms(sat, sysid, float(g.tow.iloc[0]), th, cfg, True)
        res = g.cor_pseudo.values - (np.linalg.norm(sat - th[:3], axis=1) + th[[3 + systems.index(s) for s in sysid]] + corr)
        el, _ = spp.el_az(sat, th[:3])
        recs.append(pd.DataFrame(dict(epoch=ep, sys=sysid[k], prn=g.sat_id.values[k], el=np.degrees(el[k]),
                                      snr=g.sat_snr.values[k], res=res[k], w=spp.lqlc_weight(res[k], cfg["s"]))))
    return pd.concat(recs, ignore_index=True)


def indicators(S):
    S = S.sort_values(["sys", "prn", "epoch"]).reset_index(drop=True)
    S["abs_res"] = S.res.abs()
    S["bad"] = S.abs_res > BAD
    S["el_bin"] = (S.el // 5) * 5
    S["snr_deficit"] = S.groupby(["sys", "el_bin"]).snr.transform("median") - S.snr        # (b)
    past = S.groupby(["sys", "prn"]).snr.transform(lambda x: x.shift(1).rolling(HIST, min_periods=10).median())
    S["snr_jump"] = (S.snr - past).abs()                                                    # (c), causal
    return S


def evaluate(S):
    """For each indicator: rank correlation with |res|, and how well flagging its worst 10% catches bad ones."""
    base = S.bad.mean()
    rows = []
    for name, col, worst_high in [("(a) low SNR", "snr", False), ("(b) SNR below normal for elevation", "snr_deficit", True),
                                  ("(c) sudden SNR change", "snr_jump", True)]:
        T = S.dropna(subset=[col])
        score = T[col] if worst_high else -T[col]
        rank = score.sample(frac=1, random_state=0).rank(method="first").reindex(score.index)
        flag = rank > (1 - FLAG) * len(score)                      # exactly FLAG, ties broken randomly
        prec = T.bad[flag].mean()
        rows.append(dict(indicator=name, n=len(T), spearman=score.rank().corr(T.abs_res.rank()),
                         flagged_bad_pct=100 * prec, caught_pct=100 * T.bad[flag].sum() / max(T.bad.sum(), 1),
                         lift=prec / base, w_flagged=T.w[flag].mean(), w_others=T.w[~flag].mean()))
    return pd.DataFrame(rows)


def binned(S, col, edges, label):
    b = pd.cut(S[col], edges)
    t = S.groupby(b, observed=True).agg(n=("abs_res", "size"), median_res=("abs_res", "median"),
                                        p90_res=("abs_res", lambda x: x.quantile(.9)), bad_pct=("bad", "mean"),
                                        mean_w=("w", "mean"))
    t["bad_pct"] *= 100
    return f"{label}\n{t.round(2).to_string()}\n"


def verdict(lift):
    return ("strong link" if lift >= 2 else "weak link" if lift >= 1.3 else "no useful link")


if __name__ == "__main__":
    data = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "data.csv")
    n = 1 + max([int(p.name[4:]) for p in Path(".").glob("feas*") if p.name[4:].isdigit()], default=0)
    out = Path(f"feas{n}"); out.mkdir()
    json.dump(dict(cfg=CFG, BAD=BAD, FLAG=FLAG, HIST=HIST), open(out / "config.json", "w"), indent=1)

    S = indicators(per_satellite(data, CFG))
    S.to_csv(out / "sats.csv.gz", index=False)
    E = evaluate(S)

    lines = [f"Measurements: {len(S)} (satellite x epoch) | bad = |residual| > {BAD} m: {100 * S.bad.mean():.1f}% "
             f"| median |residual| {S.abs_res.median():.2f} m\n",
             f"Each indicator flags its worst {100 * FLAG:.0f}%. lift = how much more often flagged ones are bad than "
             f"average (1 = no better than random).\n",
             E.round(2).to_string(index=False) + "\n",
             "Verdict per indicator: " + "; ".join(f"{r.indicator}: {verdict(r.lift)} (lift {r.lift:.2f})"
                                                  for r in E.itertuples()) + "\n",
             binned(S, "snr", [0, 25, 30, 35, 40, 45, 60], "|residual| by SNR (dB-Hz)"),
             binned(S, "snr_deficit", [-60, -3, 0, 3, 6, 10, 60], "|residual| by SNR deficit vs normal for elevation (dB)"),
             binned(S, "snr_jump", [-0.1, 1, 2, 4, 6, 60], "|residual| by SNR change vs own last 30 epochs (dB)"),
             binned(S, "el", [0, 15, 30, 45, 60, 90], "|residual| by elevation (deg), for context")]
    text = "\n".join(lines)
    print(text)
    (out / "feasibility.txt").write_text(text)
    print("saved to", out)
