"""
Sweep elevation mask x s for LQLC_sig_w_s only (idea 8+9: weight 1/sigma_i^2 and per-satellite scale s*sigma_i).
Uses spp_noise.py as a library. One LS run per mask is added as the reference line the sweep page needs.
The reference point for accuracy is spp_noise.REF_LLA.

Usage:  python sweep_noise.py data.csv   -> next free run<N>/ with comparison.csv, sweep.txt, sweep.html, out_*.csv
"""
import json
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
import spp_switches as spp
import spp_noise as nz

MASKS = [5, 10, 15, 20]                # elevation mask (deg)
SS = [0.3, 0.5, 1, 2, 3, 5]            # s (m) for a GPS satellite at zenith; satellite i uses s * sigma_i
METRICS = ["h_spread", "h95_spread", "v_spread", "res_rms", "mean_h"]
ACC = ["hRMS", "h95", "vRMS", "biasU"]  # only when REF_LLA is set


def configs():
    out = {}
    for m in MASKS:
        out[f"LS_m{m}"] = {**nz.BASE, "estimator": "LS", "mask_deg": m}
        for s in SS:
            out[f"LQLC_sig_w_s_m{m}_s{s}"] = {**nz.BASE, "estimator": "LQLC", "sig_w": True, "sig_s": True,
                                              "mask_deg": m, "s": s}
    return out


def load(path):                        # runs once in each worker process
    global DATA
    DATA = pd.read_csv(path)


def job(item):
    name, cfg = item
    df = nz.run(cfg, DATA)
    return name, df, spp.score(df, nz.ref_xyz())


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data.csv"
    n = 1 + max([int(p.name[3:]) for p in Path(".").glob("run*") if p.name[3:].isdigit()], default=0)
    out = Path(f"run{n}"); out.mkdir()
    cfgs = configs()
    json.dump(cfgs, open(out / "config.json", "w"), indent=1)

    rows = []
    with Pool(initializer=load, initargs=(path,)) as pool:
        for name, df, score in pool.imap_unordered(job, cfgs.items()):
            df.to_csv(out / f"out_{name}.csv", index=False)
            cfg = cfgs[name]
            rows.append(dict(config=name, estimator=cfg["estimator"], mask_deg=cfg["mask_deg"],
                             s=cfg["s"] if cfg["estimator"] == "LQLC" else None, **score))
            print("done:", name, flush=True)

    T = pd.DataFrame(rows).set_index("config").loc[list(cfgs)].reset_index()
    T.to_csv(out / "comparison.csv", index=False)

    report = [f"LQLC_sig_w_s sweep on {path} | reference: {nz.REF_LLA}\n"]
    for k in METRICS + (ACC if nz.REF_LLA is not None else []):
        P = T[T.estimator == "LQLC"].pivot(index="mask_deg", columns="s", values=k)
        P.insert(0, "LS", T[T.estimator == "LS"].set_index("mask_deg")[k])
        report.append(f"{k}  (rows: mask deg | columns: LS, then LQLC_sig_w_s by s in m)\n{P.round(2).to_string()}\n")
    text = "\n".join(report)
    print(text)
    (out / "sweep.txt").write_text(text)
    print("saved to", out)
    subprocess.run([sys.executable, str(Path(__file__).with_name("build_sweep.py")), str(out)], check=True)
