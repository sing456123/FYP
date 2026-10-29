"""Build run<N>/sweep.html from a sweep run's comparison.csv. Stdlib only.
Usage: python3 build_sweep.py [run<N>]   (default: newest run that has sweep.txt)"""
import csv, glob, json, os, sys

METRICS = ["h_spread", "h95_spread", "v_spread", "res_rms", "res_max", "mean_h"]
RUNS = sorted((p for p in glob.glob("run*") if p[3:].isdigit() and os.path.exists(f"{p}/sweep.txt")),
              key=lambda p: int(p[3:]))
RUN = sys.argv[1] if len(sys.argv) > 1 else RUNS[-1]

comp = list(csv.DictReader(open(f"{RUN}/comparison.csv")))
rows = [dict(config=r["config"], estimator=r["estimator"], mask_deg=float(r["mask_deg"]),
             s=float(r["s"]) if r["s"] else None, **{k: round(float(r[k]), 3) for k in METRICS}) for r in comp]
cfg = next(iter(json.load(open(f"{RUN}/config.json")).values()))
data = dict(run=RUN, systems=cfg["systems"], epochs=int(comp[0]["epochs"]), rows=rows)

html = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep_template.html")).read()
open(f"{RUN}/sweep.html", "w").write(html.replace("__DATA__", json.dumps(data)))
print(f"{RUN}/sweep.html")
