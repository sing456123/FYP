"""Build run<N>/viewer.html from the out_<config>.csv files listed in run<N>/comparison.csv (ECEF -> WGS-84 lat/lon/h). Stdlib only. Usage: python3 build_viewer.py [run<N>]  (default: newest run)"""
import csv, glob, json, math, os, sys

RUNS = sorted((p for p in glob.glob("run*") if p[3:].isdigit()), key=lambda p: int(p[3:]))
RUN = sys.argv[1] if len(sys.argv) > 1 else RUNS[-1]

COMP = list(csv.DictReader(open(f"{RUN}/comparison.csv")))
NAMES = [r["config"] for r in COMP]
A, E2 = 6378137.0, 6.69437999014e-3

def geodetic(x, y, z):
    lon = math.atan2(y, x); r = math.hypot(x, y); lat = math.atan2(z, r * (1 - E2))
    for _ in range(6):
        N = A / math.sqrt(1 - E2 * math.sin(lat) ** 2); h = r / math.cos(lat) - N
        lat = math.atan2(z, r * (1 - E2 * N / (N + h)))
    return math.degrees(lat), math.degrees(lon), h

cfgs, ep = [], None
for name, comp in zip(NAMES, COMP):
    rows = list(csv.DictReader(open(f"{RUN}/out_{name}.csv")))
    e = [r["epoch"] for r in rows]
    assert ep is None or e == ep, "epochs differ"
    ep = e
    g = [geodetic(float(r["x"]), float(r["y"]), float(r["z"])) for r in rows]
    cfgs.append(dict(name=name, lat=[round(p[0], 8) for p in g], lon=[round(p[1], 8) for p in g],
                     h=[round(p[2], 3) for p in g], nsat=[int(r["nsat"]) for r in rows],
                     conv=[int(r["conv"] == "True") for r in rows],
                     res=[round(float(r["res_rms"]), 2) for r in rows],
                     spread={k: round(float(comp[k]), 2)
                             for k in ("h_spread", "h95_spread", "v_spread", "res_rms", "res_max")}))

html = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates", "viewer_template.html")).read().replace("__DATA__", json.dumps(cfgs))
open(f"{RUN}/viewer.html", "w").write(html)
print(f"{RUN}/viewer.html", len(html) // 1024, "KB")
