"""Convert a sppStruct .mat (one struct per epoch) to the CSV layout spp_switches.py reads.
Usage:  python mat2csv.py "sppStruct_user copy.mat" data_midurban.csv"""
import sys
import numpy as np
import pandas as pd
import scipy.io

COLS = ["sat_id", "sat_const", "sat_x", "sat_y", "sat_z", "sat_el", "sat_az", "sat_snr", "raw_pseudo", "cor_pseudo"]

src, dst = sys.argv[1], sys.argv[2]
epochs = scipy.io.loadmat(src, squeeze_me=True, struct_as_record=False)["sppStruct"]
parts = []
for i, e in enumerate(np.atleast_1d(epochs)):
    cols = {k: np.atleast_1d(getattr(e, k)) for k in COLS}      # atleast_1d: epochs with one satellite
    if cols["sat_id"].size == 0:
        continue
    parts.append(pd.DataFrame(dict(epoch=i, week=e.week, tow=e.tow, **cols)))
df = pd.concat(parts, ignore_index=True)
df.to_csv(dst, index=False)
print(f"{dst}: {df.epoch.nunique()} epochs, {len(df)} rows")
