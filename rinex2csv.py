"""RINEX 3 NAV + OBS -> the data.csv layout read by spp_switches.py.
Usage:  python rinex2csv.py SAISA_..._NAV.rnx SAISA_..._OBS.rnx data_saisa.csv

Only GPS (sat_const 1) and Galileo (8) are written: this NAV file carries no
BeiDou ephemerides, and spp_switches.py has no QZSS/SBAS system codes.
sat_id = PRN (GPS) or PRN + 60 (Galileo), matching data.csv.
cor_pseudo = raw_pseudo + c*dts, dts = broadcast clock polynomial at signal
transmit time + relativistic term.  Sagnac/tropo/iono are left to spp_switches.py.
"""
import sys
import csv
import math
from datetime import datetime, timedelta
import numpy as np
from spp_switches import ecef_to_lla, enu_matrix

C = 299792458.0
OMEGA_E = 7.2921151467e-5
MU = {"G": 3.986005e14, "E": 3.986004418e14}
SAT_CONST = {"G": 1, "E": 8}
ID_OFFSET = {"G": 0, "E": 60}
PREF = {"G": ["C1C"], "E": ["C1X", "C5X"]}
GPS_EPOCH = datetime(1980, 1, 6)


def parse_nav(path):
    """{(sys,prn): [eph, ...]}; times as GPS seconds of week (absolute)."""
    lines = open(path).read().splitlines()
    i = next(k for k, l in enumerate(lines) if "END OF HEADER" in l)
    out, buf = {}, []

    def flush(b):
        if len(b) < 8 or not b[0][1:3].strip().isdigit() or b[0][0] not in MU:
            return
        sysid, prn = b[0][0], int(b[0][1:3])
        y, mo, d, hh, mi, ss = map(int, b[0][4:23].split())
        f = lambda ln: [float(ln[c:c + 19]) if ln[c:c + 19].strip() else math.nan
                        for c in (4, 23, 42, 61)]
        v1 = [float(b[0][c:c + 19]) for c in (23, 42, 61)]
        v2, v3, v4, v5, v6 = (f(x) for x in b[1:6])
        week = int(v6[2])
        if sysid == "E" and week < 2000:      # some writers store GST week (= GPS week - 1024)
            week += 1024
        toc_s = hh * 3600 + mi * 60 + ss
        out.setdefault((sysid, prn), []).append(dict(
            sys=sysid, toc=week * 604800.0 + toc_s, toe=week * 604800.0 + v4[0],
            af0=v1[0], af1=v1[1], af2=v1[2], e=v3[1], sqrtA=v3[3],
            crs=v2[1], dn=v2[2], m0=v2[3], cuc=v3[0], cus=v3[2], toe_s=v4[0],
            cic=v4[1], om0=v4[2], cis=v4[3], i0=v5[0], crc=v5[1], om=v5[2],
            omdot=v5[3], idot=v6[0]))

    for l in lines[i + 1:]:
        if not l.strip():
            continue
        if l[0] in MU and l[1:3].strip().isdigit():
            if buf:
                flush(buf)
            buf = [l]
        elif buf:
            buf.append(l)
            if len(buf) == 8:
                flush(buf)
                buf = []
    if buf:
        flush(buf)
    return out


def parse_obs(path):
    """[(t_gps, [(sys, prn, [obs tokens]), ...]), ...], obs types, approx xyz."""
    lines = open(path, encoding="utf-8").read().splitlines()
    i = next(k for k, l in enumerate(lines) if "END OF HEADER" in l)
    types, approx = {}, None
    for l in lines[:i]:
        if "SYS / # / OBS TYPES" in l:
            types.setdefault(l[0], []).extend(l[7:60].split())
        elif "APPROX POSITION XYZ" in l:
            approx = np.array([float(l[0:14]), float(l[14:28]), float(l[28:42])])
    out, t, sats = [], None, []
    for l in lines[i + 1:]:
        if l.startswith(">"):
            if sats:
                out.append((t, sats))
            y, mo, d, hh, mi, ss = l[1:].split()[:6]
            t = (datetime(int(y), int(mo), int(d), int(hh), int(mi))
                 + timedelta(seconds=float(ss)) - GPS_EPOCH).total_seconds()
            sats = []
        elif l.startswith("→"):          # PYRINEXCONV continuation: "→N" then 3 values
            sats[-1][2].extend(l[2:].split())
        elif l[0] in types and l[1:3].strip().isdigit():
            sats.append([l[0], int(l[1:3]), l[3:].split()])
    if sats:
        out.append((t, sats))
    return out, types, approx


def sat_pos(e, t):
    A = e["sqrtA"] ** 2
    tk = (t - e["toe"] + 302400.0) % 604800.0 - 302400.0
    n = math.sqrt(MU[e["sys"]] / A ** 3) + e["dn"]
    M = e["m0"] + n * tk
    E = M
    for _ in range(10):
        E = M + e["e"] * math.sin(E)
    v = math.atan2(math.sqrt(1 - e["e"] ** 2) * math.sin(E), math.cos(E) - e["e"])
    phi = v + e["om"]
    du = e["cus"] * math.sin(2 * phi) + e["cuc"] * math.cos(2 * phi)
    dr = e["crs"] * math.sin(2 * phi) + e["crc"] * math.cos(2 * phi)
    di = e["cis"] * math.sin(2 * phi) + e["cic"] * math.cos(2 * phi)
    u, r, inc = phi + du, A * (1 - e["e"] * math.cos(E)) + dr, e["i0"] + di + e["idot"] * tk
    xp, yp = r * math.cos(u), r * math.sin(u)
    om = e["om0"] + (e["omdot"] - OMEGA_E) * tk - OMEGA_E * e["toe_s"]
    return np.array([xp * math.cos(om) - yp * math.cos(inc) * math.sin(om),
                     xp * math.sin(om) + yp * math.cos(inc) * math.cos(om),
                     yp * math.sin(inc)]), E


def sat_clock(e, t, E):
    dt = t - e["toc"]
    return (e["af0"] + e["af1"] * dt + e["af2"] * dt * dt
            - 2 * math.sqrt(MU[e["sys"]]) / C ** 2 * e["e"] * e["sqrtA"] * math.sin(E))


def usable(e, t):
    A = e["sqrtA"] ** 2
    return 0.0 <= e["e"] < 0.1 and 2.0e7 < A < 5.5e7 and abs(t - e["toe"]) < 4 * 3600


def pick_obs(sysid, names, toks):
    for name in PREF[sysid] + [n for n in names if n.startswith("C")]:
        if name in names:
            i = names.index(name)
            if i < len(toks) and float(toks[i]) != 0:
                s = "S" + name[1:]
                j = names.index(s) if s in names else -1
                snr = float(toks[j]) if 0 <= j < len(toks) else math.nan
                return float(toks[i]), snr, name
    return None, None, None


def main(nav_path, obs_path, out_path):
    nav, (obs, types, approx) = parse_nav(nav_path), parse_obs(obs_path)
    lat, lon, _ = ecef_to_lla(approx)
    enu_mat = enu_matrix(lat, lon)
    rows, check, epoch = [], [], 0
    for t, sats in obs:
        week, tow = int(t // 604800), t % 604800
        ep = []
        for sysid, prn, toks in sats:
            cands = [e for e in nav.get((sysid, prn), []) if usable(e, t)]
            if not cands:
                continue
            e = min(cands, key=lambda e: abs(t - e["toe"]))
            rho, snr, _ = pick_obs(sysid, types[sysid], toks)
            if rho is None:
                continue
            t0 = t - rho / C
            p0, E0 = sat_pos(e, t0)
            ttx = t0 - sat_clock(e, t0, E0)
            pos, E = sat_pos(e, ttx)
            cor = rho + C * sat_clock(e, ttx, E)
            enu = enu_mat @ (pos - approx)
            el = math.degrees(math.atan2(enu[2], math.hypot(enu[0], enu[1])))
            az = math.degrees(math.atan2(enu[0], enu[1])) % 360.0
            ep.append([prn + ID_OFFSET[sysid], SAT_CONST[sysid], pos[0], pos[1], pos[2],
                       el, az, snr, rho, cor])
            if epoch < 20:
                check.append((sysid, cor - np.linalg.norm(pos - approx)))
        if ep:
            for r in sorted(ep):
                rows.append([epoch, week, f"{tow:.3f}"] + [r[0], r[1],
                            f"{r[2]:.3f}", f"{r[3]:.3f}", f"{r[4]:.3f}",
                            f"{r[5]:.6f}", f"{r[6]:.6f}", f"{r[7]:.3f}",
                            f"{r[8]:.3f}", f"{r[9]:.3f}"])
            epoch += 1

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["epoch", "week", "tow", "sat_id", "sat_const", "sat_x", "sat_y",
                    "sat_z", "sat_el", "sat_az", "sat_snr", "raw_pseudo", "cor_pseudo"])
        w.writerows(rows)

    # one runnable check: clock sign / parse errors show up as >100 km per-system scatter
    worst = 0.0
    for s in SAT_CONST:
        v = np.array([x for sy, x in check if sy == s])
        if len(v):
            worst = max(worst, np.median(np.abs(v - np.median(v))))
    print(f"{out_path}: {epoch} epochs, {len(rows)} rows, "
          f"{len({r[3] for r in rows})} sats, check scatter {worst / 1e3:.1f} km")
    assert worst < 50e3, "cor_pseudo inconsistent with satellite ranges (clock sign? parse?)"


if __name__ == "__main__":
    main(*sys.argv[1:4])
