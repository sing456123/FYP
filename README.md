# SPP: LS vs LQLC

`spp_switches.py` computes the receiver position every epoch from GNSS pseudoranges (Single Point Positioning).
It solves every epoch with both **LS** and **LQLC** on the receiver's `cor_pseudo` and compares them. Between runs
you change the satellites, elevation mask or LQLC scale in `BASE`.

```bash
.venv/bin/python spp_switches.py data.csv     # LS + LQLC → next free run<N>/ folder, viewer included
.venv/bin/python build_viewer.py run1         # rebuild the viewer for one run (newest if omitted)
```

**Input:** a CSV with `epoch, tow, sat_const, sat_x, sat_y, sat_z, cor_pseudo`

---

## 1. What a run compares

```mermaid
flowchart LR
    B["BASE<br/>systems · mask_deg · s"] --> LS["LS"]
    B --> LQ["LQLC"]
    LS --> R[/"run&lt;N&gt;/"/]
    LQ --> R
```

```mermaid
flowchart LR
    R1["run1<br/>GPS, 10°"] -.change one setting.-> R2["run2<br/>all GNSS, 15°"] -.-> R3["run3 …"]
```

Change **one** setting per run, so any difference between runs can be traced to it. Each run's `config.json`
records what was used.

> **`cor_pseudo` is already corrected.** The supervisor confirmed that the receiver applied ionosphere,
> troposphere and Sagnac (and the satellite clock) to `cor_pseudo`. So the `sagnac`, `tropo` and `iono`
> switches stay **off**: switching them on applies a correction twice. Testing showed exactly that: Sagnac moved
> the fix about 32 m east–west off the antenna site, and tropo + iono lowered the height by about 19 m.
> The correction functions are kept for a possible future `raw_pseudo` path.

---

## 2. Big picture

```mermaid
flowchart LR
    A[(data.csv)] --> B[main]
    B --> C{for each config<br/>in CONFIGS}
    C --> D["cfg = BASE + estimator"]
    D --> E["run(cfg, data)"]
    E --> F["score(df)"]
    F --> G[/"run&lt;N&gt;/out_config.csv"/]
    F --> H[/"run&lt;N&gt;/comparison.csv"/]
    C -. next config .-> C
    G --> V["build_viewer.py"] --> VH[/"run&lt;N&gt;/viewer.html"/]
```

## 3. Per-epoch solver flow

```mermaid
flowchart TD
    S([epoch g, previous fix θ0]) --> INIT["rho = cor_pseudo<br/>θ = [x,y,z from θ0, one clock per system]"]
    INIT --> LOOP

    subgraph LOOP [Iteration loop, max 2 × MAX_ITER]
        direction TB
        P{"‖θ_xyz‖ > 6000 km?<br/>(have a position yet)"}
        P -- no --> Z["corr = 0, keep all"]
        P -- yes --> IT["<b>iter_terms</b><br/>el, az from current θ<br/>Sagnac + tropo + iono → corr<br/>mask: el ≥ mask_deg"]
        Z --> Y
        IT --> Y["residual  y = rho − (range + clock + corr)<br/>Jacobian  H = [−LOS unit vector | system clock 1]"]
        Y --> W["w = keep (0 or 1)"]
        W --> RB{robust_on?}
        RB -- yes --> LW["w ×= lqlc_weight(y, s)"]
        RB -- no --> CHK
        LW --> CHK{"enough sats &<br/>every system has weight?"}
        CHK -- no --> NONE([return None → skip epoch])
        CHK -- yes --> STEP["dθ = (HᵀWH)⁻¹ HᵀW y<br/>θ ← θ + dθ"]
        STEP --> CONV{"‖dθ_xyz‖ < 1 mm?"}
        CONV -- no --> P
        CONV -- "yes, LQLC & not robust yet" --> ON["robust_on = True"] --> P
    end

    CONV -- "yes, done" --> OK([return θ, converged])
    LOOP -- "iterations used up" --> NC([return θ, not converged])
    OK --> NEXT["next epoch starts from this fix"]
```

All corrections are recomputed **every iteration**, because they depend on the current position estimate.

## 4. LS → LQLC switch-over

```mermaid
stateDiagram-v2
    [*] --> LS : epoch start
    LS --> LS : step ≥ 1 mm
    LS --> Done : converged and estimator = LS
    LS --> LQLC : converged and estimator = LQLC
    LQLC --> LQLC : step ≥ 1 mm (weights recomputed)
    LQLC --> Done : converged
    Done --> [*]
```

LS runs first because robust weights only make sense once the residuals are small. From Earth's centre every
residual is about 20,000 km and every satellite would get roughly the same weight.

---

## 5. Functions

```mermaid
flowchart TB
    main --> run --> solve_epoch
    main --> score
    solve_epoch --> iter_terms
    solve_epoch --> lqlc_weight
    iter_terms --> el_az --> ecef_to_lla
    el_az --> enu_matrix
    iter_terms --> sagnac
    iter_terms --> tropo_saastamoinen --> ecef_to_lla
    iter_terms --> iono_klobuchar --> ecef_to_lla
    score --> ecef_to_lla
    score --> enu_matrix
```

| Function | When it runs | What it does |
|---|---|---|
| `ecef_to_lla(p)` | helper | ECEF → geodetic lat/lon/height (WGS-84, 5 iterations) |
| `enu_matrix(lat, lon)` | helper | Rotation from an ECEF difference to local East/North/Up |
| `el_az(sat, rx)` | every iteration | Elevation and azimuth of each satellite from the receiver |
| `sagnac(sat, rx)` | every iteration | Earth-rotation range term (m) |
| `tropo_saastamoinen(el, rx)` | every iteration | Tropospheric delay (m). Returns 0 outside 0–10 km height (still converging) |
| `iono_klobuchar(el, az, rx, tow)` | every iteration | Ionospheric L1 delay (m), default coefficients |
| `iter_terms(...)` | every iteration | Adds up the switched-on corrections and applies the elevation mask |
| `lqlc_weight(y, s)` | robust phase | Robust weight from the residuals |
| `solve_epoch(g, θ0, cfg)` | per epoch | Gauss-Newton / IRLS solver |
| `run(cfg, data)` | per config | Loops over epochs and carries each converged fix forward |
| `score(df, ref)` | per config | Spread statistics, or error and bias when `REF_XYZ` is set |

---

## 6. Settings

`BASE` holds the settings of a run. Each config is `BASE` + an estimator.

| Setting | Values | Meaning |
|---|---|---|
| `systems` | list | `1` GPS, `4` GLONASS, `8` Galileo, `32` BeiDou. One clock unknown per system |
| `mask_deg` | degrees | Satellites below this elevation get weight 0 |
| `s` | metres | LQLC scale. Larger values behave more like LS |
| `estimator` | `LS` / `LQLC` | Set by `CONFIGS`: every run solves both |
| `sagnac` | bool | **Keep `False`**: already in `cor_pseudo` (only for raw pseudoranges) |
| `tropo` | bool | **Keep `False`**: already in `cor_pseudo` (only for raw pseudoranges) |
| `iono` | bool | **Keep `False`**: already in `cor_pseudo` (only for raw pseudoranges) |

Other settings: `REF_XYZ` (surveyed antenna ECEF, enables true-error scoring), `MAX_ITER`, `TOL`.

---

## 7. Math

### 7.1 Measurement model

For satellite $i$ in constellation $k$:

$$
\rho_i = \lVert \mathbf{r}^{s}_i - \mathbf{r} \rVert + c\,\delta t_k + S_i + T_i + I_i + \varepsilon_i
$$

Unknowns: $\boldsymbol\theta = [x,\ y,\ z,\ c\delta t_{GPS},\ \dots]$. Every constellation keeps its own time, so each
one gets its own receiver clock.

### 7.2 Linearisation

```mermaid
flowchart LR
    T["θₙ (current guess)"] --> PR["predicted<br/>range + clock + corr"]
    RHO["rho (measured)"] --> RES["y = rho − predicted"]
    PR --> RES
    T --> HM["H rows:<br/>[ −e_i | 0 … 1 … 0 ]"]
    RES --> WLS["dθ = (HᵀWH)⁻¹HᵀWy"]
    HM --> WLS
    WLS --> T2["θₙ₊₁ = θₙ + dθ"]
```

$$
y_i = \rho_i - \left(\lVert \mathbf{r}^{s}_i - \mathbf{r}_n \rVert + c\delta t_k + S_i + T_i + I_i\right),
\qquad
\mathbf{H}_i = \Big[\, -\mathbf{e}_i^\top \;\big|\; 0 \cdots 1_{(k)} \cdots 0 \,\Big],
\qquad
\mathbf{e}_i = \frac{\mathbf{r}^{s}_i - \mathbf{r}_n}{\lVert \mathbf{r}^{s}_i - \mathbf{r}_n \rVert}
$$

### 7.3 Weighted least squares step

$$
\Delta\boldsymbol\theta = (\mathbf{H}^\top \mathbf{W} \mathbf{H})^{-1}\mathbf{H}^\top \mathbf{W}\,\mathbf{y},
\qquad \mathbf{W} = \operatorname{diag}(w_i)
$$

Iterate until $\lVert \Delta\boldsymbol\theta_{xyz} \rVert < 1\ \text{mm}$.

```mermaid
flowchart LR
    M["mask<br/>el ≥ mask_deg → 1, else 0"] --> X(("×"))
    LQ["LQLC<br/>tanh(x)/x<br/>(robust phase only)"] --> X
    X --> W["w_i"]
```

### 7.4 LQLC robust estimator

Instead of $\sum r_i^2$, minimise

$$
J(\boldsymbol\theta) = \sum_i \ln\cosh\!\left(\frac{r_i}{2s}\right)
$$

| Residual size | Behaviour | Effect |
|---|---|---|
| $\lvert r\rvert \ll 2s$ | $\ln\cosh u \approx u^2/2$ | Same as least squares |
| $\lvert r\rvert \gg 2s$ | $\ln\cosh u \approx \lvert u\rvert$ | Like L1: an outlier can't dominate |

IRLS weight from the influence function $\psi(r) = \frac{1}{2s}\tanh\frac{r}{2s}$:

$$
w_i = \frac{\psi(r_i)}{r_i} \;\propto\; \frac{\tanh x_i}{x_i},
\qquad x_i = \frac{\lvert r_i\rvert}{2s}
$$

$w \to 1$ for small residuals and $w \approx 2s/\lvert r\rvert$ for large ones. The constant $1/(4s^2)$ cancels in
the normal equations.

### 7.5 Sagnac (switch, off: already in `cor_pseudo`)

While the signal travels for about 70 ms, the Earth rotates under it:

$$
S_i = \frac{\omega_e}{c}\left(x^{s}_i\,y_r - y^{s}_i\,x_r\right),
\qquad \omega_e = 7.2921151467\times10^{-5}\ \text{rad/s}
$$

### 7.6 Troposphere: Saastamoinen (switch, off: already in `cor_pseudo`)

$$
\begin{aligned}
P &= 1013.25\,(1 - 2.2557\times10^{-5} h)^{5.2568} \quad \text{hPa}\\
T &= 15 - 0.0065\,h + 273.16 \quad \text{K}\\
e &= 6.108 \times 0.7 \times \exp\!\left(\frac{17.15\,T - 4684}{T - 38.45}\right) \quad \text{(70\% humidity)}\\
\text{ZHD} &= \frac{0.0022768\,P}{1 - 0.00266\cos 2\varphi - 0.00028\,h/1000}\\
\text{ZWD} &= 0.002277\left(\frac{1255}{T} + 0.05\right)e\\
T_i &= \frac{\text{ZHD} + \text{ZWD}}{\sin(\max(el_i,\ 1^\circ))}
\end{aligned}
$$

### 7.7 Ionosphere: Klobuchar (switch, off: already in `cor_pseudo`)

```mermaid
flowchart LR
    A["el, az, receiver lat/lon"] --> B["pierce point<br/>φ_i, λ_i"]
    B --> C["geomagnetic lat φ_m"]
    C --> D["amplitude A, period P<br/>(cubic in φ_m)"]
    T["tow → local time t"] --> E
    D --> E["cosine model<br/>day bump at 14:00"]
    E --> F["× slant factor F<br/>× (f_L1/f)²"]
```

All angles are in semicircles. $E$ is elevation and $A_z$ is azimuth:

$$
\begin{aligned}
\psi &= \frac{0.0137}{E + 0.11} - 0.022 \\
\varphi_i &= \varphi_u + \psi\cos A_z \quad (\text{clamped to } \pm 0.416) \\
\lambda_i &= \lambda_u + \frac{\psi \sin A_z}{\cos \varphi_i} \\
\varphi_m &= \varphi_i + 0.064\cos(\lambda_i - 1.617) \\
t &= 43200\,\lambda_i + t_{GPS} \pmod{86400} \\
A &= \textstyle\max\!\left(0,\ \sum_{n=0}^{3}\alpha_n\varphi_m^n\right),\qquad
P = \max\!\left(72000,\ \sum_{n=0}^{3}\beta_n\varphi_m^n\right) \\
x &= \frac{2\pi(t - 50400)}{P},\qquad F = 1 + 16(0.53 - E)^3 \\
I_i &= c\,F \cdot
\begin{cases}
5\times10^{-9} + A\left(1 - \frac{x^2}{2} + \frac{x^4}{24}\right) & |x| < 1.57\\
5\times10^{-9} & \text{otherwise}
\end{cases}
\end{aligned}
$$

For non-GPS systems the L1 delay is scaled by $(f_{L1}/f)^2$, with nominal frequencies for GLONASS G1, Galileo E1
and BeiDou B1I.

### 7.8 Elevation and azimuth

```mermaid
flowchart LR
    A["r_sat − r_rx<br/>(ECEF)"] -->|"× R(lat, lon)"| B["E, N, U"]
    B --> C["el = atan2(U, √(E²+N²))"]
    B --> D["az = atan2(E, N) mod 2π"]
```

$$
\mathbf{R} =
\begin{bmatrix}
-\sin\lambda & \cos\lambda & 0\\
-\sin\varphi\cos\lambda & -\sin\varphi\sin\lambda & \cos\varphi\\
\cos\varphi\cos\lambda & \cos\varphi\sin\lambda & \sin\varphi
\end{bmatrix}
$$

---

## 8. Scoring

```mermaid
flowchart TD
    D["all fixes P (ECEF)"] --> Q{REF_XYZ set?}
    Q -- no --> S1["spread around own mean (ENU)<br/>h_spread RMS, h95, v std<br/><i>= precision only</i>"]
    Q -- yes --> S2["error vs truth (ENU)<br/>hRMS, h95, vRMS, bias E/N/U<br/><i>= accuracy</i>"]
```

Without `REF_XYZ` the table only measures **precision**: how tightly the fixes cluster. A correction mostly shifts
the mean position, which spread can't show. Watch `mean_h`, `lat` and `lon` between runs, and set `REF_XYZ` to a
surveyed antenna position to measure true accuracy.

---

## 9. Runs and viewer

Every run of `spp_switches.py` gets its own folder, so earlier results are never overwritten:

```
run1/
├── config.json        settings of every config in this run
├── comparison.csv     one summary row per config
├── out_LS.csv        one fix per epoch, LS
├── out_LQLC.csv      one fix per epoch, LQLC
└── viewer.html        built automatically at the end of the run
run2/
└── …
```

`build_viewer.py` takes the newest `run<N>/` (or the one you name), reads the config names from its
`comparison.csv`, converts every `out_<config>.csv` to lat/lon/h and fills `viewer_template.html`. The resulting
`viewer.html` shows every config on a map with tracks, per-epoch values and statistics against a chosen
reference.

---

## 10. Starting a new run

Change one setting in `BASE`, then run again. The results go into the next `run<N>/`:

```python
BASE = dict(
    ...
    systems=[1, 4, 8, 32],   # e.g. run3: all GNSS
    mask_deg=10.0,           #        back to 10°
    s=3.0,
)
```

To compare more than LS and LQLC in one run, add an entry to `CONFIGS`:

```python
CONFIGS = {
    "LS":      dict(estimator="LS"),
    "LQLC":    dict(estimator="LQLC"),
    "LQLC_s5": dict(estimator="LQLC", s=5.0),   # → out_LQLC_s5.csv
}
```
