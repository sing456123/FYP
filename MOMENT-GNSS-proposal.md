# Project Proposal: MOMENT Foundation Model as GNSS Bad-Signal Pre-Filter

**Date:** 2026-10-29  
**Context:** DGNSS / PPP positioning in urban canyon environments (NLOS, multipath)  
**Proposed Model:** `AutonLab/MOMENT-1-large` (341M params, CMU Auton Lab, ICML'24)

---

## 1. Objective

Replace or augment the existing LQLC robust estimator's satellite weighting with a **zero-shot, data-driven anomaly detector** that flags individual satellites as "bad signal" *before* they enter the WLS/Kalman update step.

Current pipeline:
```
Raw DD observation → Elevation mask → Fixed S → LQLC weighting → Position solution
```

Proposed pipeline:
```
Raw DD observation → Elevation mask → MOMENT anomaly score → Down-weight / veto → LQLC/Kalman
```

---

## 2. Why MOMENT?

| Requirement | MOMENT-1-large Status |
|---|---|
| Zero-shot (no training on GNSS data) | ✅ Pre-trained on Time-Series Pile |
| Anomaly detection via reconstruction error | ✅ `task_name="reconstruction"` |
| Causal (only uses history, no lookahead) | ✅ Mask last 8 epochs, reconstruct from prior |
| Open-source & reproducible | ✅ MIT license, Hugging Face |

**Key caveat verified by prior testing:** `MOMENT-1-small` and `MOMENT-1-base` ship with an **untrained reconstruction head** (all singular values = 1.41, bias = 0.0). Only **`MOMENT-1-large`** produces meaningful reconstruction. See `skills/moment-filter/SKILL.md` for forensic evidence.

---

## 3. Preliminary Validation (data.csv)

A quick prototype was run on a real GNSS dataset (`data.csv`, 6693 epochs, 38 satellites, GPS/GLONASS/Galileo/BeiDou).

### Method
- Per-satellite SNR time series (last 512 epochs) → z-score standardization
- MOMENT-1-large reconstructs last 8 epochs from prior 504
- Residual = mean |actual − reconstructed|

### Results: Top Outliers (90th percentile threshold = 0.368)

| Rank | Satellite | Residual | Mean SNR | Last SNR | Observation |
|---|---|---|---|---|---|
| 1 | **R51** (GLONASS) | 0.923 | 33.4 | 35 | SNR steady but low — model sees "unusual" pattern |
| 2 | **G24** (GPS) | 0.893 | 34.8 | 38 | Same as above |
| 3 | **R34** (GLONASS) | 0.553 | 38.5 | **33** | MOMENT catches SNR *drop* |
| 4 | **C133** (BeiDou) | 0.448 | 35.6 | **39** | MOMENT flags SNR *rise* too — any pattern break |

**Correlation(residual, last SNR) = 0.181** — weak, confirming MOMENT detects **pattern breaks**, not just low SNR.

**Files:**
- Script: `tmp/run_moment_gnss.py`
- Results: `tmp/moment_gnss_results.json`

---

## 4. What the Coder Should Implement

### 4.1 Core Module: `gnss_moment_filter.py`

```python
class GNSSMomentFilter:
    def __init__(self, ckpt="AutonLab/MOMENT-1-large", seq_len=512, patch_len=8):
        ...

    def score_satellite(self, snr_series: np.ndarray) -> float:
        """
        Input: 1D SNR time series (512 epochs, latest last)
        Output: reconstruction residual (higher = more anomalous)
        """
        ...

    def score_epoch(self, epoch_df: pd.DataFrame) -> dict[sat_id, float]:
        """
        Input: single epoch with all visible satellites
        Output: per-satellite anomaly scores
        Builds rolling 512-epoch buffer per satellite on the fly.
        """
        ...

    def veto_mask(self, scores: dict, threshold: float) -> dict[sat_id, bool]:
        """
        Returns True for satellites to EXCLUDE from position solution.
        """
        ...
```

### 4.2 Integration Point

Hook into the existing observation processing pipeline *after* elevation mask, *before* LQLC weight calculation:

```python
# In existing positioning solver:
visible_sats = apply_elevation_mask(all_sats, cutoff=5°)

# NEW: MOMENT veto layer
moment_scores = moment_filter.score_epoch(current_epoch_data)
vetoed = moment_filter.veto_mask(moment_scores, threshold=0.90_quantile)
clean_sats = [s for s in visible_sats if not vetoed[s]]

# Then existing LQLC on clean_sats only
weights = lqlc_weight(clean_sats)
position = weighted_least_squares(clean_sats, weights)
```

### 4.3 Calibration / Threshold Setting

Two approaches to explore:

| Approach | Method | Pros | Cons |
|---|---|---|---|
| **A. Quantile-based** | Use historical residuals' empirical distribution; flag above q90 | Simple, no labels needed | May need per-environment calibration |
| **B. SNR correlation gate** | Flag if residual > threshold AND last SNR < mean − k·std | Leverages known physics | May miss non-SNR anomalies (e.g., cycle slip) |

### 4.4 Input Data Requirements

The filter needs a **rolling 512-epoch history** per satellite. Required columns per epoch:
- `epoch` (integer time index)
- `sat_id` + `sat_const` (unique satellite identifier)
- `sat_snr` (signal quality metric)

Optional (for future extension):
- `cor_pseudo` / `raw_pseudo` (for pseudorange residual instead of SNR)
- `sat_el` (elevation, for validation)

---

## 5. Evaluation Plan

Compare three configurations on the same dataset:

| Config | Description | Expected Outcome |
|---|---|---|
| **Baseline** | Elevation mask (5°) + fixed S + LQLC only | Reference |
| **A** | Baseline + MOMENT veto (q90) + LQLC on clean sats | Fewer outliers, cleaner residuals |
| **B** | Baseline + MOMENT veto (q90) + adaptive R in KF | Best if temporal smoothing needed |

Metrics:
- Position RMSE vs ground truth (if available)
- Residual distribution: mean, std, max (should see max drop)
- Number of epochs with >1 satellite vetoed (indicator of urban vs open sky)

---

## 6. Known Risks & Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| MOMENT latency (341M params) | Real-time 10 Hz may struggle | Batch process offline first; profile on target hardware |
| MOMENT not trained on GNSS | May not generalize | Test on UrbanNav dataset with known NLOS labels |
| SNR alone insufficient | Cycle slips have high SNR | Future: switch input to DD pseudorange residual |
| False positives | Good sats vetoed → position degrade | Tune threshold empirically; start with loose gate |

---

## 7. Deliverables

1. **`gnss_moment_filter.py`** — core class with `score_satellite`, `score_epoch`, `veto_mask`
2. **`integrate_with_solver.py`** — example hook showing how to insert into existing pipeline
3. **`evaluate.py`** — run Baseline vs Config A/B on provided `data.csv`, output metrics + plots
4. **`README.md`** — setup instructions (how to install `momentfm`, which checkpoint to download)

---

## 8. Environment

```bash
# Python 3.11 recommended by MOMENT authors
pip install momentfm
# Or from source:
pip install git+https://github.com/moment-timeseries-foundation-model/moment.git
```

**Model checkpoint:** `AutonLab/MOMENT-1-large` (~341M, downloads on first use via Hugging Face)

---

## Appendix: Quick Reference — SNR Series → MOMENT Pipeline

```
Per-satellite SNR (1 epoch per row)
        |
        v
Build rolling buffer (last 512 epochs, per satellite)
        |
        v
Z-score standardize (per satellite, removes level bias between constellations)
        |
        v
MOMENT-1-large: mask last 8 epochs, reconstruct from 504
        |
        v
Residual = mean |actual − recon| over last 8 epochs
        |
        v
Threshold: empirical 90th percentile across all satellites in window
        |
        v
Flag satellite as outlier if residual > threshold
```

---

*Prepared for coder handoff. Contact: Chun Sing Tse*
