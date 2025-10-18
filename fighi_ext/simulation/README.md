# FIGHI – Simulation 

This folder provides **reproducible simulation scripts** for the core statistical properties and practical behavior of FIGHI.


## What’s included

- **Type I Error (Null Calibration)** — `type1_error.py`  
  Rejection rates at α ∈ {0.10, 0.05, 0.01} after BH correction when scanning many pairs under a pure null.
- **Power Grid** — `power_grid.py`  
  Power vs sample size N, odds ratio OR, and MAF, with BH correction among distractor pairs.
- **Robustness to LD & Misspecification** — `robustness_ld.py`  
  Power under moderate/strong LD and when the true link is probit but the test is logistic (score).
- **Scaling (Runtime)** — `scaling_runtime.py`  
  Wall-time vs number of SNPs P for a fixed number of scanned pairs M.

```text
simulation/
   ├─ common_sim.py
   ├─ type1_error.py
   ├─ power_grid.py
   ├─ robustness_ld.py
   ├─ scaling_runtime.py
   └─ README.md
```

All scripts write CSV + PNG into `simulation/out/`.

## Usage

```bash
cd notebooks/sim
python type1_error.py
python power_grid.py
python robustness_ld.py
python scaling_runtime.py
```
```text
Dependencies: numpy, pandas, scipy, matplotlib
```

Install once:

```bash
pip install numpy pandas scipy matplotlib
```


## What to report

### Type I Error:
```text
- t1e_summary.csv → empirical rejection rates across reps for each α
- t1e_qq.png → null QQ plot of raw p-values
```
Expectation: rates near nominal, QQ close to the diagonal.


### Power:

```text
- power_grid.csv, power_heatmap.png, power_curves.png
```
Expectation: power increases with N, OR, MAF; compare detectability threshold with planner inequality.

### Robustness:

```text
- robust_summary.csv, robust_bar.png
```
Expectation: graceful degradation with LD; minor loss under probit misspecification.

## Scaling:

```text
-scaling_pairs.csv, scaling_pairs.png
```
Expectation: near-linear growth with number of scanned candidates.

## Reproducibility

```text
- Deterministic RNG (NumPy Generator, fixed seeds in each script)
- Headless plotting (Agg)
- No downloads required
- Code is easily portable to SLURM arrays if needed
```

```yaml

---

## Notes

- These studies focus on **pairwise** interactions for clarity and speed. You can extend to **triples** by introducing a small set of true triples and scanning a random subset of 3-way combinations (keeping M modest).
```
