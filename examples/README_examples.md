# FIGHI Toy Examples

This folder contains a tiny, deterministic dataset and minimal commands to verify the FIGHI pipeline.

## Files
- `toy_merged.csv` — phenotype `case` + 8 SNPs (rs101..rs108) encoded 0/1/2.
- `toy.tped`, `toy.tfam`, `toy.phen` — same dataset in TPED/TFAM/PHEN form.
- `run_toy_csv.sh` — run FIGHI directly on the merged CSV.
- `run_toy_tped.sh` — full flow: TPED→named CSV→merge→FIGHI.

## Quick start
```bash
# From repo root (ensure `fighi_ext` is on PYTHONPATH or installable):
bash examples/run_toy_csv.sh

# Or full TPED flow:
bash examples/run_toy_tped.sh
```
