# Overview

This benchmark suite rigorously compares FIGHI (Fisher-Information–Guided Hyperinteraction Inference) with two established pairwise interaction approaches:

```text
- PLINK --epistasis (logistic / fast epistasis)

- MDR (Multifactor Dimensionality Reduction)
```

It evaluates runtime, peak memory, number of pairs tested, and predictive accuracy (AUC) using the same preprocessed dataset.

# Folder Structure

```bash
benchmarks/
├── benchmark_orchestrate.slurm     # Main orchestrator (runs all tools sequentially)
├── filter_csv_by_columns.py        # Keep only prescreened SNPs
├── csv_to_plink.py                 # Convert FIGHI CSV → PLINK format (.bed/.bim/.fam)
├── mdr_benchmark.py                # Run MDR pairwise scan with scikit-mdr
├── accuracy_eval.py                # Compute AUC/AP for FIGHI, PLINK, MDR
├── logs/                           # SLURM & runtime logs
└── out/                            # Output directory with timing, results, accuracy
```

# Input Data Requirements

```bash
| File            | Description                          | Format               |
| --------------- | ------------------------------------ | -------------------- |
| `{data}_merged.csv` | Merged phenotype + genotype CSV      | `IID,case,rs1,rs2,…` |
| `case` column   | Binary phenotype (0=control, 1=case) | Required             |
| SNP columns     | Dosages coded 0/1/2 (float or int)   | Required             |
```

# Dependencies

All dependencies are handled automatically within the SLURM job using a self-contained Conda environment.

Installed packages (via pip):

```nginx
numpy pandas scipy matplotlib statsmodels scikit-learn scikit-mdr
```

# External tools:

```text
- PLINK v1.9 (downloaded automatically if missing)
- Standard Linux tools: wget, /usr/bin/time, unzip
```

# Running the Benchmark

## Step 1. Submit the orchestrator job

```bash
sbatch benchmarks/benchmark_orchestrate.slurm
```

This job:

```text
- Prescreens SNPs using FIGHI (Stage 1) → defines a fair TOP_M subset
- Filters the dataset to those SNPs
- Runs FIGHI, PLINK, and MDR in sequence
- Records runtime and memory for each method (/usr/bin/time -v)
- Evaluates accuracy on an 80/20 held-out set
- Summarizes results in the terminal and writes JSON reports
```

## Key Parameters (edit inside benchmark_orchestrate.slurm)

```bash
| Parameter     | Description                            | Example |
| ------------- | -------------------------------------- | ------- |
| `TOP_M`       | SNPs retained after FIGHI prescreen    | `3000`  |
| `ROW_CHUNK`   | Row chunk size for streaming           | `50000` |
| `COL_BLOCK`   | Column block size                      | `1000`  |
| `MAX_ORDER`   | FIGHI max interaction order            | `3`     |
| `--cv_folds`  | MDR cross-validation folds             | `5`     |
| `--max_pairs` | Maximum MDR pairs (for runtime safety) | `1e7`   |
```

# Pipeline Steps in Detail

## Prescreening (FIGHI)

```bash
python run_cli.py \
  --csv CD_merged.csv --pheno case --trait binary \
  --prescreen_top_m 3000 --col_block 1000 --row_chunksize 50000 \
  --max_order 2 --write_top_cols
```

Outputs:

```text
prescreen/top_columns.txt (list of top informative SNPs by FI gain)
```

## Filtering the CSV

```bash
python filter_csv_by_columns.py \
  --csv CD_merged.csv --pheno case \
  --keep prescreen/top_columns.txt \
  --out cd_top3000.csv
```
Reduces genotype dimensionality → consistent subset for all tools.

## FIGHI Runtime

```bash
/usr/bin/time -v python run_cli.py \
  --csv cd_top3000.csv --pheno case --trait binary \
  --outdir out/fighi_run --max_order 3 --no_plots
```

Reports:

- Wall time
- Max RSS (memory) of pairs tested
- Interaction ΔI per order

Output:
out/fighi_run/fighi_results.csv

## PLINK Epistasis Runtime

### Convert CSV → PLINK:

```bash
python csv_to_plink.py \
  --csv cd_top3000.csv --pheno case --out_prefix out/plink/cd_top3000
```

### Run pairwise epistasis:

```bash
/usr/bin/time -v plink \
  --bfile out/plink/cd_top3000 \
  --fast-epistasis --threads 22 --allow-no-sex \
  --out out/plink/epi_fast
```

Output:
```text
out/plink/epi_fast.epi.cc
(columns: SNP_A, SNP_B, statistic, P)
```

### MDR Runtime

```bash
/usr/bin/time -v python mdr_benchmark.py \
  --csv cd_top3000.csv --pheno case \
  --cv_folds 5 --max_pairs 10000000 \
  --out out/mdr/mdr_pairs.csv
```

Output:
```text
out/mdr/mdr_pairs.csv
```

### Accuracy Evaluation (AUC/AP)
```text
Each model’s predictive performance is measured on a held-out 20% split.
(columns: SNP1, SNP2, CV_AUC)
```

### Output example:

```text
{
  "FIGHI": {
    "features": [["rs1"], ["rs2","rs3"]],
    "AUC_train": 0.81, "AUC_test": 0.76,
    "AP_train": 0.42, "AP_test": 0.38
  },
  "PLINK": {
    "pair": ["rsA","rsB"],
    "AUC_train": 0.64, "AUC_test": 0.61
  },
  "MDR": {
    "pair": ["rsX","rsY"],
    "AUC_train": 0.67, "AUC_test": 0.60
  }
}
```

### Outputs Summary

```text
| Tool      | Main Output                   | Timing File  | Accuracy                 |
| --------- | ----------------------------- | ------------ | ------------------------ |
| **FIGHI** | `fighi_run/fighi_results.csv` | `fighi.time` | `accuracy.json["FIGHI"]` |
| **PLINK** | `plink/epi_fast.epi.cc`       | `plink.time` | `accuracy.json["PLINK"]` |
| **MDR**   | `mdr/mdr_pairs.csv`           | `mdr.time`   | `accuracy.json["MDR"]`   |
```
```text
All times use /usr/bin/time -v
→ Elapsed (wall clock) time and Maximum resident set size.
```

## Interpreting Results

```text
- Wall clock time: total real time (benchmark’s key metric)
- Peak RSS (kB): peak memory footprint
- FIGHI vs PLINK:
  - FIGHI explores adaptive higher-order terms
  - PLINK/MDR limit to pairwise interactions
- Accuracy: ROC-AUC on hold-out; reflects predictive consistency
```

## Notes & Tips

```text
- All three methods share the same SNP subset from prescreening → ensures fairness.
- FIGHI automatically stops early when Fisher Information gain saturates — typical runtime advantage.
- For extremely large TOP_M, you can use BATCH_SIZE and HALO settings from the orchestration version of FIGHI to chunk analyses.
- To benchmark across diseases (e.g., T1D, CAD), replicate the same pipeline per phenotype.
```

## Example Summary Table (for manuscript/report)

```text
| Method | SNPs | Pairs Tested | Max Order | Wall Time (h) | Peak RAM (GB) | Test AUC |
| ------ | ---- | ------------ | --------- | ------------- | ------------- | -------- |
| FIGHI  | 3000 | ~4.5M        | 3         | 0.8           | 18            | **0.76** |
| PLINK  | 3000 | 4.5M         | 2         | 2.3           | 21            | 0.61     |
| MDR    | 3000 | 4.5M         | 2         | 4.1           | 33            | 0.60     |
```

## Extending the Benchmark

```text

To add:

```text
- New datasets → just replace CD_merged.csv
- Higher orders → set --max_order 4 (if memory allows)
```

FIGHI Benchmark Suite provides a transparent, reproducible environment to compare multi-order Fisher Information inference with traditional epistasis tools — quantifying both efficiency and biological informativeness under identical data and hardware settings. Parallel runs → convert each method block into its own SLURM script and orchestrate via afterok

