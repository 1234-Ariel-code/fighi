```yaml

---

```markdown
# Usage & Outputs

## CLI
```bash
python fighi_ext/run_cli.py \
  --csv <merged.csv> \
  --pheno <phen_col> \
  --trait {binary|linear} \
  --outdir <outdir> \
  --max_order 4 \
  [--target_or 1.3] [--read_chunksize 0]

### Important args

--max_order: hard cap on interaction order K (the search may stop earlier via information ratio).

--target_or: planner’s detectable odds ratio (guides feasible K under sample size & MAF).

--read_chunksize: if your merged CSV is huge, read in chunks (reduces peak RAM).

### Output files

#### fighi_results.csv

column	meaning

- hyperedge	tuple of rsIDs separated by `
order	interaction order K

- fi_gain	Fisher-Information gain (score-test)

- pval	optional, if permutation is enabled

- beta_hat	one-step estimate(s) for coefficients

- info	information (observed Fisher) contributing to cumulative ratio

#### fighi_feature_scores.csv

- Per-SNP totals aggregated across all retained edges:

- FI_total, FI_main, FI_interact, Rank, MAF

- Optional: Gene, Pathway (via annotate_fighi_features.py)

#### Hypergraphs

- fighi_hypergraph.gml — Gephi/Cytoscape import

- fighi_cytoscape.cyjs — direct Cytoscape session element JSON

- fighi_hypergraph.hyper — simple JSON with nodes + hyperedges

#### Summary & diagnostics

- fighi_summary.json — n_samples, n_snps, feasible K, runtime, etc.

- fighi_log.txt — run log with convergence notes & pruning summaries

- plots/ — FI distribution, interaction heatmap, FI vs order, convergence traces

#### Memory knobs

- Use merge_pheno_geno_nopandas.py to avoid pandas during the heaviest join.

- Use --read_chunksize in run_cli.py for very wide CSVs.

- Limit --max_order and/or lower the atom screening quota (see pipeline.py:screen).

- Control BLAS threads: OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS.
```
