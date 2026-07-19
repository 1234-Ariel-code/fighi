# Inputs and outputs

## Tabular input

Rows represent samples. Columns may contain:

- one phenotype column;
- an optional sample ID;
- optional numeric or categorical covariates;
- numeric genotype dosage columns.

CSV and TSV are supported, including `.csv.gz`, `.tsv.gz`, and `.tab.gz`. Feature values are expected in `[0, 2]` unless `--allow-non-genotype` is deliberately set. Common ID names are auto-detected; specify `--id-column` when ambiguous.

## PLINK preparation outputs

`prepare-plink` accepts a BED/BIM/FAM prefix, PGEN/PVAR/PSAM prefix, or VCF/BCF plus a required candidate list. Its main output is `fighi_input.tsv.gz`; genotypes count alternate alleles under PLINK 2 `--export A include-alt`. It also writes `fighi_candidates.txt`, `fighi_samples.txt`, an exact command, stdout/stderr logs and `prepare_plink_manifest.json` containing hashes and missing candidate IDs.

Phenotype/covariate tables are joined by a unique string IID. Samples without an observed phenotype are removed and counted. The wrapper refuses ambiguous duplicate candidate IDs or duplicate phenotype IIDs.

## Quality control

- Samples with missing phenotype are removed.
- Features above `max_missing` are removed.
- Remaining feature missingness is mean-imputed.
- Constant and near-constant features are removed.
- Genotype features below `min_maf` are removed.
- Numeric covariates use median imputation.
- Categorical covariates use mode imputation and one-hot encoding.
- Constant covariate encodings are removed.

Every change is counted in `fighi_summary.json`; nonfatal changes are also listed as warnings.

## Result columns

`fighi_results.csv` contains:

| Column | Meaning |
|---|---|
| `hyperedge` | Feature names separated by `|` |
| `order` | Interaction order |
| `statistic` | Efficient score-test statistic |
| `fi_gain` | Half the score statistic; ranking measure |
| `beta_hat` | One-step highest-order coefficient estimate |
| `odds_ratio` | `exp(beta_hat)` for binary traits |
| `information` | Effective information after nuisance projection |
| `p_value` | Asymptotic one-degree-of-freedom p-value |
| `q_order` | Discovery/inference per-order adjusted value used diagnostically |
| `q_value` | Final global adjusted p-value |
| `significant` | Whether `q_value <= alpha` |
| `selected_for_expansion` | Candidate used in the discovery beam |
| `stability` | Optional conditional subsampling frequency |
| `null_converged` | Logistic null convergence indicator |
| `null_iterations` | Logistic null iterations |

## Hypergraphs

Only multiplicity-controlled significant interactions are exported to networks. Each interaction is represented as a hyperedge node connected to its member feature nodes. JSON preserves native hyperedges; GML, GraphML, Cytoscape, and TSV use the incidence-node representation.

## Benchmark normalized columns

`benchmark_interactions.csv` stores `method`, sorted `hyperedge`, `order`, raw `p_value`, optional `rank_score`, internal `rank_value`, optional `effect`, source file, harmonized q-value and harmonized significance flag. Ranking-only methods have missing p/q values by design; their declared score is used only for ranking metrics and overlap.
