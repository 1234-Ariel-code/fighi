```markdown
# Quick Start

## 0) Environment
```bash
conda env create -f environment.yml
conda activate fighi
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"

### 1) From merged CSV

python fighi_ext/run_cli.py \
  --csv path/to/merged.csv \
  --pheno case \
  --trait binary \
  --outdir fighi_ext/fighi_out \
  --max_order 4

### 2) From TPED/TFAM + .phen

python fighi_ext/tped_to_named_csv.py \
  --tped data/CD_origin.tped --tfam data/CD_origin.tfam \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_filtered_named.csv

python fighi_ext/merge_pheno_geno_nopandas.py \
  --geno_csv data/CD_filtered_named.csv \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_merged.csv \
  --phen_name case --impute_mean

python fighi_ext/run_cli.py \
  --csv data/CD_merged.csv --pheno case \
  --trait binary --outdir fighi_ext/fighi_out \
  --max_order 4

### 3) Annotation & enrichment (optional)

SNP→Gene via cS2G or g:Profiler; add Gene/Pathway columns
python fighi_ext/annotate_fighi_features.py \
  --feature_csv fighi_ext/fighi_out/fighi_feature_scores.csv \
  --cs2g_dir fighi_ext/cS2G_1000GEUR \
  --gmt_files path/KEGG.gmt path/REACTOME.gmt path/GO_BP.gmt

### 4) Slurm example

See fighi_ext/fighi_job.slurm One line to set disease:

DISEASE_NAME="CD"

Saves outputs to fighi_ext/fighi_out.
```
