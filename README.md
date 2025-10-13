# FIGHI
FIGHI: Fisher-Information–Guided Hyper-interaction Inference for genome-wide epistasis discovery, summaries, and hypergraph outputs.

FIGHI/
├─ fighi_ext/                     # Python package + scripts
│  ├─ __init__.py
│  ├─ adaptive.py
│  ├─ apriori.py
│  ├─ fisher.py
│  ├─ glm.py
│  ├─ hypergraph.py
│  ├─ knockoff_perm.py
│  ├─ pipeline.py
│  ├─ report.py
│  ├─ run_cli.py                  # can run as `python -m fighi_ext.run_cli` or `python fighi_ext/run_cli.py`
│  ├─ tped_to_named_csv.py
│  ├─ merge_pheno_geno_nopandas.py
│  ├─ annotate_fighi_features.py
│  └─ utils.py
├─ examples/
│  └─ (tiny toy data / minimal commands)
├─ docs/
│  ├─ README_theory.md            # full math/theory (we drafted earlier)
│  └─ assets/
│     ├─ fighi_protocol.png       # flowchart
│     └─ fighi_cover.png          # cover art
├─ fighi_job.slurm                # SLURM job (disease-var enabled)
├─ README.md                      # user-facing quickstart
├─ LICENSE                        # MIT (if chosen)
└─ .gitignore


