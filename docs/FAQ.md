# FAQ

**Q: I see `ModuleNotFoundError: fighi_ext`.**
A: Ensure the *parent* of `fighi_ext/` is on `PYTHONPATH`:
```bash
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"
```

If you installed editable (pip install -e .), you can run python -m fighi_ext.run_cli.

Q: What does --max_order do?
A: It caps the interaction order K explored. The algorithm may stop earlier using the information ratio rule when additional orders add little information.

Q: Memory blows up during merging.
A: Use merge_pheno_geno_nopandas.py and, if needed, split inputs, then cat or awk join by IID. During modeling, use --read_chunksize.

Q: How do I annotate SNPs automatically?
A: Run annotate_fighi_features.py with --cs2g_dir (offline) or --use_gprofiler (online). See --help in the script.

Q: How to run on Slurm?
A: Use fighi_ext/examples/demo_slurm.job. It sets BLAS threads and writes into fighi_ext/fighi_out.

```yaml

---

## 7) `CONTRIBUTING.md`

```markdown
# Contributing

Thanks for considering contributions! Ways to help:
- Report issues with minimal reproductions.
- Improve docs, add examples, or tutorials.
- Optimize inner loops or add new pruning strategies.
- Integrate additional annotation sources.

## Dev setup
```bash
git clone https://github.com/<you>/FIGHI.git
cd FIGHI
conda env create -f environment.yml
conda activate fighi
pip install -e .[dev]
pre-commit install
```

Tests & style

Use ruff/flake8 for linting; pytest for tests.

Open a PR against main; CI must pass.

```yaml

---

## 8) `CODE_OF_CONDUCT.md`

(Use a standard Contributor Covenant v2.1 — omitted here for brevity.)

---

## 9) `LICENSE`

```text
MIT License

Copyright (c) 2025 ...

Permission is hereby granted, free of charge, to any person obtaining a copy...
```

