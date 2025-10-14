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


