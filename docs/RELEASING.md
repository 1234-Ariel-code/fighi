# Release checklist

1. Update the version in `pyproject.toml`, `src/fighi/__init__.py`, and `CHANGELOG.md`.
2. Run `scripts/verify_release.sh` and `python -m unittest discover -s tests -v` on supported Python versions.
3. Run a binary demo, a linear API example, null/pairwise simulations, a fake or public PLINK preparation test, a two-method benchmark, and TPED conversion.
4. Build with `python -m build`.
5. Validate with `python -m twine check dist/*`.
6. Install the wheel in a clean environment and run `fighi demo`.
7. Confirm every JSON file parses, GraphML opens, the benchmark HTML renders, and ARC scripts pass `bash -n`.
8. Review dependency and security alerts.
9. Build the Git bundle and source archive, regenerate `SHA256SUMS.txt`, tag `vX.Y.Z`, and create a GitHub release with the changelog.
10. Publish to TestPyPI, install and smoke-test, then publish the identical artifacts to PyPI.
11. Archive the release on Zenodo and add the DOI to `CITATION.cff` when available.

Never publish from an uncommitted or unreviewed working tree. Do not upload private genotype or phenotype data as release assets.
