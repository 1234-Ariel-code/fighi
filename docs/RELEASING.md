# Release checklist

1. Update the version in `pyproject.toml`, `src/fighi/__init__.py`, and `CHANGELOG.md`.
2. Run `python -m unittest discover -s tests -v` on supported Python versions.
3. Run a binary demo, a linear API example, a null simulation, and TPED conversion.
4. Build with `python -m build`.
5. Validate with `python -m twine check dist/*`.
6. Install the wheel in a clean environment and run `fighi demo`.
7. Confirm every JSON file parses and GraphML opens.
8. Review dependency and security alerts.
9. Tag `vX.Y.Z` and create a GitHub release with the changelog.
10. Publish to TestPyPI, install and smoke-test, then publish the identical artifacts to PyPI.
11. Archive the release on Zenodo and add the DOI to `CITATION.cff` when available.

Never publish from an uncommitted or unreviewed working tree. Do not upload private genotype or phenotype data as release assets.

