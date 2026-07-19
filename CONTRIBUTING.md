# Contributing

Thank you for improving FIGHI.

1. Open an issue describing the scientific or software problem.
2. Create a focused branch and include tests for changed behavior.
3. Keep statistical claims matched to implemented and validated behavior.
4. Run `make check` before opening a pull request.
5. Update documentation and the changelog when users will notice the change.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m unittest discover -s tests -v
```

Do not add identifiable participant data, controlled cohort data, credentials, or licensed reference databases. Synthetic test data must be generated from documented random seeds.

Changes to statistical methods should include the null and alternative hypotheses, assumptions, simulations, calibration checks, and expected failure modes. Performance improvements must preserve numerical results within a documented tolerance.

