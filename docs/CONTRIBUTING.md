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
git clone https://github.com/1234-Ariel-code/FIGHI.git
cd FIGHI
conda env create -f environment.yml
conda activate fighi
pip install -e .[dev]
pre-commit install
```
Tests & style

Use ruff/flake8 for linting; pytest for tests.

Open a PR against main; CI must pass.
```
