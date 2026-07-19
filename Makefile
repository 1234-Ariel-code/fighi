.PHONY: test demo verify build check clean

test:
	python -m unittest discover -s tests -v

demo:
	fighi demo --outdir build/demo --overwrite

verify:
	bash scripts/verify_release.sh

build:
	python -m build

check: test
	python -m compileall -q src
	bash -n scripts/*.sh scripts/arc/*.sh scripts/arc/*.sbatch

clean:
	python -c "from pathlib import Path; import shutil; [shutil.rmtree(p) for p in [Path('build'), Path('dist')] if p.exists()]"
