.PHONY: test demo build check clean

test:
	python -m unittest discover -s tests -v

demo:
	fighi demo --outdir build/demo --overwrite

build:
	python -m build

check: test
	python -m compileall -q src

clean:
	python -c "from pathlib import Path; import shutil; [shutil.rmtree(p) for p in [Path('build'), Path('dist')] if p.exists()]"

