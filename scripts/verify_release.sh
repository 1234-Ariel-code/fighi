#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:$PWD/src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/fighi-matplotlib}"

python -m compileall -q src
python -m unittest discover -s tests -v

work="$(mktemp -d "${TMPDIR:-/tmp}/fighi-release-XXXXXX")"
trap 'rm -rf "$work"' EXIT

python -m fighi simulate \
  --outdir "$work/simulation" \
  --samples 240 \
  --features 12 \
  --scenario pairwise \
  --effect 1.4 \
  --seed 17
python -m fighi benchmark-template --output "$work/simulation/benchmark.json"
python -m fighi benchmark-run \
  "$work/simulation/benchmark.json" \
  --outdir "$work/benchmark" \
  --strict

test -s "$work/benchmark/benchmark_report.html"
test -s "$work/benchmark/benchmark_method_summary.csv"
echo "Release verification passed"
