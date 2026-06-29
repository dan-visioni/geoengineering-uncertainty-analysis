#!/bin/bash
# Convert percent-format sources in python/src/*.py to executed .ipynb in
# j_notebooks/.  Usage: build.sh <name1> [name2 ...]   (name without .py)
set -e
ENVBIN=/opt/anaconda3/envs/geomip/bin
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p j_notebooks
for name in "$@"; do
  src="python/src/${name}.py"
  out="j_notebooks/${name}.ipynb"
  echo ">>> $name"
  "$ENVBIN/jupytext" --to ipynb "$src" -o "$out" --quiet
  "$ENVBIN/jupyter" nbconvert --to notebook --execute --inplace "$out" \
    --ExecutePreprocessor.timeout=1200 \
    --ExecutePreprocessor.kernel_name=geomip 2>&1 | tail -3
done
