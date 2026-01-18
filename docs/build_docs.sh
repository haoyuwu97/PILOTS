#!/usr/bin/env bash
set -euo pipefail

# Build the HTML manual using Sphinx.
# Output is written to docs/manual/.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python3 -m pip install -r "$ROOT_DIR/docs/sphinx/requirements.txt"

(cd "$ROOT_DIR/docs/sphinx" && make html)

echo "Built docs to: $ROOT_DIR/docs/manual/index.html"
