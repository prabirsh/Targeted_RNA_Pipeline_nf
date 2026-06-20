#!/usr/bin/env bash
# Initialise version control for the RNA panel pipeline.
# Run this ON YOUR OWN MACHINE (native filesystem) — not inside a mounted/DrvFs
# share, where git lock files fail.
#
#   cd rna-panel-pipeline
#   bash scripts/git_init.sh
set -euo pipefail

# If a broken .git stub exists (e.g. created on a network/mounted drive), remove it:
[ -d .git ] && { echo "Removing existing .git ..."; rm -rf .git; }

git init
git add -A

# Optional: track the large SureSelect BED with Git LFS (recommended).
if command -v git-lfs >/dev/null 2>&1; then
  git lfs install
  git lfs track "*.bed"
  git add .gitattributes
fi

git commit -m "feat: RNA panel pipeline v1.1.0 — panel selection, coverage QC, fusion read-names"
git tag -a v1.1.0 -m "RNA panel pipeline v1.1.0"
echo
echo "Done. Repo initialised and tagged v1.1.0."
echo "Next: git remote add origin <url> && git push -u origin main --tags"
