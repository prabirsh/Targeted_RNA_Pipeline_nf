# Keep every download, environment, cache and test run inside this repository.
# Usage: source tests/local_env.sh
# Afterwards, conda and Nextflow write only to .local-test/, which git ignores.

_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export RNA_PANEL_LOCAL="$_repo/.local-test"
mkdir -p "$RNA_PANEL_LOCAL"/{envs,pkgs,nxf,conda,work,tmp}

export CONDA_PKGS_DIRS="$RNA_PANEL_LOCAL/pkgs"
export CONDA_ENVS_DIRS="$RNA_PANEL_LOCAL/envs"
export CONDA_REGISTER_ENVS=false
export NXF_HOME="$RNA_PANEL_LOCAL/nxf"
export NXF_CONDA_CACHEDIR="$RNA_PANEL_LOCAL/conda"
export NXF_WORK="$RNA_PANEL_LOCAL/work"
export NXF_TEMP="$RNA_PANEL_LOCAL/tmp"
export TEST_DATA_DIR="$RNA_PANEL_LOCAL/data"
export PIP_CACHE_DIR="$RNA_PANEL_LOCAL/pip-cache"

unset _repo
