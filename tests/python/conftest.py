"""Shared helpers: run bin/ scripts the way the pipeline does, as subprocesses."""
import csv
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"


@pytest.fixture
def run_script():
    def _run(name, *args):
        return subprocess.run(
            [sys.executable, str(BIN / name), *map(str, args)],
            check=False,
            capture_output=True,
            text=True,
        )

    return _run


@pytest.fixture
def read_tsv():
    def _read(path):
        with open(path, newline="") as fh:
            csv.field_size_limit(sys.maxsize)
            return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))

    return _read
