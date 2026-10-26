"""Smoke test of the training CLI on a tiny subset (CPU). Run from Synthetic_domain/."""

import json
import subprocess
import sys

import pandas as pd
import pyarrow.parquet as pq
import pytest

import ml_data

DATA_FILES = ("training_set.parquet", "synthetic_validation.parquet")


def sets_are_regenerated():
    """True when both parquet sets hold the delayed-tau labels (their schema has an 'age' column)."""
    for path in DATA_FILES:
        try:
            if "age" not in pq.read_schema(path).names:
                return False
        except OSError:
            return False
    return True


@pytest.mark.skipif(not sets_are_regenerated(), reason="synthetic sets not yet regenerated on the PC")
def test_training_cli_writes_a_run(tmp_path):
    command = [sys.executable, "train_ml_model.py", "--embedding", "mlp", "--output-dir", str(tmp_path),
               "--epochs", "2", "--max-train-rows", "3000", "--max-validation-rows", "500", "--device", "cpu"]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=900)
    assert completed.returncode == 0, completed.stderr[-2000:]
    run_directory = tmp_path / "mlp"
    for name in ("model.pt", "config.json", "history.csv"):
        assert (run_directory / name).exists(), name
    history = pd.read_csv(run_directory / "history.csv")
    assert len(history) == 2
    config = json.loads((run_directory / "config.json").read_text())
    assert config["target_names"] == ml_data.TARGET_NAMES
