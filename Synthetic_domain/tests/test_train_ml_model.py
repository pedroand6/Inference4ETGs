"""Tests of the training CLI: paths per SFH model, and a smoke run on a tiny subset (CPU).

Run from Synthetic_domain/.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pyarrow.parquet as pq
import pytest

import ml_data
import ml_networks
import sfh_models
import train_ml_model

DATA_FILES = ("training_set.parquet", "synthetic_validation.parquet")
SET_COLUMNS = {"delayed": "age", "dblplaw": "turnover_lookback"}   # a column only that model's sets hold
OTHER_MODEL = {"delayed": "dblplaw", "dblplaw": "delayed"}


def sets_are_generated(sfh_model):
    """True when both sets of sfh_model exist and hold its labels (the schema has SET_COLUMNS[sfh_model])."""
    for name in DATA_FILES:
        path = sfh_models.model_file(name, sfh_model)
        try:
            if SET_COLUMNS[sfh_model] not in pq.read_schema(path).names:
                return False
        except OSError:
            return False
    return True


def test_data_files_and_run_directories_follow_the_model_suffix():
    assert train_ml_model.data_file("data", "training_set.parquet", "delayed") == os.path.join(
        "data", "training_set.parquet")
    assert train_ml_model.data_file("data", "training_set.parquet", "dblplaw") == os.path.join(
        "data", "training_set_dblplaw.parquet")
    assert train_ml_model.run_directory("ml_runs", "mlp", "delayed") == os.path.join("ml_runs", "mlp")
    assert train_ml_model.run_directory("ml_runs", "mlp", "dblplaw") == os.path.join("ml_runs", "mlp_dblplaw")


def test_sfh_model_argument_defaults_to_delayed():
    assert train_ml_model.parse_arguments([]).sfh_model == "delayed"
    assert train_ml_model.parse_arguments(["--sfh-model", "dblplaw"]).sfh_model == "dblplaw"
    with pytest.raises(SystemExit):
        train_ml_model.parse_arguments(["--sfh-model", "exponential"])


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_training_cli_writes_a_run(tmp_path, sfh_model):
    if not sets_are_generated(sfh_model):
        pytest.skip(f"the {sfh_model} sets are not yet generated on the PC")
    command = [sys.executable, "train_ml_model.py", "--sfh-model", sfh_model, "--embedding", "mlp",
               "--output-dir", str(tmp_path), "--epochs", "2", "--max-train-rows", "3000",
               "--max-validation-rows", "500", "--device", "cpu"]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=900)
    assert completed.returncode == 0, completed.stderr[-2000:]
    run_directory = train_ml_model.run_directory(str(tmp_path), "mlp", sfh_model)
    for name in ("model.pt", "config.json", "history.csv"):
        assert os.path.exists(os.path.join(run_directory, name)), name
    history = pd.read_csv(os.path.join(run_directory, "history.csv"))
    assert len(history) == 2
    with open(os.path.join(run_directory, "config.json")) as source:
        config = json.load(source)
    assert config["sfh_model"] == sfh_model
    assert config["target_names"] == ml_data.target_names(sfh_model)
    assert config["n_targets"] == len(ml_data.target_names(sfh_model))
    ml_networks.load_run(run_directory, sfh_model)
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(run_directory, OTHER_MODEL[sfh_model])
