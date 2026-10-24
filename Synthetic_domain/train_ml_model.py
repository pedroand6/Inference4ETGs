"""
Train the SFH posterior estimator.

Usage (from Synthetic_domain/):
    python train_ml_model.py --embedding mlp
    python train_ml_model.py --embedding cnn
Each run writes ml_runs/<embedding>/{model.pt, config.json, history.csv}. The training set is noise-free;
every epoch draws a fresh brightness, S-PLUS noise and observed redshift for every model (ml_data).
"""

import argparse
import os
import subprocess
import time

import numpy as np
import pandas as pd
import torch

import ml_data
import ml_networks

TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"
CATALOGUE_FILE = "dr3_corrected.csv"
MIN_PROB_GALAXY = 0.8
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}


def to_tensors(observations):
    """Keep valid rows and convert flat, banded, targets to float32 tensors (not yet standardised)."""
    valid = observations["valid"]
    tensors = {}
    for name in ("flat", "banded", "targets"):
        tensors[name] = torch.as_tensor(observations[name][valid], dtype=torch.float32)
    return tensors


def standardise(tensors, standardisers):
    """Standardised copy of a tensor dict (masked colours become 0)."""
    flat = tensors["flat"].numpy()
    banded = tensors["banded"].numpy()
    flat_mask, banded_mask = ml_data.input_masks(flat, banded)
    return {
        "flat": torch.as_tensor(standardisers["flat"].transform(flat, flat_mask), dtype=torch.float32),
        "banded": torch.as_tensor(standardisers["banded"].transform(banded, banded_mask), dtype=torch.float32),
        "targets": torch.as_tensor(standardisers["targets"].transform(tensors["targets"].numpy()),
                                   dtype=torch.float32),
    }


def git_commit():
    """Current git commit hash, or 'unknown' outside a repository (e.g. inside a job)."""
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main():
    """Parse arguments, build the data pipeline, train one embedding and save the run."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--embedding", choices=sorted(ml_networks.EMBEDDINGS), required=True)
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--output-dir", default="ml_runs")
    parser.add_argument("--epochs", type=int, default=ml_networks.DEFAULT_CONFIG["epochs"])
    parser.add_argument("--max-train-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--max-validation-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--seed", type=int, default=ml_networks.DEFAULT_CONFIG["seed"])
    arguments = parser.parse_args()

    device = arguments.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    # 1. Data
    training = pd.read_parquet(os.path.join(arguments.data_dir, TRAINING_FILE))
    if arguments.max_train_rows is not None:
        training = training.iloc[:arguments.max_train_rows]
    validation = pd.read_parquet(os.path.join(arguments.data_dir, VALIDATION_FILE))
    if arguments.max_validation_rows is not None:
        validation = validation.iloc[:arguments.max_validation_rows]
    noise_table = pd.read_csv(os.path.join(arguments.data_dir, NOISE_MODEL_FILE))
    catalogue = pd.read_csv(os.path.join(arguments.data_dir, CATALOGUE_FILE), usecols=["PROB_GAL", "r_auto", "zml"])
    catalogue = catalogue[catalogue["PROB_GAL"] > MIN_PROB_GALAXY]
    catalogue_r = catalogue["r_auto"].to_numpy()
    catalogue_z = catalogue["zml"].to_numpy()

    # 2. Standardisation statistics from one simulated epoch (seeded)
    first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                           np.random.default_rng(arguments.seed)))
    flat_mask, banded_mask = ml_data.input_masks(first_epoch["flat"].numpy(), first_epoch["banded"].numpy())
    standardisers = {
        "flat": ml_data.Standardiser.fit(first_epoch["flat"].numpy(), flat_mask),
        "banded": ml_data.Standardiser.fit(first_epoch["banded"].numpy().reshape(-1, ml_data.N_BAND_CHANNELS),
                                           banded_mask.reshape(-1, ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser.fit(first_epoch["targets"].numpy()),
    }

    # 3. Fixed validation copies: spectroscopic and photometric redshifts
    validation_sets = []
    for mode, seed in VALIDATION_SEEDS.items():
        observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed))
        validation_sets.append(standardise(to_tensors(observations), standardisers))

    def make_epoch(epoch_index):
        """Fresh noisy realisation for this epoch (seeded by run seed and epoch)."""
        rng = np.random.default_rng([arguments.seed, epoch_index])
        observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng)
        return standardise(to_tensors(observations), standardisers)

    # 4. Train and save
    config = dict(ml_networks.DEFAULT_CONFIG)
    config.update(embedding=arguments.embedding, epochs=arguments.epochs, seed=arguments.seed,
                  training_commit=git_commit(), torch_version=torch.__version__,
                  n_training_models=len(training), device=device)
    estimator = ml_networks.build_estimator(config)
    start = time.time()
    history = ml_networks.train_estimator(estimator, make_epoch, validation_sets, config, device,
                                          log=lambda message: print(message, flush=True))
    config["training_minutes"] = (time.time() - start) / 60
    output_directory = os.path.join(arguments.output_dir, arguments.embedding)
    ml_networks.save_run(output_directory, estimator, config, history, standardisers)
    print(f"saved {output_directory} after {config['training_minutes']:.1f} min", flush=True)


if __name__ == "__main__":
    main()
