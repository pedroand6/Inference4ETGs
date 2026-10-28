"""
Train the SFH posterior estimator of one SFH model.

Usage (from Synthetic_domain/):
    python train_ml_model.py                         # the MLP of the delayed-tau model
    python train_ml_model.py --sfh-model dblplaw     # the MLP of the double-power-law model
    python train_ml_model.py --embedding cnn         # kept for comparison; not part of the pipeline
Each run reads its model's sets. sfh_models.model_file names them: training_set.parquet and
synthetic_validation.parquet for delayed-tau, training_set_dblplaw.parquet and
synthetic_validation_dblplaw.parquet for the double power law.
Each run writes ml_runs/<embedding>/ (delayed-tau) or ml_runs/<embedding>_dblplaw/, with model.pt,
config.json and history.csv. config.json records the sfh_model and the ml_data.target_names the flow was
trained for.
The training set is noise-free. Every epoch draws a fresh brightness, S-PLUS noise and observed redshift for
every model (ml_data).
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
import sfh_models

TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"
CATALOGUE_FILE = "dr3_corrected.csv"
MIN_PROB_GALAXY = 0.8
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}
DEFAULT_EMBEDDING = "mlp"        # the CNN lost the 16-bin comparison, so only the MLP is trained
DEFAULT_SFH_MODEL = "delayed"    # the model of the existing run, ml_runs/mlp


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


def data_file(data_dir, name, sfh_model):
    """Path of one of an SFH model's data sets: name with the model's suffix (sfh_models.model_file), in data_dir."""
    return os.path.join(data_dir, sfh_models.model_file(name, sfh_model))


def run_directory(output_dir, embedding, sfh_model):
    """Folder of a run: <output_dir>/<embedding> for delayed-tau (as ml_runs/mlp), with the model's suffix otherwise."""
    return sfh_models.model_file(os.path.join(output_dir, embedding), sfh_model)


def parse_arguments(argv=None):
    """Command-line arguments (argv: list of strings, default sys.argv[1:])."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sfh-model", choices=sorted(sfh_models.MODELS), default=DEFAULT_SFH_MODEL)
    parser.add_argument("--embedding", choices=sorted(ml_networks.EMBEDDINGS), default=DEFAULT_EMBEDDING)
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--output-dir", default="ml_runs")
    parser.add_argument("--epochs", type=int, default=ml_networks.DEFAULT_CONFIG["epochs"])
    parser.add_argument("--max-train-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--max-validation-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--seed", type=int, default=ml_networks.DEFAULT_CONFIG["seed"])
    return parser.parse_args(argv)


def main(argv=None):
    """Parse arguments, build the data pipeline of one SFH model, train one embedding and save the run."""
    arguments = parse_arguments(argv)
    sfh_model = arguments.sfh_model

    device = arguments.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}, SFH model: {sfh_model}", flush=True)

    # 1. Data
    training = pd.read_parquet(data_file(arguments.data_dir, TRAINING_FILE, sfh_model))
    if arguments.max_train_rows is not None:
        training = training.iloc[:arguments.max_train_rows]
    validation = pd.read_parquet(data_file(arguments.data_dir, VALIDATION_FILE, sfh_model))
    if arguments.max_validation_rows is not None:
        validation = validation.iloc[:arguments.max_validation_rows]
    noise_table = pd.read_csv(os.path.join(arguments.data_dir, NOISE_MODEL_FILE))
    catalogue = pd.read_csv(os.path.join(arguments.data_dir, CATALOGUE_FILE), usecols=["PROB_GAL", "r_auto", "zml"])
    catalogue = catalogue[catalogue["PROB_GAL"] > MIN_PROB_GALAXY]
    catalogue_r = catalogue["r_auto"].to_numpy()
    catalogue_z = catalogue["zml"].to_numpy()

    # 2. Standardisation statistics from one simulated epoch (seeded)
    first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                           np.random.default_rng(arguments.seed), sfh_model))
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
        observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed), sfh_model)
        validation_sets.append(standardise(to_tensors(observations), standardisers))

    def make_epoch(epoch_index):
        """Fresh noisy realisation for this epoch (seeded by run seed and epoch)."""
        rng = np.random.default_rng([arguments.seed, epoch_index])
        observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng,
                                                     sfh_model)
        return standardise(to_tensors(observations), standardisers)

    # 4. Train and save
    config = ml_networks.model_config(sfh_model)
    config.update(embedding=arguments.embedding, epochs=arguments.epochs, seed=arguments.seed,
                  training_commit=git_commit(), torch_version=torch.__version__, n_training_models=len(training),
                  device=device)
    estimator = ml_networks.build_estimator(config)
    start = time.time()
    history = ml_networks.train_estimator(estimator, make_epoch, validation_sets, config, device,
                                          log=lambda message: print(message, flush=True))
    config["training_minutes"] = (time.time() - start) / 60
    output_directory = run_directory(arguments.output_dir, arguments.embedding, sfh_model)
    ml_networks.save_run(output_directory, estimator, config, history, standardisers)
    print(f"saved {output_directory} after {config['training_minutes']:.1f} min", flush=True)


if __name__ == "__main__":
    main()
