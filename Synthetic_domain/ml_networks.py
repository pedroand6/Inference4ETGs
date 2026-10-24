"""
Posterior estimator for the SFH targets: an embedding network (MLP or 1-D CNN) summarises the inputs,
and a zuko neural spline flow models p(targets | summary). Training maximises log p(true targets).
"""

import json
import math
import os

import numpy as np
import torch
from torch import nn
import zuko

import ml_data


# =============================================================================
# 1. Configuration
# =============================================================================
DEFAULT_CONFIG = {
    "embedding": "mlp",              # "mlp" or "cnn"
    "n_targets": ml_data.N_TARGETS,
    "summary_size": 64,              # length of the embedding output
    "hidden_size": 256,              # width of the dense layers
    "mlp_layers": 3,
    "cnn_channels": 64,
    "cnn_kernel": 3,
    "cnn_layers": 3,
    "head_layers": 2,                # dense layers after the CNN pooling
    "flow_transforms": 8,
    "flow_hidden": (256, 256),
    "learning_rate": 1e-3,
    "batch_size": 4096,
    "epochs": 200,
    "patience": 15,                  # epochs without validation improvement before stopping
    "seed": 0,
}
MODEL_FILE = "model.pt"
CONFIG_FILE = "config.json"
HISTORY_FILE = "history.csv"


# =============================================================================
# 2. Networks
# =============================================================================
def dense_stack(input_size, hidden_size, n_layers, output_size):
    """n_layers of Linear + SiLU, then a Linear to output_size."""
    layers = []
    current_size = input_size
    for _ in range(n_layers):
        layers.append(nn.Linear(current_size, hidden_size))
        layers.append(nn.SiLU())
        current_size = hidden_size
    layers.append(nn.Linear(current_size, output_size))
    return nn.Sequential(*layers)


class MLPEmbedding(nn.Module):
    """Dense network on the 37 flat inputs (colours, errors, flags, scalars)."""

    def __init__(self, config):
        """Inputs: config with hidden_size, mlp_layers, summary_size."""
        super().__init__()
        self.network = dense_stack(ml_data.N_INPUTS, config["hidden_size"], config["mlp_layers"],
                                   config["summary_size"])

    def forward(self, flat, banded):
        """Inputs: flat (batch, 37); banded is ignored. Returns: summary (batch, summary_size)."""
        return self.network(flat)


class CNNEmbedding(nn.Module):
    """1-D convolutions over the 12 bands ordered by wavelength, then a dense head with the scalars."""

    def __init__(self, config):
        """Inputs: config with cnn_channels, cnn_kernel, cnn_layers, hidden_size, head_layers, summary_size."""
        super().__init__()
        layers = []
        in_channels = ml_data.N_BAND_CHANNELS
        for _ in range(config["cnn_layers"]):
            layers.append(nn.Conv1d(in_channels, config["cnn_channels"], config["cnn_kernel"],
                                    padding=config["cnn_kernel"] // 2))
            layers.append(nn.SiLU())
            in_channels = config["cnn_channels"]
        self.convolutions = nn.Sequential(*layers)
        n_scalars = ml_data.SCALAR_INPUT_SLICE.stop - ml_data.SCALAR_INPUT_SLICE.start
        self.head = dense_stack(config["cnn_channels"] + n_scalars, config["hidden_size"],
                                config["head_layers"], config["summary_size"])

    def forward(self, flat, banded):
        """Inputs: flat (batch, 37) for the scalars; banded (batch, 12, 4). Returns: summary (batch, summary_size)."""
        channels_first = banded.transpose(1, 2)              # (batch, 4, 12): Conv1d wants channels first
        features = self.convolutions(channels_first).mean(dim=2)
        scalars = flat[:, ml_data.SCALAR_INPUT_SLICE]
        return self.head(torch.cat([features, scalars], dim=1))


EMBEDDINGS = {"mlp": MLPEmbedding, "cnn": CNNEmbedding}


class PosteriorEstimator(nn.Module):
    """Embedding + conditional neural spline flow over the standardised targets."""

    def __init__(self, config):
        """Inputs: config (see DEFAULT_CONFIG)."""
        super().__init__()
        self.embedding = EMBEDDINGS[config["embedding"]](config)
        self.flow = zuko.flows.NSF(features=config["n_targets"], context=config["summary_size"],
                                   transforms=config["flow_transforms"],
                                   hidden_features=tuple(config["flow_hidden"]))

    def log_prob(self, theta, flat, banded):
        """log p(theta | inputs) per row. Inputs: tensors (batch, n_targets), (batch, 37), (batch, 12, 4)."""
        summary = self.embedding(flat, banded)
        return self.flow(summary).log_prob(theta)

    def sample(self, flat, banded, n_samples):
        """n_samples posterior draws per row. Returns: tensor (batch, n_samples, n_targets)."""
        summary = self.embedding(flat, banded)
        draws = self.flow(summary).sample((n_samples,))       # (n_samples, batch, n_targets)
        return draws.transpose(0, 1)


def build_estimator(config):
    """A fresh PosteriorEstimator with weights seeded by config['seed']."""
    torch.manual_seed(config["seed"])
    return PosteriorEstimator(config)


# =============================================================================
# 3. Training
# =============================================================================
def mean_loss(estimator, data_set, batch_size, device):
    """Mean negative log-probability over a data set (no gradients)."""
    estimator.eval()
    total = 0.0
    n_rows = len(data_set["targets"])
    with torch.no_grad():
        for start in range(0, n_rows, batch_size):
            stop = start + batch_size
            theta = data_set["targets"][start:stop].to(device)
            flat = data_set["flat"][start:stop].to(device)
            banded = data_set["banded"][start:stop].to(device)
            total += -estimator.log_prob(theta, flat, banded).sum().item()
    return total / n_rows


def train_estimator(estimator, make_epoch, validation_sets, config, device, log=print):
    """Train with AdamW + cosine decay and early stopping; keeps the best weights.

    Inputs: estimator; make_epoch(epoch_index) -> dict of standardised tensors {flat, banded, targets}
            (a fresh noisy realisation each epoch); validation_sets: list of such dicts (fixed);
            config: batch_size, learning_rate, epochs, patience, seed; device: "cpu" or "cuda";
            log: function called with one progress line per epoch.
    Returns: history, a list of dicts {epoch, train_loss, validation_loss, learning_rate}.
    """
    estimator.to(device)
    optimiser = torch.optim.AdamW(estimator.parameters(), lr=config["learning_rate"])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, T_max=config["epochs"])
    shuffle_generator = torch.Generator().manual_seed(config["seed"])

    best_loss = math.inf
    best_state = None
    epochs_without_improvement = 0
    history = []
    for epoch_index in range(config["epochs"]):
        # 1. One pass over a fresh noisy realisation of the training set
        epoch = make_epoch(epoch_index)
        n_rows = len(epoch["targets"])
        order = torch.randperm(n_rows, generator=shuffle_generator)
        estimator.train()
        running = 0.0
        for start in range(0, n_rows, config["batch_size"]):
            batch = order[start:start + config["batch_size"]]
            theta = epoch["targets"][batch].to(device)
            flat = epoch["flat"][batch].to(device)
            banded = epoch["banded"][batch].to(device)
            loss = -estimator.log_prob(theta, flat, banded).mean()
            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            running += loss.item() * len(batch)
        scheduler.step()

        # 2. Validation loss: mean over the validation sets (e.g. spectroscopic and photometric z)
        validation_losses = []
        for validation_set in validation_sets:
            validation_losses.append(mean_loss(estimator, validation_set, config["batch_size"], device))
        validation_loss = float(np.mean(validation_losses))
        history.append({"epoch": epoch_index, "train_loss": running / n_rows,
                        "validation_loss": validation_loss, "learning_rate": scheduler.get_last_lr()[0]})
        log(f"epoch {epoch_index}: train {running / n_rows:.4f}  validation {validation_loss:.4f}")

        # 3. Early stopping on the validation loss
        if validation_loss < best_loss:
            best_loss = validation_loss
            best_state = {name: value.detach().cpu().clone() for name, value in estimator.state_dict().items()}
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= config["patience"]:
                log(f"early stop after epoch {epoch_index}: best validation {best_loss:.4f}")
                break

    estimator.load_state_dict(best_state)
    return history


# =============================================================================
# 4. Checkpoints and sampling
# =============================================================================
def save_run(directory, estimator, config, history, standardisers):
    """Write model.pt (weights), config.json (config + standardisation statistics) and history.csv."""
    os.makedirs(directory, exist_ok=True)
    torch.save(estimator.state_dict(), os.path.join(directory, MODEL_FILE))
    stored = dict(config)
    stored["flow_hidden"] = list(config["flow_hidden"])
    stored["standardisers"] = {}
    for name, standardiser in standardisers.items():
        stored["standardisers"][name] = standardiser.to_dict()
    with open(os.path.join(directory, CONFIG_FILE), "w") as output:
        json.dump(stored, output, indent=1)
    with open(os.path.join(directory, HISTORY_FILE), "w") as output:
        output.write("epoch,train_loss,validation_loss,learning_rate\n")
        for row in history:
            output.write(f"{row['epoch']},{row['train_loss']},{row['validation_loss']},{row.get('learning_rate', '')}\n")


def load_run(directory, device="cpu"):
    """Read a run written by save_run. Returns: (estimator on device, config, standardisers dict)."""
    with open(os.path.join(directory, CONFIG_FILE)) as source:
        stored = json.load(source)
    standardisers = {}
    for name, statistics in stored.pop("standardisers").items():
        standardisers[name] = ml_data.Standardiser.from_dict(statistics)
    estimator = PosteriorEstimator(stored)
    weights = torch.load(os.path.join(directory, MODEL_FILE), map_location=device, weights_only=True)
    estimator.load_state_dict(weights)
    estimator.to(device)
    estimator.eval()
    return estimator, stored, standardisers


def posterior_samples(estimator, standardisers, flat, banded, n_samples, seed, device="cpu", batch_size=512):
    """Posterior draws in physical target units for raw (unstandardised) inputs.

    Inputs: estimator; standardisers {flat, banded, targets}; flat (n, 37) and banded (n, 12, 4) as
            ml_data.colour_inputs returns them; n_samples per galaxy; seed (torch); device; batch_size (galaxies).
    Returns: array (n, n_samples, n_targets).
    """
    flat_mask, banded_mask = ml_data.input_masks(flat, banded)
    flat_standard = standardisers["flat"].transform(flat, flat_mask)
    banded_standard = standardisers["banded"].transform(banded, banded_mask)

    torch.manual_seed(seed)
    estimator.to(device)
    estimator.eval()
    pieces = []
    with torch.no_grad():
        for start in range(0, len(flat), batch_size):
            flat_batch = torch.as_tensor(flat_standard[start:start + batch_size], dtype=torch.float32, device=device)
            banded_batch = torch.as_tensor(banded_standard[start:start + batch_size], dtype=torch.float32,
                                           device=device)
            draws = estimator.sample(flat_batch, banded_batch, n_samples).cpu().numpy()
            pieces.append(standardisers["targets"].inverse_transform(draws))
    return np.concatenate(pieces, axis=0)
