"""Tests for synthetic_model.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest -v"""

import json

import numpy as np
import pandas as pd
import pytest

import synthetic_model as sm

LEGACY_FIXTURE = "tests/legacy_components.json"


# =============================================================================
# Prior
# =============================================================================
def test_components_match_legacy_notebook():
    with open(LEGACY_FIXTURE) as fixture:
        cases = json.load(fixture)

    for case in cases:
        new_components, new_params = sm.components(np.array(case["u"]))
        legacy_components = case["components"]

        assert set(new_components) == set(legacy_components)
        assert np.isclose(new_components["redshift"], legacy_components["redshift"], rtol=1e-12)
        assert np.isclose(new_components["t_bc"], legacy_components["t_bc"], rtol=1e-12)
        for key in ("Av", "eta", "n"):
            assert np.isclose(new_components["dust"][key], legacy_components["dust"][key], rtol=1e-12)
        assert np.isclose(new_components["nebular"]["logU"], legacy_components["nebular"]["logU"], rtol=1e-12)

        for name in legacy_components:
            if not name.startswith("constant"):
                continue
            for key in ("age_min", "age_max", "massformed", "metallicity"):
                assert np.isclose(new_components[name][key], legacy_components[name][key], rtol=1e-12), (name, key)

        for key, legacy_value in case["params"].items():
            assert np.isclose(new_params[key], legacy_value, rtol=1e-12), key


@pytest.mark.parametrize("branch, sfh_type_draw", [(0, 0.25), (1, 0.75), (2, 0.99)])
def test_each_branch_gives_fractions_that_sum_to_one(branch, sfh_type_draw):
    rng = np.random.default_rng(1)
    u = rng.random(len(sm.DIMS))
    u[sm.DIMS.index("sfh_type")] = sfh_type_draw

    assert sm.choose_sfh_branch(sfh_type_draw) == branch
    bagpipes_components, drawn = sm.components(u)
    fractions, metallicities, edges = sm.bins_from_components(bagpipes_components)

    assert drawn["sfh_type"] == branch
    assert np.isclose(fractions.sum(), 1.0, atol=1e-4)    # MIN_FRACTION floor adds at most 16e-6
    assert len(fractions) == len(edges) - 1


def test_old_bins_are_dropped_at_high_redshift():
    low_redshift_edges = sm.bin_edges(0.01)
    high_redshift_edges = sm.bin_edges(0.69)

    # The 14.2 Gyr bin centre is older than the Universe (13.37 Gyr at most here), so the
    # prior never has more than N_BINS - 1 bins; its labels are always padding.
    assert len(low_redshift_edges) - 1 == sm.N_BINS - 1
    assert len(high_redshift_edges) - 1 < sm.N_BINS - 1
    universe_age = sm.COSMOLOGY.age(0.69).value
    assert high_redshift_edges[-1] < universe_age


def test_metallicity_rises_from_oldest_to_youngest_bin():
    fractions = np.full(sm.N_BINS, 1 / sm.N_BINS)
    metallicities = sm.bin_metallicities(fractions, z_final=1.0)

    # Bins run youngest first, so metallicity must fall along the array
    assert np.all(np.diff(metallicities) < 0)
    assert metallicities[0] < 1.0
    assert metallicities[-1] > sm.ZINI_FRACTION * 1.0


def test_fixed_redshift_overrides_the_drawn_redshift():
    u = np.full(len(sm.DIMS), 0.5)
    bagpipes_components, drawn = sm.components(u, fixed_redshift=0.123)

    assert bagpipes_components["redshift"] == 0.123
    assert drawn["z"] == 0.123


# =============================================================================
# Labels and magnitudes
# =============================================================================
def test_bin_labels_pad_missing_old_bins():
    u = np.full(len(sm.DIMS), 0.5)
    bagpipes_components, _ = sm.components(u, fixed_redshift=0.69)
    fractions, metallicities, edges = sm.bins_from_components(bagpipes_components)
    labels = sm.bin_labels(fractions, metallicities, edges)

    n_existing = len(fractions)
    assert n_existing < sm.N_BINS
    total = 0.0
    for bin_index in range(sm.N_BINS):
        suffix = f"{bin_index:02d}"
        total += labels[f"frac_bin_{suffix}"]
        if bin_index < n_existing:
            assert np.isfinite(labels[f"logsfr_bin_{suffix}"])
            assert np.isfinite(labels[f"zmet_bin_{suffix}"])
        else:
            assert labels[f"frac_bin_{suffix}"] == 0.0
            assert np.isnan(labels[f"logsfr_bin_{suffix}"])
            assert np.isnan(labels[f"zmet_bin_{suffix}"])
    assert np.isclose(total, 1.0)
    assert list(labels) == sm.BIN_LABEL_COLUMNS


def test_mass_weighted_metallicity_of_equal_bins_is_the_mean():
    fractions = np.array([0.5, 0.5])
    metallicities = np.array([1.0, 2.0])
    assert np.isclose(sm.mass_weighted_metallicity(fractions, metallicities), 1.5)


def test_ab_magnitudes_are_nan_when_any_flux_is_not_positive():
    flux = np.full(len(sm.BANDS), 1e-17)
    assert np.all(np.isfinite(sm.ab_magnitudes(flux)))
    flux[3] = 0.0
    assert np.all(np.isnan(sm.ab_magnitudes(flux)))


def test_simulate_models_returns_labels_and_magnitudes():
    unit_vectors = np.random.default_rng(2).random((4, len(sm.DIMS)))
    models = sm.simulate_models(unit_vectors)

    assert len(models) == 4
    for column in ["z", "sfh_type", *sm.SUMMARY_LABEL_COLUMNS, *sm.BIN_LABEL_COLUMNS, *sm.BANDS]:
        assert column in models.columns, column
    fraction_columns = []
    for bin_index in range(sm.N_BINS):
        fraction_columns.append(f"frac_bin_{bin_index:02d}")
    assert np.allclose(models[fraction_columns].sum(axis=1), 1.0)
    assert np.isfinite(models[sm.BANDS].to_numpy()).all()
