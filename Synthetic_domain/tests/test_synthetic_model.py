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
