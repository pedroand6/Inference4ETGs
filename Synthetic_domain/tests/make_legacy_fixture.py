"""Freeze the prior of the original training_set.ipynb (commit c04268e) as a test fixture.

Run once from Synthetic_domain/:  ../.venv/bin/python tests/make_legacy_fixture.py
The refactored module must reproduce these components exactly.
"""

import json
import subprocess

import numpy as np

LEGACY_COMMIT = "c04268e"
NOTEBOOK_PATH = "Synthetic_domain/training_set.ipynb"
CONSTANTS_CELL = 1      # imports and prior constants
PRIOR_CELL = 3          # bin_edges, mass_fractions, components
NEW_Z_RANGE = (0.002, 0.7)   # the only prior change in the new design
N_VECTORS = 30
SEED = 123
OUTPUT_FILE = "tests/legacy_components.json"


# 1. Read the original notebook from git
notebook_text = subprocess.run(
    ["git", "show", f"{LEGACY_COMMIT}:{NOTEBOOK_PATH}"],
    capture_output=True, text=True, check=True,
).stdout
notebook = json.loads(notebook_text)

# 2. Run its constants and prior cells in a private namespace
namespace = {}
for cell_index in (CONSTANTS_CELL, PRIOR_CELL):
    source = "".join(notebook["cells"][cell_index]["source"])
    exec(source, namespace)
namespace["Z_RANGE"] = NEW_Z_RANGE

# 3. Evaluate the legacy prior on fixed unit-cube vectors
rng = np.random.default_rng(SEED)
unit_vectors = rng.random((N_VECTORS, len(namespace["DIMS"])))
cases = []
for u in unit_vectors:
    legacy_components, legacy_params = namespace["components"](u)
    cases.append({"u": u.tolist(), "components": legacy_components, "params": legacy_params})

# 4. Save (numpy floats are Python floats; anything else is cast to float)
with open(OUTPUT_FILE, "w") as output:
    json.dump(cases, output, default=float, indent=1)
print(f"{len(cases)} legacy cases written to {OUTPUT_FILE}")
