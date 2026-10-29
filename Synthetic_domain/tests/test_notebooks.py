"""Static checks of the three notebooks (they are executed only on the PC, so a broken name shows up there late).

Run from Synthetic_domain/.
"""

import importlib
import json
import re

import pytest

NOTEBOOKS = ["training_set.ipynb", "spectroscopic_reference.ipynb", "ml_model.ipynb"]
PROJECT_MODULES = ["synthetic_model", "sfh_models", "ml_data", "ml_networks", "ml_evaluation",
                   "spectral_validation", "train_ml_model"]
IMPORT_PATTERN = re.compile(r"^\s*import (\w+)(?: as (\w+))?\s*$", re.MULTILINE)


def notebook_code(name):
    """The source of every code cell of a notebook, joined into one string."""
    notebook = json.load(open(name))
    sources = []
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            sources.append("".join(cell["source"]))
    return "\n".join(sources)


@pytest.mark.parametrize("name", NOTEBOOKS)
def test_every_module_attribute_a_notebook_uses_exists(name):
    """Why: removing a module name (e.g. the delayed-tau aliases of Task 11) broke notebook 2 only on the PC."""
    code = notebook_code(name)
    missing = []
    for match in IMPORT_PATTERN.finditer(code):
        module_name = match.group(1)
        alias = match.group(2) or module_name
        if module_name not in PROJECT_MODULES:
            continue
        module = importlib.import_module(module_name)
        for attribute in sorted(set(re.findall(rf"\b{re.escape(alias)}\.(\w+)", code))):
            if not hasattr(module, attribute):
                missing.append(f"{alias}.{attribute}")
    assert missing == [], f"{name} uses names its modules do not define: {missing}"
