import mrv_tools
from os import path
import json

# Loadng and creating periodic table
periodic_table = json.load(
    open(path.join(path.split(mrv_tools.__file__)[0], "periodic_table/PeriodicTableJSON.json"), "r")
)["elements"]
periodic_table = {e["symbol"]: e for e in periodic_table}

# Transition element:
# An element whose atom has an incomplete d sub-shell,
# or which can give rise to cations with an incomplete d sub-shell.
transition_metals = [
    periodic_table[e]["symbol"] for e in periodic_table
    if periodic_table[e]["category"] == "transition metal"
]


def tm_label(atom):
    """
    Transition metal label.
    e.g. Fe(III)
    """
    roman = {
        1: "I",
        2: "II",
        3: "III",
        4: "IV",
        5: "V",
        6: "VI",
        7: "VII",
        8: "VIII",
        9: "IX",
        10: "X"
    }
    if (atom.formal_charge is None) or (atom.formal_charge == 0):
        return atom.element_type

    try:
        return f"{atom.element_type}({roman[atom.formal_charge]})"
    except KeyError:
        return f"{atom.element_type}({atom.formal_charge})"
