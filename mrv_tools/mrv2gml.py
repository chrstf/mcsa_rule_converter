from typing import List
from mrv_tools import MrvData, Atom, Bond


def atom_label(atom: Atom):
    if atom.mrv_alias is not None and atom.mrv_alias != "0":
        return atom.convert_alias()
    if atom.mrv_query_props is not None and atom.mrv_query_props != "0":
        return atom.convert_query_props()

    label = atom.element_type
    chg_label = ""

    if atom.formal_charge is not None:
        charge = atom.formal_charge
        if charge == 1:
            chg_label = "+"
        elif charge > 1:
            chg_label = f"{charge}+"
        elif charge == -1:
            chg_label = "-"
        elif charge < -1:
            chg_label = f"{(charge * -1)}-"

    return f"{label}{chg_label}"


def bond_label(bond: Bond):
    bond_type_map = {
        1.: "-",
        2.: "=",
        3.: "#",
        1.5: ":",
        "cxn:coord": ">"
    }
    if bond.order is None:
        return bond_type_map[bond.convention]

    return bond_type_map[bond.order]


def atom2gml(atom: Atom):
    assert(atom.id[0] == "a")
    return f'node [ id {atom.id[1:]} label "{atom_label(atom)}" ]'


def bond2gml(bond: Bond):
    src, tar = bond.atom_refs_2
    return f'edge [ source {src.id[1:]} target {tar.id[1:]} label "{bond_label(bond)}"]'


def convert(mrv_reactants: MrvData, mrv_products: MrvData):
    nodes_left = frozenset(map(atom2gml, mrv_reactants.flatten_atoms()))
    edges_left = frozenset(map(bond2gml, mrv_reactants.flatten_bonds()))

    nodes_right = frozenset(map(atom2gml, mrv_products.flatten_atoms()))
    edges_right = frozenset(map(bond2gml, mrv_products.flatten_bonds()))

    # Compute maximal overlapping context
    nodes_context = nodes_left.intersection(nodes_right)
    nodes_left = nodes_left - nodes_context
    nodes_right = nodes_right - nodes_context

    edges_context = edges_left.intersection(edges_right)
    edges_left = edges_left - edges_context
    edges_right = edges_right - edges_context

    # convert to GML
    out: List[str] = ["rule ["]

    span = [("left", nodes_left, edges_left),
            ("context", nodes_context, edges_context),
            ("right", nodes_right, edges_right)]

    for name, nodes, edges in span:
        out.append(f"\t{name} [")
        out.extend([f"\t\t{node}" for node in nodes])
        out.extend([f"\t\t{edge}" for edge in edges])
        out.append("\t]")

    out.append("]")

    gml_str = "\n".join(out)
    return gml_str
