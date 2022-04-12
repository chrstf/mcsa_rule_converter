"""
Marvin module
"""
from mrv_tools.tools import *
from mrv_tools.known_attribuets import *
from mrv_tools.chemistry import transition_metals, periodic_table
from bs4 import BeautifulSoup as Bs
import networkx as nx
import re
import inspect
from typing import List


def formal_charge_label(label: str, charge: int):
    """
    Create an atom gml label with charge
    :param label: atom label
    :param charge: the charge on the atom
    """
    chg_label = ""
    if charge == 1:
        chg_label = "+"
    elif charge > 1:
        chg_label = f"{charge}+"
    elif charge == -1:
        chg_label = "-"
    elif charge < -1:
        chg_label = f"{(charge * -1)}-"

    return f"{label}{chg_label}"


def nxgraph_to_gmlstr(nxgraph: nx.Graph):
    """
    Convert Networkx graph to GML string.
    :param nxgraph: Newrokx graph
    :return: GML string
    """
    gml_str = "graph [\n"
    # Nodes
    for gn in nxgraph.nodes:
        gml_str += f"\tnode [ id {gn} label \"{nxgraph.nodes[gn]['label']}\" ]\n"
    # Edges
    for ge in nxgraph.edges:
        gml_str += f"\tedge [ source {ge[0]} target {ge[1]} " \
                   f"label \"{nxgraph.edges[ge]['label']}\"]\n"
    gml_str += "]"

    return gml_str


def separate_components(gml: str, renumber_ids: bool = True, return_nx_graphs: bool = False):
    """
    Function to separate connected components
    :param gml: GML string
    :param renumber_ids: trigger for renumbering node ids for each graph
    :param return_nx_graphs: Trigger for only returning networkx graphs.
    :return: List of separated graph GML strings
    """
    graph = nx.parse_gml(gml, label="id")
    sep_gml_str = []  # list for seperated gml strings
    # evaluating whether graph is one connected component
    if not nx.is_connected(graph):
        components = nx.connected_components(graph)  # separate component IDs
        if return_nx_graphs is True:
            return components
        for c_set in components:
            edges = []
            id_dict = {}  # used for renumber_ids
            # constructing GML
            this_gml = "graph [\n"
            # nodes
            edges_visited = []
            for i, n in enumerate(c_set):
                id_dict[n] = i
                this_id = n
                if renumber_ids is True:
                    this_id = i
                this_gml += f"\tnode [ id {this_id} label \"{graph.nodes[n]['label']}\" ]\n"
                # collect edges
                for e in graph.edges:
                    edge_str = ""
                    if (n in e) and (e not in edges_visited):
                        edges.append((e[0], e[1], graph.edges[e]['label']))
                        edges_visited.append(e)

            # add edges to gml string
            for this_e in edges:
                e_source, e_target, e_label = this_e
                if renumber_ids:
                    e_source = id_dict[e_source]
                    e_target = id_dict[e_target]
                this_gml += f"\tedge [ source {e_source} target {e_target} "
                this_gml += f"label \"{e_label}\"]\n"
            this_gml += "]"
            sep_gml_str.append(this_gml)

        return sep_gml_str

    elif nx.is_connected(graph):
        return [gml]


class Atom:
    """
    A Marvin XML atom line
    """

    def __init__(self, line: Bs, the_bonds: list, mol: str, origin=None, verbose: int = 1, v_level: int = 1):
        for a in line.attrs:
            if a not in known_atom_attrs:
                warning_m(f"Unknown atom attribute: {a}", verbose=verbose)

        self.the_bonds = the_bonds
        self.mol = mol
        self.origin = origin
        self.holder = None  # A container for anything that might be usful being attached to an atom.
        self.id = line.attrs["id"]
        self.x2 = float(line.attrs["x2"])
        self.y2 = float(line.attrs["y2"])
        self.element_type = line.attrs["elementType"]
        self.content = remove_empty_lines(line.text)
        if v_level >= 2:
            if self.element_type == "*":
                warning_m(f"Atom type: {self.element_type} in atom id {self.id}", verbose=verbose)
        try:
            self.mrv_extra_label = line.attrs["mrvExtraLabel"]
        except KeyError:
            self.mrv_extra_label = None
        try:
            self.lone_pair = int(line.attrs["lonePair"])
        except KeyError:
            self.lone_pair = 0
        try:
            self.formal_charge = int(line.attrs["formalCharge"])
        except KeyError:
            self.formal_charge = 0
        try:
            self.mrv_alias = line.attrs["mrvAlias"]
            # a certain mrv format puts 0s where there should be nothing
            if self.mrv_alias == "0":
                self.mrv_alias = None
        except KeyError:
            self.mrv_alias = None
        try:
            self.h_count = line.attrs["hydrogenCount"]
        except KeyError:
            self.h_count = None
        try:
            self.sgroup_ref = line.attrs["sgroupRef"]
        except KeyError:
            self.sgroup_ref = None
        self.atom_parity = []
        for item in line.find_all("atomParity"):
            attributes = {
                "attrs": item.attrs,
                "content": remove_empty_lines(item.text)
            }
            self.atom_parity.append(attributes)

        if len(self.atom_parity) == 0:
            self.atom_parity = None
        try:
            self.r_group_ref = line.attrs["rgroupRef"]
        except KeyError:
            self.r_group_ref = None
        try:
            self.mrv_map = line.attrs["mrvMap"]
        except KeyError:
            self.mrv_map = None
        try:
            self.radical = line.attrs["radical"]
            # Radical center:
            # "monovalent" (doublet),
            # "divalent",
            # "divalent1" (singlet),
            # "divalent3" (triplet),
            # "trivalent",
            # "trivalent2" (doublet)
            # "trivalent4" (quartet)
            # In M-CSA there are only monovalent, divalent and divalent1
        except KeyError:
            self.radical = None
        try:
            self.mrv_pseudo = line.attrs["mrvPseudo"]
            # a certain mrv format puts 0s where there should be nothing
            if self.mrv_pseudo == "0":
                self.mrv_pseudo = None
        except KeyError:
            self.mrv_pseudo = None
        try:
            self.mrv_query_props = line.attrs["mrvQueryProps"]
            # a certain mrv format puts 0s where there should be nothing
            if self.mrv_query_props == "0":
                self.mrv_query_props = None
        except KeyError:
            self.mrv_query_props = None
        try:
            self.isotope = line.attrs["isotope"]
        except KeyError:
            self.isotope = None
        try:
            self.spin_multiplicity = line.attrs["spinMultiplicity"]
        except KeyError:
            self.spin_multiplicity = None
        try:
            self.ligand_order = line.attrs["ligandOrder"]
        except KeyError:
            self.ligand_order = None
        try:
            self.attachment_order = int(line.attrs["attachmentOrder"])
        except KeyError:
            self.attachment_order = None
        try:
            self.attachment_point = int(line.attrs["attachmentPoint"])
        except KeyError:
            self.attachment_point = None
        try:
            self.sgroup_attachment_point = line.attrs["sgroupAttachmentPoint"]
        except KeyError:
            self.sgroup_attachment_point = None

        self.s_e_push_status = 0

        self.electron_container = None

    def __str__(self):
        line = []
        not_these = ["the_bonds", "x2", "y2", "content"]
        for attri in self.__dict__:
            if attri in not_these:
                continue
            val = getattr(self, attri)
            if val is not None:
                if attri == "electron_container":
                    val = sum([len(ec.electrons) for ec in val])
                    attri = "electron number"
                line.append(f"{attri}: {val}")
        neighbors = self.get_connected_atoms(with_order=True)
        if neighbors is not None:
            neighbors = " ".join(map(str, [(a[0].str_short(), a[1]) for a in neighbors]))
            neighbors = f"1-neighbors: {neighbors}"
            line.append(neighbors)
        bonds = self.bonds()
        bond_cnt = 0
        coord_b_from = 0
        coord_b_to = 0
        for b in bonds:
            if b.order is not None:
                bond_cnt += b.order
            if b.convention == "cxn:coord":
                if self == b.atom_refs_2[0]:
                    coord_b_from += 1
                elif self == b.atom_refs_2[1]:
                    coord_b_to += 1
        line.append(f"num of bonds: {bond_cnt}")
        line.append(f"num of departing cxn:coord: {coord_b_from}")
        line.append(f"num of recieving cxn:coord: {coord_b_to}")
        return "\n".join(line)

    def str_short(self):
        return f"({self.id}, {self.element_type})"

    def xml(self, ind: str = "  ", mrv_extra_label: bool = True):
        """
        :return: XML style string
        """
        # example:
        # <atom id="<string>" elementType="<string>" x2="<float>" y2="<float>" mrvAlias="<string>"/>

        # building string
        attributes = []
        xml_line = f"<atom "
        known_atom_attrs_inv = invert_dict(known_atom_attrs)
        for key in self.__dict__:
            try:
                attri = known_atom_attrs_inv[key]
                if (mrv_extra_label is False) and (attri == "mrvExtraLabel"):
                    continue
                value = self.__dict__[key]
                if value is not None:
                    attributes.append(f"{attri}=\"{value}\"")
            except KeyError:
                # If key is not in the known list, do nothing.
                continue

        xml_line += " ".join(attributes)
        if (self.content != "") or self.atom_parity is not None:
            xml_line += ">\n"
            if self.atom_parity is not None:
                for thing in self.atom_parity:
                    attributes = []
                    xml_line += f"{ind}<atomParity "
                    for attri in thing["attrs"]:
                        val = thing["attrs"][attri]
                        attributes.append(f"{attri}=\"{val}\"")
                    xml_line += " ".join(attributes)
                    xml_line += ">"

                    xml_line += f"{thing['content']}"
                    xml_line += "</atomParity>\n"
            elif self.content != "":
                for line in self.content.splitlines():
                    xml_line += f"{ind}{line}\n"
            xml_line += "</atom>"
        else:
            xml_line += "/>"

        return xml_line

    def bonds(self):
        atom_bonds = []
        for b in self.the_bonds:
            if self in b.bond_atoms:
                atom_bonds.append(b)
        return atom_bonds

    def get_connected_atoms(self, with_order: bool = False):
        """
        Get connected atoms
        """
        collected = []
        if len(self.bonds()) == 0:
            return None

        for b in self.bonds():
            b: Bond
            con_atom = list(b.bond_atoms.difference({self}))[0]
            if with_order is True:
                order = b.order
                if order is None and b.convention == "cxn:coord":
                    order = "cxn:coord"
                collected.append((con_atom, order))
            else:
                collected.append(con_atom)

        return collected

    def convert_alias(self):
        """
        Create an alias label
        """
        if self.mrv_alias is None:
            return None

        label = self.element_type
        alias = self.mrv_alias
        caps = ''.join([lol.capitalize() for lol in alias])
        # checking if there is an amino acid label
        # standard amino acid pattern
        pattern = re.compile(r"([A-Z][A-Z][A-Z])(\d+)([A-Z]*)")
        p_result = pattern.search(caps)
        if p_result is not None:
            if p_result.group(1) in amino_acids:
                aa_label = p_result.group(1).capitalize()
                aa_id = p_result.group(2)
                chain_id = p_result.group(3)
                if len(chain_id) == 0:
                    chain_id = "None"
                alias = f"Amino({label},{aa_label},{aa_id},{chain_id})"
                return alias

        alias = f"Alias({label},{alias})"
        return alias

    def convert_query_props(self):
        """
        Create an alias label from query props.
        Copied from convert_alias()
        """
        if self.mrv_query_props is None:
            return None
        label = self.element_type
        alias = self.mrv_query_props
        caps = ''.join([lol.capitalize() for lol in alias])
        # checking if there is an amino acid label
        # standard amino acid pattern
        pattern = re.compile(r"([A-Z][A-Z][A-Z])(\d+)([A-Z]*)")
        p_result = pattern.search(caps)
        if p_result is not None:
            if p_result.group(1) in amino_acids:
                aa_label = p_result.group(1).capitalize()
                aa_id = p_result.group(2)
                chain_id = p_result.group(3)
                if len(chain_id) == 0:
                    chain_id = "None"
                alias = f"Amino({label},{aa_label},{aa_id},{chain_id})"
                return alias
        # Trim trailing :
        if re.search(r"\S+:", alias) is not None:
            alias = alias[:-1]
        alias = f"Alias({label},{alias})"
        return alias

    def atom_label_with_charge(self):
        """
        Create an atom label with charge
        """
        chg_label = ""
        if self.formal_charge == 1:
            chg_label = "+"
        elif self.formal_charge > 1:
            chg_label = f"{self.formal_charge}+"
        elif self.formal_charge == -1:
            chg_label = "-"
        elif self.formal_charge < -1:
            chg_label = f"{(self.formal_charge * -1)}-"

        return f"{self.element_type}{chg_label}"


def create_atom(aid: str = "a1", element: str = "H", bonds=None, mol_id: str = "m1"):
    """
    Function to create an atom
    :param aid: Atom id
    :param element: Element type
    :param bonds: The bonds of the whole molecule
    :param mol_id: The molID
    """
    if bonds is None:
        bonds = []
    atom_xml = f"<atom id=\"{aid}\" elementType=\"{element}\" x2=\"0.0\" y2=\"0.0\"/>"
    the_atom = Atom(Bs(atom_xml, "xml").find("atom"), bonds, mol_id)
    return the_atom


def create_bond_from_atom_list(atom_refs_2: list, verbose: int = 1):
    """
    Function to create a Bond class object.
    A small work around since the Bond class requires a Beautiful soup input.
    Atoms for the atomRefs2 attribute should be a list of Atom objects.
    :param atom_refs_2: A lsit of Atoms.
    :param verbose: Verbose output.
    :return: Bond object or Flase
    """
    # Some checks
    for a in atom_refs_2:
        if not isinstance(a, Atom):
            warning_m("Create bond: atom_refs_2 items must be of type Atom class.", verbose=verbose)
            return False

    # Start constructing a string
    bond_xml = f"<bond atomRefs2=\"a1 a2\" />"

    # create the bond using Beautiful soup
    the_bond = Bond(Bs(bond_xml, "xml").find("bond"), atom_refs_2[0].mol, verbose=verbose)
    # fix the atomRefs2 atoms
    the_bond.atom_refs_2 = atom_refs_2
    the_bond.bond_atoms = set(atom_refs_2)
    the_bond.order = 1
    return the_bond


def renumber_bond_ids(bonds: dict):
    for m in bonds:
        if bonds[m] is not None:
            for i, b in enumerate(bonds[m], start=1):
                b.id = f"b{i}"
    return True


class Bond:
    """
    A Marvin XML bond line
    """

    def __init__(self, line: Bs, mol: str, verbose: int = 1):
        for b in line.attrs:
            if b not in known_bond_attrs:
                warning_m(f"Unknown bond attribute: {b}", verbose=verbose)
        self.mol = mol
        self.holder = None
        self.aromatic = False
        try:
            self.id = line.attrs["id"]
        except KeyError:
            self.id = None
        try:
            self.atom_refs_2 = line.attrs["atomRefs2"].split()
        except KeyError:
            self.atom_refs_2 = None
        try:
            self.atom_refs_3 = line.attrs["atomRefs3"].split()
        except KeyError:
            self.atom_refs_3 = None
        try:
            self.atom_refs_4 = line.attrs["atomRefs4"].split()
        except KeyError:
            self.atom_refs_4 = None
        try:
            # "1", "S" (single), "2", "D" (double), "3", "T" (triple) or "A" (aromatic).
            try:
                self.order = float(line.attrs["order"])
            except ValueError:
                # aromatic case
                if line.attrs["order"] == "A":
                    self.order = 1.5
                    self.aromatic = True
        except KeyError:
            self.order = None
        try:
            self.convention = line.attrs["convention"]
        except KeyError:
            self.convention = None
        try:
            self.bond_stereo = line.attrs["bondStereo"]
        except KeyError:
            self.bond_stereo = None

        self.single_e_flow_buildup = 0  # used for creating a bond from single arrow pushes

        self.aref2_bond_stereo = None
        self.aref2_bond_stereo_aref4 = None
        if self.atom_refs_2 is not None:
            if line.find("bondStereo") is not None:
                self.aref2_bond_stereo = remove_empty_lines(line.find("bondStereo").text)
                try:
                    self.aref2_bond_stereo_aref4 = line.find("bondStereo").attrs["atomRefs4"].split()
                except KeyError:
                    self.aref2_bond_stereo_aref4 = None

        try:
            self.query_type = line.attrs["queryType"]
        except KeyError:
            self.query_type = None

        self.bond_atoms = set(self.atom_refs_2)

    def alt_str(self):
        alist = []
        b_order = {0: "-", 1: "—", 2: "═", 3: "≡"}
        try:
            v_b_order = f"{3 * b_order[self.order]}"
        except KeyError:
            v_b_order = "···"
        if self.convention == "cxn:coord":
            v_b_order = "——►"
        for fun_a in self.atom_refs_2:
            fun_a: Atom
            alist.append(f"({fun_a.id} {fun_a.element_type})")
        v_id = self.id
        if v_id is None:
            v_id = "No ID"
        return f"{v_id}: " + v_b_order.join(alist)

    def __str__(self):
        return self.xml()

    def xml(self, ind: str = " ", make_bond_stereo: bool = True):
        """
        :param ind: Indent
        :param make_bond_stereo: Create the bondStereo tag
        :return: XML style string
        """
        # example
        # <bond id="<string>" atomRefs2="<string>" order="<string>"/>

        # building string
        attributes = []

        xml_line = f"<bond "
        known_bond_attrs_inv = invert_dict(known_bond_attrs)
        for key in self.__dict__:
            try:
                attri = known_bond_attrs_inv[key]
                value = self.__dict__[key]
                if attri == "atomRefs2":
                    value = " ".join([a.id for a in self.atom_refs_2])
                if value is not None:
                    if attri == "order":
                        if value != 1.5:
                            value = int(value)
                        else:
                            value = "A"
                    attributes.append(f"{attri}=\"{value}\"")
            except KeyError:
                # If key is not in the known list, do nothing.
                continue

        xml_line += " ".join(attributes)
        if self.aref2_bond_stereo is not None and make_bond_stereo is True:
            xml_line += f">\n{ind}<bondStereo"
            if self.aref2_bond_stereo_aref4 is not None:
                xml_line += " atomRefs4=\"" + " ".join([a.id for a in self.aref2_bond_stereo_aref4]) + "\">"
            else:
                xml_line += ">"
            xml_line += self.aref2_bond_stereo
            xml_line += "</bondStereo>\n</bond>"
        else:
            xml_line += "/>"

        return xml_line


class Electron:
    """
    Marvin XML electron line
    """

    def __init__(self, line: Bs, atoms: dict, verbose: int = 1):
        for e in line.attrs:
            if e not in known_electron_attrs:
                warning_m(f"Unknown electron attribute: {e}", verbose=verbose)
        # m_id, a_id = line.attrs["atomRefs"].split(".")
        self.atom_refs = [
            lookup(atoms[lll[0]], "id", lll[1]) for lll in
            [
                lulu.split(".") for lulu in [lu for lu in line.attrs["atomRefs"].split()]
            ]
        ]
        self.dif_loc = line.attrs["difLoc"]
        self.content = line.text

    def __str__(self):
        return self.get_xml_string()

    def get_xml_string(self):
        xml_str = "<MElectron "
        attrs = []
        for attr in known_electron_attrs:
            val = getattr(self, known_electron_attrs[attr])
            if attr == "atomRefs":
                val = " ".join([f"{a.mol}.{a.id}" for a in val])
            attrs.append(f"{attr}=\"{val}\"")
        xml_str += " ".join(attrs)
        xml_str += "/>"
        return xml_str


class EContainer:
    """
    Marvin XML electron Container line
    """

    def __init__(self, line: Bs, atoms: dict, verbose: int = 1):
        for e in line.attrs:
            if e not in known_e_container_attrs:
                warning_m(f"Unknown e_container attribute: {e}", verbose=verbose)

        self.id = line.attrs["id"]
        self.occupation = line.attrs["occupation"]
        self.radical = line.attrs["radical"]
        self.content = remove_empty_lines(line.text)

        # electrons
        self.electrons = line.find_all("MElectron")
        self.atoms = None
        if self.electrons is not None:
            self.electrons = [Electron(e, atoms, verbose=verbose) for e in self.electrons]
            self.atoms = []
            for e in self.electrons:
                for e_atom in e.atom_refs:
                    self.atoms.append(e_atom)
            self.atoms = set(self.atoms)
            for e_atom in self.atoms:
                if e_atom.electron_container is None:
                    e_atom.electron_container = [self]
                else:
                    e_atom.electron_container.append(self)

    def __str__(self):
        return self.get_xml_string()

    def get_xml_string(self, ind: str = "  "):
        xml_str = "<MElectronContainer "
        attrs = []
        for attr in known_e_container_attrs:
            val = getattr(self, known_e_container_attrs[attr])
            if val is not None:
                attrs.append(f"{attr}=\"{val}\"")
        xml_str += " ".join(attrs)
        if (self.content != "") or (self.electrons is not None):
            xml_str += ">\n"
            if self.content != "":
                for item in self.content.splitlines():
                    xml_str += f"{ind}{item}\n"
            if self.electrons is not None:
                for e in self.electrons:
                    xml_str += f"{ind}{e.get_xml_string()}\n"
            xml_str += "</MElectronContainer>"
        else:
            xml_str += "/>"
        return xml_str


class AtomSetPoint:
    """
    Marvin XML MAtomSetPoint line
    """

    def __init__(self, line: Bs, atoms: dict, verbose: int = 1):
        for e in line.attrs:
            if e not in known_atom_set_point_attrs:
                warning_m(f"Unknown atom set point attribute: {e}", verbose=verbose)
        self.atom_refs_original = line.attrs["atomRefs"]
        self.atom_refs = [
            lookup(atoms[lll[0]], "id", lll[1]) for lll in
            [
                lulu.split(".") for lulu in [lu for lu in self.atom_refs_original.split()]
            ]
        ]
        try:
            self.weights = line.attrs["weights"]
        except KeyError:
            self.weights = None
        self.content = remove_empty_lines(line.text)

    def get_xml_string(self, ind: str = "  "):
        xml_str = "<MAtomSetPoint "
        attrs = []
        for attr in known_atom_set_point_attrs:
            val = getattr(self, known_atom_set_point_attrs[attr])
            if attr == "atomRefs":
                val = " ".join(f"{a.mol}.{a.id}" for a in getattr(self, "atom_refs"))
            if val is not None:
                attrs.append(f"{attr}=\"{val}\"")
        xml_str += " ".join(attrs)
        if self.content != "":
            xml_str += f"\n"
            for item in self.content.splitlines():
                xml_str += f"{ind}{item}\n"
            xml_str += "\n</MAtomSetPoint>"
        else:
            xml_str += "/>"
        return xml_str


class FlowBasePoint:
    """
    Marvin XML MEFlowBasePoint line
    """

    def __init__(self, line: Bs, atoms: dict, verbose: int = 1):
        for e in line.attrs:
            if e not in known_flow_base_point:
                warning_m(f"Unknown flow base point attribute: {e}", verbose=verbose)
        lu = line.attrs["atomRef"].split(".")
        self.atom_ref = lookup(atoms[lu[0]], "id", lu[1])
        self.content = remove_empty_lines(line.text)

    def get_xml_string(self, ind: str = "  "):
        xml_str = "<MEFlowBasePoint "
        attrs = []
        for attr in known_flow_base_point:
            val = getattr(self, known_flow_base_point[attr])
            if attr == "atomRef":
                val = f"{self.atom_ref.mol}.{self.atom_ref.id}"
            attrs.append(f"{attr}=\"{val}\"")
        xml_str += " ".join(attrs)
        if self.content != "":
            xml_str += f"\n"
            for item in self.content.splitlines():
                xml_str += f"{ind}{item}\n"
            xml_str += "\n</MEFlowBasePoint>"
        else:
            xml_str += "/>"
        return xml_str


class EFlow:
    """
    Marvin XML EFlow line.
    The defenition is as follows:
        Curved electron flow arrow.
        MEFlow is a subclass of MPolyline thus it has the same attributes, but it can only contain two points.
    """

    def __init__(self, line: Bs, atoms: dict, verbose: int = 1):
        for e in line.attrs:
            if e not in known_e_flow_attrs:
                warning_m(f"Unknown e_flow attribute: {e}", verbose=verbose)
        self.content = remove_empty_lines(line.text)
        self.id = line.attrs["id"]
        self.atoms = []  # atoms in EFlow
        self.has_transition_metal = False
        # Angle of the electron flow arrow arc.
        self.arc_angle = line.attrs["arcAngle"]
        # Describes which electroncontainer holds electron.
        try:
            self.b_e_container_index = line.attrs["baseElectronContainerIndex"]
        except KeyError:
            self.b_e_container_index = None
        # Describes which electron is base electron.
        try:
            self.b_e_index_in_container = line.attrs["baseElectronIndexInContainer"]
        except KeyError:
            self.b_e_index_in_container = None
        # Child of MEFlow. It represents an atom or atom pair (bond or incipient bond).
        self.atom_set_points = line.find_all("MAtomSetPoint")
        if self.atom_set_points is not None:
            self.atom_set_points = [AtomSetPoint(e, atoms, verbose=verbose) for e in self.atom_set_points]
            for asp in self.atom_set_points:
                self.atoms += asp.atom_refs
        else:
            self.atom_set_points = None
        # Child of MEFlow. Starting point of the electron flow arrow if the source is an atom.
        self.flow_base_point = line.find_all("MEFlowBasePoint")
        if len(self.flow_base_point) == 1:
            self.flow_base_point = FlowBasePoint(self.flow_base_point[0], atoms, verbose=verbose)
            self.atoms.append(self.flow_base_point.atom_ref)
        elif len(self.flow_base_point) > 1:
            self.flow_base_point = [FlowBasePoint(e, atoms, verbose=verbose) for e in self.flow_base_point]
        else:
            self.flow_base_point = None

        # Arrow head. Subclass of the Marvin XML MPolyline.
        # There might also be tailFlags, however, I colud not grep it in the downloaded M-CSA files. Put it anyway.
        try:
            self.head_flags = line.attrs["headFlags"]
        except KeyError:
            self.head_flags = None
        # Arrow tail. Subclass of the Marvin XML MPolyline.
        try:
            self.tail_flags = line.attrs["tailFlags"]
        except KeyError:
            self.tail_flags = None
        try:
            self.head_width = line.attrs["headWidth"]
        except KeyError:
            self.head_width = None
        try:
            self.head_length = line.attrs["headLength"]
        except KeyError:
            self.head_length = None
        try:
            self.tail_skip = line.attrs["tailSkip"]
        except KeyError:
            self.tail_skip = None
        try:
            self.head_skip = line.attrs["headSkip"]
        except KeyError:
            self.head_skip = None

        self.atoms = set(self.atoms)
        for a in self.atoms:
            if a.element_type in transition_metals:
                self.has_transition_metal = True

    def __str__(self):
        out = []
        if self.head_flags is not None:
            out.append(f"Head flags: {self.head_flags}")
        if self.flow_base_point is not None:
            if isinstance(self.flow_base_point, list):
                for base in self.flow_base_point:
                    out.append(f"Base point: {base.atom_ref.str_short()}")
            else:
                out.append(f"Base point:{self.flow_base_point.atom_ref.str_short()}")
        if self.atom_set_points is not None:
            for seti in self.atom_set_points:
                out.append(f"Atom set point: {' '.join([a.str_short() for a in seti.atom_refs])}")

        out.append(f"{len(self.atoms)} different atoms")

        return "\n".join(out)

    def get_xml_string(self, ind: str = "  "):
        xml_str = "<MEFlow "
        attrs = []
        for attr in known_e_flow_attrs:
            val = getattr(self, known_e_flow_attrs[attr])
            if val is not None:
                attrs.append(f"{attr}=\"{val}\"")
        xml_str += " ".join(attrs)
        if (self.content != "") or (self.atom_set_points is not None) or (self.flow_base_point is not None):
            xml_str += ">\n"
            if self.content != "":
                for item in self.content.splitlines():
                    xml_str += f"{ind}{item}\n"

            if self.flow_base_point is not None:
                xml_str += f"{ind}{self.flow_base_point.get_xml_string(ind=ind)}\n"
            if self.atom_set_points is not None:
                for item in self.atom_set_points:
                    xml_str += f"{ind}{item.get_xml_string(ind=ind)}\n"
            xml_str += "</MEFlow>"
        else:
            xml_str += "/>"

        return xml_str


class MrvData:
    """
    Class that holds Marvin XML data.
    """

    def __init__(self, content: Bs, origin: str = None, verbose: int = 0):
        """
        :param content: Beautiful Soup content.
        :param origin: File name or other origin of the data.
        :param verbose: verbose output.
        """
        self.explicit_h = False
        self.origin = origin
        self.verbose = verbose
        self.content = content
        try:
            self.cml_attrs = content.find("cml").attrs
        except AttributeError:
            warning_m("No \"cml\" tag in marvin cml data.",
                      verbose=verbose, verbose_threshold=2)
            self.cml_attrs = {"version": ""}
        ######################
        # Chemical Structure #
        ######################
        self.chem_struc = self.content.find("MChemicalStruct")
        if self.chem_struc is None:
            self.chem_struc = self.content
            warning_m("No \"MChemicalStruct\" tag in marvin CML data.\nJumping to \"molecule\" tag.",
                      verbose=verbose, verbose_threshold=2)
        self.molecule = {}
        self.molecule_other = {}
        self.molecule_attrs = {}
        for i, mol in enumerate(self.chem_struc.find_all("molecule"), start=1):
            # After running through openbabel molID get's lost and replaced by id
            # Some Marvin files also have an id tag instead of an molID tag
            mol_attrs: dict = mol.attrs
            if ("id" in mol_attrs) and ("molID" not in mol_attrs):
                mol_attrs["molID"] = mol_attrs["id"]
                warning_m(f"No \"molID\" in mollecule attributes. oing with \"id={mol_attrs['id']}\".",
                          verbose=verbose, verbose_threshold=2)
            elif ("id" not in mol_attrs) and ("molID" not in mol_attrs):
                mol_attrs["molID"] = f"m{i}"
                mol_attrs["id"] = f"m{i}"
                warning_m(
                    f"No \"molID\" and \"id\"in mollecule attributes.Going with \"m{i}\".",
                    verbose=verbose
                )
            self.molecule_attrs[mol.attrs["molID"]] = mol_attrs
            self.molecule[mol.attrs["molID"]] = mol
            self.molecule_other[mol.attrs["molID"]] = None
            mol_attrs: list = [key for key in mol.attrs]
            mol_attrs.remove("molID")
            if len(mol_attrs) > 0:
                self.molecule_other[mol.attrs["molID"]] = {this_attr: mol.attrs[this_attr] for this_attr in
                                                           mol_attrs}
        self.ismol = {mol: True for mol in self.molecule}
        self.mol_ids = [mol for mol in self.molecule]

        # atom and bonds
        self.atom_array = {key: arr.find("atomArray") for key, arr in self.molecule.items()}
        # check whether the molecule has atoms
        # if not it is probably not a complete molecule
        for mol in self.atom_array:
            if self.atom_array[mol] is None:
                self.ismol[mol] = False

        # Collect bonds
        # Has to be done before the atoms
        self.bond_array = {key: arr.find("bondArray") for key, arr in self.molecule.items()}
        self.bonds = {}
        for key, arr in self.bond_array.items():
            if arr is not None:
                self.bonds[key] = [Bond(a, key, verbose=verbose) for a in arr.find_all("bond")]
            else:
                self.bonds[key] = None
        # collect atoms
        self.atoms = {}
        self.has_transition_metal = False
        self.mass = {key: 0 for key in self.atom_array}
        for key, arr in self.atom_array.items():
            if arr is not None:
                self.atoms[key] = [Atom(a, self.bonds[key], key, origin=origin, verbose=verbose, v_level=1) for a in
                                   arr.find_all("atom")]
                for a in self.atoms[key]:
                    if a.element_type in transition_metals:
                        self.has_transition_metal = True
                    try:
                        self.mass[key] += periodic_table[a.element_type]["atomic_mass"]
                    except KeyError:
                        if a.element_type not in ["*", "R"]:
                            warning_m(f"Skipping unknown element \"{a.element_type}\" for mass calculation.",
                                      function_name=inspect.stack()[0][3],
                                      verbose=verbose)
            else:
                self.atoms[key] = None

        # In bonds, replace atom IDs with atoms
        self.coordination_bonds = None
        if self.bonds is not None:
            for m in self.bonds:
                if self.bonds[m] is not None:
                    replace_here = ["atom_refs_2", "atom_refs_3", "atom_refs_4", "bond_atoms"]
                    for bond in self.bonds[m]:
                        if bond.convention == "cxn:coord":
                            if self.coordination_bonds is None:
                                self.coordination_bonds = {key: 0 for key in self.mol_ids}
                            self.coordination_bonds[m] += 1
                        for attr in replace_here:
                            if getattr(bond, attr) is not None:
                                new_list = [lookup(self.atoms[m], "id", a_id) for a_id in getattr(bond, attr)]
                                if False not in new_list:
                                    if isinstance(getattr(bond, attr), set):
                                        new_list = set(new_list)
                                    setattr(bond, attr, new_list)
                                else:
                                    if self.verbose is True:
                                        warning_m(
                                            f"{attr}: Could not match one of the follwoing atom IDs "
                                            f"{' '.join(getattr(bond, attr))}"
                                            f"mol: {m}"
                                        )
                        if bond.aref2_bond_stereo_aref4 is not None:
                            new_list = [lookup(self.atoms[m], "id", a_id) for a_id in bond.aref2_bond_stereo_aref4]
                            if False not in new_list:
                                bond.aref2_bond_stereo_aref4 = new_list
                            else:
                                if self.verbose is True:
                                    warning_m(
                                        f"aref2_bond_stereo_aref4: Could not match one of the follwoing atom IDs "
                                        f"{' '.join(bond.aref2_bond_stereo_aref4)}"
                                    )
                    if None in [b.id for b in self.bonds[m]]:
                        # If there is only one None re-do all labels
                        for i, b in enumerate(self.bonds[m], start=1):
                            b.id = f"b{i}"

        # replace atomRefs in special SruSgroup molecule
        for key in self.molecule_other:
            if self.molecule_other[key] is not None:
                if ("role" in self.molecule_other[key]) and ("atomRefs" in self.molecule_other[key]):
                    if self.molecule_other[key]["role"] == "SruSgroup":
                        # CAUTION: Assuming that all atoms are defined in self.atoms["m1"]
                        atom_refs = [lookup(self.atoms["m1"], "id", a) for a in
                                     self.molecule_other[key]["atomRefs"].split()]
                        self.molecule_other[key]["atomRefs"] = atom_refs

        ##################
        # Electron Stuff #
        ##################
        self.electron_containers = self.content.find_all("MElectronContainer")
        if len(self.electron_containers) == 0:
            self.electron_containers = None
        if self.electron_containers is not None:
            self.electron_containers = [
                EContainer(line, self.atoms, verbose=verbose) for line in self.electron_containers
            ]

        self.e_flow_containers = self.content.find_all("MEFlow")
        self.head_flags = []
        self.e_flow_transition_metal = False
        if len(self.e_flow_containers) == 0:
            self.e_flow_containers = None
        if self.e_flow_containers is not None:
            self.e_flow_containers = [
                EFlow(line, self.atoms, verbose=verbose) for line in self.e_flow_containers
            ]
            for e in self.e_flow_containers:
                if e.has_transition_metal is True:
                    self.e_flow_transition_metal = True
                if e.head_flags is not None:
                    self.head_flags.append(e.head_flags)

    def e_flow_atoms(self):
        if self.e_flow_containers is None:
            return None
        atoms = []
        for ef in self.e_flow_containers:
            atoms += ef.atoms
        return set(atoms)

    def find_highest_id(self, which: str = "atoms"):
        """
        Function for finding the highest ID.
        :param which: atoms or bonds
        """
        if which not in ["atoms", "bonds"]:
            return False
        out = {key: None for key in self.ismol if self.ismol[key] is True}
        for key in out:
            ids = sorted([int(re.search(r"(\d+)", item.id).groups()[0]) for item in getattr(self, which)[key]])
            out[key] = f"{which[0]}{ids[-1]}"
        return out

    def flatten_atoms(self) -> List[Atom]:
        atoms_flat = []
        for m in self.mol_ids:
            if self.atoms[m] is not None:
                atoms_flat += self.atoms[m]
        return atoms_flat

    def flatten_bonds(self) -> List[Bond]:
        bonds_flat = []
        for m in self.mol_ids:
            if self.bonds[m] is not None:
                bonds_flat += self.bonds[m]
        return bonds_flat

    def redo_atom_ids(self, verbose: int = 0):
        """
        Function for changing atom IDs to "integers".
        Unfortunately, the IDs cannot be converted to an Int type
        since further processing in the class requires Str type.
        CAUTION: Be aware that changing atom IDs to solely integers (removing the 'a') results in an invalid Marvin file.
        A file containing atom IDs as integers will result in errors when trying to open with Marvin.
        """
        warning_m("Changing the atom IDs to integers is resulting in a corrupt Marvin file.", verbose=verbose)
        # atoms
        for key in self.atoms:
            if self.atoms[key] is not None:
                for atom in self.atoms[key]:
                    atom.id = atom.id.strip("a")

    def reset_atom_ids(self):
        """
        Function for renumbering atom IDs in order.
        """
        for key in self.atoms:
            if self.atoms[key] is not None:
                for i, atom in enumerate(self.atoms[key], start=1):
                    atom.id = f"a{i}"

    def reset_bond_ids(self):
        """
        Function for renumbering bond IDs in order.
        """
        for key in self.bonds:
            if self.bonds[key] is not None:
                for i, bond in enumerate(self.bonds[key], start=1):
                    bond.id = f"b{i}"

    def electron_contaiers_for_atom(self, the_atom: Atom):
        keep = []
        for ec in self.electron_containers:
            if ec.atoms == {the_atom}:
                keep.append(ec)
        return keep

    def next_a_id(self):
        all_ids = sorted([int(re.search(r"(\d+)", a.id).group()) for a in self.flatten_atoms()])
        return f"a{all_ids[-1] + 1}"

    def __str__(self):
        out = ""
        if self.origin is not None:
            out += f"{self.origin} | "
        out += f"{len(self.mol_ids)} molecules: "
        if len(self.mol_ids) == 1:
            out += f"{len(self.mol_ids)} molecule: "
        out += " ".join(self.mol_ids)
        out += f" | {sum([len(self.atoms[key]) for key in self.atoms if self.atoms[key] is not None]):3d} atoms"
        out += f" | {sum([len(self.bonds[key]) for key in self.bonds if self.bonds[key] is not None]):3d} bonds |"
        for key in self.mol_ids:
            out += f" {key}: {self.mass[key]:7.2f} g/mol"
        if self.e_flow_containers is None:
            out += f" | No MEFlow"
        elif self.e_flow_containers is not None:
            out += f" | {len(self.e_flow_containers)} MEFlow"

        return out

    def gml_string(self, return_ids: bool = False, renumber_ids: bool = True, cleave_coord_bonds: bool = False,
                   only_these_atoms: list = None, use_atom_ids: bool = False, no_seperate: bool = False):
        """
        Generate GML string(s).
        Seperates the connected components by using networkx.
        :param return_ids: Also return internal ID dictionary.
        :param renumber_ids: Renumber IDs in case of multiple graphs.
        :param cleave_coord_bonds: Skip coordination bonds.
        :param only_these_atoms: List of Atom class objects.
        :param use_atom_ids: Pairs with the only_these_atoms argument. Use atom IDs to identify atoms in mol.
        :param no_seperate: Do not sperate into connected components.
        :return: Dictionary containing lists. Each list item in a GML graph.
        """
        b_types = {
            1.: "-",
            2.: "=",
            3.: "#",
            1.5: ":",
            "cxn:coord": ">"
        }

        out = {}
        internal_ids = {mol: {} for mol in self.atoms}
        only_these_bonds = []
        if only_these_atoms is not None:
            if use_atom_ids is True:
                only_these_atoms = [lookup(self.atoms[aa.mol], "id", aa.id) for aa in only_these_atoms]
            for atom in only_these_atoms:
                only_these_bonds += atom.bonds()
            only_these_bonds = set(only_these_bonds)
        for mol in self.atoms:
            gml_str = "graph[\n"
            if (self.atoms[mol] is not None) and (self.ismol[mol] is True):
                for i, a in enumerate(self.atoms[mol]):
                    if only_these_atoms is not None:
                        if a not in only_these_atoms:
                            continue
                    label = a.element_type
                    if a.formal_charge is not None:
                        label = formal_charge_label(label, a.formal_charge)
                    if a.mrv_alias is not None:
                        label = a.convert_alias()

                    gml_str += f"\tnode [ id {i} label \"{label}\" ]\n"
                    internal_ids[mol][a.id] = i

                for i, b in enumerate(self.bonds[mol]):
                    if only_these_atoms is not None and len(only_these_bonds) == 0:
                        # no bonds
                        break
                    if len(only_these_bonds) != 0:
                        if b not in only_these_bonds:
                            continue
                    a1, a2 = b.atom_refs_2
                    a1 = a1.id
                    a2 = a2.id
                    b_type = b.order
                    if b_type is None:
                        if b.convention == "cxn:coord":
                            b_type = b.convention
                            if cleave_coord_bonds is True:
                                continue
                        else:
                            b_type = 1
                    try:
                        b_type = b_types[b_type]
                    except KeyError:
                        warning_m(f"Unknown bond type: {b_type}. Using \"-\" instead.", verbose=self.verbose)
                        b_type = "-"
                    gml_str += f"\tedge [ source {internal_ids[mol][a1]} target {internal_ids[mol][a2]} "
                    gml_str += f"label \"{b_type}\"]\n"

                gml_str += "]"
                if no_seperate is True:
                    out[mol] = gml_str
                elif no_seperate is False:
                    out[mol] = separate_components(gml_str, renumber_ids=renumber_ids)

        if return_ids is True:
            return [out, {mol: invert_dict(internal_ids[mol]) for mol in internal_ids}]

        return out

    def mrv_string(
            self, xml_header: bool = False, specific_mol: list = None, do_cml: bool = True,
            only_these_atoms: list = None, cleave_coord_bonds: bool = False, mrv_extra_label: bool = True,
            make_bond_stereo: bool = True
    ):
        from datetime import datetime
        now = datetime.now()
        # date and time format: dd-mm-YYYY H:M:S
        time_format = "%d-%m-%Y %H:%M:%S"

        ind = "  "  # line indent
        ind_cml = 1

        mrv_str = ""
        if xml_header is True:
            mrv_str += "<?xml version=\"1.0\"?>\n"
        # First tag
        if do_cml:
            mrv_str += "<cml"
            if len(self.cml_attrs) > 0:
                mrv_str += " "
                attrs = []
                self.cml_attrs["version"] = f"generated by Novo-Project script on {now.strftime(time_format)}"
                self.cml_attrs["group:url"] = "https://cheminf.imada.sdu.dk/"
                for attr in self.cml_attrs:
                    val = self.cml_attrs[attr]
                    attrs.append(f"{attr}=\"{val}\"")
                mrv_str += " ".join(attrs)
            mrv_str += ">\n"
        mrv_str += f"{ind * ind_cml}<MDocument>\n"
        mrv_str += f"{ind * (1 + ind_cml)}<MChemicalStruct>\n"

        # looping through molecules
        for mol in self.atoms:
            if specific_mol is not None:
                if mol not in specific_mol:
                    continue
            mrv_str += f"{ind * (2 + ind_cml)}<molecule molID=\"{mol}\""
            if self.molecule_other[mol] is not None:
                for o, o_val in self.molecule_other[mol].items():
                    if o == "atomRefs":
                        o_val = " ".join([o_a.id for o_a in self.molecule_other[mol]["atomRefs"]])
                    mrv_str += f" {o}=\"{o_val}\""
            mrv_str += ">\n"

            # atoms
            if self.atoms[mol] is not None:
                mrv_str += f"{ind * (3 + ind_cml)}<atomArray>\n"
                for atom in self.atoms[mol]:
                    if only_these_atoms is not None:
                        if atom not in only_these_atoms:
                            continue
                    for line in atom.xml(ind=ind, mrv_extra_label=mrv_extra_label).splitlines():
                        mrv_str += f"{ind * (4 + ind_cml)}" + line + "\n"
                mrv_str += f"{ind * (3 + ind_cml)}</atomArray>\n"

            # bonds
            if self.bonds[mol] is not None:
                mrv_str += f"{ind * (3 + ind_cml)}<bondArray>\n"
                for bond in self.bonds[mol]:
                    if cleave_coord_bonds is True:
                        if bond.convention == "cxn:coord":
                            continue
                    if only_these_atoms is not None:
                        skip = False
                        for a in bond.bond_atoms:
                            if a not in only_these_atoms:
                                skip = True
                                break
                        if skip is True:
                            continue
                    for line in bond.xml(ind=ind, make_bond_stereo=make_bond_stereo).splitlines():
                        mrv_str += f"{ind * (2 + ind_cml)}" + line + "\n"
                mrv_str += f"{ind * (3 + ind_cml)}</bondArray>\n"

            # closing all remaining tags
            mrv_str += f"{ind * (1 + ind_cml)}</molecule>\n"
        mrv_str += f"{ind * (1 + ind_cml)}</MChemicalStruct>\n"
        if self.electron_containers is not None:
            for entry in self.electron_containers:
                if only_these_atoms is not None:
                    skip = False
                    for e in entry.electrons:
                        for a in e.atom_refs:
                            if a not in only_these_atoms:
                                skip = True
                                break
                    if skip is True:
                        continue
                for line in entry.get_xml_string(ind=ind).splitlines():
                    mrv_str += f"{ind * (1 + ind_cml)}{line}\n"
        if self.e_flow_containers is not None:
            for entry in self.e_flow_containers:
                if only_these_atoms is not None:
                    skip = False
                    for a in entry.atoms:
                        if a not in only_these_atoms:
                            skip = True
                            break
                    if skip is True:
                        continue
                for line in entry.get_xml_string(ind=ind).splitlines():
                    mrv_str += f"{ind * (1 + ind_cml)}{line}\n"
        mrv_str += f"{ind * ind_cml}</MDocument>\n"
        if do_cml:
            mrv_str += "</cml>"

        return mrv_str

    def pdb_string(self, specific_mol: list = None):
        visited_bonds = []
        atom_id_tracker = {
            m: {atom.id: i for i, atom in enumerate(self.atoms[m], start=1)} for m in self.ismol
        }
        connect_blocks = []
        atom_blocks = []

        for m in self.ismol:
            if specific_mol is not None:
                if m not in specific_mol:
                    continue
            for i, atom in enumerate(self.atoms[m], start=1):
                pdb_line = "ATOM  "  # Record name
                pdb_line += "{:5d}".format(i)  # serial
                pdb_line += " "  # empty
                pdb_line += "{:^4s}".format(atom.element_type)  # name
                pdb_line += " "  # altLoc
                pdb_line += "{:3s}".format("UNL")  # resName
                pdb_line += " "  # empty
                pdb_line += " "  # chain
                pdb_line += "{:4d}".format(1)  # resSeq
                pdb_line += " "  # iCode
                pdb_line += " " * 3  # empty
                pdb_line += "{:8.3f}".format(float(atom.x2))  # x
                pdb_line += "{:8.3f}".format(float(atom.y2))  # y
                pdb_line += "{:8.3f}".format(.0)  # z
                pdb_line += "{:6.2f}".format(1.)  # occupancy
                pdb_line += "{:6.2f}".format(1.)  # tempFactor
                pdb_line += " " * 10  # empty
                pdb_line += "{:>2s}".format(atom.element_type)  # element
                pdb_line += " "  # charge

                atom_blocks.append(pdb_line)

                for b in atom.bonds():
                    if b.bond_atoms not in visited_bonds:
                        # There are different ways how to do the CONECT lines
                        # For simplicity we are doing each bond individually
                        bonds = list(b.bond_atoms)

                        con_line = "CONECT"  # Record name
                        con_line += "{:11d}".format(atom_id_tracker[m][bonds[0]])  # atom id 1
                        con_line += "{:11d}".format(atom_id_tracker[m][bonds[1]])  # atom id 2
                        connect_blocks.append(con_line)

                        # Double and triple bonds are added twice and thrice, respectively.
                        try:
                            order = int(b.order)
                            if 2 <= order <= 3:
                                for j in range(order - 1):
                                    connect_blocks.append(con_line)

                        except ValueError:
                            pass

                    visited_bonds.append(b.bond_atoms)

        pdb = "\n".join(atom_blocks + connect_blocks)
        pdb += "\nEND"

        return pdb

    def write_gml(self, fname: str, verbose: int = 1, cleave_coord_bonds: bool = False,
                  only_these_atoms: list = None, use_atom_ids: bool = False):
        gmls = self.gml_string(cleave_coord_bonds=cleave_coord_bonds,
                               only_these_atoms=only_these_atoms,
                               use_atom_ids=use_atom_ids)
        for i, mol in enumerate(gmls, start=1):
            for j, g in enumerate(gmls[mol], start=1):
                ext_fname = f"{fname}_{mol}_{j}.gml"
                with open(ext_fname, "w") as f:
                    f.write(g)
                message(f"Wrote XML data to \"{ext_fname}\"", verbose=verbose)

    def write_mrv(
            self, fname: str, xml_header: bool = False, verbose: int = 1, specific_mol: list = None,
            do_cml: bool = False, only_these_atoms: list = None, cleave_coord_bonds: bool = False,
            mrv_extra_label: bool = True, make_bond_stereo: bool = True
    ):
        with open(fname, "w") as f:
            f.write(
                self.mrv_string(
                    xml_header=xml_header,
                    specific_mol=specific_mol,
                    do_cml=do_cml,
                    only_these_atoms=only_these_atoms,
                    cleave_coord_bonds=cleave_coord_bonds,
                    mrv_extra_label=mrv_extra_label,
                    make_bond_stereo=make_bond_stereo
                ))
        message(f"Wrote XML data to \"{fname}\"", verbose=verbose)

    def tikz_string(self, atom_id: bool = False,
                    hidden_atoms: bool = False,
                    atom_colors: bool = False,
                    mrv_alias: bool = False):

        # elfow
        # eflow need to be done beforehand because information needs to be collected
        set_points = {}
        ef_lines = []
        if self.e_flow_containers is not None:
            for _i, ef in enumerate(self.e_flow_containers):
                # four cases
                if ef.flow_base_point is not None:
                    # in this case there should only be 1 set point
                    # case 1: 1 source, 1 target
                    if len(ef.atom_set_points[0].atom_refs) == 1:

                        ef_lines.append(f'\draw[electron jump] '
                                        f'({ef.flow_base_point.atom_ref.id}) '
                                        f'to [bend right=90, looseness=2] '
                                        f'({ef.atom_set_points[0].atom_refs[0].id});')
                    # case 2: 1 source, 2 target
                    else:
                        ef_lines.append(f'\draw[electron jump] '
                                        f'({ef.flow_base_point.atom_ref.id}) '
                                        f'to [bend right=90, looseness=2] '
                                        f'(ef{_i}.center);')
                        set_points[
                            (ef.atom_set_points[0].atom_refs[0].id, ef.atom_set_points[0].atom_refs[1].id)] = f'ef{_i}'

                else:
                    set_points[
                        ef.atom_set_points[0].atom_refs[0].id, ef.atom_set_points[0].atom_refs[1].id] = f'ef{_i}id1'
                    # case 3: 2 source, 1 target
                    if len(ef.atom_set_points[1].atom_refs) == 1:
                        ef_lines.append(f'\draw[electron jump] '
                                        f'(ef{_i}id1.center) '
                                        f'to [bend right=90, looseness=2] '
                                        f'({ef.atom_set_points[1].atom_refs[0].id});')
                    # case 4: 2 source, 2 target
                    else:
                        ef_lines.append(f'\draw[electron jump] '
                                        f'(ef{_i}id1.center) '
                                        f'to [bend right=90, looseness=2] '
                                        f'(ef{_i}id2.center);')
                        set_points[
                            ef.atom_set_points[1].atom_refs[0].id, ef.atom_set_points[1].atom_refs[1].id] = f'ef{_i}id2'
        # Start creating the string
        out = _tikzset + '\n'
        out += "\\begin{tikzpicture}\n"
        out += "\t% nodes\n"
        # atoms
        for atom in self.atoms["m1"]:
            name = formal_charge_label(atom.element_type, atom.formal_charge)
            name = "\ch{" + name + "}"
            if atom.element_type == "C":
                name = ""
                if mrv_alias is True and atom.mrv_alias is not None:
                    name = atom.mrv_alias
            if atom_id is True:
                name = atom.id
            a_option = []
            if hidden_atoms is True:
                a_option.append('hidden')
            if atom_colors is True:
                if atom.element_type == 'O':
                    a_option.append('oxygen')
                if atom.element_type == 'N':
                    a_option.append('nitrogen')
            out += f"\t\\node[{','.join(a_option)}] ({atom.id}) at ({atom.x2},{atom.y2})" + " {" + f"{name}" + "};\n"

        out += "\t% edges\n"
        b_types = {
            1.: "",
            2.: "chem double",
            3.: "triple",
            1.5: "",
            "cxn:coord": "->"
        }
        # bonds
        print(set_points)
        for bond in self.bonds["m1"]:
            b_label = b_types[bond.order]
            atoms = bond.atom_refs_2
            try:
                a1 = atoms[0].id
                a2 = atoms[1].id
                middle_node = ''
                if {a1, a2} in [{_key[0], _key[1]} for _key in set_points.keys()]:
                    try:
                        spec_node = set_points[(a1, a2)]
                        set_points.pop((a1, a2))
                    except KeyError:
                        spec_node = set_points[(a2, a1)]
                        set_points.pop((a2, a1))
                    middle_node = f' node[draw=none] ({spec_node})' + '{}'
                if mrv_alias is False:
                    if atoms[0].element_type == "C":
                        a1 = a1 + ".center"
                    if atoms[1].element_type == "C":
                        a2 = a2 + ".center"
                else:
                    if atoms[0].element_type == "C" and atoms[0].mrv_alias is None:
                        a1 = a1 + ".center"
                    if atoms[1].element_type == "C" and atoms[1].mrv_alias is None:
                        a2 = a2 + ".center"
            except AttributeError:
                a1 = atoms[0]
                a2 = atoms[1]

            b_option = []
            if hidden_atoms is True:
                b_option.append('hidden')
            b_option.append(b_label)

            out += f"\t\draw[{','.join(b_option)}] ({a1}) --{middle_node} ({a2});\n"
        # missing bonds for eflow representation
        if len(set_points) > 0:
            out += "\t% electron flow bonds\n"
            for key in set_points.keys():
                a1, a2 = key
                middle_node = f' node[draw=none] ({set_points[key]})' + '{}'
                out += f"\t\draw[hidden, dotted] ({a1}) --{middle_node} ({a2});\n"
        # write electron flow
        out += "\t% electron flow\n"
        for line in ef_lines:
            out += f'\t{line}\n'
        out += "\\end{tikzpicture}"
        return out


_tikzset = '''\\tikzset{%
	hidden/.style={
		black!25, thick
	},
	chem double/.style={
		double,
		double distance=2pt
	},
	oxygen/.style={red},
	nitrogen/.style={blue},
	sulfur/.style={yellow}
    electron jump/.style={
    -latex,
    electronjump,
    every node/.style={%
        fill=white,
        font=\scriptsize,
        ellipse,
        inner sep=0,
        outer sep=0,
        minimum size=0,
        sloped,
        anchor=center%
        }%
    }
}'''


def convert_special_format(soup: Bs):
    """
    Convert a Marvin file wehre the atomArray is in a different format.
    Here the atom tags are set as atomArray attributes.
    """
    molecule = soup.find_all("molecule")
    for mol in molecule:
        atom_array = mol.find("atomArray")
        atom_array_sperated = []
        if atom_array is None:
            continue
        for i, attr in enumerate(atom_array.attrs):
            list_data = atom_array.attrs[attr].split()
            if i == 0:
                atom_array_sperated = [{} for x in list_data]
            if len(atom_array_sperated) != len(list_data):
                error_m(f"Corrupt Marvin file.\n"
                        f"Number of items do not match:\n"
                        f"{len(atom_array_sperated)} vs {len(list_data)}")
            for j, item in enumerate(list_data):
                atom_array_sperated[j].update({attr: item})

        atom_classes = []
        for atom in atom_array_sperated:
            new_atom = create_atom(atom["atomID"], element=atom["elementType"])
            for attr in atom:
                if attr in ["atomID", "elementType"]:
                    continue
                setattr(new_atom, known_atom_attrs[attr], atom[attr])

            atom_classes.append(new_atom)

        atom_soup = Bs("<atomArray></atomArray>", "xml")
        newarray = atom_soup.atomArray
        for atom in atom_classes:
            newarray.append(Bs(atom.xml(), "xml"))
        mol.find("atomArray").replace_with(atom_soup)
    return soup


def read_from_file(fname: str, origin: str = None, verbose: int = 0):
    """
    Read from Marvin file.
    :param fname: File name
    :param origin: Can be used to set the origin (e.g., filename) of the marvin data
    :param verbose: Verbose output
    :return: MrvData class object
    """
    soup = Bs(open(fname, "r"), "xml")
    # check format of the marvin file
    # Read special mrv format where atomArray is in a different format
    if len(soup.find("atomArray").attrs) != 0:
        soup = convert_special_format(soup)

    return MrvData(soup, origin=origin, verbose=verbose)


def read_from_string(a_string: str, origin: str = None, verbose: int = 0, verbose_level: int = 1):
    """
    Read from Marvin string.
    :param a_string: marvin as string
    :param origin: Can be used to set the origin (e.g., filename) of the marvin data
    :param verbose: Verbose output
    :param verbose_level: Verbosity level
    :return: MrvData class object
    """
    soup = Bs(a_string, "xml")
    # check format of the marvin file
    # Read special mrv format where atomArray is in a different format
    if len(soup.find("atomArray").attrs) != 0:
        soup = convert_special_format(soup)
    return MrvData(soup, origin=origin, verbose=verbose)


def initialize_mrv_data() -> MrvData:
    dummy = read_from_string("""<cml>
  <MDocument>
    <MChemicalStruct>
      <molecule molID="m1">
        <atomArray>
        </atomArray>
        <bondArray>
        </bondArray>
      </molecule>
    </MChemicalStruct>
  </MDocument>
</cml>""")
    return dummy


def get_coordination_bonds(molecule: MrvData):
    """
    Get coordination bonds.
    :param molecule: a molecule
    :return: A dictionary containing lists of coorindation bonds
    """
    out = {}
    for mol in molecule.bonds:
        if molecule.ismol[mol] is True:
            if mol not in out:
                out.update({mol: []})
            for b in molecule.bonds[mol]:
                if b.convention == "cxn:coord":
                    out[mol].append(b)
    return out


def set_coordination_bonds(molecule: MrvData, coord_bonds: dict):
    """
    Set coordination bonds.
    :param molecule: The Marvin file data.
    :param coord_bonds: The output from get_coordination_bonds()
    """
    for mol in molecule.bonds:
        if molecule.ismol[mol] is True:
            for b in molecule.bonds[mol]:
                if b.bond_atoms in coord_bonds[mol]:
                    b.order = None
                    b.convention = "cxn:coord"


def separate_components_molecule(
        mol: MrvData, return_atoms_bonds: bool = False, cleave_coord_bonds: bool = False
):
    """
    Separate a MrvData molecule into its connected components.
    :param mol: MrvData object.
    :param return_atoms_bonds: only retrun a atoms and bonds (dict).
    :param cleave_coord_bonds: Ignore coordination bonds and hence cleave them.
    :return: List of MrvData objects
    """
    gml, ids = mol.gml_string(renumber_ids=False, return_ids=True, cleave_coord_bonds=cleave_coord_bonds)

    conected_comp = {key: [] for key in mol.mol_ids}
    mol_components = []

    node_pattern = re.compile("node \[ id (\d+) label \"(.+)\" \]")
    edge_pattern = re.compile("edge \[ source (\d+) target (\d+) label \"(.+)\"\]")

    for key in gml:
        if len(gml[key]) == 1:
            # If there is only one connected component
            mol_components.append(mol)
            continue

        for g in gml[key]:
            collected_atoms = []
            collected_bonds = []
            for line in g.splitlines():
                node_result = re.search(node_pattern, line)
                edge_result = re.search(edge_pattern, line)
                if node_result is not None:
                    n_id, n_label = node_result.groups()
                    atom = lookup(mol.atoms[key], "id", ids[key][int(n_id)])
                    if atom is not False:
                        collected_atoms.append(atom)
                if edge_result is not None:
                    e_id1, e_id2, e_label = edge_result.groups()
                    e_a1 = lookup(mol.atoms[key], "id", ids[key][int(e_id1)])
                    e_a2 = lookup(mol.atoms[key], "id", ids[key][int(e_id2)])
                    if (e_a1 is not False) and (e_a2 is not False):
                        bond = lookup_bond(mol.bonds[key], {e_a1, e_a2})
                        if bond is not False:
                            collected_bonds.append(bond)

            # Creating new MrvData
            mrv_str = mol.mrv_string(only_these_atoms=collected_atoms)
            new_mol = read_from_string(mrv_str)
            # Renumber the atoms
            for akey in new_mol.atoms:
                try:
                    for ai, aa in enumerate(new_mol.atoms[akey], start=1):
                        aa.id = f"a{ai}"
                except TypeError as err:
                    continue
            mol_components.append(new_mol)
            conected_comp[key].append({"atoms": collected_atoms, "bonds": collected_bonds})

    if return_atoms_bonds is True:
        return conected_comp

    return mol_components


def get_electron_configuration(atom: Atom, verbose: int = 0):
    """
    Function to return Dictionary with electron configuration.
    :param atom: an Atom class object.
    :param verbose: verbose output.
    :return: Dictionary with electron configuration information.
    """
    f = re.compile("(\d)f(\d+)")
    d = re.compile("(\d)d(\d+)")
    s = re.compile("(\d)s(\d)")
    p = re.compile("(\d)p(\d)")
    try:
        e_conf = periodic_table[atom.element_type]["electron_configuration"]
        e_conf_sem = periodic_table[atom.element_type]["electron_configuration_semantic"]
        group = periodic_table[atom.element_type]["xpos"]
    except KeyError:
        if atom.element_type not in ["R", "*"]:
            warning_m(f"get_electron_configuration(): Unknown element: '{atom.element_type}'", verbose=verbose)
        return None
    e_configurations = {
        "electron_configuration": e_conf,
        "electron_configuration_semantic": e_conf_sem,
        "f": re.findall(f, e_conf),
        "d": re.findall(d, e_conf),
        "s": re.findall(s, e_conf),
        "p": re.findall(p, e_conf),
        "f_out": re.findall(f, e_conf_sem),
        "d_out": re.findall(d, e_conf_sem),
        "s_out": re.findall(s, e_conf_sem),
        "p_out": re.findall(p, e_conf_sem),
        "total": 0,
        "total_valence": 0
    }
    for conf in ["f", "d", "s", "p"]:
        if len(e_configurations[conf]) > 0:
            e_configurations[conf] = [(int(item[0]), int(item[1])) for item in e_configurations[conf]]
            e_configurations["total"] += sum([item[1] for item in e_configurations[conf]])
        if len(e_configurations[f"{conf}_out"]) > 0:
            e_configurations[f"{conf}_out"] = [(int(item[0]), int(item[1])) for item in e_configurations[f"{conf}_out"]]
            if (1 <= group <= 3) or (13 <= group <= 18):
                e_configurations["total_valence"] += sum(
                    [int(item[1]) for item in e_configurations[f"{conf}_out"] if conf in ["s", "p"]])
            else:
                e_configurations["total_valence"] += sum([int(item[1]) for item in e_configurations[f"{conf}_out"]])

    return e_configurations


def calculate_formal_charge(atom: Atom, verbose: int = 0):
    """
    Function to calculate the formal charge of an atom in a molecule.
    FC = V - N - (B/2)
    where FC is the formal charge, V the valence electrons (ground state), N the number
    electrons not engaged in a bond and B the number of electrons shared in bonds.
    :param atom: Atom class object
    :param verbose: verbose output.
    :return: formal charge (int)
    """
    v = get_electron_configuration(atom, verbose=verbose)
    if v is None:
        return None
    v = v["total_valence"]
    n = 0  # number of e not engaged in bonds

    # look for MElectronContainers
    # CAUTION: This method is not accurate enough
    # some atoms have lone pairs but no corresponding electrons
    # counting lone pairs is better
    if atom.lone_pair is not None:
        n += atom.lone_pair * 2
    if atom.radical == "monovalent":
        n += 1
    if (atom.radical == "divalent") or (atom.radical == "divalent1"):
        n += 2
    if atom.radical == "trivalent":
        n += 3
    if verbose is True:
        if n != 0:
            if (atom.lone_pair is None) and (atom.radical != "monovalent"):
                warning_m(
                    f"{atom.id} {atom.element_type}: No lone pairs or radical but {n} free electron(s).",
                    function_name=inspect.stack()[0][3]
                )
    b = 0

    # check for coordination bonds
    for ab in atom.bonds():
        if ab.order is not None:
            b += ab.order

    # Assuming that always two electrons are shared in a bond
    # The number of bonds equals the number of electrons shared in bonds devided by 2
    # print(f"{v - n - b} = {v} - {n} - {b}")
    return v - n - b


def atom_formal_charge_check(atom: Atom, verbose: int = 0):
    """
    Check whether an atom formal charge is OK.
    FC = V - N - (B/2)
    where FC is the formal charge, V the valence electrons (ground state), N the number
    electrons not engaged in a bond and B the number of bonds.
    :param atom: Atom class object
    :param verbose: verbose output.
    :return: True False
    """

    # no unknown elements
    if atom.element_type not in periodic_table:
        return None
    # no aliases (amino acid labels do not get protonated)
    if atom.mrv_alias is not None:
        return None
    # no pseudos (e.g., ubiquitin), they are just placeholders
    if atom.mrv_pseudo is not None:
        return None
    # no query properties (e.g., A etc)
    if atom.mrv_query_props is not None:
        return None
    # no transition metals. Too tricky at the moment
    if atom.element_type in transition_metals:
        return None

    # calculate formal charge
    calc_fc = calculate_formal_charge(atom, verbose=verbose)
    if calc_fc is None:
        return None

    # compare to existing formal charge
    # Does match
    if calc_fc == atom.formal_charge:
        return True

    # Does not match
    elif calc_fc != atom.formal_charge:
        # look for "Any" bonds
        for b in atom.bonds():
            if b.query_type == "Any":
                return None
        return False


def mass_from_gml(gml_string: str, verbose: int = 1):
    """
    Calculate molecular mass from GML file.
    Might not be accurate due to inaccurate atom labels.
    :param gml_string: GML string.
    :param verbose: verbose output.
    """
    tot_mass = 0
    # look for Alias labels
    if re.search(re.compile("Alias\(([A-Z][a-z]*),.*\)"), gml_string) is not None:
        gml_string = replace_special_labels_in_gml(gml_string, "Alias")
    if re.search(re.compile("Amino\(([A-Z][a-z]*),.*\)"), gml_string) is not None:
        gml_string = replace_special_labels_in_gml(gml_string, "Amino")
    pattern = re.compile("node \[\s*id \d+ label \"([A-Za-z]+).*\"\s*\]")
    result = re.findall(pattern, gml_string)
    for atom in result:
        try:
            mass = periodic_table[atom]["atomic_mass"]
            tot_mass += mass
        except KeyError as err:
            if atom != "R":
                warning_m(f"Skipping unknown element \"{atom}\" for mass calculation.",
                          function_name=inspect.stack()[0][3],
                          verbose=verbose)

    return tot_mass


def replace_special_labels_in_gml(gml_string: str, which: str):
    """
    Function for replacing Amino or Alias labels with their respective atom labels.
    :param gml_string: GML string with line breaks.
    :param which: Amino or Alias
    """

    pattern = {
        "Amino": re.compile("Amino\(([A-Z][a-z]*),.*\)"),
        "Alias": re.compile("Alias\(([A-Z][a-z]*),.*\)")
    }
    if which not in pattern:
        return gml_string

    gml_string = gml_string.splitlines()
    for li, line in enumerate(gml_string):
        if re.search(pattern[which], line) is not None:
            gml_string[li] = re.sub(pattern[which], r"\1", line)
    gml_string = "\n".join(gml_string)

    return gml_string


def decompose_amino_label(label: str):
    """
    Function to extract Amino() label components.
    :param label: Amino(atom, amino acid, reisdue id, chain id)
    :return: tuple of (name, resn, resi, chain)
    """
    amino_pattern = re.compile("Amino\(([A-Z][a-z]*),\s*([A-Za-z]+),\s*(\d+),\s*([A-Za-z]+)\)")
    result = re.search(amino_pattern, label)
    if result is None:
        return None

    return result.groups()


def decompose_alias_label(label: str):
    """
    Function to extract Alias() label components.
    :param label: Alias(atom, alias)
    :return: tuple of (name, alias)
    """
    alias_pattern = re.compile("Alias\(([A-Z][a-z]*),\s*(\S+)\)")
    result = re.search(alias_pattern, label)
    if result is None:
        return None

    return result.groups()
