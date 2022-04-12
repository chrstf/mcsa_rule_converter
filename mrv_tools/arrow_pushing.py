from mrv_tools import *
from mrv_tools.tools import warning_m
from mrv_tools.chemistry import *


def arrow_pushing(mol: MrvData, do_transition_metals: bool = False, do_single_electron: bool = False,
                  exclude_fc_mismatches: bool = False, verbose: int = 0, do_coord_bonds: bool = False):
    """
    Apply arrow pushes
    :param mol: The molecule
    :param do_transition_metals: Include transition metals
    :param do_single_electron: Include single barbed arrows
    :param do_coord_bonds: Include coordination bonds
    :param exclude_fc_mismatches: Exclude reactions where e-flow atoms have a wrong formal charge
    :param verbose: Increase verbosity
    """
    e_flow_done = []  # for keeping track of EFlow containers
    single_e_head_flags = ["2", "4"]
    abord_message = f"{mol.origin}: Skipping arrow pushes."

    # Check whether something is to be done
    if mol.e_flow_containers is None:
        message(f"No arrows to push.\n{mol.__str__()}\n{mol.origin}", verbose=verbose, verbose_level_threshold=2)
        return None

    ef_atoms = mol.e_flow_atoms()

    # Formal charge check
    if exclude_fc_mismatches is True:
        for atom in ef_atoms:
            if atom_formal_charge_check(atom) is False:
                warning_m("--\nAtom formal charge porbably worng. \n"
                          f"{mol.origin} \n"
                          f"{atom.__str__()}\n"
                          f"{abord_message}\n--", verbose=verbose, verbose_threshold=2)
                return False

    # Check for transition metal in MEFlow
    if mol.e_flow_transition_metal is not False:
        message(f"{mol.origin} contains at least one transition metal.", verbose=verbose, verbose_level_threshold=2)
        for atom in ef_atoms:
            if atom.element_type in transition_metals:
                message(f"Transition metal ({tm_label(atom)}) among electron flow atoms.", verbose=verbose,
                        verbose_level_threshold=2)
                if do_transition_metals is False:
                    warning_m(abord_message, verbose=verbose, verbose_threshold=2)
                    return None

    # Check for single e pushes aka special arrow heads
    if len(mol.head_flags) > 0:
        message(f"Arrow heads: {', '.join(map(str, mol.head_flags))}", c="BLUE", verbose=verbose,
                verbose_level_threshold=2)
        if do_single_electron is False:
            warning_m("Detected single barbed arrow(s).\n"
                      f"{abord_message}", verbose=verbose, verbose_threshold=2)
            return None
        # Collect only single barbed arrow pushes
        if evaluate_single_barbed(mol, verbose=verbose) is False:
            return False

    # Check for coordination bonds being part of the arrow pushes
    if mol.coordination_bonds is not None and do_coord_bonds is False:
        for ef in mol.e_flow_containers:
            for _atomset in ef.atom_set_points:
                _atoms = _atomset.atom_refs
                if len(_atoms) == 2:
                    coord_check = lookup_bond(mol.flatten_bonds(), set(_atoms))
                    if coord_check is False:
                        continue
                    if coord_check.convention is None:
                        continue
                    if "cxn:coord" in coord_check.convention:
                        warning_m("Detected a coordination bond among atom set points.\n"
                                  f"{abord_message}", verbose=verbose, verbose_threshold=2)
                        return False

    # Full arrow head pushes
    # Loop through the MEFlow containers
    for ef in mol.e_flow_containers:
        # 1) atom set to atom set
        # 2) atom set to single atom
        # 3) single atom to atom set
        # 4) signle atom to signle atom
        # 5) single e: atom set to atom set
        # 6) single e: atom set to single atom
        # 7) single e: single atom to atom set
        # 8) single e: signle atom to signle atom
        rxn = 0
        if ef.flow_base_point is not None:
            # there should be only one atom_set_point
            if len(ef.atom_set_points[0].atom_refs) == 2:
                # single atom to atom set
                if ef.head_flags in single_e_head_flags:
                    rxn = 7
                    if base_to_set_s_e(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)
                else:
                    rxn = 3
                    if base_to_set(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)
            elif len(ef.atom_set_points[0].atom_refs) == 1:
                # signle atom to signle atom
                if ef.head_flags in single_e_head_flags:
                    rxn = 8
                    if single_to_single_s_e(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)
                else:
                    rxn = 4
                    if single_to_single(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)

        elif ef.flow_base_point is None:
            # there should be two atom_set_points
            # The frist should be a pair
            if len(ef.atom_set_points[1].atom_refs) == 1:
                # atom set to single atom
                if ef.head_flags in single_e_head_flags:
                    rxn = 6
                    if set_to_single_s_e(ef, mol, verbose=verbose):
                        e_flow_done.append(ef)
                else:
                    rxn = 2
                    if set_to_single(ef, mol, verbose=verbose):
                        e_flow_done.append(ef)
            elif len(ef.atom_set_points[1].atom_refs) == 2:
                # atom set to atom set
                if ef.head_flags in single_e_head_flags:
                    rxn = 5
                    if set_to_set_s_e(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)
                else:
                    rxn = 1
                    if set_to_set(ef, mol, verbose=verbose) is True:
                        e_flow_done.append(ef)
        message(f"Reaction: {rxn}", verbose=verbose, verbose_level_threshold=3)

    # check if all e flows are done:
    ef_diff = len(mol.e_flow_containers) - len(e_flow_done)
    if ef_diff != 0:
        if ef_diff == 1:
            warning_m(f"{ef_diff} arrow push was not successful.", verbose=verbose, verbose_threshold=2)
        else:
            warning_m(f"{ef_diff} arrow pushes were not successful.", verbose=verbose, verbose_threshold=2)
        return False

    # renumber bonds
    mol.reset_bond_ids()

    # Set all negative lone pairs to 0
    for atom in mol.e_flow_atoms():
        if atom.lone_pair < 0:
            atom.lone_pair = 0

    # remove e_flow containers
    mol.e_flow_containers = [efc for efc in mol.e_flow_containers if efc not in e_flow_done]
    # return keep_this

    # Formal charge check
    # exclude_fc_mismatches = False
    if exclude_fc_mismatches is True:
        for atom in ef_atoms:
            if atom_formal_charge_check(atom) is False:
                warning_m("--\nAtom formal charge porbably worng after arrow pushes. \n"
                          f"{mol.origin} \n"
                          f"{atom.__str__()}\n"
                          f"{abord_message}\n--", verbose=verbose, verbose_threshold=2)
                return False

    return True


def base_to_set(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Electron flow from single to bond. No special head flags.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 3
    base_point: Atom = e_flow_container.flow_base_point.atom_ref
    # In case of a flow base point there is only one set point
    set_points = e_flow_container.atom_set_points[0]
    s_point_atoms = set_points.atom_refs

    # select atom where the electron pair from the base is going to
    # target_atom: list = [a for a in s_point_atoms if a.id != base_point.id]
    target_atom: set = set(s_point_atoms) - {base_point}
    if len(target_atom) > 1:
        error_m("Multiple target atoms")
    target_atom: Atom = list(target_atom)[0]
    # verbose graphics
    message(f"base-to-set arrow push\n"
            f"┌——┐\n"
            f"│  ▼\n"
            f"{base_point.element_type}----{target_atom.element_type}",
            verbose=verbose, verbose_level_threshold=3)

    # Bond
    increase_bond(s_point_atoms, mol, verbose=verbose)

    # Atoms
    update_atom_e_loss(base_point, verbose=verbose)
    update_atom_e_gain(target_atom, verbose=verbose)

    return True


def set_to_single(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Electron flow from a bond to a single atom. No special head flags.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # reaction 2
    arrow_head = e_flow_container.head_flags
    sp1: AtomSetPoint = e_flow_container.atom_set_points[0]
    sp2: AtomSetPoint = e_flow_container.atom_set_points[1]
    the_atoms = set(sp1.atom_refs).union(set(sp2.atom_refs))
    sp1_atoms: list = sp1.atom_refs
    sp2_atom: Atom = sp2.atom_refs[0]  # There is only one atom
    # print(sp2_atom.element_type, sp2_atom.formal_charge, calculate_formal_charge(sp2_atom))

    if len(the_atoms) == 2:
        # This means that the target atom is in the bond
        # One atom looses an electron pair

        # verbose graphics
        v_a1 = list(set(sp1_atoms).difference({sp2_atom}))[0].element_type
        v_a2 = list(set(sp1_atoms).intersection({sp2_atom}))[0].element_type
        message(f"set-to-single (2 atoms) push\n"
                f"{len(v_a1) * ' '} ┌—┐\n"
                f"{len(v_a1) * ' '} │ ▼\n"
                f"{v_a1}———{v_a2}",
                verbose=verbose, verbose_level_threshold=3)

        leaving_atom: Atom = list(the_atoms - {sp2_atom})[0]

        # get bond information
        sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))

        # Modifying bond
        decrease_bond(sp1_b, mol, verbose=verbose)

        # update sp2
        update_atom_e_gain(sp2_atom, verbose=verbose)

        # update leaving atom
        update_atom_e_loss(leaving_atom, verbose=verbose)

        return True

    elif len(the_atoms) == 3:
        # This means that the electron pair leaving the bond will go to a third atom.

        # verbose graphics
        v_a1 = sp1_atoms[0].element_type
        v_a2 = sp1_atoms[1].element_type
        v_a3 = sp2_atom.element_type
        message(f"set-to-single (3 atoms) push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}————┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}    ▼\n"
                f"{v_a1}———{v_a2}   {v_a3}",
                verbose=verbose, verbose_level_threshold=3)

        sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))

        # Modifying bond
        decrease_bond(sp1_b, mol, verbose=verbose)

        # update atoms
        # source atoms
        for atom in sp1_atoms:
            update_atom_e_loss(atom, verbose=verbose)
        # target atom
        update_atom_e_gain(sp2_atom, verbose=verbose)
        return True


def set_to_set(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Electron flow between two atom pairs. No special head flags.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 1
    # Collect atom information
    sp1: AtomSetPoint = e_flow_container.atom_set_points[0]
    sp2: AtomSetPoint = e_flow_container.atom_set_points[1]
    sp1_atoms: list = sp1.atom_refs
    sp2_atoms: list = sp2.atom_refs
    case = len(set(sp1_atoms).intersection(set(sp2_atoms)))
    # case 0: no common atoms
    # case 1: one atom is common between the two sets

    # Collect bond information
    sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))
    sp2_b: Bond = lookup_bond(mol.bonds[sp2_atoms[0].mol], set(sp2_atoms))

    # check if sp2_b exists
    v_b2 = "---"
    if sp2_b is False:
        v_b2 = "———"

    # verbose graphics
    # case 0
    if case == 0:
        v_a1 = sp1_atoms[0].element_type
        v_a2 = sp1_atoms[1].element_type
        v_a3 = sp2_atoms[0].element_type
        v_a4 = sp2_atoms[1].element_type
        message(f"set-to-set push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}————{(len(v_a3) - 1) * '-'}——┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}    {(len(v_a3) - 1) * '-'}  ▼\n"
                f"{v_a1}———{v_a2}   {v_a3}{v_b2}{v_a4}",
                verbose=verbose, verbose_level_threshold=3)
    # case 1
    if case == 1:
        v_a2 = set(sp1_atoms).intersection(set(sp2_atoms))
        v_a1 = list(set(sp1_atoms).difference(set(v_a2)))[0].element_type
        v_a4 = list(set(sp2_atoms).difference(set(v_a2)))[0].element_type
        v_a2 = list(v_a2)[0].element_type
        v_a3 = v_a2
        message(f"set-to-set push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}——┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}  ▼\n"
                f"{v_a1}———{v_a2}{v_b2}{v_a4}",
                verbose=verbose, verbose_level_threshold=3)
    # sp1_b modification
    decrease_bond(sp1_b, mol, verbose=verbose)
    # sp2_b modification
    increase_bond(sp2_atoms, mol, verbose=verbose)

    # update atoms
    # sp1, e loss
    for atom in sp1_atoms:
        update_atom_e_loss(atom, verbose=verbose)
    # sp2, e gain
    for atom in sp2_atoms:
        update_atom_e_gain(atom, verbose=verbose)
    return True


def single_to_single(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Electon flow between two atoms. No special flags.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 4

    warning_m(f"atom-to-atom arrow push in {mol.origin}. "
              "Verify if this is a valid move.\n", verbose=verbose, verbose_threshold=2)
    e_source: Atom = e_flow_container.flow_base_point.atom_ref
    e_target: Atom = e_flow_container.atom_set_points[0].atom_refs[0]
    # verbose graphics
    v_a1 = e_source.element_type
    v_a2 = e_target.element_type
    message(f"single-to-single push\n"
            f"┌{(len(v_a1) - 1) * '—'}———┐\n"
            f"│{(len(v_a1) - 1) * ' '}   ▼\n"
            f"{v_a1}   {v_a2}",
            verbose=verbose, verbose_level_threshold=3)

    e_source.formal_charge += 1
    e_target.formal_charge -= 1

    # TODO: finish
    return True


def evaluate_single_barbed(mol: MrvData, verbose: int = 0):
    sb_flows = []
    # Collect EFlow containers with single barbed arrows
    for ef in mol.e_flow_containers:
        if ef.head_flags in ["2", "4"]:
            sb_flows.append(ef)

    # Gather cases
    cases = []
    for ef in sb_flows:
        # cases
        # (1,1) atom to atom
        # (1,2) atom to bond
        # (2,1) bond to atom
        # (2,2) bond to bond

        # check for start
        origin = None
        target = None
        if ef.flow_base_point is not None:
            # e comes from one atom
            origin = [ef.flow_base_point.atom_ref]
            # There should only be one set point
            target = ef.atom_set_points[0].atom_refs
        else:
            # e comes from a bond
            origin = ef.atom_set_points[0].atom_refs
            target = ef.atom_set_points[1].atom_refs

        # case = len(target), len(origin)
        cases.append([origin, target])

    # Check cases
    missing_arrow_msg = "Missing single-barbed arrow(s)"
    # If (2, 1) or (2, 2) this means that a bond is being split up.
    # Therefore, we need matching origin bonds.
    from_bonds = [c for c in cases if len(c[0]) == 2]
    # This should be pair
    # if len(from_bonds) % 2 != 0:
    #     warning_m(f"{missing_arrow_msg}: bond-to-single or bond-to-bond")
    #     return False
    # check for matching bond
    origin_bond = [set(c[0]) for c in from_bonds]
    if len(origin_bond) > 0:
        for b in origin_bond:
            if origin_bond.count(b) != 2:
                warning_m(f"{missing_arrow_msg}: "
                          f"Missing pair for origin bond {' '.join([a.str_short() for a in b])}",
                          verbose=verbose, verbose_threshold=2)
                return False

    # If (1, 2) or (2, 2) this means that a bond is being created.
    # We need matching target bonds
    to_bonds = [c for c in cases if len(c[1]) == 2]
    # This should be pair
    # if len(to_bonds) % 2 != 0:
    #     warning_m(f"{missing_arrow_msg}: single-to-bond or bond-to-bond")
    #     return False
    # check for matching bond
    target_bond = [set(c[1]) for c in to_bonds]
    if len(to_bonds) > 0:
        for b in target_bond:
            if target_bond.count(b) != 2:
                warning_m(f"{missing_arrow_msg}: "
                          f"Missing pair for target bond {' '.join([a.str_short() for a in b])}",
                          verbose=verbose, verbose_threshold=2)
                return False

    return True


def base_to_set_s_e(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Single electron flow from atom to bond.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 7
    base_point: Atom = e_flow_container.flow_base_point.atom_ref
    # In case of a flow base point there is only one set point
    set_points = e_flow_container.atom_set_points[0]
    s_point_atoms = set_points.atom_refs

    # select atom where the electron is going to
    target_atom: set = set(s_point_atoms) - {base_point}
    if len(target_atom) > 1:
        warning_m("Multiple target atoms", verbose=verbose, verbose_threshold=2)
        return False
    target_atom: Atom = list(target_atom)[0]
    # verbose graphics
    message(f"single electron: base-to-set arrow push\n"
            f"┌——┐\n"
            f"│  │/\n"
            f"{base_point.element_type}----{target_atom.element_type}",
            verbose=verbose, verbose_level_threshold=3)

    # Set new base atom radical
    update_atom_s_e_loss(base_point, verbose=verbose)
    # modify for bond
    if increase_bond_s_e(s_point_atoms, mol, verbose=verbose) is False:
        return False

    return True


def set_to_set_s_e(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Single electron flow between two atom pairs.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 5
    # Collect atom information
    sp1: AtomSetPoint = e_flow_container.atom_set_points[0]
    sp2: AtomSetPoint = e_flow_container.atom_set_points[1]
    sp1_atoms: list = sp1.atom_refs
    sp2_atoms: list = sp2.atom_refs
    case = len(set(sp1_atoms).intersection(set(sp2_atoms)))
    # case 1: one atom is common between the two sets
    # case 0: no common atoms

    # Collect bond information
    sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))
    sp2_b: Bond = lookup_bond(mol.bonds[sp2_atoms[0].mol], set(sp2_atoms))

    # check if sp2_b exists
    v_b2 = "---"
    if sp2_b is True:
        v_b2 = "———"
        if sp2_b.order == 2:
            v_b2 = "═══"

    # verbose graphics
    # case 0
    if case == 0:
        v_a1 = sp1_atoms[0].element_type
        v_a2 = sp1_atoms[1].element_type
        v_a3 = sp2_atoms[0].element_type
        v_a4 = sp2_atoms[1].element_type
        message(f"set-to-set push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}————{(len(v_a3) - 1) * '-'}——┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}    {(len(v_a3) - 1) * '-'}  │/\n"
                f"{v_a1}———{v_a2}   {v_a3}{v_b2}{v_a4}",
                verbose=verbose, verbose_level_threshold=3)
    # case 1
    if case == 1:
        v_a2 = set(sp1_atoms).intersection(set(sp2_atoms))
        v_a1 = list(set(sp1_atoms).difference(set(v_a2)))[0].element_type
        v_a4 = list(set(sp2_atoms).difference(set(v_a2)))[0].element_type
        v_a2 = list(v_a2)[0].element_type
        message(f"set-to-set push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}——┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}  │/\n"
                f"{v_a1}———{v_a2}{v_b2}{v_a4}",
                verbose=verbose, verbose_level_threshold=3)

    # sp2_b modification
    sp2_b_mod = increase_bond_s_e(sp2_atoms, mol, verbose=verbose)
    if sp2_b_mod is False:
        return False

    # sp1_b modification
    sp1_b_mod = decrease_bond_s_e(sp1_b, mol, verbose=verbose)

    # radical
    # lone pairs
    if case == 1:
        # nothing should change for the atoms at this point.
        # Charges, radical etc will be take care by the matching functions.
        # e.g. set_to_single_s_e() or base_to_set_s_e()
        return True
    if case == 0:
        warning_m("Unknown case: single electron transfer with 4 different atoms.",
                  function_name=inspect.stack()[0][3],
                  verbose=verbose, c="RED", verbose_threshold=2)
        return False


def set_to_single_s_e(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Single electron flow from a bond to a single atom
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # reaction 6
    arrow_head = e_flow_container.head_flags
    sp1: AtomSetPoint = e_flow_container.atom_set_points[0]
    sp2: AtomSetPoint = e_flow_container.atom_set_points[1]
    the_atoms = set(sp1.atom_refs).union(set(sp2.atom_refs))
    sp1_atoms: list = sp1.atom_refs
    sp2_atom: Atom = sp2.atom_refs[0]  # There is only one atom

    if len(the_atoms) == 2:
        # This means that the target atom is in the bond

        # verbose graphics
        v_a1 = list(set(sp1_atoms).difference({sp2_atom}))[0].element_type
        v_a2 = list(set(sp1_atoms).intersection({sp2_atom}))[0].element_type
        message(f"set-to-single (2 atoms) push\n"
                f"{len(v_a1) * ' '} ┌—┐\n"
                f"{len(v_a1) * ' '} │ │/\n"
                f"{v_a1}———{v_a2}",
                verbose=verbose, verbose_level_threshold=3)

        leaving_atom: Atom = list(the_atoms - {sp2_atom})[0]

        # get bond information
        sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))

        # sp1_b modification
        decrease_bond_s_e(sp1_b, mol, verbose=verbose)

        # sp2 adjustment
        update_atom_s_e_gain(sp2_atom, verbose=verbose)
        return True

    elif len(the_atoms) == 3:
        v_a1 = sp1_atoms[0].element_type
        v_a2 = sp1_atoms[1].element_type
        v_a3 = sp2_atom.element_type
        message(f"set-to-single (3 atoms) push\n"
                f"{len(v_a1) * ' '} ┌—{(len(v_a2) - 1) * '-'}————┐\n"
                f"{len(v_a1) * ' '} │ {(len(v_a2) - 1) * ' '}    │/\n"
                f"{v_a1}———{v_a2}   {v_a3}",
                verbose=verbose, verbose_level_threshold=3)
        sp1_b: Bond = lookup_bond(mol.bonds[sp1_atoms[0].mol], set(sp1_atoms))
        decrease_bond_s_e(sp1_b, mol, verbose=verbose)
        update_atom_s_e_gain(sp2_atom, verbose=verbose)

        warning_m("Unknown case: Single electron transfer "
                  "from bond to single with three different atoms",
                  verbose=verbose, verbose_threshold=2,
                  function_name=inspect.stack()[0][3], c="RED")
        return False


def single_to_single_s_e(e_flow_container: EFlow, mol: MrvData, verbose: int = 0):
    """
    Single electon flow between two atoms.
    :param e_flow_container: Electron flow container.
    :param mol: The Molecule.
    :param verbose: Verbose.
    :return: True False
    """
    # Reaction 8
    e_source: Atom = e_flow_container.flow_base_point.atom_ref
    e_target: Atom = e_flow_container.atom_set_points[0].atom_refs[0]
    # verbose graphics
    v_a1 = e_source.element_type
    v_a2 = e_target.element_type
    message(f"single-to-single push\n"
            f"┌{(len(v_a1) - 1) * '—'}———┐\n"
            f"│{(len(v_a1) - 1) * ' '}   │/\n"
            f"{v_a1}   {v_a2}",
            verbose=verbose, verbose_level_threshold=3)

    update_atom_s_e_loss(e_source, verbose=verbose)
    update_atom_s_e_gain(e_target, verbose=verbose)

    return True


def remove_bond(bond: Bond, mol: MrvData, verbose: int = 0):
    this_bond: int = lookup_bond(mol.bonds[bond.mol], bond.bond_atoms, return_index=True)
    mol.bonds[bond.mol].pop(this_bond)
    message(f"Removed bond: {bond.alt_str()}", verbose=verbose, verbose_level_threshold=3)
    return True


def decrease_bond(bond: Bond, mol: MrvData, verbose: int = 0):
    """
    Function for decreasing a bond.
    """
    # Modifying bond
    if (bond.order == 0) or (bond.order is None):
        remove_bond(bond, mol, verbose=verbose)
        return True
    bond.order -= 1
    if bond.order == 0:
        remove_bond(bond, mol, verbose=verbose)
    return True


def increase_bond(bond_atoms: list, mol: MrvData, verbose: int = 0):
    """
    Function for increasing the bond order.
    """
    # Check whether target bond exists
    bond_exists = lookup_bond(mol.bonds[bond_atoms[0].mol], set(bond_atoms))
    if bond_exists is False:
        # Create a new bond
        new_bond = create_bond_from_atom_list(bond_atoms)
        if new_bond is not False:
            mol.bonds[bond_atoms[0].mol].append(new_bond)
            message(f"Created bond: {new_bond.alt_str()}", verbose=verbose,
                    verbose_level_threshold=3)
            return True
        warning_m("Failed to create bond.", verbose=verbose, verbose_threshold=3)
        return False
    bond_exists.order += 1
    return True


def decrease_bond_s_e(bond: Bond, mol: MrvData, verbose: int = 0):
    """
    Function for decreasing the bond order for single arrow pushes
    """
    if np.sign(bond.single_e_flow_buildup) <= 0:
        if bond.single_e_flow_buildup % 2 == 0:
            bond.single_e_flow_buildup -= 1
            return None
        bond.order -= 1
        if (bond.order == 0) or (bond.order is None):
            remove_bond(bond, mol, verbose=verbose)
        return True
    return False


def increase_bond_s_e(bond_atoms: list, mol: MrvData, verbose: int = 0):
    """
    Function for increasing the bond order for single arrow pushes
    """
    # check if bond exists
    this_bond: Bond = lookup_bond(mol.flatten_bonds(), set(bond_atoms))
    if this_bond is False:
        # create bond
        new_bond: Bond = create_bond_from_atom_list(bond_atoms, verbose=verbose)
        new_bond.single_e_flow_buildup = 1
        if new_bond is not False:
            mol.bonds[bond_atoms[0].mol].append(new_bond)
            message(f"Created bond: {new_bond.alt_str()}", verbose=verbose, verbose_level_threshold=2)
            return True
        warning_m("Failed to create bond.", verbose=verbose, verbose_threshold=2)
        return False

    # modify bond
    if this_bond.single_e_flow_buildup % 2 != 0:
        # the bond is in the process of being changed.
        # Finalise it.
        this_bond.single_e_flow_buildup += 1
        return True
    # This is the first step in bond order change
    this_bond.single_e_flow_buildup += 1
    this_bond.order += 1

    return None


def update_atom_e_loss(atom: Atom, verbose: int = 0):
    # Increase foraml charge of base atom +1
    atom.formal_charge += 1
    # Decrease number of Lone pairs
    atom.lone_pair -= 1
    return True


def update_atom_e_gain(atom: Atom, verbose: int = 0):
    # Decrease formal charge on target atom
    atom.formal_charge -= 1
    # Increase number of lone pairs of sp2 by one
    # Finding group in the periodic table
    # CAUTION: This is probably only valid for the 2nd period. Maybe the 3rd.
    group = periodic_table[atom.element_type]["xpos"]
    if group >= 14:
        atom.lone_pair += 1
    return True


def update_atom_s_e_loss(atom: Atom, verbose: int = 0):
    if atom.radical == "monovalent":
        atom.radical = None
    elif atom.radical == "divalent":
        atom.radical = "monovalent"
    elif atom.radical == "trivalent":
        atom.radical = "divalent"
    elif atom.radical is None:
        # Adjust lone pairs
        if atom.lone_pair == 0:
            warning_m("No free electrons to shuffle.", verbose=verbose, verbose_threshold=2)
            return False
        atom.lone_pair -= 1
        atom.radical = "monovalent"
        # TODO: Formal charge
    return True


def update_atom_s_e_gain(atom: Atom, verbose: int = 0):
    if atom.radical is None or atom.radical == "0":
        atom.radical = "monovalent"
    elif atom.radical == "monovalent":
        atom.radical = "divalent"
    elif atom.radical == "divalent":
        atom.radical = "trivalent"
    return True
