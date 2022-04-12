import os
import shutil
import subprocess
import multiprocessing as mp
import re
import mrv_tools
from mrv_tools.tools import message


def tool_exists(name):
    """Check whether `name` is on PATH and marked as executable."""
    return shutil.which(name) is not None


def collect_file_names(path):
    fnames = []
    for _, _, files in os.walk(path):
        fnames.extend(files)
    return fnames


def _add_h_molconvert(in_file, out_file):
    cmd = ["molconvert", "mrv:H", in_file, "-o", out_file]
    subprocess.call(cmd)


def _add_h_patch(in_file, out_file, patch):
    message(in_file)
    imp_mol = mrv_tools.read_from_file(in_file)

    # if there is a <molecule molID="m2"/> tag things can get messy
    # sice it usually comes after the molID="m1".
    # Therefore, Marvin renumbers the IDs of m2 after adding Hydrogen
    # atoms between m1 and m2.
    # A renumbering is needed
    not_m1 = [k for k in imp_mol.ismol if (imp_mol.ismol[k] is True) and (k != "m1")]
    update_these = None
    if len(not_m1) > 0:
        update_these = {}
        for m_key in not_m1:
            update_these.update({m_key: {"atoms": imp_mol.atoms[m_key]}})
            if imp_mol.bonds[m_key] is not None:
                update_these[m_key].update({"bonds": imp_mol.bonds[m_key]})

    for atom in patch["atoms"]:
        imp_mol.atoms[atom.mol].append(atom)
    for bond in patch["bonds"]:
        imp_mol.bonds[bond.mol].append(bond)

    if update_these is not None:
        # update atom and bond IDs
        num_pattern = re.compile(r"(\d+)")
        for m_key in update_these:
            prev_key = imp_mol.mol_ids[imp_mol.mol_ids.index(m_key) - 1]
            highest_atm = int(num_pattern.search(imp_mol.find_highest_id()[prev_key]).groups()[0])
            highest_bnd = int(num_pattern.search(imp_mol.find_highest_id(which="bonds")[prev_key]).groups()[0])
            for i, _ in enumerate(imp_mol.atoms[m_key]):
                imp_mol.atoms[m_key][i].id = f"a{highest_atm + i + 1}"
            for i, _ in enumerate(imp_mol.bonds[m_key]):
                imp_mol.bonds[m_key][i].id = f"a{highest_bnd + i + 1}"

    imp_mol.write_mrv(out_file, verbose=False)


def convert_with_molconvert(implicit_dir: str, explicit_dir: str, verbose: int, njobs: int):
    message("Adding explicit hydrogens with molconvert.", verbose=verbose)
    message("\"Stereochemical descriptor mismatch after 2D coordinate generation\" "
            "messages can be ignored.")
    file_names = collect_file_names(implicit_dir)
    pool_args = [(os.path.join(implicit_dir, fn), os.path.join(explicit_dir, fn)) for fn in file_names]
    with mp.Pool(processes=njobs) as pool:
        pool.starmap(_add_h_molconvert, pool_args)


def convert(implicit_dir: str, explicit_dir: str, verbose: int, njobs: int):
    message("Convert marvin files with implicit hydrogens to explicit.", verbose=verbose)

    if tool_exists("molconvert"):
        message("Found molconvert on system and will use it to add hydrogens.", verbose=verbose)
        convert_with_molconvert(implicit_dir, explicit_dir, verbose, njobs)
        return True
    else:
        message("did not find molconvert on system. using static patch instead.", verbose=verbose)
        return False


