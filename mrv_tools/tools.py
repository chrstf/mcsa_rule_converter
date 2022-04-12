import sys
import numpy as np

amino_acids = {
  "VAL": "V",
  "ILE": "I",
  "LEU": "L",
  "GLU": "E",
  "GLN": "Q",
  "ASP": "D",
  "ASN": "N",
  "HIS": "H",
  "TRP": "W",
  "PHE": "F",
  "TYR": "Y",
  "ARG": "R",
  "LYS": "K",
  "SER": "S",
  "THR": "T",
  "MET": "M",
  "ALA": "A",
  "GLY": "G",
  "PRO": "P",
  "CYS": "C"
}


def invert_dict(a_dict: dict):
    try:
        return {item[1]: item[0] for item in a_dict.items()}
    except TypeError:
        return {str(item[1]): str(item[0]) for item in a_dict.items()}


class BColors:
    # ANSI codes
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    BLINK = '\033[5m'

    ENDC = '\033[0m'
    NORMAL = '\033[0m'


def color(acolor: str, string: str):
    return getattr(BColors, acolor) + string + BColors.ENDC


def message(print_this: str, c: str = "NORMAL", verbose: int = 1, end: str = "\n", function_name: str = None,
            lead_symbol: str = "#", verbose_level_threshold: int = 1):
    """
    A function for priting messages.
    :param print_this: Message to print
    :param c: print color
    :param verbose: verbose
    :param verbose_level_threshold: This hurdle needs to be overcome in oder to get verbose
    :param end: end argument of the print() function
    :param function_name: For referring to a function name
    :param lead_symbol: Symbol prepended to the message
    """

    if verbose == 0:
        return None

    if verbose < verbose_level_threshold:
        return None

    if isinstance(print_this, list) is True:
        print_this = " ".join(map(str, print_this))
    elif isinstance(print_this, str) is not True:
        print_this = str(print_this)

    if function_name is not None:
        print_this = f"{function_name}(): {print_this}"

    lines = print_this.splitlines()
    if len(lines) == 1:
        print(f"{lead_symbol} {color(c, print_this)}", end=end)
    elif len(lines) > 1:
        for line_id, line in enumerate(lines):
            if function_name is not None:
                if line_id == 0:
                    print(f"{lead_symbol} {color(c, line)}")
                else:
                    print(f"{lead_symbol} {color(c, ' ' * (len(function_name) + 4) + line)}")
            else:
                print(f"{lead_symbol} {color(c, line)}")


def warning_m(print_this: str, c: str = "YELLOW", verbose: int = 0, function_name: str = None,
              verbose_threshold: int = 1):
    message(print_this,
            c=c,
            lead_symbol="@",
            verbose=verbose,
            function_name=function_name,
            verbose_level_threshold=verbose_threshold)


def error_m(print_this: str, c: str = "RED", function_name: str = None):
    message(print_this, c=c, lead_symbol="!", verbose=True, function_name=function_name)
    print(color(c, "Aborting!"))
    sys.exit()


def remove_empty_lines(a_string: str):
    """
    Function for removing empty lines
    :param a_string: A string
    :return: A string without empty lines
    """
    if len(a_string.splitlines()) > 0:
        new_content = []
        for item in a_string.splitlines():
            if len(item) > 0:
                item = item.strip()
                if len(item) > 0:
                    new_content.append(item)

        if len(new_content) > 0:
            a_string = "\n".join(new_content)
        else:
            a_string = ""
    return a_string


def lookup(obj_list: list, attr: str, value: str, return_index: bool = False):
    """
    Function for looking up an object in a list based on a given attribute value.
    :param obj_list: list of objects
    :param attr: object attriubte
    :param value: objec.attribute value we are looking for
    :param return_index: only return the list index
    :return: looked up object or False
    """
    for i, item in enumerate(obj_list):
        item_val = getattr(item, attr)
        if item_val == value:
            if return_index is True:
                return i
            elif return_index is False:
                return item

    return False


def lookup_bond(bonds: list, bond: set, return_index: bool = False):
    """
    Function for looking up a bond.
    :param bonds: a list of bond. (e.g., MrvData.bonds["m1"])
    :param bond: a set of atoms of the Atom class.
    :param return_index: only return the list index.
    :return: looked up object or False
    """
    for i, b in enumerate(bonds):
        if bond == b.bond_atoms:
            if return_index is True:
                return i
            return b

    return False


def progress_bar(current: int, list_len: int, columns: int, lead: str = "", verbose: int = True):
    """
    Progress bar.
    :param current: current list index.
    :param list_len: Length of the list.
    :param columns: Number of terminal columns.
    :param lead: Leading string.
    :param verbose: Verbose output.
    """
    symbol = "█"
    bar_end = "\r"
    if verbose is False:
        return None
    lead += f" {((float(current) + 1) / float(list_len)) * 100:5.1f} %"
    to_substract = len(lead)
    bar_range = np.linspace(0, columns - (to_substract+1), list_len, endpoint=True)
    bar = f"{lead} {symbol * (int(bar_range[current])-1)}"
    if (current + 1) == list_len:
        bar_end = "\n"
    print(f"{bar}{' ' * (columns - len(bar) - 1)}|", end=bar_end)
    return None


def flatten_list(alist: list):
    new_list = []
    for item in alist:
        new_list += item
    return new_list
