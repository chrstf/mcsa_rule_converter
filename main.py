"""
Main script.
An example.
"""

import argparse
import mrv_tools as mt
import mrv_tools.mrv2gml
from mrv_tools.arrow_pushing import *
from os import walk, path, remove, makedirs
import re
import json
import multiprocessing as mp
import mrv_tools.mcsa.download
import mrv_tools.mcsa.implicit2explicit
from datetime import datetime


def remove_existing_data(this_folder: str):
    for _, _, ff in walk(this_folder):
        for ffi in ff:
            remove(path.join(this_folder, ffi))
            message(
                f"Removed {path.join(this_folder, ffi)}",
                c="YELLOW",
                verbose=args.verbose,
                verbose_level_threshold=3
            )


def process_mrv_file(fname: str):
    if path.split(fname)[-1] in blacklist:
        warning_m(f"{fname} blacklisted", verbose=args.verbose)
        return None
    mcsa_filename_match = mrv_fname_pattern.search(fname)
    entry, mechanism, step = None, None, None
    if mcsa_filename_match is not None:
        if mcsa_filename_match.group() in args.blacklist:
            message(f"{fname} is blacklisted. Skipping.", c="YELLOW", verbose=args.verbose, lead_symbol="!",
                    verbose_level_threshold=2)
            return None
        entry, mechanism, step = mrv_fname_pattern.search(fname).groups()
        message(f"entry: {entry:>4s}, proposal: {mechanism:1s}, step: {step:>2s}", c="BLUE", verbose=args.verbose,
                verbose_level_threshold=2)
    else:
        message(f"non M-CSA file {fname}", c="BLUE", verbose=args.verbose, verbose_level_threshold=2)

    mrv_reactants = mt.read_from_file(fname,
                                      origin=fname,
                                      verbose=args.verbose)
    mrv_products = MrvData(mrv_reactants.content, origin=fname)
    push_result = arrow_pushing(mrv_products,
                                do_transition_metals=False,
                                do_single_electron=False,  # is not correctly implemented yet
                                do_coord_bonds=False,
                                exclude_fc_mismatches=False,
                                verbose=args.verbose)
    if push_result is not True and args.exclude_failed is True:
        warning_m(f"{mrv_reactants.origin} failed", verbose=args.verbose, verbose_threshold=1)
        return None
    if not args.marvin_products:
        gml_rule = mt.mrv2gml.convert(mrv_reactants, mrv_products)
        if entry is not None:
            return {
                "entry": int(entry),
                "proposal": int(mechanism),
                "step": int(step),
                "gml": gml_rule
            }
        else:
            return {
                "gml": gml_rule
            }
    if not path.isdir(args.output):
        makedirs(args.output)
    # else:
    #     remove_existing_data(args.output)
    mrv_products.write_mrv(path.join(args.output, mrv_fname_pattern.search(fname)[0]), verbose=args.verbose)
    return None


def handle_input_download():
    entries = mrv_tools.mcsa.download.request_entry_data(args.verbose, args.njobs)
    if not path.exists(args.input):
        makedirs(args.input)

    with open(path.join(args.input, "entries.json"), "w") as out_file:
        json.dump(entries, out_file, indent=True)

    implicit_dir = path.join(args.input, "implicit")
    if not path.exists(implicit_dir):
        makedirs(implicit_dir)
    else:
        remove_existing_data(implicit_dir)

    mrv_tools.mcsa.download.request_marvin_files(args.input, entries, args.verbose, args.njobs)


def handle_input_explicit():
    explicit_dir = path.join(args.input, "explicit")
    if not path.exists(explicit_dir):
        makedirs(explicit_dir)
    else:
        remove_existing_data(explicit_dir)

    _result = mrv_tools.mcsa.implicit2explicit.convert(
        path.join(args.input, "implicit"),
        explicit_dir,
        args.verbose,
        args.njobs
    )
    return _result


class TimePoint:
    def __init__(self, _now):
        self.now = _now
        self.short_date = f'{_now.year:4d}{_now.month:02d}{_now.day:02d}'
        self.long = f'{_now.year:4d}{_now.month:02d}{_now.day:02d}{_now.hour:2d}{_now.minute:02d}'


def save_args(_loc: str = 'arguments.json'):
    # Saving arguments and timestamp
    global args
    global now
    _loc = f'{now.long}_{_loc}'
    _out: dict = args.__dict__
    _out.update({'timestamp': now.long})
    with open(_loc, 'w') as _f:
        json.dump(_out, _f, indent=True)
    message(f'Wrote arguments to \"{_loc}\"', verbose=args.verbose, verbose_level_threshold=3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Rule pipeline",
        usage="ipython main.py -- -devj 4",
        description="Script for generating GML rules from the M-CSA database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", type=str, default="steps",
                        help="Root folder of the steps. Based on the arguments "
                             "implicit, explicit and products folders are being created within this root folder.")
    parser.add_argument("-o", "--output", type=str, default="gml_rules.json",
                        help="Where the rules are saved. Default: gml_rules.json")
    parser.add_argument("-d", "--download", action="store_true",
                        help="If set, input files will be downloaded from MCSA DB. "
                             "In this case, the input argument works as the root "
                             "destination folder for the downloaded files.")
    parser.add_argument("-e", "--explicit", action="store_true",
                        help="Add explicit hydrogen atoms to the downloaded marvin files."
                             "the marvin files must be found in an 'implicit/' folder within the input arg")
    parser.add_argument("--marvin_products", action="store_true",
                        help="Output the products of each arrow pushing diagram as marvin files. "
                             "If this flag is set, the output must be a folder.")
    parser.add_argument("-j", "--njobs", type=int, default=4, help="number of processes to use.")
    parser.add_argument("--blacklist", type=str, default="blacklist.json",
                        help="List of steps that have issues and should be skipped. Default: data/blacklist.json")
    parser.add_argument("--specific_files", nargs="+", help="Specify files that should be processed")
    parser.add_argument("--specific_files_json", type=str, help="Specify files that should be processed (JSON list)")
    parser.add_argument("--exclude_failed", action="store_true",
                        help="Exclude failed arrow pushes. Aka, do not write empty GML rules")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase output verbosity by setting this flag multiple times.")
    args = parser.parse_args()

    now = TimePoint(datetime.now())

    if args.output == 'gml_rules.json':
        args.output = f'{now.long}_{args.output}'

    save_args()

    if args.download:
        handle_input_download()

    if args.download or args.explicit:
        if handle_input_explicit() is False:
            exit()

    if args.marvin_products is False:
        if not args.output.endswith(".json"):
            args.output = f"{args.output}.json"
    if args.marvin_products is True:
        if args.output == "gml_rules.json":
            args.output = "steps/products"

    steps_path = path.join(args.input, "explicit")
    mrv_fname_pattern = re.compile(r"(\d+)-(\d+)-(\d+).mrv")
    blacklist = json.load(open(args.blacklist, "r"))
    # Collecting file names
    fnames = []
    if args.specific_files_json is not None:
        with open(args.specific_files_json) as thisf:
            args.specific_files = json.load(thisf)
    if args.specific_files is None:
        for _, _, files in walk(steps_path):
            fnames += [path.join(steps_path, f) for f in files]
    else:
        fnames = args.specific_files
    # fnames = fnames[:20]
    message("Processing marvin files.", verbose=args.verbose)
    with mp.Pool(processes=args.njobs) as pool:
        out = pool.map(process_mrv_file, fnames)
    out = [item for item in out if item is not None]
    out.sort(key=lambda step: (step["entry"], step["proposal"], step["step"]) if "entry" in step else (-1, -1, -1))

    if not args.marvin_products:
        with open(args.output, "w") as f:
            json.dump(out, f, indent=2)
        message(f"Wrote {len(out)} rules to \"{args.output}\"", verbose=args.verbose, c="GREEN", lead_symbol=">")
        message(f"({len(out)} out of {len(fnames)} input files)", verbose=args.verbose, c="GREEN", lead_symbol=">")
