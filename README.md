# Generate Rules from the M-CSA Database

This workflow downloads Marvin files from the M-CSA database and will at some point shape them into crude rules.

## First thing

Update the periodic table submodule.

```bash
git submodule update --init --recursive
```

## Dependencies

### Setting up a conda environment

```bash
conda env create -f environment.yml
```

#### Python Packages

- bs4
- networkx
- numpy
- requests

### Molecule File Converter
Molecule File Converter, version 20.20.0, (C) 1999-2020 ChemAxon Ltd.

This tool is used to add explicit hydrogen atoms to the molecules.

Must be accessible via `molconvert` in your `$PATH`.

## main.py script

Execute the following to download the enzyme mechanism steps from the M-CSA database,
add explicit Hydrogen atoms and create GML rules.
```commandline
ipython main.py -- -devj 4
```

For listing all available command line options execute the following.

```commandline
ipython main.py -- --help
```

## Output

The generated GML rules are stored in a json file named `YYYYMMDDHHMM_gml_rules.json`
(YYYYMMDDHHMM indicating the time stamp).

The command line arguments are stored in `YYYYMMDDHHMM_arguments.json`.
