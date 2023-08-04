# mcule analysis

This is a simple RDKit-based script to find instances of a given functional group in a list of SMILES strings.

Here, it's used to extract molecules containing given functional groups from the Mcule building blocks database.
(Functional groups are defined in ``matches.txt``.)

## Usage

To run these scripts for yourself:

```
unzip mcule_bb.smi.zip
mkdir output
python find_matches.py
```

This may take a little while.

## Visualization

See ``visualize.ipynb`` for examples of how to plot the output ``.smi`` files.

*Corin Wagen, 2023*
