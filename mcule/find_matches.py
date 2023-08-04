import re, tqdm, sys, os, time
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import rdBase
from rdkit.Chem import rdSubstructLibrary

def filter_smiles_list(input_filename, output_dir, match_filename, exclude=False):
    """
    Takes a file of SMILES strings and returns a file of SMILES strings that do or don't contain given substructures.
    """
    start = time.time()

    # create containers
    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()

    # extract molecules and add to containers - slow for large datasets (duh)
    print("loading input smiles")
    num_start = 0
    supplier = Chem.SmilesMolSupplier(input_filename, smilesColumn=0, nameColumn=0, titleLine=False, sanitize=True)
    for mol in tqdm.tqdm(supplier, position=0, leave=True, ncols=100):
        try:
            mol = Chem.AddHs(mol)
            fp = Chem.PatternFingerprint(mol)
            mols.AddSmiles(Chem.MolToSmiles(mol))
            fps.AddFingerprint(fp)
            num_start += 1
        except Exception as e:
            pass
    print(f"{num_start} molecules loaded")

    # create library and perform search
    print("performing substructure search")
    library = rdSubstructLibrary.SubstructLibrary(mols,fps)

    with open(match_filename, "r") as matches:
        for line in matches.readlines():
            match_name, match_smarts = line.split()

            substructure_mol = Chem.MolFromSmarts(match_smarts)
            match_indices = library.GetMatches(substructure_mol, maxResults=num_start)
            print(f"{len(match_indices)} matches found")

            # figure out which indices we want to keep
            match_indices = set(match_indices)
            total_indices = set(range(num_start))
            desired_indices = total_indices.intersection(match_indices)
            if exclude:
                desired_indices = total_indices - match_indices

            # turn indices into smiles (and save molecules)
            print("building result list")
            result_smiles = []
            result_mols = []
            for index in desired_indices:
                result_smiles.append(Chem.MolToSmiles(library.GetMol(index)))
                result_mols.append(library.GetMol(index))

            # write the smiles to output file
            print(f"writing {len(result_smiles)} results to {output_dir}/{match_name}.smi")
            with open(f"{output_dir}/{match_name}.smi", "w") as output:
                for smiles in tqdm.tqdm(result_smiles, position=0, leave=True, ncols=100):
                    output.write(f"{smiles}\n")

    end = time.time()
    elapsed = end - start
    print(f"total run time: {elapsed:.2f} s")

filter_smiles_list("mcule_bb.smi", "output", "matches.txt")
