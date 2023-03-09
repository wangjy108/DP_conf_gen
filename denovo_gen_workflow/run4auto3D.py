import os
import sys
import pandas as pd
from rdkit import Chem


def run(path_dir, sdf_tag):
    _smi = []
    name_list = []
    for ff in os.listdir(path_dir):
        if sdf_tag in ff:
            mol = [mm for mm in  Chem.SDMolSupplier(f"{path_dir}/{ff}", removeHs=False) if mm][0]
            _name = ff.split(".")[0]
            _smi.append(Chem.MolToSmiles(mol))
            name_list.append(_name)

    df = pd.DataFrame({"smi":_smi, "name": name_list})
    df.to_csv(os.path.join(path_dir,"input.smi"), sep='\t', index=None,header=None)
    return 
