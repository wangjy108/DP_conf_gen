import os
import sys
import pandas as pd
from rdkit import Chem
import argparse


def genSMI(sdf_tag):
    _smi = []
    name_list = []
    path_dir= os.getcwd()
    for ff in os.listdir(path_dir):
        if sdf_tag in ff:
            mol = [mm for mm in  Chem.SDMolSupplier(f"{path_dir}/{ff}", removeHs=False) if mm][0]
            _name = ff.split(".")[0]
            _smi.append(Chem.MolToSmiles(mol))
            name_list.append(_name)

    df = pd.DataFrame({"smi":_smi, "name": name_list})
    df.to_csv(os.path.join(path_dir,"input.smi"), sep='\t', index=None,header=None)
    return 

def prepare_auto3D_para(input_smi:str):
    work_dir=os.getcwd()
    full_input = os.path.join(work_dir, input_smi)
    with open("/home/jywang/tool_file/DP_conf_gen/denovo_gen_workflow/auto3D_param/parameters.yaml", 'r+') as ff:
        param_content = [cc for cc in ff.readlines()]
    
    with open(os.path.join(work_dir, "parameters.yaml"), "w+") as cc:
        for pp in param_content:
            if pp.startswith("path"):
                new_pp = f"path: {full_input}"
                cc.write(new_pp)
            else:
                cc.write(pp)
    
    return 

def run(sdf_tag):
    genSMI(sdf_tag)
    prepare_auto3D_para("input.smi")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='genSMI and prepare auto3D run')
    parser.add_argument('--sdf_tag', type=str, required=True, help='input serial sdf tag')
    #parser.add_argument('--input_smi', type=str, required=True, help='input smi file name')
    args = parser.parse_args()

    run(args.sdf_tag)





