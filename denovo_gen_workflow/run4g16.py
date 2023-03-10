import os
import sys
import pandas as pd
from rdkit import Chem
import argparse

## gen gaussian input from auto3D gen
## defualt method B3LYP/6-311G* em=GD3BJ opt freq scrf(SMD,solvent=water)
##                B2PLYPD3/def2TZVP scrf(SMD,solvent=water) geom=allcheck
## defualt run platform is dp-bohrium, should config to login

def prepare_opt_gaussian(path, prefix):
    os.system(f"obabel -isdf {os.path.join(path, prefix + '.sdf')} -O {os.path.join(path, prefix + '.xyz')} -h")

    with open(f"{os.path.join(path, prefix+'.xyz')}", "r") as f_xyz:
        atom_xyz = [line.strip() for line in f_xyz.readlines() if len(line.strip()) > 0][2:]

    with open(f"{os.path.join(path, prefix+'.gjf')}", "w+") as c:
        c.write(f"%chk={prefix}.chk\n")
        c.write("%mem=20GB\n")
        c.write("%nproc=16\n")
    
        cmd = f"# B3LYP/6-311G* em=GD3BJ opt freq scrf(SMD,solvent=water) nosymm \n"
        c.write(cmd)
        c.write("\n")
        c.write(f"opt {prefix}\n")
        c.write("\n")
        c.write(f"0 1\n")

        for line in atom_xyz:
            c.write(line + "\n")

        c.write("\n")
        ## link section
        c.write("--link1-- \n")
        c.write(f"%chk={prefix}.chk\n")
        c.write("%mem=10GB\n")
        c.write("%nproc=16\n")
        c.write("# B2PLYPD3/def2TZVP scrf(SMD,solvent=water) nosymm geom=allcheck \n")
        c.write("\n")
        c.write("\n")
        

    return 

def run(input_sdf:str):
    work_dir = os.getcwd()
    genMol = [mm for mm in Chem.SDMolSupplier(os.path.join(work_dir, input_sdf)) if mm]
    for mm in genMol:
        cc = Chem.SDWriter(os.path.join(work_dir, f"{mm.GetProp('ID')}.sdf"))
        cc.write(mm)
        cc.close()
    
    g16_set = [nn.split('.')[0] for nn in os.listdir(work_dir) if 'sdf' in nn and 'input_out' not in nn]

    for prefix in g16_set:
        prepare_opt_gaussian(work_dir, prefix)

    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process g16 run')
    parser.add_argument('--input_sdf', type=str, required=True, help='input auto3D gen sdf')
    args = parser.parse_args()

    run(args.input_sdf)
