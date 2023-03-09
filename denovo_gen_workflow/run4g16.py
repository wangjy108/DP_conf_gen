import os
import sys
import pandas as pd
from rdkit import Chem

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
        c.write("%mem=50GB\n")
        c.write("%nproc=32\n")
    
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
        c.write("%mem=50GB\n")
        c.write("%nproc=32\n")
        c.write("# B2PLYPD3/def2TZVP scrf(SMD,solvent=water) nosymm geom=allcheck \n")
        c.write("\n")
        c.write("\n")
        

    return 


