#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

import argparse
import logging
from rdkit import Chem
import os
import sys
import pandas as pd
from util.ConfGenbyMM import ConfGen
from util.Cluster import cluster
from util.OptbySQM import System as sysopt
from util.SPcalc import System as syssp


def main(input_sdf:str):
    ## navigation
    main_dir = os.getcwd()
    _name = input_sdf.split(".")[0]
    work_dir = os.path.join(main_dir, f"{_name}")
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    os.chdir(work_dir)
    os.system(f"mv ../{input_sdf} ./")

    ## get_smi
    try:
        this_smi = Chem.MolToSmiles([mm for mm in Chem.SDMolSupplier(input_sdf) if mm][0])
    except Exception as e:
        return 
    
    with open("_input.smi", "w+") as write_smi:
        write_smi.write(f"{this_smi}\n")
    
    ## run sampling
    ConfGen(input_smi_file="_input.smi", method="MMFF94").run()

    ## run align 
    cluster(inputSDF_fileName="SAVE.sdf", save_n=20).run()

    ## run xtb opt
    collect_from_opt = sysopt(input_sdf="FILTER.sdf").run_process()
    cc = Chem.SDWriter("opted.sdf")
    for mm in collect_from_opt:
        cc.write(mm)
    cc.close()
    logging.info("Finish xtb opt")

    ## run final SP
    sorted_input_pose, input_mol = syssp(input_sdf=input_sdf).run_pyscf()
    logging.info("Get input pose single point energy")
    stable_pose, _ = syssp(input_sdf="opted.sdf", charge_method="read").run_pyscf()
    logging.info("Get theoretical stable pose single point energy")

    #input_energy = float(input_pose[0][0][0])
    stable_energy = float(stable_pose[0][0])

    cc = Chem.SDWriter(os.path.join(work_dir, f"{_name}_withEneTag.sdf"))
    for each in sorted_input_pose:
        get_energy = float(each[0])
        diff = abs(stable_energy - get_energy) * 627.51
        get_idx = int(each[1])
        this_mol = input_mol[get_idx]
        this_mol.SetProp("Energy_dft", str(diff))
        cc.write(this_mol)
    cc.close()

    os.system(f"obabel -isdf opted.sdf -O opt_{_name}.pdb")
    os.system(f"rm -f _input.smi SAVE.sdf FILTER.sdf opted.sdf")
    #logging.info(f"energy diff for {input_sdf} is : {diff}")
    os.chdir(main_dir)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sample naive strain energy workflow, \
                                    $WORKDIR has same name with input sdf, \
                                    save energy labeled sdf in [**_withEneTag.sdf], \
                                    and stable pose in [opt_***.sdf]')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    args = parser.parse_args()
    main(args.input_sdf)




    


