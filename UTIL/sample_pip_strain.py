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

def shift_sdf(self, original_sdf, update_sdf, save_prefix):
    with open(original_sdf, "r+") as f1:
        ori_content = [ff for ff in f1.readlines()]

    rdmolobj_ori = [mm for mm in Chem.SDMolSupplier(original_sdf, removeHs=False) if mm][0]
    
    with open(update_sdf, "r+") as f2:
        upt_content = [ff for ff in f2.readlines()]
    
    rdmolobj_upt = [mm for mm in Chem.SDMolSupplier(update_sdf, removeHs=False) if mm][0]

    ## header should stop at 
    xyz_upt = rdmolobj_upt.GetConformer().GetPositions()
    xyz_ori = rdmolobj_ori.GetConformer().GetPositions()
    ## xyz_upt[-1][0]
    
    upper_end_idx = [ii for ii in range(len(upt_content)) if "END" in upt_content[ii]][0]
    
    middle_start_idx = [ii for ii in range(len(ori_content)) if ori_content[ii].strip().startswith(str(xyz_ori[-1][0]))][0] + 1
    middle_end_idx = [ii for ii in range(len(ori_content)) if "END" in ori_content[ii]][0]

    header_replace_idx_upt = [ii for ii in range(len(upt_content)) \
                            if upt_content[ii].strip().startswith(str(xyz_upt[0][0]))][0] -1
    header_replace_idx_ori = [ii for ii in range(len(ori_content)) \
                            if ori_content[ii].strip().startswith(str(xyz_ori[0][0]))][0] -1
    
    upper = upt_content[:upper_end_idx]
    upper[header_replace_idx_upt] = ori_content[header_replace_idx_ori]

    assemble_content = upper \
                    + ori_content[middle_start_idx:middle_end_idx] \
                    + upt_content[upper_end_idx:]
    
    with open(f"{save_prefix}.sdf", "w+") as cc:
        for line in assemble_content:
            cc.write(line)
    
    return 

def main(input_sdf:str, if_csv:bool):
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
    logging.info("Start geom optimization")
    _ = sysopt(input_sdf="FILTER.sdf").run_process()

    ## run final SP
    logging.info("Start Single point energy calc")

    sorted_input_pose, input_mol = syssp(input_sdf=input_sdf).run_pyscf()
    logging.info("Get input pose single point energy")
    stable_pose, _ = syssp(input_sdf="_OPT.sdf", charge_method="read").run_pyscf()
    logging.info("Get theoretical stable pose single point energy")

    #input_energy = float(input_pose[0][0][0])
    stable_energy = float(stable_pose[0][0])

    track_energy = []
    track_mol_label = []

    cc = Chem.SDWriter(os.path.join(work_dir, f"{_name}_withEneTag.sdf"))
    for each in sorted_input_pose:
        get_energy = float(each[0])
        diff = abs(stable_energy - get_energy) * 627.51
        get_idx = int(each[1])
        this_mol = input_mol[get_idx]
        this_mol.SetProp("Energy_dft", f"{diff:.3f}")

        track_energy.append(f"{diff:.3f}")
        track_mol_label.append(f"{_name}_{get_idx}")

        cc.write(this_mol)
    cc.close()

    if if_csv:
        df = pd.DataFrame({"mol_label":track_mol_label, \
                           "Strain_ene(kcal/mol)": track_energy})
        df.to_csv(f"StrainEne_{_name}.csv", index=None)
        logging.info(f"Save strain energy in {os.getcwd()}/StrainEne_{_name}.csv")

    os.system(f"mv _OPT.sdf stable_{_name}.sdf")
    os.system(f"rm -f _input.smi SAVE.sdf FILTER.sdf")
    logging.info(f"Strain energy for input is labeled in {_name}_withEneTag.sdf \
                  Tag name is [Energy_dft], unit is [kcal/mol]")
    
    os.chdir(main_dir)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sample naive strain energy workflow, \
                                    $WORKDIR has same name with input sdf, \
                                    save energy labeled sdf in [**_withEneTag.sdf], \
                                    and stable pose in [opt_***.sdf]')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    parser.add_argument('--if_csv', type=bool, default=True, 
                        help="if save csv to record calc result")
    args = parser.parse_args()
    main(args.input_sdf, args.if_csv)




    


