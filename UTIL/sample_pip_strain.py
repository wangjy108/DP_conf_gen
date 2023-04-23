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
from util.Align import Align as align
from util.ConfRelaxbySQM import System as MDsample


class main():
    def __init__(self, **args):
        try:
            self.db_name = args["input_sdf"]
        except Exception as e:
            self.db_name = None

        try:
            self.method = args["method"]
        except Exception as e:
            self.method = 'denovo'
        
        try:
            self.if_csv = args["if_csv"]
        except Exception as e:
            self.if_csv = True
        else:
            if not isinstance(self.if_csv, bool):
                self.if_csv = True
        
        try:
            self.N_gen_conformer = args["N_gen_conformer"]
        except Exception as e:
            self.N_gen_conformer = 50
        else:
            try:
                self.N_gen_conformer = int(self.N_gen_conformer)
            except Exception as e:
                self.N_gen_conformer = 50

    def pip_denovo(self):
        ## get_smi
        try:
            this_smi = Chem.MolToSmiles([mm for mm in Chem.SDMolSupplier(self.db_name) if mm][0])
        except Exception as e:
            logging.info("Bad input sdf, check and run again")
            return 
        logging.info("Run strain calc in de-novo mode")

        _name = self.db_name.split(".")[0]

        if not self.method == "sample":
            with open("_input.smi", "w+") as write_smi:
                write_smi.write(f"{this_smi}\n")
            ## run sampling
            try:
                ConfGen(input_smi_file="_input.smi", method="MMFF94").run()
            except Exception as e:
                logging.info("Initial MM conf gen failed, use MD sampling instead")
                MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer).run()
            
            if (not os.path.isfile("SAVE.sdf")) or (not os.path.getsize("SAVE.sdf")):
                logging.info("Initial MM conf gen failed, use MD sampling instead")
                MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer).run()
        else:
            logging.info("use MD relax for sampling")
            MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer).run()

        ## run align 
        cluster(inputSDF_fileName="SAVE.sdf", save_n=self.N_gen_conformer).run()

        ## run xtb opt
        logging.info("Start geom optimization")
        _ = sysopt(input_sdf="FILTER.sdf").run_process()

        ## run final SP
        logging.info("Start Single point energy calc")

        sorted_input_pose, input_mol = syssp(input_sdf=self.db_name).run_pyscf()
        logging.info("Get input pose single point energy")
        stable_pose, stable_mol = syssp(input_sdf="_OPT.sdf", charge_method="read").run_pyscf()
        logging.info("Get theoretical stable pose single point energy")

        #input_energy = float(input_pose[0][0][0])
        stable_energy = float(stable_pose[0][0])

        track_energy = []
        track_mol_label = []

        cc = Chem.SDWriter(os.path.join(os.getcwd(), f"{_name}_withEneTag.sdf"))
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

        if self.if_csv:
            df = pd.DataFrame({"mol_label":track_mol_label, \
                            "Strain_ene(kcal/mol)": track_energy})
            df.to_csv(f"StrainEne_{_name}.csv", index=None)
            logging.info(f"Save strain energy in {os.getcwd()}/StrainEne_{_name}.csv")

        #os.system(f"mv _OPT.sdf stable_{_name}.sdf")

        ## should be update
        _stable_mol = [mm for mm in Chem.SDMolSupplier("_OPT.sdf", removeHs=False) if mm][0]
        _ref_mol = [mm for mm in Chem.SDMolSupplier(self.db_name, removeHs=False) if mm][0]
        get_aligned_opt_pose = align(SearchMolObj=_stable_mol, RefMolObj=_ref_mol, method="crippen3D").run()
        cc_opt = Chem.SDWriter(f"stable_{_name}.sdf")
        cc_opt.write(get_aligned_opt_pose)
        cc.close()

        os.system("rm -f _input.smi SAVE.sdf FILTER.sdf _OPT.sdf")
        logging.info(f"Strain energy for input is labeled in {_name}_withEneTag.sdf \
                    Tag name is [Energy_dft], unit is [kcal/mol]")

        return
    
    def pip_local(self):
        ## get mol to do opt
        logging.info("Run strain calc in local mode")

        logging.info("Start geom optimization")
        _ = sysopt(input_sdf=self.db_name).run_process()

        _name = self.db_name.split(".")[0]

        ## run final SP
        logging.info("Start Single point energy calc")

        sorted_input_pose, input_mol = syssp(input_sdf=self.db_name).run_pyscf()
        logging.info("Get input pose single point energy")
        stable_pose, stable_mol = syssp(input_sdf="_OPT.sdf", charge_method="read").run_pyscf()
        logging.info("Get theoretical stable pose single point energy")

        #input_energy = float(input_pose[0][0][0])
        stable_energy = float(stable_pose[0][0])

        track_energy = []
        track_mol_label = []

        cc = Chem.SDWriter(os.path.join(os.getcwd(), f"{_name}_withEneTag.sdf"))
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

        if self.if_csv:
            df = pd.DataFrame({"mol_label":track_mol_label, \
                            "Strain_ene(kcal/mol)": track_energy})
            df.to_csv(f"StrainEne_{_name}.csv", index=None)
            logging.info(f"Save strain energy in {os.getcwd()}/StrainEne_{_name}.csv")

        os.system(f"mv _OPT.sdf stable_{_name}.sdf")

        os.system("rm -f *.smi SAVE.sdf FILTER.sdf _*.sdf")
        logging.info(f"Strain energy for input is labeled in {_name}_withEneTag.sdf \
                    Tag name is [Energy_dft], unit is [kcal/mol]")

        return
    
    def run(self):
        try:
            _name = self.db_name.split(".")[0]
        except Exception as e:
            logging.info("Bad input sdf, check and run again")
            return
        
        main_dir = os.getcwd()
        work_dir = os.path.join(main_dir, f"{_name}")
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        os.chdir(work_dir)
        os.system(f"mv ../{self.db_name} ./")

        run_dict = {"denovo": self.pip_denovo, \
                    "local": self.pip_local, \
                    "sample": self.pip_denovo}
        try:
            run_dict[self.method]
        except Exception as e:
            logging.info(f"Wrong method setting, please choose from [{[kk for kk in run_dict.keys()]}]")
            return None
        
        run_dict[self.method]()
        os.chdir(main_dir)

        return 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sample naive strain energy workflow, \
                                    $WORKDIR has same name with input sdf, \
                                    save energy labeled sdf in [**_withEneTag.sdf], \
                                    and stable pose in [opt_***.sdf]')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    parser.add_argument('--method', type=str, default='denovo',
                        help='strain method, available from ["denovo","local","sample"], \
                        default [denovo]')
    parser.add_argument('--N_gen_conformer', type=int, default=50, 
                        help='available for [denovo] mode, define N conformers in sampling stage, \
                        default 50')
    parser.add_argument('--if_csv', type=bool, default=True, 
                        help="if save csv to record calc result")
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, \
         method=args.method, \
         N_gen_conformer=args.N_gen_conformer, \
         if_csv=args.if_csv).run()




    
##

