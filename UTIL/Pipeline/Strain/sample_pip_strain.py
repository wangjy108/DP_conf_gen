#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

import argparse
import logging
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
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
            db_name = args["input_sdf"]
        except Exception as e:
            db_name = None
        
        self.main_dir = os.getcwd()
        
        ## initial
        _dic_rotableBond = {0: 100, 1: 300}
        try:
            self.get_mol = [mm for mm in Chem.SDMolSupplier(db_name) if mm][0]
        except Exception as e:
            self.get_mol = None
            rotable_bond_index = 0
        else:
            rotable_bond = rdMolDescriptors.CalcNumRotatableBonds(self.get_mol)
            _def_func = lambda x: 0 if max(5, x) == 5 else 1
            rotable_bond_index = _def_func(rotable_bond)

        try:
            self.method = args["method"]
        except Exception as e:
            self.method = 'sample'
        
        try:
            self.N_gen_conformer = args["N_gen_conformer"]
        except Exception as e:
            self.N_gen_conformer = _dic_rotableBond[rotable_bond_index]    
        else:
            try:
                self.N_gen_conformer = int(self.N_gen_conformer)
            except Exception as e:
                self.N_gen_conformer = _dic_rotableBond[rotable_bond_index]
        
        try:
            self.charge = args["define_charge"]
        except Exception as e:
            self.charge = None
        else:
            try:
                self.charge = int(self.charge)
            except Exception as e:
                self.charge = None
        
        #print(f"gen conformer: {self.N_gen_conformer}")
        
        self.HA_constrain = args["use_constrain"]
        self.if_csv = args["if_csv"]

        ## perform initial_opt by SQM
        if db_name:
            self.prefix = ".".join(db_name.split(".")[:-1])
            work_dir = os.path.join(self.main_dir, f"{self.prefix}")
            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
            os.chdir(work_dir)
            os.system(f"mv {self.main_dir}/{db_name} ./")

            if self.charge:
                _, _ = sysopt(input_sdf=db_name, HA_constrain=True, charge_method="define", define_charge=self.charge).run_process()
            else:
                _, read_charge = sysopt(input_sdf=db_name, HA_constrain=True).run_process()
                self.charge = read_charge[0]

            ## rename _OPT.sdf for other run 
            if os.path.isfile("_OPT.sdf"):
                #_prefix = db_name.split(".")[0]
                self.db_name = f"{self.prefix}.initial_opt.sdf"
                os.system(f"mv _OPT.sdf {self.db_name}")
            else:
                self.db_name = None
                logging.info("Failed at initial optimization for input")
        else:
            self.db_name = None
            logging.info("Bad input sdf file")
        

    def pip_denovo(self):
        ## get_smi
        try:
            this_smi = Chem.MolToSmiles(self.get_mol)
        except Exception as e:
            logging.info("check and run again")
            return 
        logging.info("Run strain calc in de-novo gen mode")
        
        if not self.method == "sample":
            with open("_input.smi", "w+") as write_smi:
                write_smi.write(f"{this_smi}\n")
            ## run sampling
            try:
                ConfGen(input_smi_file="_input.smi", method="MMFF94").run()
            except Exception as e:
                logging.info("Initial MM conf gen failed, use MD sampling instead")
                MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer, define_charge=self.charge).run()
            
            if (not os.path.isfile("SAVE.sdf")) or (not os.path.getsize("SAVE.sdf")):
                logging.info("Initial MM conf gen failed, use MD sampling instead")
                MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer, define_charge=self.charge).run()
        else:
            logging.info("use MD relax for sampling")
            MDsample(input_sdf=self.db_name, save_frame=self.N_gen_conformer, define_charge=self.charge).run()

        ## run align 
        cluster(inputSDF_fileName="SAVE.sdf", save_n=self.N_gen_conformer).run()

        ## run xtb opt
        logging.info("Start geom optimization")
        _, _ = sysopt(input_sdf="FILTER.sdf", HA_constrain=self.HA_constrain, define_charge=self.charge).run_process()

        ## run final SP
        logging.info("Start Single point energy calc")

        sorted_input_pose, input_mol = syssp(input_sdf=self.db_name, charge_method="define", define_charge=self.charge).run_pyscf()
        logging.info("Get input pose single point energy")
        stable_pose, stable_mol = syssp(input_sdf="_OPT.sdf", charge_method="define", define_charge=self.charge).run_pyscf()
        logging.info("Get theoretical stable pose single point energy")

        #input_energy = float(input_pose[0][0][0])
        stable_energy = float(stable_pose[0][0])

        track_energy = []
        track_mol_label = []

        cc = Chem.SDWriter(os.path.join(os.getcwd(), f"{self.prefix}_withEneTag.sdf"))
        for each in sorted_input_pose:
            get_energy = float(each[0])
            diff = abs(stable_energy - get_energy) * 627.51
            get_idx = int(each[1])
            this_mol = input_mol[get_idx]
            this_mol.SetProp("Energy_dft", f"{diff:.3f}")

            track_energy.append(f"{diff:.3f}")
            track_mol_label.append(f"{self.prefix}_{get_idx}")

            cc.write(this_mol)
        cc.close()

        if self.if_csv:
            df = pd.DataFrame({"mol_label":track_mol_label, \
                            "Strain_ene(kcal/mol)": track_energy})
            df.to_csv(f"StrainEne_{self.prefix}.csv", index=None)
            logging.info(f"Save strain energy in {os.getcwd()}/StrainEne_{self.prefix}.csv")

        #os.system(f"mv _OPT.sdf stable_{_name}.sdf")

        ## should be update
        _stable_mol = [mm for mm in Chem.SDMolSupplier("_OPT.sdf", removeHs=False) if mm][0]
        _ref_mol = [mm for mm in Chem.SDMolSupplier(self.db_name, removeHs=False) if mm][0]
        get_aligned_opt_pose = align(SearchMolObj=_stable_mol, RefMolObj=_ref_mol, method="crippen3D").run()
        cc_opt = Chem.SDWriter(f"stable_{self.prefix}.sdf")
        cc_opt.write(get_aligned_opt_pose)
        cc.close()

        os.system("rm -f _input.smi SAVE.sdf FILTER.sdf _OPT.sdf")
        logging.info(f"Strain energy for input is labeled in {self.prefix}_withEneTag.sdf \
                    Tag name is [Energy_dft], unit is [kcal/mol]")

        return
    
    def pip_local(self):
        ## get mol to do opt
        logging.info("Run strain calc in local mode")

        logging.info("Start geom optimization")
        _, _ = sysopt(input_sdf=self.db_name, HA_constrain=False, define_charge=self.charge).run_process()

        ## run final SP
        logging.info("Start Single point energy calc")

        sorted_input_pose, input_mol = syssp(input_sdf=self.db_name, charge_method="define", define_charge=self.charge).run_pyscf()
        logging.info("Get input pose single point energy")
        stable_pose, stable_mol = syssp(input_sdf="_OPT.sdf", charge_method="define", define_charge=self.charge).run_pyscf()
        logging.info("Get theoretical stable pose single point energy")

        #input_energy = float(input_pose[0][0][0])
        stable_energy = float(stable_pose[0][0])

        track_energy = []
        track_mol_label = []

        cc = Chem.SDWriter(os.path.join(os.getcwd(), f"{self.prefix}_withEneTag.sdf"))
        for each in sorted_input_pose:
            get_energy = float(each[0])
            diff = abs(stable_energy - get_energy) * 627.51
            get_idx = int(each[1])
            this_mol = input_mol[get_idx]
            this_mol.SetProp("Energy_dft", f"{diff:.3f}")

            track_energy.append(f"{diff:.3f}")
            track_mol_label.append(f"{self.prefix}_{get_idx}")

            cc.write(this_mol)
        cc.close()

        if self.if_csv:
            df = pd.DataFrame({"mol_label":track_mol_label, \
                            "Strain_ene(kcal/mol)": track_energy})
            df.to_csv(f"StrainEne_{self.prefix}.csv", index=None)
            logging.info(f"Save strain energy in {os.getcwd()}/StrainEne_{self.prefix}.csv")

        os.system(f"mv _OPT.sdf stable_{self.prefix}.sdf")

        os.system("rm -f *.smi SAVE.sdf FILTER.sdf _*.sdf")
        logging.info(f"Strain energy for input is labeled in {self.prefix}_withEneTag.sdf \
                    Tag name is [Energy_dft], unit is [kcal/mol]")

        return
    
    def run(self):
        if not self.db_name:
            logging.info("check and run again")
            return 
        
        run_dict = {"denovo": self.pip_denovo, \
                    "local": self.pip_local, \
                    "sample": self.pip_denovo}
        try:
            run_dict[self.method]
        except Exception as e:
            logging.info(f"Wrong method setting, please choose from [{[kk for kk in run_dict.keys()]}]")
            return None
        
        run_dict[self.method]()
        os.chdir(self.main_dir)
        
        return 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sample naive strain energy workflow, \
                                    $WORKDIR has same name with input sdf, \
                                    save energy labeled sdf in [**_withEneTag.sdf], \
                                    and stable pose in [opt_***.sdf]')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    parser.add_argument('--method', type=str, default="sample",
                        help='strain method, available from ["denovo","local","sample"], \
                        default [sample]')
    parser.add_argument('--N_gen_conformer', type=int, 
                        help='available for [denovo and sample] mode, define N conformers in sampling stage, \
                        if not defined, 100 if rotable bond <=5, else 300')
    parser.add_argument('--define_charge', type=str, default=None,
                        help='define charge for input, use you sense that rdkit may fail to get corret bond order, default None')
    parser.add_argument('--not_save_csv', default=True, action='store_false', \
                        help="adding this option will not save final result in csv file")
    parser.add_argument('--no_constrain', default=True, action='store_false', \
                        help="adding this option will turnoff constrain on heavy atoms when performing geometry optimization")
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, \
         method=args.method, \
         N_gen_conformer=args.N_gen_conformer, \
         define_charge=args.define_charge, \
         if_csv=args.not_save_csv,\
         use_constrain=args.no_constrain).run()




    
##

