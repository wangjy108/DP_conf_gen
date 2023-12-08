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
import numpy as np
import shutil
from util.ConfRelaxbySQM import System as MDsample
from util.Cluster import cluster


class main():
    def __init__(self, **args):
        try:
            self.db_name = args["input_sdf"]
        except Exception as e:
            self.db_name = None
        
        self.main_dir = os.getcwd()
        
        ## initial
        _dic_rotableBond = {0: 100, 1: 300}
        #_dic_HA_index = {}
        try:
            self.get_mol = [mm for mm in Chem.SDMolSupplier(self.db_name) if mm][0]
        except Exception as e:
            self.get_mol = None
            rotable_bond_index = 0
        else:
            rotable_bond = rdMolDescriptors.CalcNumRotatableBonds(self.get_mol)
            _def_func = lambda x: 0 if max(5, x) == 5 else 1
            rotable_bond_index = _def_func(rotable_bond)
        
        
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
        
        ##parameter for sample & cluster
        try:
            self.k = args["cluster_n"]
        except Exception as e:
            self.k = None
        else:
            if not isinstance(self.k, int):
                self.k = None

        try:
            self.distance_cutoff = args["rmsd_cutoff_cluster"]
        except Exception as e:
            self.distance_cutoff = 0.85
        else:
            if not isinstance(self.distance_cutoff, float):
                self.distance_cutoff = 0.85
    
    
    def step2_sampling(self, input_file_name):
        sampled_file_name = None

        if not self.charge:
            MDsample(input_sdf=input_file_name, 
                    save_frame=self.N_gen_conformer).run()
        else:
            MDsample(input_sdf=input_file_name, 
                    save_frame=self.N_gen_conformer, 
                    define_charge=self.charge).run()

        if os.path.isfile("SAVE.sdf") and os.path.getsize("SAVE.sdf"):
            logging.info("Samping success")
            sampled_file_name = "SAVE.sdf"
        
        else:
            logging.info("Failed at sampling")

        return sampled_file_name
    
    def step3_cluster(self, input_sdf_name):

        filtered_mol = cluster(input_sdf=input_sdf_name,
                cluster_n=self.k,
                rmsd_cutoff_cluster=self.distance_cutoff,
                if_write_cluster_mol=False).run()
        
        #filtered_file_name = [cc for cc in os.listdir() if ("FILTER_" in cc) and (cc.endswith(".sdf"))]
        if not filtered_mol:
            logging.info("Failed at cluster")
            return None
        
        cc = Chem.SDWriter("FILTER.sdf")
        for each in filtered_mol:
            cc.write(each)
        cc.close()

        return "FILTER.sdf"
    
    def run(self):
        if not self.db_name:
            logging.info("check and run again")
            return
        
        self.prefix = ".".join(self.db_name.split(".")[:-1])
        work_dir = os.path.join(self.main_dir, f"{self.prefix}")
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)

        os.chdir(work_dir)
        os.system(f"mv {self.main_dir}/{self.db_name} ./")
        
        logging.info("Step1: Geometrical Sampling ....")
        self.sampled_file = self.step2_sampling(self.db_name)
        if not self.sampled_file:
            logging.info("Terminated at Step1.")
            return 
        
        logging.info("Step2: Clustering ...")
        self.filter_name = self.step3_cluster(self.sampled_file)
        if not self.filter_name:
            logging.info("Terminated at Step2.")
            return
        
        logging.info("Done")
        
        #os.system(f"rm -f SAVE.sdf")

        os.chdir(self.main_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='conf sample')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file, docking pose(s) for single mol')
    parser.add_argument('--define_charge', type=str, default=None,
                        help='define charge for input, use you sense that rdkit may fail to get corret bond order, default None')
    parser.add_argument('--cluster_n', type=int, 
                        help='set static n for clustering after sampling')
    parser.add_argument('--rmsd_cutoff', type=float,
                        help='rmsd cutoff for reduce redundancy, default 0.85')
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, \
         define_charge=args.define_charge, \
         clsuter_n=args.cluster_n, \
         rmsd_cutoff_cluster=args.rmsd_cutoff).run()

    
##

