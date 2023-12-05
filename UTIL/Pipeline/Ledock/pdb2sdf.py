import os
import logging
import pandas as pd
import argparse
import configparser
import math
import json
import subprocess
from rdkit import Chem

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self, **args):
        self.work_dir = args["abs_path"]
        
        os.chdir(self.work_dir)

        self.pdb_file_set = [cc for cc in os.listdir("./") if cc.endswith(".pdb")]
    
    def col_and_seperate(self, pdb):
        _dic = {}

        with open(pdb, 'r+') as f:
            content = [ff for ff in f.readlines()]
        
        end_idx = [-1] + [idx for idx, line in enumerate(content) if "END" in line][:-1]
        
        for ii, idx in enumerate(end_idx):
            try:
                get_each = content[end_idx[ii]+1:end_idx[ii+1]+1]
            except Exception as e:
                get_each = content[end_idx[ii]+1:]
            
            header_idx = [idx for idx in range(len(get_each)) if \
                       get_each[idx].startswith("ATOM") and \
                       len(get_each[idx].split()) > 6][0] - 1
            header = get_each[header_idx]
            
            _dic.setdefault(header, self.reframe(get_each))
                
        return _dic
    
    def reframe(self, original_block):
        update = []
        for ii, line in enumerate(original_block):
            if len(line.split()) > 6 and line.startswith("ATOM"):
                if len(line[11:17].strip()) > 1:
                    new_line = line[:12] + line[11:17].strip()[0] + line[11:17].strip()[1].lower() + " " + line[15:]
                else:
                    new_line = line
            else:
                new_line = line
            
            update.append(new_line)

        return update
    
    def pdb2sdf(self, kk_header, vv_pdb_block):
    
        name = kk_header.split()[0]
        docking_score = kk_header.strip().split(":")[-1].strip()

        with open(f"_TEMP_{name}.pdb", "w+") as cc:
            for ii, line in enumerate(vv_pdb_block):
                if line.startswith("ATOM") or line.startswith("END"):
                    cc.write(line)
        
        cmd = f"obabel -ipdb _TEMP_{name}.pdb -O _TEMP_{name}.sdf"
        
        try:
            p = subprocess.run(cmd.split(), timeout=20, check=True, stdout=subprocess.PIPE)
        except subprocess.TimeoutExpired:
            logging.info("Timeout with obabel format trans")
            mol = None
        else:
            try:
                mol = [cc for cc in Chem.SDMolSupplier(f"_TEMP_{name}.sdf", removeHs=False) if cc]
            except Exception as e:
                return None 
            if not mol:
                return None
            else:
                mol = mol[0]
                mol.SetProp("_Name", name)
                mol.SetProp("DockingScore", docking_score)

                os.system(f"rm -f _TEMP_{name}.sdf")
                os.system(f"rm -f _TEMP_{name}.pdb")  
        
        return mol
    
    def run(self):
        if not self.pdb_file_set:
            logging.info("No available pdb file to process, terminate")
            return 
        
        _dic_assemble = {}
        for ii, pdb in enumerate(self.pdb_file_set):
            get_docked_pose = self.col_and_seperate(pdb)
            _dic_assemble.update(get_docked_pose)
        
        cc = Chem.SDWriter(f"Assembled.sdf")
        
        for kk, vv in _dic_assemble.items():
            try:
                get_mol = self.pdb2sdf(kk, vv)
            except Exception as e:
                get_mol = None
            
            if not get_mol:
                continue
            else:
                cc.write(get_mol)
        
        cc.close()

        logging.info("Save assembled sdf in Assemble.sdf")
        logging.info("left out *.pdb (if any), are not included")

        return 
    
            


            


                
                

            
                    

        




        
        
