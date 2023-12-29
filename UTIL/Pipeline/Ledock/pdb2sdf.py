import os
import logging
import pandas as pd
import argparse
import configparser
import math
import json
import subprocess
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import numpy as np

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self, **args):
        parser = argparse.ArgumentParser(description="ledock pdb2sdf, only works with lbg submited type")

        parser.add_argument("--config", help="path of config file", default=False)
        args = parser.parse_args()

        config = configparser.ConfigParser()
        config.read(args.config)

        general = config["general"]

        self.main_dir = general["work_dir"]
        
        self.work_dir = os.path.join(self.main_dir, "DockingPose")
        self.input_db = os.path.join(self.main_dir, general["ligand_db_name"])

        #os.chdir(self.work_dir)

        #self.pdb_file_set = [cc for cc in os.listdir("./") if cc.endswith(".pdb")]
    
    def read_sdf(self):
        rdmol_obj_dic = {}
        #_name_idx_pair = {}
        #_label = self.input_db.split(".")[0]

        with open(self.input_db, "r+") as f:
            sdf_content = [ff for ff in f.readlines()]
        
        end_idx = [-1] + [idx for idx, line in enumerate(sdf_content) if "$$$$" in line][:-1]

        for ii, idx in enumerate(end_idx):
            try:
                get_mol = sdf_content[end_idx[ii]+1:end_idx[ii+1]+1]
            except Exception as e:
                get_mol = sdf_content[end_idx[ii]+1:]
            
            mol_block = []
            #get_mol_real_name = get_mol[0].strip()
            for j, line in enumerate(get_mol):
                mol_block.append(line)
            
            #try:
            #    rdmol_obj = Chem.MolFromMolBlock(mol_block, removeHs=False)
            #except Exception as e:
            #    rdmol_obj = None
            
            #try:
            #    get_real_name = rdmol_obj.GetProp("_Name")
            #except Exception as e:
            #    get_real_name = ""

            get_real_name = mol_block[0].strip()
            
            #if rdmol_obj and get_real_name:
            rdmol_obj_dic.setdefault(get_real_name, mol_block)
            #_name_idx_pair.setdefault(ii, get_real_name)
                
        return rdmol_obj_dic
    
    def pdb2xyz(self, pdb_atom_block):
        xyz = np.zeros((len(pdb_atom_block), 3))
        atom = []

        special_type = ["BR", "CL"]

        _func = lambda x: f"{x[0]}{x[1].lower()}" if x in special_type else f"{x[0]}"


        for idx, line in enumerate(pdb_atom_block):
            atom.append(_func(line[11:17].strip()))
            this_xyz = np.array([float(cc.strip()) for cc in line[26:].split() if cc.strip()])
            xyz[idx] = this_xyz
        
        df = pd.DataFrame({"atom": atom, \
                           "x": xyz[:, 0], \
                           "y": xyz[:, 1], \
                           "z": xyz[:, 2]})
        
        xyz_block = f"{xyz.shape[0]}\nDocked Pose\n"

        for idx, row in df.iterrows():
            xyz_block += f"{row['atom']:<3}{row['x']:>15.3f}{row['y']:>15.3f}{row['z']:>15.3f}\n"
        
        mol = rdmolfiles.MolFromXYZBlock(xyz_block)

        return mol
    
    def shift_sdf(self,
                  ori_content,
                  upt_content,
                  upt_moj):

        upt_content[3] = ori_content[3]
        
        #shape_xyz = int(upt_content[3].strip().split()[0])
        shape_xyz = upt_moj.GetConformer().GetPositions().shape[0]

        upper_save = upt_content[:3+shape_xyz+1]

        ori_end_idx = [cc for cc in range(len(ori_content)) if ori_content[cc].strip().endswith("END")][0]

        middile_save = ori_content[3+shape_xyz+1:ori_end_idx]

        up_end_part_start_idx = [cc for cc in range(len(upt_content)) if upt_content[cc].strip().endswith("END")][0]
        end_save = upt_content[up_end_part_start_idx:]

        assemble_content = upper_save + middile_save + end_save
        
        return assemble_content
    
    def read_input_ligand(self, pdb):

        _dic = {}
        #_dic_score = {}

        ligand_name = "_".join(pdb.split("/")[0].split(".")[0].split("_")[1:])

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
            docking_score = float(header.split(":")[-1].strip())

            pdb_atom_block = [cc for cc in get_each if cc.startswith("ATOM") and len(cc.strip().split()) > 6]

            get_mol = self.pdb2xyz(pdb_atom_block)
            
            try:
                get_mol.SetProp("DockingScore", f"{docking_score}")
                get_mol.SetProp("_Name", f"{ligand_name}:dock_{ii}")
            except Exception as e:
                logging.info(f"{ligand_name} failed")
                return None

            _dic.setdefault(f"{ligand_name}:{ii}", get_mol)
            #_dic_score.setdefult(f"{ligand_name}:{ii}", docking_score)
                
        return _dic

    def run(self):
        os.chdir(self.work_dir)
        issued = []
        assembled = []
        pdb_file_list = [cc for cc in os.listdir("./") if cc.endswith(".pdb")]

        if not pdb_file_list:
            logging.info("Nothing to process, terminate")
            return 
        
        input_sdf_rdmol_molcontent = self.read_sdf()

        for idx, pdb_file in enumerate(pdb_file_list):
            try:
                get_mol = self.read_input_ligand(pdb_file)
            except Exception as e:
                issued.append(pdb_file)
                continue
            
            if not get_mol:
                issued.append(pdb_file)
                continue

            get_name = [kk for kk in get_mol.keys()][0].split(":")[0]

            try:
                ori_content = input_sdf_rdmol_molcontent[get_name]
            except Exception as e:
                issued.append(pdb_file)
                continue
            else:
                for kk, vv in get_mol.items():
                    upt_obj = vv

                    cc = Chem.SDWriter(f"TEMP_{kk}.sdf")
                    cc.write(upt_obj)
                    cc.close()

                    with open(f"TEMP_{kk}.sdf", "r+") as f:
                        upt_content = [cc for cc in f.readlines()]
                    
                    os.system(f"rm -f TEMP_{kk}.sdf")

                    new_mol_content = self.shift_sdf(ori_content, upt_content, upt_obj)

                    assembled += new_mol_content
        
        with open(os.path.join(self.main_dir, "DockingPose.sdf"), "w+") as cc:
            for line in assembled:
                cc.write(line)
        if issued:
            df = pd.DataFrame({"LigandName": issued})
            df.to_csv(f"{os.path.join(self.main_dir,'Issued.LIST')}", index=None)
            logging.info(f"Save issued ligand name in Issued.LIST in {self.main_dir}")

        
        logging.info(f"Save DockingPose.sdf in {self.main_dir}")
        logging.info(f"Done")

        os.chdir(self.main_dir)

        return 

if __name__ == "__main__":
    main()





            


                
                

            
                    

        




        
        
