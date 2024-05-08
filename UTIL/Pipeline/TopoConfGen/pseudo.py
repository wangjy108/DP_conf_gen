import os
import logging
import subprocess
import pandas as pd
from rdkit import Chem

os.environ["PYTHONPATH"] = "/opt/scripts"

from gen import gen_former

class confGen_topo():
    def __init__(self,**args):
        self.input_smi = args["input_csv"]
        ## consist with uniqsar, comma seperated csv, SMILES for smiles, NAME for ligand label
        self.work_dir = os.getcwd()
        self.gen_EnergyWindow = args["energy_window"]
        self.gen_RmsdCutoff = args["rmsd_cutoff"]
        self.dock_SaveN = args["n_pose_per_conformer"]
        
    def read_in_process(self):
        df = pd.read_csv(self.input_smi, header=0)
        if not "SMILES" in df.columns:
            logging.info("Input smiles should be under [SMILES] header, abort")
            return None
        
        if not "NAME" in df.columns:
            saved_colunms = df[["SMILES"]]
        else:
            saved_colunms = df[["SMILES", "NAME"]]
        
        if not saved_colunms.shape[0]:
            logging.info("Nothing to do, abort")
            return None
        saver = os.path.join(self.work_dir, "_TEMP_input.smi")
        saved_colunms.to_csv(saver, header=None, index=None, sep="\t")

        return "_TEMP_input.smi"

    def NLdock(self, ligand_mol2):
        command = f"NLDock /root/topoI_NLdock/topo_receptor.mol2 {ligand_mol2} -site /root/topoI_NLdock/bds.pdb -write_max {self.dock_SaveN}"
        (_, _) = subprocess.getstatusoutput(command)

        saver = os.path.join(self.work_dir, "nldock.mol2")

        return saver

    def reshape(self, docked_mol2):
        _dic = {}

        with open(docked_mol2, "r") as f:
            content = [cc for cc in f.readlines()]
        
        indicator = [idx for idx in range(len(content)) if content[idx].startswith("########## Ligand Number")]

        i = 0
        while i < len(indicator):
            try:
                _raw_block = content[indicator[i]:indicator[i+1]]
            except Exception as e:
                _raw_block = content[indicator[i]:]
            
            mol_block_s_idx = [idx for idx, line in enumerate(_raw_block) if line.startswith("@<TRIPOS>MOLECULE")]
            #mol_block_e_idx = [idx for idx, line in enumerate(_raw_block) if line.startswith("@<TRIPOS>SUBSTRUCTURE")]

            mol_block = _raw_block[mol_block_s_idx[0]:]

            get_ligand_name = [line for line in _raw_block if line.startswith("########## Ligand Name")][0].split(":")[-1].strip()
            get_ligand_idx = [line for line in _raw_block if line.startswith("########## Orient Number")][0].split(":")[-1].strip()
            get_docking_score = [line for line in _raw_block if line.startswith("########## ITScore")][0].split(":")[-1].strip()

            with open(f"{get_ligand_name}_{get_ligand_idx}.TEMP.mol2", "w+") as cc:
                for each in mol_block:
                    cc.write(each)
            
            cmd = f"obabel -imol2 {get_ligand_name}_{get_ligand_idx}.TEMP.mol2 -O {get_ligand_name}_{get_ligand_idx}.TEMP.sdf"

            _, _ = subprocess.getstatusoutput(cmd)

            if os.path.isfile(f"{get_ligand_name}_{get_ligand_idx}.TEMP.sdf"):
                try:
                    mol = [cc for cc in Chem.SDMolSupplier(f"{get_ligand_name}_{get_ligand_idx}.TEMP.sdf", removeHs=False) if cc]
                except Exception as e:
                    mol = None

            if mol:
                mol[0].SetProp("NLDockingScore", get_docking_score)
                os.system(f"rm -f *.TEMP.*")
            
            try:
                _dic[get_ligand_name]
            except Exception as e:
                _dic.setdefault(get_ligand_name, [])
            
            _dic[get_ligand_name].append(mol[0])

            i += 1
        
        saver = os.path.join(self.work_dir, "output.sdf")
        writer = Chem.SDWriter(saver)

        for kk, vv in _dic.items():
            for idx, mm in enumerate(vv):
                Name = f"{kk}:dock_{idx}"
                mm.SetProp("_Name", Name)
                writer.write(mm)
        
        writer.close()

        with open(saver, "r+") as ff:
            saver_content = [ll for ll in ff.readlines()]

        return saver_content
    
    def run(self):
        logging.info("Pre-process: readin ----->")
        c_input_smi = self.read_in_process()
        if not c_input_smi:
            return 
        logging.info("Step 1: LigPrep ----->")
        gen_former(input_file=c_input_smi,
                    energy_window=self.gen_EnergyWindow,
                    rmsd_cutoff=self.gen_RmsdCutoff)
        logging.info("Step 2: Constrain sampling from docking ----->")
        try:
            ligand = [os.path.join(self.work_dir, cc) for cc in os.listdir(self.work_dir) if cc.startswith("Gen") and cc.endswith(".sdf")]
        except Exception as e:
            return 
        else:
            ligand_sdf = ligand[0]
            
        ligand_mol2 = os.path.join(self.work_dir, "_TEMP_GEN.mol2")
        cmd = f"obabel -isdf {ligand_sdf} -O {ligand_mol2}"
        _, _ = subprocess.getstatusoutput(cmd)
        if not os.path.isfile(ligand_mol2):
            return 
        docked_mol2 = self.NLdock(ligand_mol2)
        if not os.path.isfile(docked_mol2):
            logging.info("Failed at docking, abort")
            return

        saved_sdf = self.reshape(docked_mol2)
        if saved_sdf:
            logging.info("Done")
            if os.path.isfile("ERROR.smi"):
                with open("ERROR.smi", "r+") as ee:
                    saved_error = [cc for cc in ee.readlines()]
            else:
                saved_error = None
                

            rm_cmd = "rm -f *.sdf *TEMP* nldock* *.smi"
            _, _ = subprocess.getstatusoutput(rm_cmd)

        return saved_sdf, saved_error
