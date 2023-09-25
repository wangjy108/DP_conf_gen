import os
import logging
import pandas as pd
import argparse
import configparser
import subprocess
from multiprocessing import cpu_count
import math
from joblib import Parallel, delayed

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self):
        parser = argparse.ArgumentParser(description="ledock")

        parser.add_argument("--config", help="path of config file", default=False)
        args = parser.parse_args()

        config = configparser.ConfigParser()
        config.read(args.config)

        general = config["general"]
    
        #work_dir = general["work_dir"]
        
        #os.chdir(work_dir)
        
        self.receptor_name = general["receptor_name"]

        try:
            self.ligand_db_name = general["ligand_db_name"]
        except Exception as e:
            self.ligand_db_name = None
        
        if not self.ligand_db_name:
            try:
                self.ligand_file_type = general["ligand_file_type"]
            except Exception as e:
                self.ligand_file_type = None

        docking = config["docking"]

        try:
            self.get_center_x = docking.getfloat("center_x")
        except Exception as e:
            self.get_center_x = None
        
        try:
            self.get_center_y = docking.getfloat("center_y")
        except Exception as e:
            self.get_center_y = None
        
        try:
            self.get_center_z = docking.getfloat("center_z")
        except Exception as e:
            self.get_center_z = None
        
        try:
            self.get_size_x = docking.getfloat("size_x")
        except Exception as e:
            self.get_size_x = None
        
        try:
            self.get_size_y = docking.getfloat("size_y")
        except Exception as e:
            self.get_size_y = None
        
        try:
            self.get_size_z = docking.getfloat("size_z")
        except Exception as e:
            self.get_size_z = None

        try:
            self.rmsd_cutoff = docking.getfloat("rmsd_cutoff")
        except Exception as e:
            self.rmsd_cutoff = 1.0
        
        try:
            self.save_n_poses = docking.getint("max_n_poses")
        except Exception as e:
            self.save_n_poses = 20
        
        #_filter = config["filter"]
        #try:
        #    self.save_n = _filter.getint("save_n")
        #except Exception as e:
        #    self.save_n = self.save_n_poses

    def transform_db(self, db_name):
        db_type = db_name.split("/")[-1].split(".")[-1]
        db_prefix = db_name.split("/")[-1].split(".")[0]

        cmd_split = f"obabel  -i{db_type} {db_name} -O {db_prefix}_.mol2 -m"
        (_, _) = subprocess.getstatusoutput(cmd_split)

        try:
            validation_set = [cc for cc in os.listdir("./") if f"{db_prefix}_" in cc and "mol2" in cc]
        except Exception as e:
            logging.info("fail to get splited ligand, please check and run again")
            return None
        
        if validation_set:
            sorted_validate_set = sorted(validation_set, 
                                        key=lambda x: int(x.split(".")[0].split("_")[-1]))
            
            n_thread = cpu_count()
            #n_thread = 16

            run_ligand_list = []
            
            if n_thread > 1:
                while True:
                    n_in_thread = math.ceil(len(sorted_validate_set)/n_thread)
                    if math.ceil(len(sorted_validate_set)/n_in_thread) == n_thread:
                        break
                    else:
                        n_thread -= 1
            else:
                n_in_thread = math.ceil(len(sorted_validate_set)/n_thread)
            
            for idx in range(n_thread):
                _this = list(sorted_validate_set)[idx*n_in_thread: (idx+1)*n_in_thread]
                df = pd.DataFrame({"name": _this})
                df.to_csv(f"ligands.list.{idx}", header=None, index=None)
                run_ligand_list.append(f"ligands.list.{idx}")

            return run_ligand_list
        else:
            return None

    def get_seperate_db(self, file_type):
        validation_set = [cc for cc in os.listdir("./") if file_type in cc and 
                        "rec.pdb" not in cc]
        if not validation_set:
            logging.info("no available ligand, please check and run again")
            return None
        
        if "mol2" not in file_type:
            new_validation_set = []
            for idx, mm in enumerate(validation_set):
                cmd = f"obabel -i{file_type.strip('.')} {mm} -O {mm.split('.')[0]}.mol2"
                (_, _) = subprocess.getstatusoutput(cmd)
                if os.path.isfile(f"{mm.split('.')[0]}.mol2") and os.path.getsize(f"{mm.split('.')[0]}.mol2"):
                    new_validation_set.append(f"{mm.split('.')[0]}.mol2")
                    os.system(f"rm -f {mm}")
        else:
            new_validation_set = validation_set
        
        if not new_validation_set:
            logging.info("failed when ligand type transfer")
            return None
        
        sorted_validate_set = sorted(new_validation_set, 
                                        key=lambda x: int(x.split(".")[0].split("_")[-1]))
        n_thread = cpu_count()

        run_ligand_list = []
        
        if n_thread > 1:
            while True:
                n_in_thread = math.ceil(len(sorted_validate_set)/n_thread)
                if math.ceil(len(sorted_validate_set)/n_in_thread) == n_thread:
                    break
                else:
                    n_thread -= 1
        else:
            n_in_thread = math.ceil(len(sorted_validate_set)/n_thread)
        
        for idx in range(n_thread):
            _this = list(sorted_validate_set)[idx*n_in_thread: (idx+1)*n_in_thread]
            df = pd.DataFrame({"name": _this})
            df.to_csv(f"ligands.list.{idx}", header=None, index=None)
            run_ligand_list.append(f"ligands.list.{idx}")

        return run_ligand_list


    def transform_xyz(self, center_x, center_y, center_z, size_x, size_y, size_z):
        try:
            center_x * center_y * center_z * size_x * size_y * size_z
        except Exception as e:
            return None
        
        xmin = f"{center_x - size_x/2:.3f}"
        xmax = f"{center_x + size_x/2:.3f}"
        ymin = f"{center_y - size_y/2:.3f}"
        ymax = f"{center_y + size_y/2:.3f}"
        zmin = f"{center_z - size_z/2:.3f}"
        zmax = f"{center_z + size_z/2:.3f}"

        return [[xmin, xmax], [ymin,ymax],[zmin,zmax]]


    def get_docking_input(self):
        termination_flag = 0
    
        try:
            naive_input_pdb = [cc for cc in os.listdir("./") if self.receptor_name in cc]
        except Exception as e:
            naive_input_pdb = None
        if not naive_input_pdb:
            logging.info("No input receptor file, check again")
            return termination_flag
        cmd_lepro = f"lepro {naive_input_pdb[0]}"
        (_, _) = subprocess.getstatusoutput(cmd_lepro)

        if os.path.isfile("dock.in"):
            with open("dock.in", "r+") as f:
                content_original_dockin = [ff.strip() for ff in f.readlines() if ff]
            
            _dic = {"Receptor": "pro.pdb\n",
                    "RMSD": f"{self.rmsd_cutoff}\n",
                    "Binding pocket": "",
                    "Number of binding poses": f"{self.save_n_poses}\n",
                    "Ligands list": ""
                    }
            
            for idx, ff in enumerate(content_original_dockin):
                if "Binding pocket" in ff:
                    docking_box = [xx.split() for xx in content_original_dockin[idx+1:idx+4]]
            
            try:
                validate = float(docking_box[0][0])
            except Exception as e:
                get_dock_box = self.transform_xyz(self.get_center_x, 
                                             self.get_center_y, 
                                             self.get_center_z, 
                                             self.get_size_x,
                                             self.get_size_y,
                                             self.get_size_z)
            else:
                get_dock_box = docking_box
            
            if not get_dock_box:
                logging.info("Assign box information or include ligand in input_receptor_pdb")
                os.system("rm -f pro.pdb dock.in")
                return termination_flag
            
            for idx, line in enumerate(get_dock_box):
                _dic["Binding pocket"] += f"{line[0]} {line[1]}\n"

            ## db
            if self.ligand_db_name:
                get_ligand_list_name = self.transform_db(self.ligand_db_name)
            else:
                if not self.ligand_file_type:
                    logging.info("no ligand indication, please check and run again")
                    os.system("rm -f pro.pdb dock.in")
                    return termination_flag
                    
                
                get_ligand_list_name = self.get_seperate_db(self.ligand_file_type)
            
            if not get_ligand_list_name:
                os.system("rm -f pro.pdb dock.in")
                return termination_flag
            
            for idx, LigName in enumerate(get_ligand_list_name):
                _dic["Ligands list"] = LigName
                with open(f"dock.{idx}.in", "w+") as cc:
                    for kk, vv in _dic.items():
                        cc.write(f"{kk}\n")
                        cc.write(vv)
                        cc.write("\n")
            
            logging.info("normal preparison, then run for docking")
            os.system("rm -f dock.in")
            termination_flag = len(get_ligand_list_name)
            return termination_flag


    def run_ledock(self, _prefix_idx):
        cmd = f"ledock dock.{_prefix_idx}.in"
        (_, _) = subprocess.getstatusoutput(cmd)

        issued = []

        dic_energy = {}

        ## get and process result
        with open(f"ligands.list.{_prefix_idx}", "r+") as ff:
            ligand_label_list = [ss.strip() for ss in ff.readlines() if ss]
        
        for idx, ligand in enumerate(ligand_label_list):
            ligand_prefix = ligand.split(".")[0]
            try:
                get_docked_result = [cc for cc in os.listdir("./") if ligand_prefix in cc and cc.endswith(".dok")]
            except Exception as e:
                issued.append(ligand_prefix)
            else:
                if get_docked_result:
                    cmd_split = f"ledock -spli {get_docked_result[0]}"
                    (_, _) = subprocess.getstatusoutput(cmd_split)
                    try:
                        get_all_pose = [cc for cc in os.listdir("./") if f"{ligand_prefix}_dock" in cc \
                                                                        and cc.endswith(".pdb") \
                                                                        and os.path.getsize(cc)]
                    except Exception as e:
                        issued.append(ligand_prefix)
                    else:
                        #sorted_all_pose = sorted(get_all_pose, key=lambda x: int(x.split(".")[0].split("_")[-1].strip("dock")))
                        _dic = {}
                        #saved = list(sorted_all_pose)[:self.save_n]
                        saved = get_all_pose
                        #saved_in_sdf = []
                        for idx, each in enumerate(saved):
                            with open(each,"r+") as ff:
                                content = [line.strip() for line in ff.readlines() if line]
                            
                            energy_line_idx = [idx for idx in range(len(content)) if content[idx].startswith("ATOM")][0] - 1

                            get_energy = float(content[energy_line_idx].split()[-2])

                            new_line = f"{each.split('.')[0]} Score: {get_energy}"

                            content[energy_line_idx] = new_line

                            #trans_cmd = f"obabel -ipdb {each} -O {each.split('.')[0]}.sdf"
                            #(_,_) = subprocess.getstatusoutput(trans_cmd)

                            _dic.setdefault(f"{each.split('.')[0]}", [str(get_energy), content])

                            #if os.path.isfile(f"{each.split('.')[0]}.sdf") and os.path.getsize(f"{each.split('.')[0]}.sdf"):
                            #    try:
                            #        mol = [cc for cc in Chem.SDMolSupplier(f"{each.split('.')[0]}.sdf", removeHs=False) if cc]
                            #    except Exception as e:
                            #        mol = None
                            #    if mol:
                            #        this = mol[0]
                            #        this.SetProp("_Name", f"{each.split('.')[0]}")
                            #        this.SetProp("DockingScore", str(get_energy))
                            #        saved_in_sdf.append(this)
                        
                        
                        os.system(f"rm -f {ligand_prefix}.dok")

                        with open(f"./DockingPose/docked_{ligand_prefix}.pdb", "w+") as c:
                            for kk, vv in _dic.items():
                                for line in vv[-1]:
                                    c.write(f"{line}\n")
                                
                                dic_energy.update({kk: vv[0]})

                        #os.system(f"touch ./DockingPose/docked_{ligand_prefix}.pdb")
                        #for each in saved:
                        #    os.system(f"cat {each} >> ./DockingPose/docked_{ligand_prefix}.pdb")
                        
                        os.system(f"rm -f {ligand_prefix}_dock*.pdb")
                        
                        #dic_energy.update(_dic)

        return dic_energy, issued

    def run(self):
        termination_flag = self.get_docking_input()
        #termination_flag = 16
        
        if termination_flag:
            if not os.path.exists("DockingPose"):
                os.mkdir("DockingPose")
                
            col = Parallel(n_jobs=termination_flag)(delayed(self.run_ledock)(i) for i in range(termination_flag))
        else:
            logging.info("Terminated")
            return 

        energy = {}
        issued_col = []

        for each in col:
            energy.update(each[0])
            if each[1]:
                issued_col += each[1]

        
        if issued_col:
            error = pd.DataFrame({"LigName": issued_col})
            error.to_csv("ERROR", header=None, index=None)
            logging.info("Ligand(s) failed is(are) saved in ERROR")
        
        df = pd.DataFrame({"LigName_PoseID": [kk for kk in energy.keys()],
                           "DockingScore": [vv for  vv in energy.values()]})
        df.to_csv("RESULT.dat", sep="\t", index=None)

        logging.info("Normal termination, docking results saved in RESULT.csv")
        logging.info("docking poses for each ligand saved in ./DockingPose")

        #if os.path.exists("Input"):
        #    os.mkdir("Input")
        ## delete temp file
        if self.ligand_db_name:
            _prefix = self.ligand_db_name.split(".")[0]
            need_to_delect = [cc for cc in os.listdir("./") if cc.startswith(_prefix) and cc.endswith(".mol2")]
            for each in need_to_delect:
                os.remove(each)
        
        ## tmp file
        if not os.path.exists("Setting"):
            os.mkdir("Setting")
        os.system("mv dock.*.in ./Setting")
        os.system("mv ligands.list.* ./Setting")
        os.system("mv pro.pdb ./Setting")

        return 
    

if __name__ == "__main__":
    main().run()
