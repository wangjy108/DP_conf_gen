import os
import logging
import pandas as pd
import argparse
import configparser
import subprocess
import math
import json
import subprocess
import shutil

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self):
        parser = argparse.ArgumentParser(description="ledock")

        parser.add_argument("--config", help="path of config file", default=False)
        parser.add_argument("--run_type", help="submit or col_result", default=False)
        args = parser.parse_args()

        config = configparser.ConfigParser()
        config.read(args.config)

        lbg_config = config["lbg_config"]

        try:
            self.lbg_project = lbg_config["lbg_project"]
        except Exception as e:
            self.lbg_project = None
        
        try:
            self.lbg_machine_type = lbg_config["lbg_machine_type"]
        except Exception as e:
            self.lbg_machine_type = "c16_m128_cpu"
        if not self.lbg_machine_type:
            self.lbg_machine_type = "c16_m128_cpu"

        try:
            self.max_N = lbg_config.getint("max_N")
        except Exception as e:
            self.max_N = 20000

        
        self.general = config["general"]

        self.main_work_dir = self.general["work_dir"]

        self.docking = config["docking"]

        self.filter = config["filter"]

        try:
            self.receptor_name = self.general["receptor_name"]
        except Exception as e:
            self.receptor_name = None
        else:
            if not self.receptor_name:
                self.receptor_name = None

        try:
            self.ligand_db_name = self.general["ligand_db_name"]
        except Exception as e:
            self.ligand_db_name = None
        
        if not self.ligand_db_name:
            try:
                self.ligand_file_type = self.general["ligand_file_type"]
            except Exception as e:
                self.ligand_file_type = None
        
        os.chdir(self.main_work_dir)

        self.run_type = args.run_type
    
    def col_from_seperate_ligand(self, ligand_file_type):
        track_path = []
        track_name = {}
        
        ligand_list = [cc for cc in os.listdir("./") if cc.endswith(ligand_file_type)]
        if ligand_list:
            group = math.ceil(len(ligand_list) / self.max_N)

            for i in range(group):
                content_in_this_group = ligand_list[i*group:(i+1)*group]

                if not os.path.exists(f"{i}"):
                    os.mkdir(f"{i}")
            
                save_dir = os.path.join(self.main_work_dir, f"{i}")

                for each in content_in_this_group:
                    os.system(f"mv {each} {save_dir}")
                
                track_path.append(f"{i}")

                for each in content_in_this_group:
                    track_name.setdefault(each.split(".")[0], each.split(".")[0])

        else:
            logging.info(f"No ligand in work dir by type {ligand_file_type}, abort")
            return None
        
        return track_path, track_name
    
    def splite_from_db(self):
        _dic = {"mol2": self.split_from_mol2,
                "sdf": self.split_from_sdf}
        
        get_db_type = self.ligand_db_name.split(".")[-1]

        try:
            _dic[get_db_type]
        except Exception as e:
            logging.info(f"DB type {get_db_type} not recognized, available for mol2 and sdf, abort")
            return None
        
        track_path, track_name = _dic[get_db_type](self.ligand_db_name)
        
        return track_path, track_name
    
    def split_from_mol2(self, db_name):
        track_path = {}
        track_name = {}
        _label = db_name.split(".")[0]
        _type = db_name.split(".")[-1]

        with open(db_name, "r+") as cc:
            content = [ff for ff in cc.readlines()]
        
        start_idx = [idx for idx in range(len(content)) if content[idx].startswith("@<TRIPOS>MOLECULE")]

        group = math.ceil(len(start_idx)/self.max_N)

        for i in range(group):
            try:
                content_in_this_group = content[start_idx[i*self.max_N]:start_idx[(i+1)*self.max_N]]
            except Exception as e:
                content_in_this_group = content[start_idx[i*self.max_N]:]
            
            if not os.path.exists(f"{i}"):
                os.mkdir(f"{i}")
            
            save_dir = os.path.join(os.path.join(self.main_work_dir, f"{i}"),f"{_label}_{i}.{_type}")

            with open(save_dir, "w+") as cc:
                for each_line in content_in_this_group:
                    cc.write(each_line)
            
            track_path.setdefault(f"{i}", f"{_label}_{i}.{_type}")

            group_start_idx = [idx for idx in range(len(content_in_this_group)) if content_in_this_group[idx].startswith("@<TRIPOS>MOLECULE")]

            for idx, line_idx in enumerate(group_start_idx):
                get_real_name = content_in_this_group[line_idx+1].strip()
                if not get_real_name:
                    get_real_name = f"{_label}_{i*self.max_N + i}"
                track_name.setdefault(f"{_label}_{i}_{idx+1}", get_real_name)

        return track_path, track_name
    
    def split_from_sdf(self, db_name):
        track_path = {}
        _label = db_name.split(".")[0]
        _type = db_name.split(".")[-1]

        track_name = {}

        with open(db_name, "r+") as f:
            content = [ff for ff in f.readlines()]
        
        end_idx = [-1] + [idx for idx, line in enumerate(content) if "$$$$" in line][:-1]

        group = math.ceil(len(end_idx)/self.max_N)

        for i in range(group):
            try:
                content_in_this_group = content[end_idx[i*self.max_N]+1:end_idx[(i+1)*self.max_N]+1]
            except Exception as e:
                content_in_this_group = content[end_idx[i*self.max_N]+1:]
            
            if not os.path.exists(f"{i}"):
                os.mkdir(f"{i}")
            
            save_dir = os.path.join(os.path.join(self.main_work_dir, f"{i}"),f"{_label}_{i}.{_type}")

            with open(save_dir, "w+") as cc:
                for each_line in content_in_this_group:
                    cc.write(each_line)
            
            track_path.setdefault(f"{i}", f"{_label}_{i}.{_type}")

            group_end_idx = [-1] + [idx for idx, line in enumerate(content_in_this_group) if "$$$$" in line][:-1]

            for idx, line_idx in enumerate(group_end_idx):
                get_real_name = content_in_this_group[line_idx+1].strip()
                if not get_real_name:
                    get_real_name = f"{_label}_{i*self.max_N + i}"
                track_name.setdefault(f"{_label}_{i}_{idx+1}", get_real_name)
    
        return track_path, track_name
    
    def json_config(self, _path):
        save_file = os.path.join(_path, "input.json")
        
        lbg_json = {
            "job_name": "ledock",
            "command": "python run_ledock.py --config config.in",
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.lbg_project}" ,
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 128,
            "machine_type": f"{self.lbg_machine_type}",
            "job_type": "container",
            "image_name": "registry.dp.tech/dptech/prod-1364/ledock:public0.0.2"
        }

        with open(save_file, "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
    
    def config_file(self,
                    _path,
                    config_dic):
        
        save_file = os.path.join(_path, "config.in")
        with open(save_file, "w+") as cc:
            for kk, vv in config_dic.items():
                cc.write(f"[{kk}]\n")
                for kkk, vvv in vv.items():
                    cc.write(f"{kkk}={vvv}\n")
                cc.write("\n")
        
        return
    
    def col_result(self, csv, _name):

        logging.info("Collecting results ....")
        if not os.path.exists("DockingPose"):
            os.mkdir("DockingPose")

        df = pd.read_csv(csv, header=0)
        df_name = pd.read_csv(_name, header=0)
        name_label = {}

        list_result = []
        should_not_delet = []

        for _idx, _row in df_name.iterrows():
            name_label.setdefault(_row["LigandSepName"], _row["LigandRealName"])

        for idx, row in df.iterrows():
            os.chdir(row["local_path"])
            cmd = f"lbg job download {row['JOB_ID']}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isdir(str(row['JOB_ID'])):
                (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")
                
                if os.path.isfile("RESULT.dat") and os.path.getsize("RESULT.dat"):
                    _result = pd.read_csv("RESULT.dat", sep='\t', header=0)

                    _result["LigPose"] = _result["LigName_PoseID"].apply(\
                        lambda x: f"{name_label['_'.join(x.split('_')[:-1])]}_pose{x[-3:]}")
                    
                    _result = _result[["LigPose", "DockingScore"]]
                    
                    list_result.append(_result)
                else:
                    should_not_delet.append(row["local_path"])
                
                if os.path.exists("DockingPose"):
                    os.chdir("DockingPose")
                    for each in os.listdir("./"):
                        _need_change = "_".join(each.strip().split(".")[0].split("_")[1:])
                        get_real_name = name_label[_need_change]
                        os.system(f"mv {each} docked_{get_real_name}.pdb")

                    os.system(f"mv ./* {os.path.join(self.main_work_dir, 'DockingPose')}")
                    os.chdir(row["local_path"])
                else:
                    should_not_delet.append(row["local_path"])

            else:
                logging.info(f"Fail to collect result in {row['local_path']}")
            
            os.chdir(self.main_work_dir)
        
        if list_result:
            df_final = pd.concat(list_result)
            df_final.to_csv("RESULT.dat", index=None)
        
        for _path in list(df["local_path"]):
            if _path not in should_not_delet:
                shutil.rmtree(_path)

        
        logging.info("Finish collect")
        
        return 
    
    def submit(self, submission_track:list):
        logging.info("Run submit ....")

        submission_cmd = "lbg job submit -i input.json -p ./ > TRACK_JOB"

        col_result = {}

        for _path in submission_track:
            os.chdir(_path)
            (_, _) = subprocess.getstatusoutput(submission_cmd)

            if os.path.isfile("TRACK_JOB") and os.path.getsize("TRACK_JOB"):
                with open("TRACK_JOB", "r+") as f:
                    content = [ll.strip() for ll in f.readlines() if ll]
                
                if content:
                    try:
                        JOB_ID = content[-1].split()[-1]
                    except Exception as e:
                        JOB_ID = None
                else:
                    JOB_ID = None
                
                if JOB_ID:
                    col_result.setdefault(_path, JOB_ID)
                else:
                    logging.info(f"Submission failed at {_path}, check in folder for more information")
            else:
                logging.info(f"Submission failed at {_path}, check in folder for more information")

            os.chdir(self.main_work_dir)
        
        df_track = pd.DataFrame({"local_path": [kk for kk in col_result.keys()],
                                 "JOB_ID": [vv for vv in col_result.values()]})
        
        df_track.to_csv("submission_track_list", index=None)

        logging.info("Finish submission")

        return 

    def prepare(self):
        submission_track = []
        if not self.lbg_project:
            logging.info("Please set lbg project id, abort")
            return None
        if not self.receptor_name:
            logging.info("No receptor defined, abort")
            return None
        if not os.path.isfile(self.receptor_name):
            logging.info("Can't find defined receptor file, abort")
            return 
        
        receptor = os.path.join(self.main_work_dir, self.receptor_name)

        ## define input ligand type and rewrite path
        if self.ligand_db_name:
            track_path, track_name = self.splite_from_db()
            if not track_path:
                return None
            logging.info("Prepare for cloud run distribution")
            for _path, ligand_db_name in track_path.items():
                config_dic = {}
                self.general["work_dir"] = os.path.join(self.main_work_dir, _path)
                submission_track.append(os.path.join(self.main_work_dir, _path))
                self.general["ligand_db_name"] = ligand_db_name

                config_dic.setdefault("general", self.general)
                config_dic.setdefault("docking", self.docking)
                config_dic.setdefault("filter", self.filter)

                self.config_file(_path=_path, config_dic=config_dic)
                self.json_config(_path=_path)
                os.system(f"cp /opt/scripts/run_ledock.py {os.path.join(self.main_work_dir, _path)}")
                os.system(f"cp {receptor} {os.path.join(self.main_work_dir, _path)}")
        else:
            if self.ligand_file_type:
                track_path, track_name = self.col_from_seperate_ligand(self.ligand_file_type)
                if not track_path:
                    return None
                for _path in track_path:
                    config_dic = {}
                    self.general["work_dir"] = os.path.join(self.main_work_dir, _path)
                    submission_track.append(os.path.join(self.main_work_dir, _path))

                    config_dic.setdefault("general", self.general)
                    config_dic.setdefault("docking", self.docking)
                    config_dic.setdefault("filter", self.filter)

                    self.config_file(_path=_path, config_dic=config_dic)
                    self.json_config(_path=_path)
                    os.system(f"cp /opt/scripts/run_ledock.py {os.path.join(self.main_work_dir, _path)}")
                    os.system(f"cp {receptor} {os.path.join(self.main_work_dir, _path)}")
            else:
                logging.info("Assign at least one ligand related parameter in config.in, abort")
                return None
        
        df_name = pd.DataFrame({"LigandSepName": [kk for kk in track_name.keys()],
                                "LigandRealName": [vv for vv in track_name.values()]})
        df_name.to_csv("_TEMP_name", index=None)

        
        return submission_track
    
    def run(self):
        if "submit" in self.run_type:
            track = self.prepare()
            self.submit(track)
        elif "col_result" in self.run_type:
            os.chdir(self.main_work_dir)
            self.col_result(csv="submission_track_list", _name="_TEMP_name")
            os.system("rm -f _TEMP_name submission_track_list")
        else:
            logging.info(f"Not available run type: {self.run_type}")
        
        return 


if __name__ == "__main__":
    main().run()
    
        
        





        

        
        

            





                    
                    



    





