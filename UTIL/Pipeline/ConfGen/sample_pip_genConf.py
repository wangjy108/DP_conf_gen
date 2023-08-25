import argparse
import configparser
import json
import os
import math
import logging
import shutil
import subprocess
import pandas as pd

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


class main():
    def __init__(self):
        parser = argparse.ArgumentParser(description="DPDH_pipeline_for_confGen(sub-vision from strain)")

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
            self.lbg_image_name = lbg_config["lbg_image_name"]
        except Exception as e:
            self.lbg_image_name = None
        
        try:
            self.max_N = int(lbg_config["max_N"])
        except Exception as e:
            self.max_N = 10000
        
        general = config["general"]
        try:
            self.input_smi_db = general["input_smi_db"]
        except Exception as e:
            self.input_smi_db = None
        
        try:
            self.energy_window = general.getfloat("energy_window")
        except Exception as e:
            self.energy_window = 5.0
        
        try:
            self.rmsd_cutoff = general.getfloat("rmsd_cutoff")
        except Exception as e:
            self.rmsd_cutoff = 1.0
        
        self.main_dir = os.getcwd()
        self.run_type = args.run_type

        self._prefix = self.input_smi_db.split(".")[0]
        self._tail = self.input_smi_db.split(".")[-1]

    def prepare_db(self):
        if not self.input_smi_db:
            return None
        
        with open(self.input_smi_db) as f:
            try:
                content = [ff for ff in f.readlines() if ff]
            except Exception as e:
                return None
            if not content:
                return None
        
        batch = math.ceil(len(content)/self.max_N)

        track_path = {}

        for idx in range(batch):
            save_line = content[idx*self.max_N: (idx+1)*self.max_N]
            work_dir = os.path.join(self.main_dir, str(idx))
            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
            save_file = os.path.join(work_dir, f"{self._prefix}_{idx}.{self._tail}")
            with open(save_file, "w+") as cc:
                for each in save_line:
                    cc.write(each)
            track_path.setdefault(work_dir, f"{self._prefix}_{idx}.{self._tail}")
                
        return track_path
    
    def json_config(self, _path, cmd):
        save_file = os.path.join(_path, "input.json")
        
        lbg_json = {
            "job_name": "confGen",
            "command": f"{cmd}",
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.lbg_project}" ,
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 128,
            "machine_type": f"{self.lbg_machine_type}",
            "job_type": "container",
            "image_name": f"{self.lbg_image_name}"
        }

        with open(save_file, "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
    
    def submit(self):
        if not (self.lbg_image_name and self.lbg_project):
            logging.info("Please assign lbg related parameter in config.in, abort")
            return None
        
        logging.info("Preparing ...")
        
        track_path = self.prepare_db()
        if not track_path:
            logging.info("Failed with input smile DB, abort")
            return None
        
        track_run_done = {}

        logging.info("Run submit ....")
        
        for _path, _db_name in track_path.items():
            os.chdir(_path)
            each_cmd = f"python gen.py --input_file {_db_name} --energy_window {self.energy_window} --rmsd_cutoff {self.rmsd_cutoff}"
            
            self.json_config(_path=_path, cmd=each_cmd)

            cmd_file = os.path.join("/opt/scripts", "gen.py")
            if not os.path.isfile(cmd_file):
                cmd_file = os.path.join(self.main_dir, "gen.py")
            
            if not os.path.isfile(cmd_file):
                logging.info("No available command file, abort")
                return None
            
            cp_cmd = f"cp {cmd_file} {_path}"
            (_, _) = subprocess.getstatusoutput(cp_cmd)

            if os.path.isfile("input.json") and os.path.isfile("gen.py"):
                submission_cmd = "lbg job submit -i input.json -p ./ > TRACK_JOB"
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
                        track_run_done.setdefault(_path, JOB_ID)
                else:
                    logging.info(f"Submission failed at {_path}, check in folder for more information")
                    
            else:
                logging.info(f"Submission failed at {_path}, check in folder for more information")
    
            
        df_track = pd.DataFrame({"local_path": [kk for kk in track_run_done.keys()],
                                "JOB_ID": [vv for vv in track_run_done.values()]})
        
        df_track.to_csv(os.path.join(self.main_dir,"submission_track_list"), index=None)

        #logging.info("Finish submission")

        os.chdir(self.main_dir)

        return "submission_track_list"
    
    def collect(self):
        os.chdir(self.main_dir)
        logging.info("Collecting results ....")

        _failed_download = {}
        _success_download = {}

        if not os.path.isfile("submission_track_list"):
            logging.info("No tracking information to collect from, abort")
            return 
        
        df = pd.read_csv("submission_track_list", header=0)

        for idx, row in df.iterrows():
            _path = row["local_path"]
            JOB_ID = row["JOB_ID"]

            os.chdir(_path)

            download_cmd = f"lbg job download {JOB_ID}"
            (_, _) = subprocess.getstatusoutput(download_cmd)

            if os.path.isdir(str(JOB_ID)):
                logging.info(f"Collect JOB {JOB_ID}")
                (_, _) = subprocess.getstatusoutput(f"mv {JOB_ID}/* ./")
                assembled_sdf = [cc for cc in os.listdir("./") if ".sdf" in cc and "Gen" in cc]
            
            else:
                logging.info(f"Failed collect JOB {JOB_ID}")
                assembled_sdf = None
            
            if assembled_sdf:
                get_sdf = os.path.join(_path, assembled_sdf[0])
                _success_download.setdefault(_path, get_sdf)
                shutil.rmtree(str(JOB_ID))
            else:
                _failed_download.setdefault(_path, JOB_ID)
        
        if _success_download:
            os.chdir(self.main_dir)
            
            root_save = os.path.join(self.main_dir, f"{self._prefix}.sdf")

            if not os.path.isfile(root_save):
                os.system(f"touch {root_save}")
            
            for _path, sdf in _success_download.items():
                _cmd = f"cat {sdf} >> {root_save}"
                (_, _) = subprocess.getstatusoutput(_cmd)
                shutil.rmtree(_path)
        else:
            logging.info("Can not collect result files")
        
        if _failed_download:
            logging.info("Failed for collecting results and corresponding JOB IDs are saved in ERROR")
            df_error = pd.DataFrame({"Path": [kk for kk in _failed_download.keys()],
                                     "JOB ID": [vv for vv in _failed_download.values()]})
            df_error.to_csv(os.path.join(self.main_dir, "ERROR"), index=None)
            return 
        else:
            logging.info("Collecting results done")
        
        return 
    
    def run(self):
        if "submit" in self.run_type:
            track = self.submit()
            if track:
                logging.info("Submission done")

        elif "col_result" in self.run_type:
            self.collect()
            os.system("rm -f submission_track_list gen.py")
        else:
            logging.info(f"Not available run type: {self.run_type}")
        
        return
            

if __name__ == "__main__":
    main().run()










                    



