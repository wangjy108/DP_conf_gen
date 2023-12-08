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
        parser = argparse.ArgumentParser(description="Conf Sample")

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
            self.max_N = lbg_config.getint("max_N")
        except Exception as e:
            self.max_N = 100
        
        self.general = config["general"]

        self.work_dir = self.general["work_dir"]
        self.main_dir = os.getcwd()

        
        try:
            self.rmsd_cutoff = float(self.general["rmsd_cutoff"])
        except Exception as e:
            self.rmsd_cutoff = 1.5

        
        #os.chdir(self.main_work_dir)

        self.run_type = args.run_type
    
    def col_from_seperate_ligand(self):
         
        os.chdir(self.work_dir)

        track_path = []
        ligand_list = [cc for cc in os.listdir("./") if cc.endswith(".sdf")]
        if ligand_list:
            group = math.ceil(len(ligand_list) / self.max_N)

            for i in range(group):
                content_in_this_group = ligand_list[i*self.max_N:(i+1)*self.max_N]

                if not os.path.exists(f"{i}"):
                    os.mkdir(f"{i}")
            
                save_dir = os.path.join(self.work_dir, f"{i}")

                for each in content_in_this_group:
                    os.system(f"mv {each} {save_dir}")
                
                track_path.append(f"{i}")

        else:
            logging.info(f"No ligand in work dir, abort")
            return None
        
        return track_path
    
    def setup_local(self, _path):
        file_name = os.path.join(_path, "run.sh")
        
        with open(file_name, "w+") as cc:
            cc.write("#!/bin/sh \n")
            cc.write("\n")
            cc.write("runfile=`ls *.sdf`\n")
            cc.write("\n")
            cc.write("for each in ${runfile}\n")
            cc.write("do\n")
            cc.write(f"\tpython ConfSample.py --input_sdf ${{each}} --rmsd_cutoff {self.rmsd_cutoff}\n")
            cc.write("done\n")
            cc.write("\n")
        #logging.info("Prepared for local run, run by 'nohup bash run.sh &' to start running")
        return
        
    def json_config(self, _path):
        save_file = os.path.join(_path, "input.json")
        
        lbg_json = {
            "job_name": "confGen",
            "command": "bash ./run.sh",
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
    
    def col_result(self, csv):

        logging.info("Collecting results ....")

        df = pd.read_csv(csv, header=0)

        for idx, row in df.iterrows():
            os.chdir(row["local_path"])
            cmd = f"lbg job download {row['JOB_ID']}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isdir(str(row['JOB_ID'])):
                (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")
                (_, _) = subprocess.getstatusoutput("rm -f ST* TRACK* input.json *.sh *.py")
                (_, _) = subprocess.getstatusoutput(f"mv * {self.work_dir}")

                os.chdir(self.main_dir)
                shutil.rmtree(row["local_path"]) 
            else:
                logging.info(f"Fail to collect result in {row['local_path']}")
                os.chdir(self.main_dir)

        logging.info("Finishe collect")
        
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

            os.chdir(self.main_dir)
        
        df_track = pd.DataFrame({"local_path": [kk for kk in col_result.keys()],
                                 "JOB_ID": [vv for vv in col_result.values()]})
        
        df_track.to_csv("submission_track_list", index=None)

        logging.info("Finish submission")

        return 

    def prepare(self):
        os.chdir(self.main_dir)
        submission_track = []
        if not self.lbg_project:
            logging.info("Please set lbg project id, abort")
            return None
        if not self.lbg_image_name:
            logging.info("No available image name assigned, abort")
            return None

        ## define input ligand type and rewrite path
        tracker = self.col_from_seperate_ligand()
        if not tracker:
            logging.info("Nothing to submit, abort")
            return None
        for _path in tracker:
            submission_track.append(os.path.join(self.work_dir, _path))

            self.setup_local(_path=_path)
            self.json_config(_path=_path)
            os.system(f"cp {os.path.join(self.main_dir, 'ConfSample.py')} {_path}")
            
        return submission_track
    
    def run(self):
        if "submit" in self.run_type:
            track = self.prepare()
            if not track:
                return 
            self.submit(track)
        elif "col_result" in self.run_type:
            os.chdir(self.main_dir)
            self.col_result(csv="submission_track_list")
            os.system("rm -f submission_track_list")
        else:
            logging.info(f"Not available run type: {self.run_type}")
        
        return 


if __name__ == "__main__":
    main().run()
    
        
        





        

        
        

            





                    
                    



    





