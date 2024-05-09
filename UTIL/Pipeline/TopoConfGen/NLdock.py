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
from pathlib import Path


logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self):
        parser = argparse.ArgumentParser(description="NLdock")

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
            self.max_N = lbg_config.getint("max_N")
        except Exception as e:
            self.max_N = 20

        
        self.general = config["general"]

        try:
            self.input_csv = self.general["input_csv"]
        except Exception as e:
            self.input_csv = None

        self.handler = config["handler"]
        
        self.work_dir = os.getcwd()

        self.run_type = args.run_type
    
    def json_config(self, _path, cmd):
        save_file = os.path.join(_path, "input.json")
        
        lbg_json = {
            "job_name": "topo conf gen",
            "command": cmd,
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.lbg_project}" ,
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 64,
            "machine_type": "c16_m64_cpu",
            "job_type": "container",
            "image_name": "registry.dp.tech/dptech/prod-1364/confgen:topoI0.3"
        }

        with open(save_file, "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
    
    def check_input(self):

        submission_track = []

        if not self.input_csv:
            logging.info("Bad input csv, abort")
            return None

        if not os.path.isfile(os.path.join(self.work_dir, "main.py")):
            logging.info("missing main.py for processing, abort")
            return None
        
        if not self.lbg_project:
            logging.info("Please set lbg project id, abort")
            return None
        
        df = pd.read_csv(self.input_csv, header=0).reset_index()

        group = math.ceil(df.shape[0] / self.max_N)

        i = 0

        while i < group:
            if not os.path.exists(os.path.join(self.work_dir, str(i))):
                os.mkdir(os.path.join(self.work_dir, str(i)))
            os.chdir(os.path.join(self.work_dir, str(i)))
            current_df = df.iloc[i*self.max_N:(i+1)*self.max_N, :]

            current_df.to_csv("input.csv", index=None)

            des = os.path.join(os.path.join(self.work_dir, str(i)), "main.py")

            shutil.copy2(os.path.join(self.work_dir, "main.py"), des)

            cmd = f"python main.py --input_csv input.csv "

            if self.handler:
                for kk, vv in self.handler.items():
                    if vv:
                        cmd += f"--{kk} {vv} "
            
            self.json_config(_path=os.path.join(self.work_dir, str(i)),
                             cmd=cmd)
            
            submission_track.append(os.path.join(self.work_dir, str(i)))
            
            i += 1
            os.chdir(self.work_dir)

        logging.info("Well prepared")
        return submission_track
    
    
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

            os.chdir(self.work_dir)
        
        df_track = pd.DataFrame({"local_path": [kk for kk in col_result.keys()],
                                 "JOB_ID": [vv for vv in col_result.values()]})
        
        df_track.to_csv("submission_track_list", index=None)

        logging.info("Finish submission")

        return

    def col_result(self):

        logging.info("Collecting results ....")

        if not os.path.isfile(os.path.join(self.work_dir, "submission_track_list")):
            logging.info("Nothing to collect, abort")
            return 

        df = pd.read_csv("submission_track_list", header=0)
        sdf_contents = []
        error_message = []

        should_not_delete = []

        for idx, row in df.iterrows():
            os.chdir(row["local_path"])
            cmd = f"lbg job download {row['JOB_ID']}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isdir(str(row['JOB_ID'])):
                (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")

                if os.path.exists(os.path.join(row["local_path"], "output.sdf")):
                    with open(os.path.join(row["local_path"], "output.sdf"), "r+") as f:
                        content = [ll for ll in f.readlines()]
                    
                    sdf_contents += content
                
                if os.path.exists(os.path.join(row["local_path"], "ERROR.smi")):
                    with open(os.path.join(row["local_path"], "ERROR.smi"), "r+") as e:
                        error = [ll for ll in e.readlines()]

                    error_message += error

            else:
                logging.info(f"Fail to collect result in {row['local_path']}")
                should_not_delete.append(row['local_path'])
            
            os.chdir(self.work_dir)
        
        for _path in list(df["local_path"]):
            if _path not in should_not_delete:
                shutil.rmtree(_path)
        
        logging.info("Finish collect")

        output_dir = Path(self.work_dir)

        output_dir.mkdir(parents=True, exist_ok=True)
        with (output_dir / "output.sdf").open("w+") as f:
            for each in sdf_contents:
                f.write(each)

        if error_message:
            with (output_dir / "ERROR.smi").open("w+") as ff:
                for ee in error_message:
                    ff.write(ee)
        
        logging.info("Check ./output.sdf for final results and ERROR.smi for errored ones if exists")
        
        return
    
    def run(self):
        if "submit" in self.run_type:
            track = self.check_input()
            self.submit(track)
        elif "col_result" in self.run_type:
            os.chdir(self.work_dir)
            self.col_result()
            os.system("rm -f submission_track_list")
        else:
            logging.info(f"Not available run type: {self.run_type}")
        
        return 

if __name__ == "__main__":
    main().run()
    


