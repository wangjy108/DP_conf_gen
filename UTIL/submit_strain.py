import os
import json
import math
import subprocess
import argparse

def setup_lbg(projectID, _path):
    lbg_json = {
        "job_name": "strain",
        "command": "bash run.sh",
        "log_file": "tmp_log",
        "backward_files": [],
        "program_id": projectID,
        "platform": "ali",
        "job_group_id": "",
        "disk_size": 128,
        "machine_type": "c16_m128_cpu",
        "job_type": "container",
        "image_name": "registry.dp.tech/dptech/prod-1364/strain:run0.0.6"
    }
    
    file_name = os.path.join(_path, "input.json")
    with open(file_name, "w+") as cc:
        json.dump(lbg_json, cc, indent=4)
    return

def setup_local(_path):
    file_name = os.path.join(_path, "run.sh")
    
    with open(file_name, "w+") as cc:
        cc.write("#!/bin/sh \n")
        cc.write("\n")
        cc.write("runfile=`ls *.sdf`\n")
        cc.write("\n")
        cc.write("for each in ${runfile}\n")
        cc.write("do\n")
        cc.write("\tpython sample_pip_strain.py --input_sdf ${each}\n")
        cc.write("done\n")
        cc.write("\n")
    #logging.info("Prepared for local run, run by 'nohup bash run.sh &' to start running")
    return

def run(N_in_group, projectID):
    main_dir = os.getcwd()
    _all = [sdf for sdf in os.listdir(main_dir) if ".sdf" in sdf]

    n_group = math.ceil(len(_all) / N_in_group)

    for i in range(n_group):
        _work_dir = os.path.join(main_dir, str(i))
        if not os.path.exists(_work_dir):
            os.mkdir(_work_dir)
        os.chdir(_work_dir)
        print(f"At {i}th path")
        get_set = _all[i*N_in_group:(i+1)*N_in_group]

        _flag = 0

        for each in get_set:
            original_file = os.path.join(main_dir, each)
            target_file = os.path.join(_work_dir, each)
            cmd = f"mv {original_file} {target_file}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isfile(target_file):
                _flag += 1
        
        if _flag == len(get_set):
            setup_lbg(projectID=projectID, _path=_work_dir)
            setup_local(_path=_work_dir)
            cmd_cp = f"cp {os.path.join(main_dir, 'sample_pip_strain.py')} ./"
            (_, _) = subprocess.getstatusoutput(cmd_cp)

            cmd_submit = f"lbg job submit -i input.json -p ./ > TRACK_JOB"
            (_,_) = subprocess.getstatusoutput(cmd_submit)
        
        else:
            print(f"{i} failed")

        os.chdir(main_dir)
    return 

if __name__ == "__main__":
    #_path = os.getcwd()
    parser = argparse.ArgumentParser(description='serial submission of strain calc')
    parser.add_argument('--projectID', type=int, required=True, 
                        help='lbg project id')
    parser.add_argument('--N_in_group', type=int, required=True, 
                        help='N mol in each serial run')
    args = parser.parse_args()

    run(args.N_in_group, args.projectID)
    


    