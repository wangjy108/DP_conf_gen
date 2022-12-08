import pandas as pd
import os
import random
import math
import numpy as np
import time
import argparse
import logging
#from logging.config import fileConfig
#from joblib import Parallel, delayed
#fileConfig('logging_config.ini')
logging.basicConfig(format='[%(levelname)s] %(message)s', \
                    filename='out.STDERROR',\
                    filemode='w', \
                    level=logging.INFO)

def get_charge(name_label):
    fchk = f"{name_label}.fchk"
    os.system(f"Multiwfn {fchk} < resp.in")
    return f"{name_label}.chg"

def get_raw_mol2(name_label):
    os.system(f"antechamber -i {name_label}.log -fi gout -o {name_label}.mol2 -fo mol2 \
               -at sybyl -dr no -pf y")
    return f"{name_label}.mol2"

def get_final_mol2(name_label):
    try:
        raw_mol2 = get_raw_mol2(name_label)
    except Exception as e:
        print(f"check ambertool for raw mol2 generation with {name_label}")
        return False

    try:
        chg = get_charge(name_label)
    except Exception as e:
        print(f"check multiwfn for charge generation with {name_label}")
        return False

    with open(chg, "r") as ff:
        charge = [float(line.strip()[38:]) for line in ff.readlines() if len(line)>0]

    with open(raw_mol2, "r") as ff:
        info = [line for line in ff.readlines()]

    if min(len(info), len(charge)) == 0:
        print(f"something wrong with {name_label}")
        return False

    idx = [ii for ii in range(len(info)) if info[ii].startswith("@")]

    with open(f"{name_label}.2.mol2", "w+") as cc:
        header = info[idx[0]:idx[1]+1]

        for hh in header:
            cc.write(hh)

        ori_atom = [ll[:73] for ll in info[idx[1]+1:idx[2]]]

        for i in range(len(ori_atom)):
            each_atom = ori_atom[i] + str(charge[i])
            cc.write(each_atom + "\n")

        other = info[idx[2]:]
        for oo in other:
            cc.write(oo)

    return True


if __name__ == "__main__":

    if not os.path.exists("resp.in"):
        resp_in = ['7', '18', '1', ' ', 'y', '0', '0', 'q']
        with open("resp.in", "w+") as cc:
            for ii in resp_in:
                cc.write(ii + '\n')

    get_input_name = [cc.strip().split(".")[0] for cc in os.listdir() if 'fchk' in cc]
    for name in get_input_name:
        #print(name)
        try:
            flag = get_final_mol2(name)
        except Exception as e:
            logging.warning(f"ERROR for {name.strip().split()[0]}")
            continue
        else:
            if not flag:
                logging.warning(f"ERROR for {name.strip().split()[0]}")
                continue










###
