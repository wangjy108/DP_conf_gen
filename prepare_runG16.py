import pandas as pd
from rdkit import Chem
import os
import random
import shutil
from rdkit.Chem import Draw
import math
import numpy as np
import time
import argparse

def define_gemopt(initial_gjf, prefix):
    with open(initial_gjf, "r") as f:
        line = [line.strip() for line in f.readlines()]

    space_idx = [i for i in range(len(line)) if len(line[i])==0]

    out_put = os.path.join(os.getcwd(), f"{prefix}.gjf")
    with open(out_put, "w+") as c:
        chk = f"%chk={prefix}.chk"
        c.write(chk + "\n")
        mem = "%mem=50GB"
        c.write(mem + "\n")
        nproc = "%nprocshared=32"
        c.write(nproc + "\n")
        cmd = "# B3LYP/6-311G* em=GD3BJ opt scrf(SMD,solvent=water)"
        c.write(cmd + "\n")
        c.write("\n")
        c.write("opt run\n")
        c.write("\n")

        save = line[space_idx[1]+1:space_idx[2]]
        for ll in save:
            c.write(ll + "\n")

        c.write("\n")
        c.write("\n")
        c.write("\n")
    return print(f"opt job prepared at {out_put}")

def define_sp(prefix):
    out_put = os.path.join(os.getcwd(), f"sp_{prefix}.gjf")
    with open(out_put, "w+") as c:
        chk = f"%oldchk={prefix}.chk"
        c.write(chk + "\n")
        mem = "%mem=10GB"
        c.write(mem + "\n")
        nproc = "%nprocshared=14"
        c.write(nproc + "\n")
        cmd = "# B2PLYPD3/def2TZVP geom=allcheck"
        c.write(cmd + "\n")

        c.write("\n")
        c.write("\n")
        c.write("\n")
    return print(f"High level single point calculation job prepared at {out_put}")


def define_gas(prefix):
    out_put = os.path.join(os.getcwd(), f"gas_{prefix}.gjf")
    with open(out_put, "w+") as c:
        chk = f"%oldchk={prefix}.chk"
        c.write(chk + "\n")
        mem = "%mem=10GB"
        c.write(mem + "\n")
        nproc = "%nprocshared=14"
        c.write(nproc + "\n")
        cmd = "# M052X/6-31G* geom=allcheck"
        c.write(cmd + "\n")

        c.write("\n")
        c.write("\n")
        c.write("\n")
    return print(f"Gas calculation job prepared at {out_put}")

def define_sol(prefix):
    out_put = os.path.join(os.getcwd(), f"solv_{prefix}.gjf")
    with open(out_put, "w+") as c:
        chk = f"%oldchk={prefix}.chk"
        c.write(chk + "\n")
        mem = "%mem=10GB"
        c.write(mem + "\n")
        nproc = "%nprocshared=14"
        c.write(nproc + "\n")
        cmd = "# M052X/6-31G* scrf=SMD geom=allcheck"
        c.write(cmd + "\n")

        c.write("\n")
        c.write("\n")
        c.write("\n")
    return print(f"Solv calculation job prepared at {out_put}")

def run_gaussian(gjf):
    if not os.path.isfile("run_g16.sh"):
        os.system("cp /data1/zhengdan/jywang/tool_file/qm/run_g16.sh ./")
    os.system(f"sbatch run_g16.sh {gjf}")

def run_prepre(gjf, prefix):
    define_gemopt(gjf, prefix)
    define_sp(prefix)
    define_gas(prefix)
    define_sol(prefix)

def run_opt(prefix):
    sub = f"{prefix}.gjf"
    run_gaussian(sub)

def run_prepre_opt(gjf, prefix):
    define_gemopt(gjf, prefix)

def run_sp(prefix):
    sub_group = [file_name for file_name in os.listdir() if f"_{prefix}.gjf" in file_name]
    for ff in sub_group:
        run_gaussian(ff)


def run_mod(gjf, prefix, mod):
    if mod == "prepare":
        run_prepre(gjf, prefix)
    elif mod == "prepare_opt":
        run_prepre_opt(gjf, prefix)
    elif mod == "run_opt":
        run_opt(prefix)
    elif mod == "run_other":
        run_sp(prefix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prepare/run system for pka calculation based on gaussian16')
    parser.add_argument('--gjf', type=str, required=True, \
                        help='initial gjf file that define the geometry')
    parser.add_argument('--prefix', type=str, required=True, \
                        help='system prefix')
    parser.add_argument('--mod', type=str, required=True, \
                        help='executive mode: {prepare} to generate .gjf file for calculation, \
                                              {run_opt} to run initial geometry optimization, \
                                              {run_other} to run subsequent solvation free energy calculation, \
                                              {prepare_opt} to prepare for a single optimization run')

    args = parser.parse_args()
    run_mod(args.gjf, args.prefix,args.mod)


"""
    parser = argparse.ArgumentParser(description='prepare system for pka calculation based on gaussian16')
    parser.add_argument('--gjf', type=str, required=True, \
                        help='initial gjf file that define the geometry')
    parser.add_argument('--prefix', type=str, required=True, \
                        help='system prefix')
    args = parser.parse_args()
    run_prepre(args.gjf, args.prefix)
"""



    ## find charge + mul && xyz
