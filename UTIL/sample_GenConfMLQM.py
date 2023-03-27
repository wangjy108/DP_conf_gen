#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

from rdkit import Chem
import argparse
import logging
import pandas as pd
from util.ConfGenbyMLQM import *

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

def main():
    parser = argparse.ArgumentParser(description='Gen conf from smi with auto3D')
    parser.add_argument('--InputSMI', type=str, required=True, 
                        help='input smile file name')
    parser.add_argument('--k', type=int, default=None, 
                        help='gen and save top k-th conformers for each smile, default None')
    parser.add_argument('--window', type=float, default=20.0, 
                        help="energy window for each smile, default 20.0")
    parser.add_argument('--enumerate_tautomer', type=bool, default=False, 
                        help="if enumerate tautomers for each input, \
                        deault False")
    parser.add_argument('--mpi_np', type=int, default=8, 
                        help="CPU cores for the isomer generation step, defaul 8")
    parser.add_argument('--optimizing_engine', type=str, default="AIMNET", 
                        help="opt enginee(ANI model), default AIMNET, availabel from \
                        [ANI2x, AIMNET])")
    parser.add_argument('--use_gpu', type=bool, default=False, 
                        help="if use gpu, default False")
    parser.add_argument('--threshold', type=float, default=0.5, 
                        help="RMSD cutoff to define SAME conformer, default 0.5")
    

    args = parser.parse_args()

    ConfGen(input_smi_file=args.InputSMI, \
            k=args.k, \
            window=args.window, \
            enumerate_tautomer=args.enumerate_tautomer, \
            mpi_np=args.mpi_np, 
            optimizing_engine=args.optimizing_engine, \
            use_gpu=args.use_gpu, \
            threshold=args.threshold).run()

    return 

if __name__ == '__main__':
    main()