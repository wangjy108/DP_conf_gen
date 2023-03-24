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
from util.ConfGenbyMM import *

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description='Gen conf from smi with MM method')
    parser.add_argument('--InputSMI', type=str, required=True, 
                        help='input smile file name')
    parser.add_argument('--method', type=str, default="MMFF94", 
                        help='gen method, default is MMFF94')
    parser.add_argument('--genConfNum', type=int, default=20, 
                        help="gen conformer N for each input smi, default 20")
    parser.add_argument('--saveConfNum', type=int, default=20, 
                        help="save conformer N from generated, default 20, \
                        should be <= genConfNum")
    parser.add_argument('--saveFileName', type=str, default="SAVE.sdf", 
                        help="saved sdf file name, default is SAVE.sdf")

    args = parser.parse_args()

    ConfGen(input_smi_file=args.InputSMI, \
            method=args.method, \
            genConfNum=args.genConfNum, \
            saveConfNum=args.saveConfNum, \
            fileName=args.saveFileName).run()

    return 
    

if __name__ == '__main__':
    main()
    




