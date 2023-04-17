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
from util.ConfRelaxbySQM import System as MDgen

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

def main():
    parser = argparse.ArgumentParser(description='sampling conf with xtb MD')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file name')
    parser.add_argument('--run_time', type=int, default=50, 
                        help='md run time, int, default == 50')
    parser.add_argument('--run_temperature', type=int, default=400, 
                        help='md run temperature, int, default 400')
    parser.add_argument('--save_frame', type=int, default=100, 
                        help='save N frame, int, default 100')
    
    args = parser.parse_args()
    MDgen(input_sdf=args.input_sdf, 
          run_time=args.run_time,
          run_temperature=args.run_temperature, 
          save_frame=args.save_frame).run()
    
    return 


if __name__ == '__main__':
    main()



    