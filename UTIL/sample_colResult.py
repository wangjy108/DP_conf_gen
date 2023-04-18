#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

import argparse
import logging
from rdkit import Chem
import os
import sys
import pandas as pd

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

def run(prefix="StrainFinal"):
    available = {}
    should_have = []
    df_col = []

    for a, b, c in os.walk("./"):
        if b:
            #print(f"upper level:{a}")
            #print(f"second level: {b}")
            for each in b:
                for aa, cc, cc in os.walk(os.path.join(a, each)):
                    if [ff for ff in cc if "csv" in ff]:
                        each_csv = [ff for ff in cc if "csv" in ff][0]
                        available.setdefault(each, (os.path.join(os.path.join(a, each), each_csv)))
                    if [ff for ff in cc if "sdf" in ff]:
                        should_have.append(each)
    #print(available)
    

    for _dir, _csv in available.items():
        df = pd.read_csv(_csv, header=0)
        df_col.append(df)

        df_all = pd.concat(df_col)

        df_all.to_csv(f"{prefix}.csv", index=False)

    careful = ""
    for each in should_have:
        if each not in [kk for kk in available.keys()]:
            careful += str(each) + ","
    
    logging.info(f"Final strain result for this dataset saved in {prefix}.csv")
    if careful:
        logging.info(f"Questionable mol as {careful} do not have strain result, check each file for further information")
    
    return 

run()
