#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

from rdkit import Chem
from util.CalRMSD import *
import argparse
import logging
import pandas as pd

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description='calc rmsd, maybe should do align first')
    parser.add_argument('--InputSDFmol', type=str, required=True, 
                        help='input search sdf file N >= 1')
    parser.add_argument('--InputSDFref', type=str, required=True, 
                        help='input reference sdf file N == 1')
    parser.add_argument('--saveValue', type=bool, default=False, 
                        help="if save calculated rmsd in file, default False, \
                        if set True, calc RMSD value will be saved in RMSD.dat")
    parser.add_argument('--rmsd_method', type=str, default="selfWhole", 
                        help="rmsd method, default is selfWhole, could be \
                        [crippen3D, selfWhole]")
    args = parser.parse_args()

    search_mol = [mm for mm in Chem.SDMolSupplier(args.InputSDFmol, removeHs = False) if mm]
    ref_mol = [mm for mm in Chem.SDMolSupplier(args.InputSDFref, removeHs = False) if mm][0]

    rmsd_record = {}
    for ii in range(len(search_mol)):
        try:
            get_rmsd = RMSD(rdmolobj_mol=search_mol[ii], rdmolobj_ref=ref_mol, method=args.rmsd_method).run()
        except Exception as e:
            logging.info(f"Failed with {ii}th mol in {args.InputSDFmol}, mol name is {search_mol[ii].GetProp('_Name')}")
        else:
            rmsd_record.setdefault(ii, get_rmsd)
    
    df = pd.DataFrame({"mol_idx": [cc for cc in rmsd_record.keys()], \
                           "calc rmsd": [rr for rr in rmsd_record.values()]})
            
    if args.saveValue:
        df.to_csv("RMSD.dat", index=None)
    else:
        logging.info(f"calculated rmsd according show: {[rr for rr in rmsd_record.values()]}")

    return 
    

if __name__ == '__main__':
    main()
    




