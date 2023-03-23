#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

from rdkit import Chem
from util.Align import *
import argparse
import logging

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description='single align')
    parser.add_argument('--InputSDFmol', type=str, required=True, 
                        help='input search sdf file')
    parser.add_argument('--InputSDFref', type=str, required=True, 
                        help='input reference sdf file')
    parser.add_argument('--saveFilename', type=str, default="ALIGNED.sdf", 
                        help="saved aligned sdf file name, default is ALIGNED.sdf")
    parser.add_argument('--align_method', type=str, default="crippen3D", 
                        help="align method, default is crippen3D, could be \
                        [crippen3D, LSalignFlex], latter should have extra application installed")
    args = parser.parse_args()

    search_mol = [mm for mm in Chem.SDMolSupplier(args.InputSDFmol, removeHs = False) if mm]
    ref_mol = [mm for mm in Chem.SDMolSupplier(args.InputSDFref, removeHs = False) if mm][0]

    cc = Chem.SDWriter(args.saveFilename)
    for ii in range(len(search_mol)):
        #args_set = {"SearchMolObj": search_mol[ii], \
        #            "RefMolObj": ref_mol, \
        #            "method": args.align_method}
        try:
            aligned_each = Align(SearchMolObj=search_mol[ii], RefMolObj=ref_mol, method=args.align_method).run()
        except Exception as e:
            logging.info(f"Failed with {ii}th mol in {args.InputSDFmol}, mol name is {search_mol[ii].GetProp('_Name')}")
        else:
            if aligned_each:
                cc.write(aligned_each)
            else:
                logging.info(f"Failed with {ii}th mol in {args.InputSDFmol}, mol name is {search_mol[ii].GetProp('_Name')}")
    cc.close()

    logging.info(f"Save aligned mol in {args.saveFilename}")

    return 
    

if __name__ == '__main__':
    main()
    




