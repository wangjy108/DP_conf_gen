#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

from util.Cluster import *
import argparse
import logging

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description='align, cluster and save')
    parser.add_argument('--InputSDFmol', type=str, required=True, 
                        help='input sdf file')
    parser.add_argument('--save_N', type=int, required=True, 
                        help='save n conformer for each mol')
    parser.add_argument('--method', type=str, default="RMSD", 
                        help="rdkit property name tag, should be numeric energy related one")
    parser.add_argument('--name_tag', type=str, default="_Name", 
                        help="rdkit property name tag used to collect conformer for each mol, \
                        if input rdkit file do not have _Name tag, this option should be assigned\
                        with availble property tag")
    parser.add_argument("--align_method", type=str, default="crippen3D", 
                        help="applied align method, default crippen3D\
                        could be from [crippen3D, LSalignFlex], latter should have extra application installed")
    parser.add_argument('--rmsd_method', type=str, default="selfWhole", 
                        help="rmsd method, default is selfWhole, could be \
                        [crippen3D, selfWhole]")
    parser.add_argument('--reference', type=str, default=None, 
                        help="reference file name, if not define, \
                        will set the 1st one from conformer set as reference")
    args = parser.parse_args()

    cluster(inputSDF_fileName=args.InputSDFmol, \
            save_n=args.save_N, \
            method=args.method, \
            name_tag=args.name_tag, \
            align_method=args.align_method, \
            rmsd_method=args.rmsd_method, \
            reference=args.reference).run()

    return 
    

if __name__ == '__main__':
    main()
    




