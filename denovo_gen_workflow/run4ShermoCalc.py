import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
import numpy as np
import scipy.spatial
import argparse


def rmsd_self_whole(rdmolobj_mol, rdmolobj_ref):
    mol_xyz = rdmolobj_mol.GetConformer().GetPositions()
    ref_xyz = rdmolobj_ref.GetConformer().GetPositions()
    
    mol_atomIdx = np.array([a.GetAtomicNum() for a in rdmolobj_mol.GetAtoms()])
    ref_atomIdx = np.array([a.GetAtomicNum() for a in rdmolobj_ref.GetAtoms()])
    
    if mol_xyz.shape[0] > ref_xyz.shape[0]:
        search = ref_xyz
        search_atomIdx = ref_atomIdx
        target = mol_xyz
        target_atomIdx = mol_atomIdx
    else:
        search = mol_xyz
        search_atomIdx = mol_atomIdx
        target = ref_xyz
        target_atomIdx = ref_atomIdx
    
    match_target = []
    match_atomIdx = []
    
    for i in range(search.shape[0]):
        this_xyz = search[i]
        dis = scipy.spatial.distance.cdist(this_xyz.reshape(1, -1), target)
        minP = np.argmin(dis)
        match_target.append(target[minP])
        match_atomIdx.append(target_atomIdx[minP])
        np.delete(target, minP, axis=0)
        np.delete(target_atomIdx, minP)
    
    np_match_target = np.array(match_target).reshape(search.shape[0], 3)
    np_match_atomIdx = np.array(match_atomIdx)
    ## calc naive rmsd
    #naive_rmsd = np.power(sum((np.power(np.sum((search - np_match_target)**2, axis=1), 0.5))**2)/search.shape[0], 0.5)
    naive_rmsd = np.power(sum((np.power(np.sum((search - np_match_target)**2, axis=1), 0.5) \
                        * np.exp(np_match_atomIdx - search_atomIdx))**2)/search.shape[0], 0.5)
    
    return naive_rmsd


def align_by_3Dcrippen(rdmolobj_mol, rdmolobj_ref):
    crippen_contrib_mol = rdMolDescriptors._CalcCrippenContribs(rdmolobj_mol)
    crippen_contrib_ref = rdMolDescriptors._CalcCrippenContribs(rdmolobj_ref)


    crippen3D = rdMolAlign.GetCrippenO3A(rdmolobj_mol, rdmolobj_ref,                                          crippen_contrib_mol, crippen_contrib_ref)
    crippen3D.Align()
    #map = crippen3D.Matches()
    rmsd = rmsd_self_whole(rdmolobj_mol, rdmolobj_ref)

    return rmsd, rdmolobj_mol
    #return rdmolobj_mol

def run():
    ## step1 assemble results from g16 calc
    #mol = [ff for ff in os.listdir("./")]
    mol = []
    for each in os.listdir("./"):
        for a, b, c in os.walk(each):
            if len(c) and "TRACK_JOB" in c:
                mol.append(each)
    
    real_name = set([nn.split("_")[0] for nn in mol])

    for nnmm in real_name:
        os.system(f"mkdir {nnmm}")
        os.system(f"mv {nnmm}_*/*.log {nnmm}")
        os.system(f"rm -rf {nnmm}_*")
        ## processing shermo

        log_file_name = [f"{nnmm}/{fff}" for fff in os.listdir(nnmm) if '.log' in fff]
        energy_tag = []
        for ll in log_file_name:
            with open(ll, 'r+') as muchff:
                energy_tag.append([tt.split()[4] for tt in muchff if 'E(RB2PLYPD3)' in tt][0])
        
        df = pd.DataFrame({"name":log_file_name, "ene": energy_tag})
        df.to_csv("list.txt", header=None, index=None, sep=';')
        
        os.system(f"Shermo list.txt | grep 'Boltzmann weight' | awk '{{print $9}}' > {nnmm}/Dis.dat")
        collect = []

        df1 = pd.read_csv(f"list.txt", sep=';', header=None)
        df2 = pd.read_csv(f"{nnmm}/Dis.dat", header=None)
        df = pd.concat([df1, df2], axis=1)
        df.columns = ["name", "energy", "partition"]
        df.sort_values(by="partition", ascending=False, inplace=True)
        for name in df["name"].to_list():
            os.system(f"obabel -ig16 {name} -O {nnmm}/{name.split('/')[-1].split('.')[0]}.sdf")
            collect.append([mm for mm in Chem.SDMolSupplier(f"{nnmm}/{name.split('/')[-1].split('.')[0]}.sdf") if mm][0])
            os.system(f"rm -f {nnmm}/*.sdf")
        
        partition_label = df["partition"].to_list()
        name_label = df["name"].to_list()

        os.system(f"mv list.txt {nnmm}/")

        for ii in range(len(collect)):
            rmsd, aligned_mol  = align_by_3Dcrippen(collect[ii], collect[0])
            save_name = f"{name_label[ii].split('/')[-1].split('.')[0].split('_')[0]}_{partition_label[ii]}_{rmsd}.sdf"
            cc = Chem.SDWriter(f"{nnmm}/{save_name}")
            cc.write(aligned_mol)
            cc.close()
    
    return 


if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description='processing final results')
    #parser.add_argument('--sdf_tag', type=str, required=True, help='input serial sdf tag')
    #parser.add_argument('--input_smi', type=str, required=True, help='input smi file name')
    #args = parser.parse_args()

    run()