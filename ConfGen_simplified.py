#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

from rdkit import rdBase, Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import numpy as np
import os
import argparse
import copy
import scipy.spatial

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


    crippen3D = rdMolAlign.GetCrippenO3A(rdmolobj_mol, rdmolobj_ref,crippen_contrib_mol, crippen_contrib_ref)
    crippen3D.Align()
    #map = crippen3D.Matches()
    #rmsd = rdMolAlign.AlignMol(rdmolobj_mol, rdmolobj_ref, atomMap=map)
    rmsd = rmsd_self_whole(rdmolobj_mol, rdmolobj_ref)

    return rmsd, rdmolobj_mol


def gen_conf_mm(inputSmi:list, genConfNum:int, saveConfNum:int, **arg):
    """
    gen3D mol by rdkit
    optimize by mm force filed, default mmff94, return with conformation set obj
    =====
    input: inputfile name;
    genConfNum: number of conformation generate for each smi;
    saveConfNum: number of conformation to save;
    args: 
        //fileName: saved sdf file name
    """

    # GetConformer(id=???)

    #wrong_smi = []

    try:
        save_file_name = os.path.join(os.getcwd(), arg["fileName"])
    except Exception as e:
        save_file_name = os.path.join(os.getcwd(), "SAVE.sdf")

    cc = Chem.SDWriter(save_file_name)

    for i in range(len(inputSmi)):
        mol = Chem.MolFromSmiles(inputSmi[i])
        #mol = input_mol[i]
        m3d = Chem.AddHs(mol)

        try:
            nGenConfs = AllChem.EmbedMultipleConfs(m3d,numConfs=genConfNum, numThreads=0)
        except Exception as e:
            #wrong_smi.append(inputSmi[i])
            continue
        else:
            if len(nGenConfs) == 0:
                #wrong_smi.append(inputSmi[i])
                continue

        res = AllChem.MMFFOptimizeMoleculeConfs(m3d, numThreads=0)

        if len(res) == 0:
            #wrong_smi.append(inputSmi[i])
            continue
        elif min([rr[0] for rr in res]) < 0:
            m3d.SetProp("cSMILES", inputSmi[i])
            m3d.SetProp("conf_idx", "0")
            m3d.SetProp("MM_energy", "0.0")
            cc.write(m3d, confId=0)
            continue

        if len([(i, res[i][-1]) for i in range(len(res)) if res[i][0] == 0]) < saveConfNum:
            stable_conf = [(i, res[i][-1]) for i in range(len(res)) if res[i][0] == 0] + \
                           sorted([(i, res[i][-1]) for i in range(len(res)) if res[i][0] == 1], \
                           key=lambda x: x[-1], reverse=False)[:saveConfNum-len([(i, res[i][-1]) \
                           for i in range(len(res)) if res[i][0] == 0])]
        else:
            stable_conf = sorted([(i, res[i][-1]) for i in range(len(res))],\
                          key=lambda x: x[-1], reverse=False)[:saveConfNum]


        for ii in stable_conf:
            m3d.SetProp("cSMILES", inputSmi[i])
            m3d.SetProp("conf_idx", str(ii[0]))
            m3d.SetProp("MM_energy", f"{ii[1]:.7f}")
            cc.write(m3d, confId=ii[0])

    cc.close()

    try:
        Chem.SDMolSupplier(save_file_name)
    except Exception as e:
        os.system(f"rm -f {save_file_name}")
        return inputSmi

    sdf = set([mm.GetProp("cSMILES") for mm in Chem.SDMolSupplier(save_file_name) if mm])

    wrong_smi = set([ii for ii in inputSmi if ii not in sdf])

    return list(wrong_smi)


def save_by_cluster(molsdf_set:str, ref_single, n_cluster:int, **arg):
    """
    args:
        //proTag: property tag to differentiate conformer
        //eneTag: energy property name in sdf file to be used in cluster section, 
                  defult to MM_energy
    """

    try:
        pro_tag = arg["proTag"]
    except Exception as e:
        pro_tag = "_Name"

    try:
        ene_tag = arg["eneTag"]
    except Exception as e:
        ene_tag = "MM_energy"

    #ref = [mm for mm in Chem.SDMolSupplier(ref_single, removeHs=False)][0]
    ref = ref_single ## rdmolobj
    mol_set = [mm for mm in Chem.SDMolSupplier(molsdf_set) if mm != None]
    ## extrac file Name
    mol_name = molsdf_set.strip().split("/")[-1]

    avail_prop = [nn for nn in mol_set[0].GetPropNames()]

    try:
        mol_set[0].GetProp(pro_tag)
    except Exception as e:
        #logging.info(f"Choose name identifier from properName: {avail_prop}")
        return

    try:
        mol_set[0].GetProp(ene_tag)
    except Exception as e:
        #logging.info(f"Choose energy identifier from properName: {avail_prop}")
        return 

    mol_collect = {}
    for i in range(len(mol_set)):
        mm_name = mol_set[i].GetProp(pro_tag)
        if not mm_name in mol_collect.keys():
            mol_collect.setdefault(mm_name, [])
        mol_collect[mm_name].append(mol_set[i])

    save = []

    for k, v in mol_collect.items():
        #logging.info(f"Working with {k}")
        if len(v) > n_cluster and n_cluster > 1:
            rmsd_mol = sorted([align_by_3Dcrippen(v[i], ref) for i in range(len(v))], key=lambda x:x[0])
            rmsd = [cc[0] for cc in rmsd_mol]
            #sorted(rmsd_mol, key=lambda x:x[0])
            bins = np.linspace(min(rmsd), max(rmsd)+10e-6, n_cluster+1)
            ## cut1 = bins[1]
            ## cut2 = bins[2]
            #collect ={0:[rmsd_mol[0][1]], 1:[], 2:[rmsd_mol[-1][1]]}
            collect = {}
            for i in range(n_cluster):
                collect.setdefault(i, [])

            collect[0].append(rmsd_mol[0][1])
            #collect[n_cluster-1].append(rmsd_mol[-1][1])

            i = 1
            j = 1
            while i <  len(v):
                item = rmsd_mol[i]
                if item[0] < bins[j]:
                    collect[j-1].append(item[1])
                else:
                    j += 1
                    collect[j-1].append(item[1])
                i += 1

            #save += [cc[0] for cc in collect.values() if len(cc) >=1]
            #count = [len(pp) for pp in collect.values()]
            #sel_from = list(collect.values())
            col_img = copy.deepcopy(collect)
            for kk, vv in col_img.items():
                if len(vv) == 0:
                    del collect[kk]

            re_collect = sorted(list(collect.values()), key=lambda x:len(x), reverse=True)

            for cc in re_collect:
                cc.sort(key=lambda x:float(x.GetProp(ene_tag)), reverse=True)

            current_save = []

            while True:
                sel_from_which = n_cluster - len(current_save)
                if sel_from_which > 0:
                    for ccc in re_collect[:sel_from_which]:
                        if ccc:
                            current_save.append(ccc.pop())
                else:
                    break
            save += current_save
            #logging.info(f"collect as {[cc.GetProp('E_tot') for cc in current_save]}")

        else:
            rmsd_mol = sorted([align_by_3Dcrippen(v[i], ref) for i in range(len(v))], key=lambda x:x[0])
            #save += [align_by_3Dcrippen(v[i], ref)[1] for i in range(len(v))]
            save += [cc[-1] for cc in rmsd_mol][:n_cluster]
            #print(rmsd_mol)
                

    cc = Chem.SDWriter(f"FILTERED_{mol_name}")
    for ii in save:
        cc.write(ii)

    cc.close()

    #logging.info(f"Save in need molecules in FILTERED_{mol_name}")
    return


def run(InputSmiFile:str, Ref3DFile:str, ConfGenNum:int, SaveConfNumEach:int):
    ##read-in reference
    RefFileType = Ref3DFile.split(".")[-1]
    if RefFileType == 'sdf':
        refMolObj = [mm for mm in Chem.SDMolSupplier(Ref3DFile)][0]
    elif RefFileType == 'mol2':
        refMolObj = Chem.MolFromMol2File(Ref3DFile)
    else:
        raise Exception("Wrong input mol type, use ['.sdf', '.mol2'] type as input")
        return 
    
    ## get smiles info
    get_df = pd.read_csv(InputSmiFile, header=None, sep="\\s+")
    searchMolSmi = get_df.iloc[:,0].to_list()
    nameTag = get_df.iloc[:,-1].to_list()

    pariwise_smile_name = {}

    for ii in range(len(searchMolSmi)):
        pariwise_smile_name.setdefault(searchMolSmi[ii], nameTag[ii])
    
    ## gen conf from smiles
    ## default method: MMFF94, default mataching property: _cSMILES
    LeftSmi = gen_conf_mm(searchMolSmi, genConfNum=ConfGenNum, saveConfNum=20)
    ## this will generate a file named SAVE.sdf with all generated 3D conformer
    
    ## align, ranking and save
    save_by_cluster(molsdf_set="SAVE.sdf", ref_single = refMolObj, n_cluster = SaveConfNumEach, proTag='cSMILES')
    ## all available 3D conformations (1 for each simle) are saved in FILTERED_SAVE.sdf

    ## extrac setting:
    before_finale = [mm for mm in Chem.SDMolSupplier("FILTER_SAVE.sdf", removeHs=False)]

    finale = Chem.SDWriter("FINALE.sdf")
    for each in before_finale:
        getSMI = each.GetProp("cSMILES")
        _name = pariwise_smile_name[getSMI]
        each.SetProp("_Name", _name)
        finale.write(each)
    finale.close()

    os.system("rm -f SAVE.sdf FILTER_SAVE.sdf")
    print("FINALE.sdf save 3D conformation information")
    
    
    if len(LeftSmi) != 0:
        df = pd.DataFrame({"":LeftSmi})
        df.to_csv("Qestionable.smi", index=None, header=None)
        print("SMILES in question are save in [Questionable.smi] upon check")
        return 
    else:
        return 
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='generate 3D conformation based on MMFF94(rdkit), align by 3Dcrippen,\
                                                  and save according to rmsd')
    parser.add_argument('--InputSmiFile', type=str, required=True, help='input smile file name, single column, no header')
    parser.add_argument('--Ref3DFile', type=str, required=True, help='input 3D reference file name, sdf should be better')
    parser.add_argument('--ConfGenNum', type=int, default=20, help='initial conformer generation number from rdkit MMFF94, default==20')
    parser.add_argument('--SaveConfNumEach', type=int, default=1, help='finale saved conformation number, default==1')

    args = parser.parse_args()

    run(args.InputSmiFile, args.Ref3DFile, args.ConfGenNum, args.SaveConfNumEach)

    #run("testSMILES.txt", "R003007.sdf", 20, 1)




