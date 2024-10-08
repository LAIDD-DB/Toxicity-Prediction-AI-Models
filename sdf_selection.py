# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:22:14 2022

@author: Gil
"""
import os
import numpy as np
from rdkit import Chem
import numpy as np

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

import io
import pandas as pd
from rdkit import Chem

sdf_files = Chem.SDMolSupplier('BindingDB_herg_ic50.sdf')

mols = [mol for mol in sdf_files]

print (mols[0].GetPropsAsDict())

ic50_dict = {}
smi_dict = {}
check_unconvertibles = []

for i in range(len(mols)):
    # "PubChem CID of Ligand" 속성이 있는지 확인
    if "PubChem CID of Ligand" in mols[i].GetPropsAsDict():
        # 분자가 존재하는지 확인
        if mols[i]:
            # IC50 값과 PubChem CID, SMILES 문자열을 가져옵니다.
            ic50 = mols[i].GetProp("IC50 (nM)")
            cid = mols[i].GetProp("PubChem CID of Ligand")
            smi = Chem.MolToSmiles(mols[i])
            
            # IC50 값이 float으로 변환 가능한지 확인
            if isfloat(ic50):
                # IC50 값을 float으로 변환하고 pIC50 값을 계산
                ic50_data = float(ic50)
                pic50 = -np.log10(ic50_data)
                
                # 이미 존재하는 SMILES 문자열인 경우, IC50 및 pIC50 값을 추가
                if smi in smi_dict:
                    ic50_dict[smi_dict[smi]]['ic50'].append(ic50_data)
                    ic50_dict[smi_dict[smi]]['pic50'].append(pic50)
                else:
                    # 새로운 SMILES 문자열인 경우, 딕셔너리에 추가
                    smi_dict[smi] = i
                    mols[i].SetProp('pIC50 (nM)', str(pic50))
                    ic50_dict[i] = {'ic50': [ic50_data], 'pic50': [pic50], 'cid': cid, 'smi': smi, 'mol': mols[i]}
            else:
                # IC50 값을 float으로 변환할 수 없는 경우 리스트에 추가
                check_unconvertibles.append(ic50)
    else:
        # PubChem CID 속성이 없는 분자는 출력
        print(mols[i], " doesn't have PubChem CID of Ligand")

def check_error(target_list, cutoff):
    mean_v = np.mean(target_list)
    max_v = np.max(target_list)
    min_v = np.min(target_list)
    
    diff_max_mean = abs(max_v - mean_v)
    diff_mean_min = abs(mean_v - min_v)
    if max(diff_max_mean, diff_mean_min) < cutoff:
        return True
    else:
        return False
        
remove_i_list = []
for i, subdict in ic50_dict.items():
    ic50_list = subdict['ic50']
    pic50_list = subdict['pic50']
    if len(pic50_list) > 1:
        if check_error(pic50_list, 0.5):
            ic50_mean = np.mean(ic50_list)
            subdict['ic50']=[ic50_mean]
            subdict['mol'].SetProp('IC50 (nM)',str(ic50_mean))

            pic50_mean = np.mean(pic50_list)
            subdict['pic50']=[pic50_mean]
            subdict['mol'].SetProp('pIC50 (nM)',str(pic50_mean))
        else:
            remove_i_list.append(i)

for i in remove_i_list:
    ic50_dict.pop(i)

sdf_writer = Chem.SDWriter('selected_pic50.sdf')
for i, subdict in ic50_dict.items():
    sdf_writer.write(subdict['mol'])
sdf_writer.close()
