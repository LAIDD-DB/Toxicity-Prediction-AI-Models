# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:11:46 2022

@author: Gil

Document for pubchempy
https://pubchempy.readthedocs.io/en/latest/

Pubchem fingerprint (how to convert fingerprint to binary string)
https://github.com/mcs07/PubChemPy/blob/master/examples/Chemical%20fingerprints%20and%20similarity.ipynb
"""

import os
import pandas as pd
import pubchempy as pcp
#Compound name-based search
wdir = r'your/directory'
excel = 'dilirank_curate.xlsx'

df = pd.read_excel(os.path.join(wdir, excel))
df = df[df['vDILIConcern'].isin(['vMost-DILI-Concern','vNo-DILI-Concern'])]
df = df.reset_index(drop=True)
#you can search compounds by name or cas number.
#pcp.get_compounds('!input your cas number here!', 'name')
#result = pcp.get_compounds(df.loc[0,'Compound Name'], 'name')

#downloaded data can be checked in dictionary format. By checking all the details in the data, we can decide what to be stored in the dataframe.
#pcp_dict = result[0].to_dict()

for index, row in df.iterrows():
    if index%100 == 0:
        print (index, 'downloading...')
    result = pcp.get_compounds(row['Compound Name'], 'name')
    if result:
        df.loc[index, 'smi'] = result[0].canonical_smiles
        df.loc[index, 'h_bond_acc'] = result[0].h_bond_acceptor_count
        df.loc[index, 'h_bond_don'] = result[0].h_bond_donor_count
        df.loc[index, 'heavy'] = result[0].heavy_atom_count
        df.loc[index, 'rotatable_bond'] = result[0].rotatable_bond_count
        df.loc[index, 'tpsa'] = result[0].tpsa  #TPSA: Topological polar surface area
        df.loc[index, 'mw'] = result[0].molecular_weight
        df.loc[index, 'xlogp'] = result[0].xlogp
        df.loc[index, 'fp'] = result[0].fingerprint
        df.loc[index, 'cactvs_fp'] = result[0].cactvs_fingerprint

df = df[~df['smi'].isna()]
df.to_csv(os.path.join(wdir, excel.replace('.xlsx','_pcp.csv')), index=False)

#Convert fingerprint to dataframe and save it for model development.
fp_arr = []
cactvs_fp_arr = []
for index, row in df.iterrows():
    fp_arr.append(list(bin(int(row['fp'], 16))[2:]))
    cactvs_fp_arr.append(list(row['cactvs_fp']))

df_fp = pd.DataFrame(fp_arr)
df_fp = pd.concat([df[['Compound Name','vDILIConcern']],df_fp],axis=1)
df_fp.to_csv(os.path.join(wdir, 'dilirank_fp.csv'),index=False)

df_cactvs = pd.DataFrame(cactvs_fp_arr)
df_cactvs = pd.concat([df[['Compound Name','vDILIConcern']],df_cactvs],axis=1)
df_cactvs.to_csv(os.path.join(wdir, 'dilirank_cactvs.csv'),index=False)