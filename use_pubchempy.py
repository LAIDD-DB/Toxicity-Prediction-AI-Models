# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:11:46 2022

@author: Gil

Document for pubchempy
https://pubchempy.readthedocs.io/en/latest/

Pubchem fingerprint (how to convert fingerprint to binary string)
https://github.com/mcs07/PubChemPy/blob/master/examples/Chemical%20fingerprints%20and%20similarity.ipynb
"""
#you can search compounds by name or cas number.
#pcp.get_compounds('!input your cas number here!', 'name')
#result = pcp.get_compounds(df.loc[0,'Compound Name'], 'name')

#downloaded data can be checked in dictionary format. By checking all the details in the data, we can decide what to be stored in the dataframe.
#pcp_dict = result[0].to_dict()
#전체 df에서 수행 시 서버 문제(PUGREST.ServerBusy)가 발생하여, 이에따라 100개의 데이터에서만 테스트 수행하였으며, 필요 시 df_test -> df로 변경하여 수행
import pubchempy as pcp  # PubChem 데이터를 가져오기 위한 라이브러리
import pandas as pd  # 데이터 조작 및 분석을 위한 라이브러리

# 엑셀 파일을 읽어와 데이터프레임으로 변환
df = pd.read_excel('dilirank_curate.xlsx')

# 특정 'vDILIConcern' 값들만 필터링하여 데이터프레임을 업데이트
df = df[df['vDILIConcern'].isin(['vMost-DILI-Concern', 'vNo-DILI-Concern'])]

# 인덱스를 초기화하여 데이터프레임을 업데이트
df = df.reset_index(drop=True)

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
df.to_csv('dilirank_curate_pcp.csv', index=False)

fp_arr = []
cactvs_fp_arr = []
for index, row in df.iterrows():
    try:
        # 'fp' 값이 'NaN'이 아닌 경우에만 처리
        if pd.notna(row['fp']):
            # 16진수로 변환하여 리스트에 추가
            fp_arr.append(list(bin(int(row['fp'], 16))[2:]))
        if pd.notna(row['cactvs_fp']):
            cactvs_fp_arr.append(list(row['cactvs_fp']))
    except ValueError as e:
        print(f"Error converting fp at index {index}: {row['fp']}")

df_fp = pd.DataFrame(fp_arr)
df_fp = pd.concat([df[['Compound Name','vDILIConcern']],df_fp],axis=1)
df_fp.to_csv('dilirank_fp.csv',index=False)

df_cactvs = pd.DataFrame(cactvs_fp_arr)
df_cactvs = pd.concat([df[['Compound Name','vDILIConcern']],df_cactvs],axis=1)
df_cactvs.to_csv('dilirank_cactvs.csv',index=False)
