# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:14:56 2022

@author: Gil
"""
import os
import pandas as pd

"""
Error occured since pd.read_csv assumes that csv file is separated by ','.
If delimiter is not ',', it should be specifically inputted to 'sep' in pd.read_csv
"""
#df = pd.read_csv(os.path.join(wdir, csv))

# CSV 파일을 읽어 데이터프레임(df)으로 로드합니다.
# '../CHEMBL240_Toxicity_Data.csv' 파일을 ';' 구분자로 읽어들입니다.
df = pd.read_csv('../CHEMBL240_Toxicity_Data.csv', sep=';')

#In order to check if the data is correctly loaded.
print ('how to use head:', df.head())
print ('how to use tail:', df.tail())

#To check column names of csv files
print ('Check the list of column names in the csv file:', df.columns)
print ('*'*40)

#To check the number of data
print ('Statistics of Data Validity Comment\n:',df['Data Validity Comment'].value_counts())
print ('*'*40)
print ('Statistics of Standard Relation\n',df['Standard Relation'].value_counts())
print ('*'*40)

"""
Error! because in the column the value is '=', not just =.
Therefore matching pattern should be '=', not =
"""
#df = df[df['Standard Relation']=='=']
df = df[df['Standard Relation']=='\'=\'']
df = df[df['Data Validity Comment'].isna()]

print ('Statistics of unit\n',df['Standard Units'].value_counts())
print ('*'*40)
#Separate data whose unit is ug/mL and nM.
df_diffunit = df[df['Standard Units']=='ug.mL-1']
print (df_diffunit.shape)
df_nM = df[df['Standard Units']!='ug.mL-1']
print (df_nM.shape)
print ('*'*40)

"""
Calculation fail since one of values were considered as str.
"""
#converted = 10**(-6)*df_diffunit['Standard Value']/df_diffunit['Molecular Weight']
#print (df_diffunit['Standard Value'].to_list(), type(df_diffunit['Standard Value'].to_list()[0]))
#print (df_diffunit['Molecular Weight'].to_list(), type(df_diffunit['Molecular Weight'].to_list()[0]))

converted = 10**(6)*df_diffunit['Standard Value']/df_diffunit['Molecular Weight'].astype('float')
print ('Maximum IC50 value:',df['Standard Value'].max())
print ('Minimum IC50 value:',df['Standard Value'].min())

df_diffunit['Standard Value']=converted
df_diffunit['Standard Units']='nM'
print (df_diffunit[['Standard Value','Standard Units']])
print ('*'*40)
#Combine two dataframe after conversion.
df = pd.concat([df_nM, df_diffunit])

#Let's check duplicated structures
print ('ChEMBL ID duplicated data:', df.duplicated(subset='Molecule ChEMBL ID',keep=False).sum())
print ('Smiles duplicated data:', df.duplicated(subset='Smiles', keep=False).sum())

#Simple way to identify same structures.
from rdkit import Chem
"""
#Different smiles codes from one identical molecular structure.
smi_list = ['C1=NC2=C(N1)C(=S)N=CN2',
            'n1c[nH]c2c1NC=NC2=S',
            '[H]N1C=NC2=C1C(=S)N=C([H])N2[H]']
for smi in smi_list:
    m = Chem.MolFromSmiles(smi)
    print (Chem.MolToSmiles(m))
"""

for index, row in df.iterrows():
    m = Chem.MolFromSmiles(row['Smiles'])
    df.loc[index, 'rdkit_smi']=Chem.MolToSmiles(m)
print ('standardized smiles code duplicated data:', df.duplicated(subset='rdkit_smi',keep=False).sum())

#Check if duplicated compounds are identical.
print ('Duplicated compounds between ID and Smiles',(df[df.duplicated(subset='Smiles',keep=False)].index==df[df.duplicated(subset='Molecule ChEMBL ID',keep=False)].index).sum())
print ('Duplicated compounds between ID and rdkit_smi',(df[df.duplicated(subset='rdkit_smi',keep=False)].index==df[df.duplicated(subset='Molecule ChEMBL ID',keep=False)].index).sum())
print ('Duplicated compounds between Smiles and rdkit_smi',(df[df.duplicated(subset='rdkit_smi',keep=False)].index==df[df.duplicated(subset='Smiles',keep=False)].index).sum())

import numpy as np
#pIC50 converstion
df['pIC50'] = -np.log10(df['Standard Value'])

#check experimental errors between duplicated compounds.
df_dup = df[df.duplicated(subset='Molecule ChEMBL ID',keep=False)]
df_single = df[~df.duplicated(subset='Molecule ChEMBL ID')]

dup_smi_list = list(set(df_dup['Smiles'].to_list()))
print (df_dup[df_dup['Smiles']==dup_smi_list[0]])
print (df_dup.loc[df_dup['Smiles']==dup_smi_list[0],'Standard Value'])
print (df_dup.loc[df_dup['Smiles']==dup_smi_list[0],'pIC50'])

df_list = []
for smi in dup_smi_list:
    #Take dataframe with identical smiles code.
    temp_df = df_dup[df_dup['Smiles']==smi]
    
    #Calculate mean and 
    mean_pIC50 = temp_df['pIC50'].mean()
    max_pIC50 = temp_df['pIC50'].max()
    min_pIC50 = temp_df['pIC50'].min()
    
    dif_mean_max = abs(max_pIC50 - mean_pIC50)
    dif_mean_min = abs(mean_pIC50 - min_pIC50)
    if max(dif_mean_max, dif_mean_min) < 0.5:
        temp_df['Standard Value']=mean_pIC50
        #Collect dataframe.
        df_list.append(temp_df.iloc[0])

#Concatenate dataframes in the list.        
df_dup_final = pd.concat(df_list, axis=1).T

df_final = pd.concat([df_dup_final, df_single])
#Final data
print ('Final data size:', df_final.shape)
df_final.to_csv(os.path.join(wdir, 'hERG_preprocessed.csv'),index=False)
