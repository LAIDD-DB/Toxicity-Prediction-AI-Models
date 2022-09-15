# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:03:28 2022

@author: Gil
"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

wdir = r'directory/of/the/csv/file'
csv = 'dilirank_curate_pcp.csv'
df = pd.read_csv(os.path.join(wdir, csv))

#If molecule size is too big, it fails to optimize 3D structures.
df = df[df['mw']<1000]
df = df.reset_index(drop=True)
df.to_csv(os.path.join(wdir, 'dilirank_mw_lessthan_1000.csv'),index=False)
"""
First perform 3D structure optimization.
"""

#If directory is not created, make one.
#Directory for 3D structures.
dir_3d = '3d'
if not os.path.exists(os.path.join(wdir, dir_3d)):
    os.mkdir(os.path.join(wdir, dir_3d))
print ('Start 3D structure optimization.')    
for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['smi'])
    mol_H = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_H)
    AllChem.MMFFOptimizeMolecule(mol_H)
    print (Chem.MolToMolBlock(mol_H), file=open(os.path.join(wdir,dir_3d,str(index)+'.mol'),'w+'))
#index was used as a file name in order to correctly and easily assign qm descriptors to the source molecules.

"""
Convert mol file to xyz file using openbabel.
XYZ file only contains xyz coordinates of each atom.
xyz file is quite similar to mop file except that it doesn't have calculation keywords.
Therefore, it is easy to prepare mop file by inserting the keywords at the first line of xyz file.
"""
    
import openbabel

#Directory for xyz files.
dir_xyz = 'xyz'
if not os.path.exists(os.path.join(wdir, dir_xyz)):
    os.mkdir(os.path.join(wdir, dir_xyz))

mol_files = os.listdir(os.path.join(wdir, dir_3d))
print ('Converting mol files to xyz files.')
mol_to_xyz = openbabel.OBConversion()
mol_to_xyz.SetInAndOutFormats('mol','xyz')

mol = openbabel.OBMol()
for m in mol_files:
    mol_to_xyz.ReadFile(mol, os.path.join(wdir, dir_3d, m))
    mol_to_xyz.WriteFile(mol, os.path.join(wdir, dir_xyz, m.replace('mol','xyz')))
mol_to_xyz.CloseOutFile()

"""
Generate mop file from xyz file
"""
dir_mop = 'mop'
if not os.path.exists(os.path.join(wdir, dir_mop)):
    os.mkdir(os.path.join(wdir, dir_mop))
print ('Generating xyz files to mop files.')    
mopac_keywords = 'SUPER STATIC'+'\n'
xyz_files = os.listdir(os.path.join(wdir, dir_xyz))
for xyz in xyz_files:
    xyz_text = open(os.path.join(wdir, dir_xyz, xyz)).readlines()
    xyz_text[0] = mopac_keywords
    xyz_text[1] = xyz.replace('.xyz','')+'\n\n' #This line addees name of the molecule. Here we use file name.
    
    mop = xyz.replace('xyz','mop')
    mop_file = open(os.path.join(wdir, dir_mop, mop),'w')
    for line in xyz_text:
        mop_file.write(line)
    mop_file.close()

"""
Make batch file for multiple mop files.
"""
print ('Make batch file to calculate multiple molecules in MOPAC')
mop_files = os.listdir(os.path.join(wdir, dir_mop))
mop_batch = open(os.path.join(wdir, dir_mop, 'mop_batch.bat'),'w')

#call activate 'environment where you installed your mopac'.
#call activate 'base' is used since mopac is installed in base environment in my personal anaconda setting.
#Please check which anaconda environment you installed mopac.

mop_batch.write('call activate base\n') 
for m in mop_files:
    if m.endswith('mop'):
        mop_batch.write('mopac '+m+'\n')
mop_batch.close()
