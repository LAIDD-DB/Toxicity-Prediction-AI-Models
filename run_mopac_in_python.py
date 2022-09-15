# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 18:02:12 2022

@author: Gil

If you developed the model with MOPAC...

When you developed the model, it is possible that another researchers want to check the result.
In order to use the model conveniently, it is better to prepare a code to run the model with the given smiles code or the structure.
"""
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import datetime
import openbabel

def qm_parsing(outfile):
    qm_desc = {}
    for each_line in outfile:
        if 'FINAL HEAT OF FORMATION' in each_line:
            #print(each_line.split())
            hof = each_line.split()
            qm_desc['hof']=hof[5]
        elif 'COSMO AREA' in each_line:
            pass
        elif 'COSMO VOLUME' in each_line:
            pass
        elif 'HOMO LUMO ENERGIES' in each_line:
            #print (each_line.split())
            mo = each_line.split()
            qm_desc['homo']=mo[-2]
            qm_desc['lumo']=mo[-1]
        elif 'Mulliken electronegativity' in each_line:
            pass
        elif 'Parr & Pople absolute hardness' in each_line:
            pass
        elif 'Dn(r)' in each_line: #values next to Total: after Dn(r), De(r) q(r)-Z(r)
            pass
        elif 'Average Polarizability from' in each_line: #values under the text
            pass

    return qm_desc


#Below is the input from users.
smi = 'CC(=O)NC1=CC=C(C=C1)O'

rightnow = datetime.datetime.now().strftime('%y%m%d%H%M%S')
wdir = r'directory/to/save/files/'
if not os.path.exists(wdir):
    os.makedirs(wdir)

#From smiles to mol file.
mol = Chem.MolFromSmiles(smi)
mol_H = Chem.AddHs(mol)

AllChem.EmbedMolecule(mol_H)
AllChem.MMFFOptimizeMolecule(mol_H)
print (Chem.MolToMolBlock(mol_H), file=open(os.path.join(wdir,rightnow+'.mol'),'w+'))

#From mol file to xyz file.
mol_to_xyz = openbabel.OBConversion()
mol_to_xyz.SetInAndOutFormats('mol','xyz')
openbabel_mol = openbabel.OBMol()
mol_to_xyz.ReadFile(openbabel_mol, os.path.join(wdir, rightnow+'.mol'))
mol_to_xyz.WriteFile(openbabel_mol, os.path.join(wdir, rightnow+'.xyz'))
mol_to_xyz.CloseOutFile()

#From xyz file to mop file.
xyz_text = open(os.path.join(wdir, rightnow+'.xyz')).readlines()
xyz_text[0] = 'SUPER STATIC'+'\n'
xyz_text[1] = rightnow+'\n\n' #This line addees name of the molecule. Here we use file name.
    
mop_file = open(os.path.join(wdir, rightnow+'.mop'),'w')
for line in xyz_text:
    mop_file.write(line)
mop_file.close()

subprocess.run(['mopac', os.path.join(wdir, rightnow+'.mop')])

#Parsing out file.
out_file = open(os.path.join(wdir, rightnow+'.out')).readlines()
qm_desc = qm_parsing(out_file)
print (qm_desc)

"""
Next step! (For the model developed with the normalized descriptors!)

1) Normalize descriptors based on the model training set as below.

      desc. - min(desc.)
 ---------------------------
   max(desc.) - min(desc.)
   
2) Select features used in the model.
3) Load the model
4) Input the descriptors to the model to get prediction outcome.
"""