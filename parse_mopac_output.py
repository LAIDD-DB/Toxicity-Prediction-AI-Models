# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:09:05 2022

@author: Gil
"""

import os
import pandas as pd

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
    
wdir = '/content/Toxicity-Prediction-AI-Models/'
csv = 'dilirank_mw_lessthan_1000.csv'
df = pd.read_csv(os.path.join(wdir, csv))

for index, row in df.iterrows():
    out_dir = os.path.join(wdir, 'mop', str(index)+'.out')
    if os.path.exists(out_dir):
        out_file = open(out_dir).readlines()
        qm_desc = qm_parsing(out_file)
        for key, value in qm_desc.items():
            df.loc[index, key]=value

df.to_csv(os.path.join(wdir,'dilirank_qm.csv'))
