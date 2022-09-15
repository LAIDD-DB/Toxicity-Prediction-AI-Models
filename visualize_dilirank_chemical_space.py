# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:56:55 2022

@author: Gil
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

wdir = r'directory_of_input_file'
csv = 'dilirank_curate_pcp.csv'

df = pd.read_csv(os.path.join(wdir, csv))

"""
Basic for a scatter plot (chemical space visualization)
"""

plt.scatter(df['mw'],df['xlogp'])
plt.title('Chemical space of DILIrank')
plt.xlabel('Molecular weight')
plt.ylabel('logP')
plt.savefig(os.path.join(wdir, 'dilirank_chemspace.png'))
plt.close()

"""
Chemical space comparison between most and no DILI data.
But let's visualize drugs whose mw is less than 1000.
"""
df = df[df['mw']<1000]
df_most = df[df['vDILIConcern']=='vMost-DILI-Concern']
df_no = df[df['vDILIConcern']=='vNo-DILI-Concern']

plt.scatter(df_most['mw'],df_most['xlogp'],label='Most-DILI')
plt.scatter(df_no['mw'],df_no['xlogp'],label='No-DILI')
plt.title('Chemical space b.w. Most and No DILI')
plt.xlabel('Molecular weight')
plt.ylabel('logP')
plt.legend(loc=0)
plt.savefig(os.path.join(wdir, 'dili_most_no_chemspace.png'))
plt.close()

"""
Let's use bar plot

Normalization is essential in order to match scale of each descriptor.

Example)
Range of MW is from tens to hundreds
but scale of logP is one-digit number.
Thus, they can't be fairly compared.
"""

norm_desc = df[['h_bond_acc','h_bond_don','heavy','rotatable_bond','tpsa','mw','xlogp']]
norm_desc = (norm_desc - norm_desc.min())/(norm_desc.max()-norm_desc.min())

x = ['H-bond ACC', 'H-bond DON', 'HEAVY AT.', 'Rot. Bond', 'TPSA', 'MW', 'XlogP']
descriptors = [norm_desc['h_bond_acc'].mean(),
               norm_desc['h_bond_don'].mean(),
               norm_desc['heavy'].mean(),
               norm_desc['rotatable_bond'].mean(),
               norm_desc['tpsa'].mean(),
               norm_desc['mw'].mean(),
               norm_desc['xlogp'].mean()]
plt.bar(range(len(descriptors)), descriptors)
plt.title('Descriptor mean values in dilirank')
plt.xlabel('Descriptors')
plt.ylabel('Mean values')
plt.xticks(range(len(x)), x, rotation=-45, ha='left', rotation_mode='anchor')
plt.savefig(os.path.join(wdir, 'dili_descriptors.png'),bbox_inches='tight')
plt.close()

"""
It's useful to compare descriptor distribution between labels or training and test.
"""
import numpy as np 
x = ['H-bond ACC', 'H-bond DON', 'HEAVY AT.', 'Rot. Bond', 'TPSA', 'MW', 'XlogP']
x_loc = np.arange(len(x))

#Index of data is used to subset normalized descriptors for df_most and df_no.
norm_most = norm_desc.loc[df_most.index]
norm_no = norm_desc.loc[df_no.index]

desc_most = [norm_most['h_bond_acc'].mean(),
             norm_most['h_bond_don'].mean(),
             norm_most['heavy'].mean(),
             norm_most['rotatable_bond'].mean(),
             norm_most['tpsa'].mean(),
             norm_most['mw'].mean(),
             norm_most['xlogp'].mean()]
desc_no = [norm_no['h_bond_acc'].mean(),
           norm_no['h_bond_don'].mean(),
           norm_no['heavy'].mean(),
           norm_no['rotatable_bond'].mean(),
           norm_no['tpsa'].mean(),
           norm_no['mw'].mean(),
           norm_no['xlogp'].mean()]

plt.bar(x_loc-0.2, desc_most, width=0.4, label='Most-DILI')
plt.bar(x_loc+0.2, desc_no, width=0.4, label='No-DILI')
plt.title('Descriptor mean values in dilirank')
plt.xlabel('Descriptors')
plt.ylabel('Mean values')
plt.xticks(range(len(x)), x, rotation=-45, ha='left', rotation_mode='anchor')
plt.legend(loc=0)
plt.savefig(os.path.join(wdir, 'dili_most_no_desc.png'),bbox_inches='tight')
plt.close()

"""
Let's compare chemical space between train and test.

IMPORTANT: 
If classification data is randomly separated, labels between train and test can be unbalanced.
Therefore, label balance should be checked after data set split.
'stratify' option considers the number of labels in data split.

*ref: https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
"""
from sklearn.model_selection import train_test_split

#df_train, df_test = train_test_split(df, test_size=0.2)
df_train, df_test = train_test_split(df, test_size=0.2, stratify=df['vDILIConcern'])
print ('Label ratio (training set)')
print (df_train['vDILIConcern'].value_counts()/df_train.shape[0])
print (df_test['vDILIConcern'].value_counts()/df_test.shape[0])

#Modify codes for 'dili_most_no_chemspace.png' to visualize df_train and df_test data on a scatter plot with mw and xlogp.
#Modify codes for 'dili_most_no_desc.png' to compare descriptor distribution between df_train and df_test.

###Start your code.