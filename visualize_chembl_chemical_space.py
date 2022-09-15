# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 14:52:28 2022

@author: Gil
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

wdir = r'directory_of_input_file'
csv = 'hERG_preprocessed.csv'

df = pd.read_csv(os.path.join(wdir, csv))

"""
First visualize endpoint distribution.
Compare IC50 and pIC50.
"""

plt.hist(df['Standard Value'])
plt.title('hERG IC50 histogram.')
plt.xlabel('IC50')
plt.savefig(os.path.join(wdir, 'hERG_ic50_hist.png'))
plt.close()

plt.hist(df['pIC50'])
plt.title('hERG pIC50 histogram.')
plt.xlabel('pIC50')
plt.savefig(os.path.join(wdir, 'hERG_pic50_hist.png'))
plt.close()

"""
Chemical space distribution with endpoint.

In order to prepare a figure for the publication, specific size information can be given to adjust the size and resoultion of the figure.
"""
#AlogP columns has a noise.
df = df[df['AlogP']!='None']

plt.figure(figsize=[3,3],dpi=300)
plt.scatter(df['Molecular Weight'],df['AlogP'].astype('float'), c=df['pIC50'], s=5, edgecolors='k', linewidths=.1)
plt.title('hERG chemical space', fontsize=9)
plt.xlabel('MW', fontsize=7)
plt.ylabel('LogP', fontsize=7)
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
cbar = plt.colorbar()
cbar.set_label('pIC50', size=7)
cbar.ax.tick_params(labelsize=5)
plt.savefig(os.path.join(wdir, 'hERG_chemspace_pic50.png'))
plt.close()

"""
Let's compare train and test set distribution.
"""
import numpy as np
from sklearn.model_selection import train_test_split

df_train, df_test = train_test_split(df, test_size=0.2)

bins = np.linspace(-5.5, 1.5, 50)

plt.hist(df_train['pIC50'], bins, alpha=0.5, label='Train')
plt.hist(df_test['pIC50'], bins, alpha=0.5, label='Test')
plt.title('hERG pIC50: training vs test')
plt.legend(loc=0)
plt.savefig(os.path.join(wdir, 'herg_pic50_train_test.png'))
plt.close()

"""
Multi dimensional data can be used in chemical space visualization.
"""
from sklearn.decomposition import PCA

fp_csv = 'dilirank_cactvs.csv'
df_fp = pd.read_csv(os.path.join(wdir, fp_csv))

x_fp = df_fp.drop(['Compound Name', 'vDILIConcern'],axis=1)
pca = PCA(n_components=2)
x_pca = pca.fit_transform(x_fp)
df_pca = pd.DataFrame(x_pca)
dili_pca = pd.concat([df_fp[['Compound Name', 'vDILIConcern']],df_pca],axis=1)

"""
Try to compare chemical space of training and test set on the scatter plot.
Training and test set can be marked differently by shape or color.

Markers: https://matplotlib.org/stable/api/markers_api.html
Colors: https://matplotlib.org/stable/gallery/color/named_colors.html

*Use PCA to visualize chemical space.
You can use any other descriptors in PCA to visualize chemical space.
Try to use  QM descriptors or 2D molecular descriptors.
"""

###star your code.