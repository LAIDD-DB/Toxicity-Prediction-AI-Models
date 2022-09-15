# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:44:43 2022

@author: Gil

*Description on scikit-learn SVM.
https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
"""
import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import svm

wdir = r'directory/for/training/and/test/set'
train = 'dilirank_train.csv'
test = 'dilirank_test.csv'

df_train = pd.read_csv(os.path.join(wdir, train))
y_train = df_train['vDILIConcern']
x_train = df_train.drop(['Compound Name','vDILIConcern'],axis=1)

df_test = pd.read_csv(os.path.join(wdir, test))
y_test = df_test['vDILIConcern']
x_test = df_test.drop(['Compound Name','vDILIConcern'],axis=1)

#Chemical space between x_train and x_test?? 

dili_model = svm.SVC()
dili_model.fit(x_train, y_train)

print ('Train accuracy:', dili_model.score(x_train, y_train))
print ('Test accuracy:', dili_model.score(x_test, y_test))

"""
Cross validation and grid search.
"""

print ('Apply k-fold cross validation with grid search')
import numpy as np
from sklearn.model_selection import KFold
kf = KFold(n_splits=3)
c_list = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
cv_acc_list = []
test_acc_list = []
for c in c_list:
    dili_classifier = svm.SVC(C=c)
    acc_list = []
    
    for tr_i, val_i in kf.split(df_train):
        x_train_cv = x_train.iloc[tr_i]
        y_train_cv = y_train.iloc[tr_i]

        dili_classifier.fit(x_train_cv, y_train_cv)
        
        x_val_cv = x_train.iloc[val_i]
        y_val_cv = y_train.iloc[val_i]
        
        acc_list.append(dili_classifier.score(x_val_cv, y_val_cv))
    cv_acc_list.append(np.mean(acc_list))
    test_acc_list.append(dili_classifier.score(x_test, y_test))
    print ('Model training with C:',c,'Cross valiation:',cv_acc_list[-1],'Test:',test_acc_list[-1])
    

"""
Save and load the model.
"""

from joblib import dump
dili_model = svm.SVC(C=10)
dili_model.fit(x_train, y_train)
dump(dili_model, os.path.join(wdir, 'dili_model.joblib'))
print ('Saved model accuracy (train):', dili_model.score(x_train, y_train))
print ('Saved model accuracy (test):', dili_model.score(x_test, y_test))

from joblib import load

loaded_model = load(os.path.join(wdir, 'dili_model.joblib'))
print ('Loaded model accuracy (train):', loaded_model.score(x_train, y_train))
print ('Loaded model accuracy (test):', loaded_model.score(x_test, y_test))
