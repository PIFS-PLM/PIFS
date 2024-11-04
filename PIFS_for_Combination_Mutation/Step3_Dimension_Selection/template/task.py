#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:32:19 2022

@author: user
"""
''' Use the simplest LogisticRegression for conserved and active site residues feature'''

import embedding as Emb
import time,random, os, torch, warnings
import numpy as np
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
warnings.simplefilter("ignore")
# t1 = time.time()

SEED = 7
random.seed(SEED)
np.random.seed(SEED)
os.environ['PYTHONHASHSEED'] = str(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)

id_, seqs_embd_cls, seqs_embd_class = Emb.get_embed_dataset()
data_fold, test_data = Emb.classification_dataset(id_, seqs_embd_cls, seqs_embd_class)

def input_x12(input_x):
    x1_list, x2_list = [],[]
    for x1,x2 in input_x:
        x1, x2 = x1.mean(), x2.mean()
        x1_list.append(x1*1)
        x2_list.append(x2*1)
    return np.column_stack((x1_list,x2_list))

plot_acc, plot_bacc, plot_roc_auc = [],[],[]
plot_acc_val, plot_bacc_val, plot_roc_auc_val = [],[],[]
plot_acc_test, plot_bacc_test, plot_roc_auc_test = [],[],[]

tol_list = [1e-2,1e-3,1e-4,1e-5]
C_list = [0.001,0.01,0.1,0.3,0.5,1.0,10.0,100.0,1000.0]
max_iter_list = [100,1000,10000]
for tol_ in tol_list:
    for C_ in C_list:
        for max_iter_ in max_iter_list:
            acc_mean, val_acc_mean, test_acc_mean = [], [], []
            bacc_mean, val_bacc_mean, test_bacc_mean= [], [], []
            roc_auc_mean, val_roc_auc_mean, test_roc_auc_mean = [], [], []
            for val_keys_x, train_keys_x, val_keys_y, train_keys_y in data_fold:
                model = LogisticRegression(penalty='l2', dual=True, tol=tol_, C=C_, max_iter=max_iter_,
                                           fit_intercept=True, intercept_scaling=1, class_weight='balanced',
                                           solver='liblinear', n_jobs=1)
                model.fit(input_x12(train_keys_x), train_keys_y)
                ### train ###
                test_keys_y, test_y_pred = train_keys_y, model._predict_proba_lr(input_x12(train_keys_x))
                y_pred_ = np.argmax(test_y_pred, axis=1)
                y_pred_auc = test_y_pred[:,1] 
                acc = accuracy_score(train_keys_y, y_pred_)
                bacc = balanced_accuracy_score(train_keys_y, y_pred_)
                roc_auc = roc_auc_score(train_keys_y, y_pred_auc)
                acc_mean.append(acc)
                bacc_mean.append(bacc)
                roc_auc_mean.append(roc_auc)
                ### val ###
                test_keys_y, test_y_pred = val_keys_y, model._predict_proba_lr(input_x12(val_keys_x))
                test_y_pred_ = np.argmax(test_y_pred, axis=1)
                test_y_pred_auc = test_y_pred[:,1] 
                test_acc = accuracy_score(test_keys_y, test_y_pred_)
                test_bacc = balanced_accuracy_score(test_keys_y, test_y_pred_)
                test_roc_auc = roc_auc_score(test_keys_y, test_y_pred_auc)
                val_acc_mean.append(test_acc)
                val_bacc_mean.append(test_bacc)
                val_roc_auc_mean.append(test_roc_auc)
                ### test ###
                test_keys_y, test_y_pred = test_data[1], model._predict_proba_lr(input_x12(test_data[0]))
                test_y_pred_ = np.argmax(test_y_pred, axis=1)
                test_y_pred_auc = test_y_pred[:,1] 
                test_acc = accuracy_score(test_keys_y, test_y_pred_)
                test_bacc = balanced_accuracy_score(test_keys_y, test_y_pred_)
                test_roc_auc = roc_auc_score(test_keys_y, test_y_pred_auc)
                test_acc_mean.append(test_acc)
                test_bacc_mean.append(test_bacc)
                test_roc_auc_mean.append(test_roc_auc)

            plot_acc.append(np.array(acc_mean).mean())
            plot_bacc.append(np.array(bacc_mean).mean())
            plot_roc_auc.append(np.array(roc_auc_mean).mean())
            plot_acc_val.append(np.array(val_acc_mean).mean())
            plot_bacc_val.append(np.array(val_bacc_mean).mean())
            plot_roc_auc_val.append(np.array(val_roc_auc_mean).mean())
            plot_acc_test.append(np.array(test_acc_mean).mean())
            plot_bacc_test.append(np.array(test_bacc_mean).mean())
            plot_roc_auc_test.append(np.array(test_roc_auc_mean).mean())

index_ = plot_bacc.index(max(plot_bacc))
print('train_acc','train_bacc','train_roc_auc',
      'val_acc','val_bacc','val_roc_auc',
      'test_acc','test_bacc','test_roc_auc')
print(plot_acc[index_],plot_bacc[index_],plot_roc_auc[index_],
      plot_acc_val[index_],plot_bacc_val[index_],plot_roc_auc_val[index_],
      plot_acc_test[index_],plot_bacc_test[index_],plot_roc_auc_test[index_])

