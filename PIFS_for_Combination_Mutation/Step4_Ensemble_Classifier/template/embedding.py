#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 20:40:13 2022

@author: user
"""

from sklearn.model_selection import KFold
import pandas as pd
import time, random, torch, re, warnings
warnings.simplefilter("ignore")
t1 = time.time()

MUT12 = [290, 311, 315, 316, 385, 434, 454, 457, 468, 476, 553, 558]
INPUT_CSV = '../../../Step1_Feature_Extraction/ldtps5_train.csv'
INPUT_NPY = '../../../Step1_Feature_Extraction/train_esm.pth'
INPUT_CSV_T = '../../../Step1_Feature_Extraction/ldtps5_test.csv'
INPUT_NPY_T = '../../../Step1_Feature_Extraction/test_esm.pth'

def load_dataset(path):
    df = pd.read_csv(path,names=['id','seq','class'],skiprows=1)
    id_ = list(df['id'])
    df['seq_fixed'] = ["".join(seq.split()) for seq in df['seq']]
    df['seq_fixed'] = [re.sub(r"[UZOB]", "X", seq) for seq in df['seq_fixed']]
    seqs = [ list(seq) for seq in df['seq_fixed']]
    class_ = list(df['class'])
    return id_, seqs, class_

def get_key_residues_num(path):
    df = pd.read_csv(path)
    num_list = list(df['AA_num'])
    return num_list

def embed_dataset(npy):
    dim_x1 = [777] #[n for n in range(1280)]
    dim_x2 = [777] #[n for n in range(1280)]
    active_site5 = [1]
    conservation147 = [1]
    conservation147 = [i for i in conservation147 if i not in MUT12]
    target_x1 = [nn-1 for nn in sorted(set(conservation147))]
    target_x2 = [nn-1 for nn in sorted(set(active_site5))]

    wt = torch.load(INPUT_NPY)[0]
    wt_delta_1280 = []
    for i in range(1280):
        x_emb = wt[:,[i]]
        wt_delta_1280.append(max(x_emb)-min(x_emb))
    wt_delta_1280 = torch.tensor(wt_delta_1280)
    embs = [(tensor-wt)/wt_delta_1280*100 for tensor in torch.load(npy)]
    inputs_embedding_x1 = [tensor[target_x1,dim_x1] for tensor in embs]
    inputs_embedding_x2 = [tensor[target_x2,dim_x2] for tensor in embs]

    inputs_embedding = [[inputs_embedding_x1[i],inputs_embedding_x2[i]] for i in range(len(inputs_embedding_x1))]

    return inputs_embedding

def get_embed_dataset():

    id_, seqs, class_ = load_dataset(INPUT_CSV)
    seqs_embd = embed_dataset(INPUT_NPY)

    return id_, seqs_embd, class_
def get_embed_dataset_test():

    id_, seqs, class_ = load_dataset(INPUT_CSV_T)
    seqs_embd = embed_dataset(INPUT_NPY_T)

    return id_, seqs_embd, class_, seqs
def classification_dataset(id_, seqs, class_):
    random.seed(0)
    test_id = id_[-10:]
    
    dic_id_xy = {}
    for id_i in range(len(id_)):
        dic_id_xy[id_[id_i]] = [seqs[id_i],class_[id_i]]
    random.shuffle(id_)
    kf = KFold(n_splits=5, shuffle=True,random_state=0)
    train_val_id = [id_tv for id_tv in id_ if id_tv not in test_id]
    # print(len(test_id))
    # print(len(train_val_id))

    test_x, test_y = [dic_id_xy[_id][0] for _id in test_id], [dic_id_xy[_id][1] for _id in test_id]
    test_data = [test_x, test_y]
    id_kfold = list(kf.split(train_val_id))
    data_kfold = []
    for xx, yy in id_kfold:
        train_keys_x, train_keys_y= [dic_id_xy[train_val_id[x]][0] for x in xx], [dic_id_xy[train_val_id[x]][1] for x in xx]
        val_keys_x, val_keys_y= [dic_id_xy[train_val_id[y]][0] for y in yy], [dic_id_xy[train_val_id[y]][1] for y in yy]
        
        data_kfold.append([val_keys_x, train_keys_x, val_keys_y, train_keys_y])

    return data_kfold, test_data
