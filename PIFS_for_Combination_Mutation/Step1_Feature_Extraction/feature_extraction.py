#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 09:46:12 2023

@author: user
"""

import torch, esm, time
import pandas as pd

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

def get_seq(path):
    df = pd.read_csv(path, names=['id', 'seq', 'class'], skiprows=1)
    id_ = list(df['id'])
    seqs = list(df['seq'])
    data = []
    for i, seq in enumerate(seqs):
        data.append((id_[i], seq))
    return data
def extract_features(data_set, output):
    # Prepare data
    sequence_representations = []
    for i_num in range(len(data_set)):
        tt1 = time.time()
        data = data_set[i_num:1 + i_num]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        # Extract per-residue representations (on CPU)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33]

        # Generate per-sequence representations via averaging
        # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.

        for i, tokens_len in enumerate(batch_lens):
            # sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))
            sequence_representations.append(token_representations[i, 1 : tokens_len - 1])
        tt2 = time.time()
        print(data[0][0],'time: ',tt2-tt1)
    torch.save(sequence_representations,output)

if __name__ == '__main__':
    t1 = time.time()

    train_data = get_seq('ldtps5_train.csv')
    extract_features(train_data, 'train_esm.pth')

    test_data = get_seq('ldtps5_test.csv')
    extract_features(test_data, 'test_esm.pth')
    t2 = time.time()
    print('Time: '+str(t2-t1).split('.')[0]+' s')

