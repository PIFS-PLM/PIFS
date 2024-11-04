#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 09:46:12 2023

@author: user
"""

import torch, os, sys
import numpy as np
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from Step1_Feature_Extraction.feature_extraction import extract_features

PATH = os.getcwd()
FASTA = 'LdTPS5.fasta'
ALA_SCAN_ESM = 'AlaScan_esm.pth'
FOLD_LIST = [2.5, 1.5 ,1.0 ,0.4] # output conditions: conserved_residues decreases by 90% each time [len(conserved)/len(active_site)]
''' Get LdTPS5 sequence info '''
def get_fasta(fasta_name):
    with open(fasta_name,'r') as fr:
        lines = [line.strip() for line in fr.readlines()[1:]]
        fasta = ''.join(lines)
    return fasta

''' Get LdTPS5 mutation sequence info '''
def get_mut_seq(dic_fasta, mut):
    org_seq = dic_fasta
    if mut == 'WT':
        new_seq = org_seq
    else:
        mut_list = mut.split(',')
        mut_dic = {}
        for num_mut in mut_list:
            mut_dic[int(num_mut[1:-1])] = [num_mut[0], num_mut[-1]]
        new_seq = ''
        for aa_num in range(len(org_seq)):
            if aa_num + 1 in mut_dic:
                if org_seq[aa_num] == mut_dic[aa_num + 1][0]:
                    new_seq += mut_dic[aa_num + 1][1]
                else:
                    print('Error for ' + mut_dic[aa_num + 1][0] + str(aa_num + 1) + mut_dic[aa_num + 1][1])
            else:
                new_seq += org_seq[aa_num]
    return new_seq

''' Input_seqs for feature extraction (based on alaScan) '''
def get_esm_input_seqs(dic_fasta):
    lines = [['mut_WT', dic_fasta]]
    for i, res in enumerate(dic_fasta):
        if res == 'A':
            mut_res = 'A' + str(i + 1) + 'G'
            seqs = get_mut_seq(dic_fasta, mut_res)
        else:
            mut_res = res + str(i + 1) + 'A'
            seqs = get_mut_seq(dic_fasta, mut_res)
        line = ['mut_' + mut_res,  # id
                seqs,  # seq
                ]
        lines.append(line)
    return lines

''' Get seqence key residues'''
def get_key_residues_num(path):
    df = pd.read_csv(path)
    num_list = list(df['AA_num'])
    return num_list

''' Conserved residues reduced '''
def conserved_residues_reduced(input, output):
    ### get key residues num
    active_site5 = sorted(get_key_residues_num('../Active_Site_Residues/active_site_residues.csv'))
    conserved_147 = sorted(get_key_residues_num('../Conserved_Residues/conserved_residues.csv'))
    ### get alascan feature
    seq_emb_esm = torch.load(input)
    ### parameter ranges for each dimension in WT
    wt_delta_1280 = []
    for i in range(1280):
        x_emb = seq_emb_esm[0][:,[i]]
        wt_delta_1280.append(max(x_emb)-min(x_emb))
    wt_delta_1280 = torch.tensor(wt_delta_1280)
    ### conserved_residues decreases by 90% each time
    iter_n = 0
    fold_output_num = 0
    reduced_info = [['model0',','.join([str(i) for i in conserved_147])]]
    print('model0',
          'conserved_res_num:', len(conserved_147),
          'fold:', len(conserved_147) / len(active_site5))
    while len(conserved_147) > len(active_site5)/3:
        target_list = conserved_147
        sorted_list = [nn - 1 for nn in sorted(set(target_list))]
        ### output conditions (FOLD_LIST)
        if fold_output_num < len(FOLD_LIST):
            if len(conserved_147)/len(active_site5) < FOLD_LIST[fold_output_num]:
                reduced_info.append(['model'+str(fold_output_num+1),','.join([str(i) for i in conserved_147])])
                print('model'+str(fold_output_num+1),
                      'conserved_res_num:', len(sorted_list),
                      'fold:', len(sorted_list)/len(active_site5))
                fold_output_num += 1
        ### key_residues : active site and conserved residues
        sorted_list_AS = [nn - 1 for nn in sorted(set(target_list + active_site5))]
        ### find the dimension with the largest change
        v_top10 = {}
        for i in range(1, len(seq_emb_esm)):
            delte_all = seq_emb_esm[i] - seq_emb_esm[0]  # len(fasta),1280
            delte_select = delte_all[sorted_list_AS, :] / wt_delta_1280  # len(active_site+conserved),1280
            norm = np.linalg.norm(delte_select, axis=0)  # 1,1280
            delte_float = [float(xx) for xx in norm]
            top_10_index = sorted(range(len(delte_float)), key=lambda i: abs(delte_float[i]), reverse=True)[:10]
            for i_num, index in enumerate(top_10_index):
                if max(delte_float) != 0:
                    if index not in v_top10:
                        v_top10[index] = 1
                    else:
                        index_n = v_top10[index]
                        v_top10[index] = index_n + 1
        v_top10_list = sorted(v_top10.items(), key=lambda x: x[1], reverse=True)
        top1_dim = v_top10_list[0][0]
        # ranking of key residues in the dimension with the largest change
        update = []
        for i in range(1, len(seq_emb_esm)):
            delte_all = seq_emb_esm[i] - seq_emb_esm[0]
            delte = delte_all[sorted_list_AS, top1_dim] / wt_delta_1280[sorted_list_AS]
            update.append(delte)
        update_t = torch.stack(update)
        update_t = torch.abs(update_t).mean(0)
        dic_update = {}
        for i in range(len(sorted_list)):
            dic_update[sorted_list[i]] = update_t[i]
        dic_update_list = sorted(dic_update.items(), key=lambda x: x[1], reverse=True)
        dic_update_list_top90 = []
        for top90 in dic_update_list:
            if top90[0] + 1 in conserved_147:
                if len(dic_update_list_top90) < int(len(sorted_list) * 0.9):
                    dic_update_list_top90.append(top90[0] + 1)
        conserved_147 = sorted(dic_update_list_top90)
        iter_n += 1
    df = pd.DataFrame(reduced_info, columns=['model', 'AA_num'])
    df.to_csv(output, index=False)

if __name__ == '__main__':
    ### AlaScan feature extraction
    tar_fasta = get_fasta(FASTA)
    data_seqs = get_esm_input_seqs(tar_fasta)
    extract_features(data_seqs, ALA_SCAN_ESM)
    ### Conserved residues reduced
    conserved_residues_reduced(ALA_SCAN_ESM, 'residues_reduced.csv')

