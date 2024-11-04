#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 20:32:24 2023

@author: user
"""
import time
from collections import Counter
import pandas as pd

SEQS_NUM = 1019
TARGET_SEQ_NAME = 'J7LMP2'
TARGET_SEQ_NUM = 492 # J7LMP2 - 492
LINES = open('MSA_for_TPS1019.fasta.promals.aln','r').readlines()
'''Motif: RxR-DDxxD-NSE/DTE [R280,D317,D318,D321,D461,D462,E476]'''
MOTIF_RES_NUM = [280,317,318,321,461,462,469]

'''Output the conserved residues, base on MOTIF_RES_NUM'''
def get_conserved_residues():
    seq_num = 0
    motif_min_ratio = 1.0
    conserved_residues_info = []
    for i,line in enumerate(LINES):
        if line.startswith(TARGET_SEQ_NAME):
            tar_seq = line.split(' ')[-1].strip()
            for seq_i, seq_AA in enumerate(tar_seq):
                if seq_AA != '-':
                    seq_num+=1
                    dic_aa = []
                    n=0
                    for seq_x in LINES[i-TARGET_SEQ_NUM+1:i+SEQS_NUM-TARGET_SEQ_NUM+1]:
                        seq_x_AA = seq_x.split(' ')[-1][seq_i]
                        if seq_x_AA == seq_AA:
                            n+=1
                        dic_aa.append(seq_x_AA)
                    counter = dict(Counter(dic_aa))
                    del counter['-']
                    counter_sort = sorted(counter.items(), key=lambda x:x[1],reverse=True)
                    info = [seq_num, seq_AA, float(n/SEQS_NUM),
                            counter_sort[0][0],  str(counter_sort[0][1]/SEQS_NUM)[:17],
                            counter_sort[1][0],  str(counter_sort[1][1]/SEQS_NUM)[:17],
                            float(str(counter_sort[0][1]/SEQS_NUM)[:17])+ float(str(counter_sort[1][1]/SEQS_NUM)[:17])]
                    conserved_residues_info.append(info)
                    if info[0] in MOTIF_RES_NUM:
                        if info[-1] < motif_min_ratio:
                            motif_min_ratio = info[-1]
    conserved_residues_info = sorted(conserved_residues_info, key=lambda x: x[-1],reverse=True)
    conserved_residues_info = [line for line in conserved_residues_info if not line[-1]<motif_min_ratio]
    df = pd.DataFrame(conserved_residues_info,
                      columns=['AA_num', 'AA', 'AA_conserved_ratio',
                               'MSA_top1_AA', 'MSA_top1_AA_ratio',
                               'MSA_top2_AA', 'MSA_top2_AA_ratio',
                               'MSA_top1+2_AA_ratio'])
    df.to_csv('conserved_residues.csv', index=False)

if __name__ == '__main__':
    t1 = time.time()
    get_conserved_residues()
    t2 = time.time()
    print('Time:',str(t2-t1),'s')
