#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:16:43 2022

@author: user
"""
import subprocess, os, time
import pandas as pd

PATH = os.getcwd()
DIM_LIST = [str(i) for i in range(1280)]

def get_key_residues_num(path):
    df = pd.read_csv(path)
    num_list = list(df['AA_num'])
    return num_list

###Submit Linux command~
def linux_cmd(command):
    cmd = command
    p = subprocess.Popen(cmd,shell=True)
    p.wait()

def task1280():
    active_site5 = get_key_residues_num('../Step2_Residues_Selection/Active_Site_Residues/active_site_residues.csv')
    conserved_147 = get_key_residues_num('../Step2_Residues_Selection/Conserved_Residues_Reduced/residues_reduced.csv')
    reduced_dic = {}
    linux_cmd('mkdir task')
    for i_num, reduced in enumerate(conserved_147):
        reduced_dic['model' + str(i_num)] = reduced
    for key in reduced_dic:
        for dim in DIM_LIST:
            name = key+'_dim'+str(dim)
            print(name)
            linux_cmd('cp -r template task/'+name)
            lines = open('task/'+name+'/embedding.py','r').readlines()
            with open('task/'+name+'/embedding.py','w') as fw:
                for line in lines:
                    if line.startswith("    dim_x1 = "):
                        fw.write("    dim_x1 = ["+dim+"]\n")
                    elif line.startswith("    dim_x2 = "):
                        fw.write("    dim_x2 = ["+dim+"]\n")
                    elif line.startswith("    active_site5 = [1"):
                        fw.write("    active_site5 = ["+','.join([str(i) for i in active_site5])+"]\n")
                    elif line.startswith("    conservation147 = [1"):
                        fw.write("    conservation147 = ["+reduced_dic[key]+"]\n")
                    else:
                        fw.write(line)
def find_best_dim(input):
    ### Due to the imbalance of data, more attention is paid to train_bAcc
    ### Choose the top50 models that perform best at train_bAcc
    input.sort(key=lambda x: x[2], reverse=True)
    top50_models = input[:50]
    ### Top10 models in each indicator are counted, and the corresponding dimensions are also record
    dic_dim = {}
    for indicator in range(1, len(top50_models[0])):
        top50_models.sort(key=lambda x: x[indicator], reverse=True)
        top10_model = top50_models[:10]
        for model in top10_model:
            dim = model[0].split('dim')[1]
            if dim not in dic_dim:
                dic_dim[dim] = 1
            else:
                dim_ = dic_dim[dim] + 1
                dic_dim[dim] = dim_
    top1_dim = sorted(dic_dim.items(), key=lambda x: x[1], reverse=True)[0]
    print('Maybe you best dimension is '+ top1_dim[0])
    return top1_dim[0]

def task_run():
    list_run = sorted([i for i in os.listdir(PATH+'/task')])
    for num,file in enumerate(list_run):
        print(num,file,'Finished: ','%.2f'%((num/len(list_run))*100),'%')
        os.chdir(PATH+'/task/'+file)
        try:
            linux_cmd('python task.py > result.txt')
        except:
            print('Some error in ',file)
    csv_lines = []
    for num,file in enumerate(list_run):
        try:
            lines = open(PATH+'/task/'+file+'/result.txt').readlines()[1].strip().split(' ')
            csv_lines.append([file]+lines)
        except:
            print('Some error in '+file+'/result.txt')
    csv_lines.sort(key=lambda x: x[2], reverse=True)
    best_dim = find_best_dim(csv_lines)
    df = pd.DataFrame(csv_lines, columns=['model','train_acc', 'train_bacc', 'train_roc_auc',
                                          'val_acc', 'val_bacc', 'val_roc_auc',
                                          'test_acc', 'test_bacc', 'test_roc_auc'])
    df.to_csv(PATH+'/Best_dim_'+best_dim+'.csv', index=False)

if __name__ == '__main__':
    t1 = time.time()
    task1280()
    task_run()
    t2 = time.time()
    print('Total time taken ',t2-t1,'seconds')
