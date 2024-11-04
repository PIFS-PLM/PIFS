#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:16:43 2022

@author: user
"""
import subprocess, os, time
import pandas as pd
from collections import OrderedDict


PATH = os.getcwd()

def get_key_residues_num(path):
    df = pd.read_csv(path)
    num_list = list(df['AA_num'])
    return num_list

###Submit Linux command~
def linux_cmd(command):
    cmd = command
    p = subprocess.Popen(cmd,shell=True)
    p.wait()

def task_best_dim():
    active_site5 = get_key_residues_num('../Step2_Residues_Selection/Active_Site_Residues/active_site_residues.csv')
    conserved_147 = get_key_residues_num('../Step2_Residues_Selection/Conserved_Residues_Reduced/residues_reduced.csv')
    best_dim = [i for i in os.listdir('../Step3_Dimension_Selection') if i.endswith('csv')]
    reduced_dic = {}
    linux_cmd('mkdir task')
    for i_num, reduced in enumerate(conserved_147):
        ### Base on Step3_Dimension_Selection result, we only use the middle three reduced conserved_147 list
        if i_num in [1,2,3]:
            reduced_dic['model' + str(i_num)] = reduced
    for key in reduced_dic:
        for dim in best_dim:
            for cv_num in range(10):
                dim = dim.split('_')[-1].split('.')[0]
                name = key+'_dim'+str(dim)+'_cv'+str(cv_num)
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
                        elif line.startswith("    kf = KFold"):
                            fw.write("    kf = KFold(n_splits=5, shuffle=True,random_state=0)\n")
                            fw.write("\n")
                            fw.write("    test_n = " + str(cv_num) + "\n")
                            fw.write("    kf_test = KFold(n_splits=10, shuffle=True,random_state=0)\n")
                            fw.write("    test_id_kf = list(kf_test.split(id_))\n")
                            fw.write("    test_id = [id_[ii] for ii in list(test_id_kf[test_n][1])]\n")
                            fw.write("\n")
                        else:
                            fw.write(line)

def task_run():
    list_run = sorted([i for i in os.listdir(PATH+'/task')])
    for num,file in enumerate(list_run):
        print(num, file, 'Finished: ', '%.2f' % ((num / len(list_run)) * 100), '%')
        os.chdir(PATH+'/task/'+file)
        try:
            linux_cmd('python task.py')
        except:
            print('Some error in ',file)
    pd_data = OrderedDict()
    for num,file in enumerate(list_run):
        # print(num, file, 'Finished: ', '%.2f' % ((num / len(list_run)) * 100), '%')
        csv_list = sorted([i for i in os.listdir(PATH + '/task/'+ file) if i.endswith('csv')])
        for csv in csv_list:
            df = pd.read_csv(PATH + '/task/'+ file+'/'+csv)
            mut_id_list = list(df['mut_id'])
            pred_list = list(df['pred_proba'])
            assert len(mut_id_list) == len(pred_list)
            for i in range(len(mut_id_list)):
                if mut_id_list[i] not in pd_data:
                    pd_data[mut_id_list[i]] = pred_list[i]
                else:
                    sum_num = pd_data[mut_id_list[i]] + pred_list[i]
                    pd_data[mut_id_list[i]] = sum_num

    pd_data_sort = sorted(pd_data.items(), key=lambda x: x[1], reverse=True)
    df_merge = pd.DataFrame(pd_data_sort, columns=['Mutant','Voting'])
    df_merge.to_csv(PATH+'/Voting_result.csv', index=False)


if __name__ == '__main__':
    t1 = time.time()
    task_best_dim()
    task_run()
    t2 = time.time()
    print('Total time taken ',t2-t1,'seconds')
