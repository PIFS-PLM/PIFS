#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 20:10:01 2020

@author: Lin Guo
"""
import os,subprocess,multiprocessing
from schrodinger import structure


SCHRODINGER_PATH = '/opt/schrodinger2017-2' ### make sure you schrodinger pathway!
LIGAND = 'FPP' ### make sure you ligand name
NPROCS = 15 ### multiprocessing
PATH = os.getcwd()
LIST = [i for i in os.listdir(PATH) if i.endswith('.maegz')]

def linux_cmd(command):
    '''    Submits command to Linux    '''
    cmd = command
    p = subprocess.Popen(cmd,shell=True)
    p.wait()

def get_ligand_st(st):
    '''    Get ligand structure    '''
    ct=st.copy()
    list_=[]
    for res in ct.residue:
        if res.pdbres.strip() != LIGAND:
            for i in res.atom:
                list_.append(i.index)
    ct.deleteAtoms(list_)
    return ct

def get_protein_st(st):
    '''    Get protein structure    '''
    ct=st.copy()
    list_=[]
    for res in ct.residue:
        if res.pdbres.strip() == LIGAND:
            for i in res.atom:
                list_.append(i.index)
    ct.deleteAtoms(list_)
    return ct

def write_st(st,file_name):
    ''' Output of the input file for the calculation of MMGBSA '''
    ct1 = get_protein_st(st)
    ct2 = get_ligand_st(st)
    writer = structure.StructureWriter(file_name,overwrite=False)    
    writer.append(ct1)
    writer.append(ct2)

def sub_minimized(file):
    '''    Subtask for minimized    '''
    os.chdir(PATH+'/'+file)
    linux_cmd(SCHRODINGER_PATH+'/run minimized_st.py')

def sub_mmgbsa(file,mm):
    '''    Subtask for MMGBSA    '''
    os.chdir(PATH+'/'+file+'/mmgbsa')
    linux_cmd(SCHRODINGER_PATH+'/prime_mmgbsa '+PATH+'/'+file+'/mmgbsa/'+mm+' -WAIT')
    
def pool_minimized(list_file):
    '''    Task_pool for minimized    '''
    pool1 =  multiprocessing.Pool(processes=NPROCS)
    for file in LIST:
        pool1.apply_async(sub_minimized,(file,))
    pool1.close()
    pool1.join()
    
def split_pv(list_file):
    ''' Output of all input file for the calculation of MMGBSA '''
    for file in list_file:
        print('Split protein and ligand structure: '+file)
        os.chdir(PATH+'/'+file)
        linux_cmd('mkdir mmgbsa')
        num=0
        for st in structure.StructureReader(PATH+'/'+file+'/traj.maegz'):
            ct1 =st.copy()
            write_st(ct1,PATH+'/'+file+'/mmgbsa/pv-'+str(num)+'.maegz')
        
def pool_mmgbsa(list_file):
    '''    Task_pool for MMGBSA    '''
    pool2 =  multiprocessing.Pool(processes=NPROCS)
    for file in list_file:
        mmgbsa_list = [mf for mf in os.listdir(PATH+'/'+file+'/mmgbsa') if mf.endswith('.maegz')]
        print('Run MMGBSA: '+file)
        for mm in mmgbsa_list:
            pool2.apply_async(sub_mmgbsa,(file,mm,))
        # break
    pool2.close()
    pool2.join()

def get_mmgbsa_score_txt(files):
    '''    Read and write MMGBSA_score    '''
    for file in files:
        mmgbsa_file = [ms for ms in os.listdir(PATH+'/'+file+'/mmgbsa') if ms.endswith('-out.maegz')][0]
        n=0
        for st in structure.StructureReader(PATH+'/'+file+'/mmgbsa/'+mmgbsa_file):
            n+=1
            if n==2:
                mmgbsa_score = st.property['r_psp_MMGBSA_dG_Bind']
        with open(PATH+'/'+file+'/mmgbsa/mmgbsa.txt','w') as f:
            line = '\t'.join([file, str(mmgbsa_score)]) + '\n'
            f.write(line)

if __name__ == '__main__':
    print('#######Start task to run minimized......#######')
    pool_minimized(LIST)
    print('#######Start task to split protein and ligand structure......#######')
    split_pv(LIST)
    print('#######Start task to run MMGBSA......#######')
    pool_mmgbsa(LIST)
    print('#######Start task to get MMGBSA_score......#######')
    get_mmgbsa_score_txt(LIST)
    os.chdir(PATH)
    linux_cmd('mkdir traj_files')
    linux_cmd('mv *maegz traj_files/.')
