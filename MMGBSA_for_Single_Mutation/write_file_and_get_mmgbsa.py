#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:54:42 2020

@author: Lin Guo
"""
import os, subprocess,time
from schrodinger import structure
from collections import OrderedDict
import pandas as pd

SCHRODINGER_PATH = '/opt/schrodinger2017-2' ### make sure you schrodinger pathway!
PATH = os.getcwd()
STANDARD_AMINO_ACIDS = ['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly',
                        'His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
                        'Thr','Trp','Tyr','Val']
STANDARD_AMINO_ACIDS_UPPER = [i.upper() for i in STANDARD_AMINO_ACIDS]


'''Returns a dictionary containing information about the list of amino acids in the structure.'''
def get_dic():
    dic_wt_aa=OrderedDict()
    for st in structure.StructureReader('template_MUT19AA_/MUT_file/template.maegz'):
        for res in st.residue:
            if res.pdbres.strip() in STANDARD_AMINO_ACIDS_UPPER:
                dic_wt_aa[res.resnum]=res.pdbres.strip()
            else:
                if res.pdbres.strip() in ['HID','HIE','HIP']:
                    dic_wt_aa[res.resnum] = 'HIS'
                elif res.pdbres.strip() == 'ARN':
                    dic_wt_aa[res.resnum] = 'ARG'
                elif res.pdbres.strip() == 'ASH':
                    dic_wt_aa[res.resnum] = 'ASP'
                elif res.pdbres.strip() == 'GLH':
                    dic_wt_aa[res.resnum] = 'GLU'
                elif res.pdbres.strip() in ['LYN','LYP']:
                    dic_wt_aa[res.resnum] = 'LYS'
                elif res.pdbres.strip() in ['CYN','CYM']:
                    dic_wt_aa[res.resnum] = 'CYS'
                elif res.pdbres.strip() in ['NME','ACE']:
                    pass
                else:
                    print("#######Info: make sure it's not an amino acid --- " +str(res.resnum)+res.pdbres.strip()+"#######")
                
    print("#######Info: number of amino acids in the structure: "+str(len(dic_wt_aa))+"#######")
    print("#######Info: dictionary of amino acids information #######")
    print(dic_wt_aa)
    return dic_wt_aa

'''    Submits command to Linux    '''
def linux_cmd(command):
    cmd = command
    p = subprocess.Popen(cmd,shell=True)
    p.wait()

'''    Get MMGBSA result    '''
def get_result_csv(wt_aa_dic):
    muts_dic = OrderedDict()
    muts_aa = [i for i in os.listdir(PATH + '/traj/traj_files') if i.endswith('.maegz')]
    for mut in muts_aa:
        with open(PATH + '/traj/traj_files/' + mut + '/mmgbsa/mmgbsa.txt', 'r') as f:
            mmgbsa = f.readlines()
            assert len(mmgbsa) == 1
            mmgbsa_ = mmgbsa[0].split('\t')
            assert len(mmgbsa_) == 2
            muts_dic[mut[:-6]] = float(mmgbsa_[1].strip())

    wt_mmgbsa = muts_dic['template']
    mmgbsa_lines = []
    for wt_aa in wt_aa_dic:
        line = [wt_aa]
        for AA in STANDARD_AMINO_ACIDS_UPPER:
            if AA != wt_aa_dic[wt_aa]:
                line.append(muts_dic[str(wt_aa)+AA])
            else:
                line.append(wt_mmgbsa)
        mmgbsa_lines.append(line)
    df = pd.DataFrame(mmgbsa_lines,
                      columns=['AA_num']+STANDARD_AMINO_ACIDS_UPPER)
    df.to_csv('mmgbsa.csv', index=False)

    delta_mmgbsa_lines = []
    for wt_aa in wt_aa_dic:
        line = [wt_aa]
        for AA in STANDARD_AMINO_ACIDS_UPPER:
            if AA != wt_aa_dic[wt_aa]:
                line.append(muts_dic[str(wt_aa)+AA]-wt_mmgbsa)
            else:
                line.append(wt_mmgbsa-wt_mmgbsa)
        line.append(len([mm for mm in line[1:] if mm < 0]))
        delta_mmgbsa_lines.append(line)
    delta_mmgbsa_lines = sorted(delta_mmgbsa_lines, key=lambda x: x[-1],reverse=True)
    df = pd.DataFrame(delta_mmgbsa_lines,
                      columns=['AA_num']+STANDARD_AMINO_ACIDS_UPPER+['Number of beneficial mutations'])
    df.to_csv('delta_mmgbsa.csv', index=False)

def keeprun():
    wt_aa_dic = get_dic()
    '''Get single mutation files'''
    for key in wt_aa_dic:
        # if 316<key<322: ### you can restrict the mutated sites for testing here
            file_name = str(key)+'_'+wt_aa_dic[key]+'_MUT19AA'
            linux_cmd('cp -r template_MUT19AA_ '+file_name)
            linux_cmd('mkdir '+ file_name + '/result')
            line1 = open('template_mutA.py').readlines()
            with open(PATH+'/'+file_name+'/MUT_file/mutA.py','w') as fp:
                for s in line1:
                    fp.write(s.replace('AA_NUMBER',str(key)).replace('AA_NAME',wt_aa_dic[key]))
            # break

    list_mut_files = [i1 for i1 in os.listdir(PATH) if (i1.endswith('_MUT19AA'))]
    for mut_file in list_mut_files:
        print('Info: '+mut_file+' is in process......')
        path_ = PATH+'/'+mut_file
        os.chdir(path_+'/MUT_file')
        linux_cmd(SCHRODINGER_PATH+'/run mutA.py')
        txt_list = [i for i in os.listdir(path_+'/MUT_file') if (i.endswith('.maegz'))]
        for i in txt_list:
            os.chdir(path_+'/result')
            linux_cmd('mkdir '+ i)
            os.chdir(path_)
            linux_cmd('cp '+path_+'/MUT_file/'+i+' '+path_+'/result/'+i+'/template.maegz')
            linux_cmd('cp minimized_st.py '+path_+'/result/'+i+'/')
        os.chdir(PATH)

    for f in list_mut_files:
        linux_cmd('cp -r '+PATH+'/'+f+'/result/* '+PATH+'/traj/')
        linux_cmd('rm -r '+PATH+'/traj/template.maegz')
    linux_cmd('cp -r '+PATH+'/'+list_mut_files[0]+'/result/template.maegz '+PATH+'/traj/')
    linux_cmd('rm -r '+PATH+'/*_MUT19AA/result')
    linux_cmd('mkdir mut_files')
    linux_cmd('mv *_MUT19AA mut_files/.')


    '''Processing structure, run MMGBSA and output the results'''
    os.chdir(PATH+'/traj')
    linux_cmd(SCHRODINGER_PATH+'/run task_mmgbsa.py')
    get_result_csv(wt_aa_dic)



if __name__ == '__main__':
    t1 = time.time()
    keeprun()
    t2 = time.time()
    print('Time: '+str(t2-t1).split('.')[0]+' s')
