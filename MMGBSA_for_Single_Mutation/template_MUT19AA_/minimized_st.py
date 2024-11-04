#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
1) move the terminal res along a direction vector
2) minimize local structure with FF minimizer

Created on 2019-4-16

@author: Boxue Tian, PhD
 Jacobson lab, UCSF
Contact: boxue.tian@gmail.com
'''
import decimal,math,multiprocessing,os,time
from schrodinger import structure
import schrodinger.structutils.minimize as minimize
from schrodinger.structutils.transform import translate_structure,rotate_structure,get_angle_between_vectors
import numpy as np

PATH = os.getcwd()
DIELECTRIC = decimal.Decimal(78.3553)
STEP_SIZE = decimal.Decimal(0.05)
FORCE_CONSTANT = decimal.Decimal(1000.0)
STANDARD_AMINO_ACIDS = ('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly',
                        'His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
                        'Thr','Trp','Tyr','Val',
                        # we also include some variations of proton states
                        'Hid','Hie','Hip','Arn',
                        'Ash','Glh','Lyn','Lyp','Cyn','Cym','Nme','Ace')
STANDARD_AMINO_ACIDS_UPPER = tuple([i.upper() for i in STANDARD_AMINO_ACIDS])

TRAJ = 'traj.maegz'

def findAtomFromCt(ct,atom_index):
    '''
    find the atom in ct, save in ct_atom
    '''
    for atom in ct.atom:
        if atom.index == atom_index:
            ct_atom = atom
            break
    return ct_atom  

def findBackboneNitrogen(ct,res):
    '''
    backbone nitrogen atom object
    '''
    N = [i.index for i in res.atom if i.pdbname.strip()=='CA'][0]
    return N

def findResNitrogenByResNum(ct,resnum):
    for res in ct.residue:
        if res.resnum == resnum:            
            break
    return findBackboneNitrogen(ct,res)   

def findNTermNitrogen(ct):
    '''
    return N-term nitrogen atom object
    '''
    for res in ct.residue:
        N_term = res
        break
    return findBackboneNitrogen(ct,N_term)

def findCTermNitrogen(ct):
    '''
    return C-term nitrogen atom object
    '''
    for res in ct.residue:
        if res.pdbres.strip() in STANDARD_AMINO_ACIDS_UPPER:
            C_term = res   
    return findBackboneNitrogen(ct,C_term)  

def writeStructure(st,output_maegz):
    writer = structure.StructureWriter(output_maegz,overwrite=False)    
    writer.append(st)    

def readFirstStructure(maegz):
    for st in structure.StructureReader(maegz):
        ct = st
        break
    return ct

def minimizeLocalStructure(st,frozen_atom_index_list=[],constraint_atom_number_list=[]):
    '''
    atom in frozen atom list is frozen
    atom in constraint atom list is constraint with a force constant
    '''
    ct = st.copy()
    min = minimize.Minimizer(struct=ct,ffld_version=16,dielectric_constant=DIELECTRIC) # opls3
    min.deleteAllRestraints()
    for atom in ct.atom:
        if atom.index in constraint_atom_number_list:
            min.addPosRestraint(atom.index,FORCE_CONSTANT)
        if atom.index in frozen_atom_index_list:
            min.addPosFrozen(atom.index)        
    min.minimize()
    while min.min_converged != 1:
        min.minimize()
    minimized_ct = min.getStructure()
    final_energy = min.min_energy
    return minimized_ct,final_energy 

def moveNTermToZero(ct):
    '''
    Move NTerm atoms to zero
    '''
    for res in ct.residue:
        N_term = res        
        break   
    N_atom_index = findResNitrogenByResNum(ct,N_term.resnum)    
    N_atom = findAtomFromCt(ct,N_atom_index)
    x,y,z = N_atom.xyz   
    translate_structure(ct,-x,-y,-z)
    return ct 

def fitCTermToX(ct):
    '''
    rotate ct, make sure center of mass of C-term is on the x-axis y = 0 and z = 0 
    ''' 
    Cterm_nitrogen_index = findCTermNitrogen(ct)
    Cterm_nitrogen = findAtomFromCt(ct,Cterm_nitrogen_index)
    x,y,z = Cterm_nitrogen.xyz
    v1,v2 = np.array((0.0,y,z)),np.array((0.0,1.0,0.0))
    angle_x = 0.5*math.pi - get_angle_between_vectors(v1,v2) 
    if Cterm_nitrogen.y > 0:
        rotate_structure(ct,angle_x,0.0,0.0) 
    else:
        rotate_structure(ct,-angle_x,0.0,0.0)       
    x,y,z = Cterm_nitrogen.xyz              
    v3,v4 = np.array((x,0.0,z)),np.array((0.0,0.0,1.0))
    angle_y = 0.5*math.pi - get_angle_between_vectors(v3,v4)   
    if Cterm_nitrogen.x > 0:
        rotate_structure(ct,0.0,angle_y,0.0) 
    else:
        rotate_structure(ct,0.0,-angle_y,0.0)                      
    return ct  


if __name__ == '__main__':
    t1 = time.time()

    initial_structure = readFirstStructure('template.maegz')
    pos_restraint_atoms = []
    pos_restraint_atoms.append(findNTermNitrogen(initial_structure))
    pos_restraint_atoms.append(findCTermNitrogen(initial_structure))   
    minimized_structure, ref_energy = minimizeLocalStructure(initial_structure,[],pos_restraint_atoms)
    minimized_structure = moveNTermToZero(minimized_structure)
    minimized_structure = fitCTermToX(minimized_structure)    
    writeStructure(minimized_structure,TRAJ)     

    t2 = time.time()
    print('Minimization of Mutations '+PATH.split('/')[-1]+' is complete!    time: '+str(t2-t1).split('.')[0]+' s')

