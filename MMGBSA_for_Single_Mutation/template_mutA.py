#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 19:08:14 2020

@author: Lin Guo
"""

from schrodinger import structure
from schrodinger.structutils import build
from schrodinger.protein import rotamers

NUM = AA_NUMBER
AAA = 'AA_NAME'
STANDARD_AMINO_ACIDS = ('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly',
                        'His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
                        'Thr','Trp','Tyr','Val')
STANDARD_AMINO_ACIDS_UPPER = tuple([i.upper() for i in STANDARD_AMINO_ACIDS])
def mutation(num, AA, ST):
    for res in ST.residue:
        if res.resnum == num:
            residue = res
    # Here, mut_type is the three-letter code for the new residue.
            build.mutate(ST, residue.atom[1], AA)
            if AA != 'ALA' and AA != 'GLY' and AA != 'PRO':
                rotamer_lib = rotamers.Rotamers( ST, residue.atom[1] )
                for rotamer in rotamer_lib.rotamers:
                    rotamer.apply()
    # Operate on the new structure - compute energy, analyze or write
    # the structure.
Binding = False
for i in range(len(STANDARD_AMINO_ACIDS_UPPER)):
    if STANDARD_AMINO_ACIDS_UPPER[i] != AAA:
        writer1 = structure.StructureWriter(str(NUM) + STANDARD_AMINO_ACIDS_UPPER[i] + '.maegz')
        for st in structure.StructureReader('template.maegz'):
            st_copy = st.copy()
            try:
                mutation(NUM, STANDARD_AMINO_ACIDS_UPPER[i] ,st)
                st._setTitle(str(NUM) + STANDARD_AMINO_ACIDS_UPPER[i])
                writer1.append(st)
            except:
                st_copy._setTitle(str(NUM) + STANDARD_AMINO_ACIDS_UPPER[i])
                writer1.append(st_copy)
                Binding = True
if Binding == True:
    print('#######Info: maybe '+AAA+str(NUM)+' is the binding site, so it cannot be mutated (returning WT output) #######')
else:
    pass
