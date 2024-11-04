#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:54:42 2020

@author: Lin Guo
"""
from schrodinger import structure
import pandas as pd

STANDARD_AMINO_ACIDS = ['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly',
                        'His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
                        'Thr','Trp','Tyr','Val']
STANDARD_AMINO_ACIDS_UPPER = [i.upper() for i in STANDARD_AMINO_ACIDS]

def get_active_site_residues():
    active_sites_info = []
    for st in structure.StructureReader('Active_site_residues_within_5_Å_of_the_ligand_in_pocket.maegz'):
        for res in st.residue:
            if res.pdbres.strip() in STANDARD_AMINO_ACIDS_UPPER:
                active_sites_info.append([res.resnum, res.pdbres])
    df = pd.DataFrame(active_sites_info, columns=['AA_num','AA'])
    df.to_csv('active_site_residues.csv', index=False)

if __name__ == '__main__':
    get_active_site_residues()

