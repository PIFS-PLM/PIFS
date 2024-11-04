#!/bin/sh

cd MMGBSA_for_Single_Mutation 
python write_file_and_get_mmgbsa.py

cd ../PIFS_for_Combination_Mutation/Step1_Feature_Extraction
python feature_extraction.py

cd ../Step2_Residues_Selection/Active_Site_Residues
python get_active_site_residues.py

cd ../Conserved_Residues
python get_conserved_residues.py

cd ../Conserved_Residues_Reduced
python conserved_residues_reduced.py

cd ../../Step3_Dimension_Selection
python dimension_selection.py

cd ../Step4_Ensemble_Classifier
python ensemble_classifier.py
