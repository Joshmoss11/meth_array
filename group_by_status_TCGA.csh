#!/bin/csh
set parent_dir = '/mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA/betas'
set by_status_dir = ${parent_dir}/by_status
mkdir $by_status_dir
mkdir ${by_status_dir}/Primary_Tumor
mkdir ${by_status_dir}/Solid_Tissue_Normal
mkdir ${by_status_dir}/Metastatic
cp ${parent_dir}/by_cancer/*/*Primary_Tumor.csv ${by_status_dir}/Primary_Tumor
cp ${parent_dir}/by_cancer/*/*Solid_Tissue_Normal.csv ${by_status_dir}/Solid_Tissue_Normal
cp ${parent_dir}/by_cancer/*/*Metastatic.csv ${by_status_dir}/Metastatic
