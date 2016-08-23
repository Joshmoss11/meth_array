#!/bin/csh

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA
cd $parent_path
mkdir betas/by_cancer
mkdir betas/by_cancer/group_by_cancer_out

foreach a (`cat cancer_queue`)
	echo Working on $a
	sbatch --output=betas/by_cancer/group_by_cancer_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/group_by_cancer_type_TCGA.csh $parent_path $a
end
