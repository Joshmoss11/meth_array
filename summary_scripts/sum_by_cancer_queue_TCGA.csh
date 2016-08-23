#!/bin/csh

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA
cd $parent_path
mkdir summary
mkdir summary/by_cancer
mkdir summary/by_cancer/group_by_cancer_out

foreach a (`cat cancer_queue`)
	if (-e summary/by_cancer/${a} != 1) then
		echo Working on $a
		mkdir summary/by_cancer/${a}
		setenv input_betas_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA/betas/by_cancer/${a}
		setenv output_betas_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA/summary/by_cancer/${a}

	        sbatch --output=summary/by_cancer/group_by_cancer_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_cancer_TCGA.csh $input_betas_path $output_betas_path
	endif
end

