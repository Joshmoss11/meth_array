#!/bin/csh

if ($#argv != 1) then
    echo "Usage: $0 <data source (GEO/lab_data)>"
    exit 0
endif

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/${1}
cd $parent_path
mkdir summary
mkdir summary/by_tissue
mkdir summary/by_tissue/group_by_tissue_out

setenv working_path ${parent_path}/betas/by_tissue

foreach a (`cat tissue_queue`) 
        echo Working on $a
	setenv input_betas_file ${working_path}/${a}.csv
	setenv output_betas_file ${parent_path}/summary/by_tissue/${a}.csv

        sbatch --output=summary/by_tissue/group_by_tissue_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_tissue.csh ${input_betas_file} ${output_betas_file}
end

