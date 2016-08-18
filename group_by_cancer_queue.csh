#!/bin/csh

if ($#argv != 1) then
    echo "Usage: $0 < data source (GEO/TCGA/lab_data)#>"
        exit 0
endif

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/$1
cd $parent_path
mkdir by_cancer
mkdir by_cancer/group_process_array_out

foreach a (`cat cancer_queue`)
	echo Working on $a
	sbatch --output=by_cancer/group_process_array_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/group_by_cancer_type.csh $parent_path $a
end
