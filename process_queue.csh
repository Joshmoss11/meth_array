#!/bin/tcsh
if ($#argv != 1) then
    echo "Usage: $0 <data source (GEO/TCGA/lab_data)>"
    exit 0
endif
setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/$1
cd $parent_path
mkdir betas
mkdir betas/by_array
mkdir process_array_out
foreach a (`cat queue`)
	if ( -e ${parent_path}/betas/by_array/${a}.csv != 1) then
		echo Working on $a
		sbatch --output=process_array_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/process_array.csh $parent_path $a
	endif
end
