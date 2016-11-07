#!/bin/tcsh

if ($#argv != 1) then
    echo "Usage: $0 <data source (GEO/lab_data)>"
    exit 0
endif


setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/${1}
cd $parent_path
mkdir betas/by_tissue
mkdir betas/by_tissue/group_by_tissue_out

foreach a ("`cat tissue_queue`")
	if ( -e ${parent_path}/betas/by_tissue/${a}.csv != 1) then
		echo working on $a
		sbatch --output=betas/by_tissue/group_by_tissue_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/group_by_tissue_type.csh $parent_path $a
	endif
end

