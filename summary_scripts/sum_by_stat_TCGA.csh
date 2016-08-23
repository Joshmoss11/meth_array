#!/bin/csh

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/TCGA
cd $parent_path
mkdir summary/by_statistics
mkdir summary/by_statistics/sum_by_stat_out

setenv input_sum_dir ${parent_path}/summary/by_cancer
setenv output_stat_dir ${parent_path}/summary/by_statistics

foreach a (`cat sample_type_queue`)
	if (-e summary/by_statistics/${a} != 1) then
		echo working on $a
		mkdir summary/by_statistics/${a}
		setenv output_stat_dir ${parent_path}/summary/by_statistics/${a}
		sbatch --output=summary/by_statistics/sum_by_stat_out/${a}.out /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_statistics_TCGA.csh ${a} ${input_sum_dir} ${output_stat_dir}
	endif
end
