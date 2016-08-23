#!/bin/csh

setenv parent_path /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/lab_data
cd $parent_path
mkdir summary/by_statistics

setenv input_sum_dir ${parent_path}/summary/by_tissue
setenv output_stat_dir ${parent_path}/summary/by_statistics

        sbatch --output=summary/by_statistics/sum_by_stat.out /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_statistics_lab_data.csh ${input_sum_dir} ${output_stat_dir}
