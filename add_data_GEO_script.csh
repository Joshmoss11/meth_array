

# parent_path = /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/GEO

# 1. TO_DO - update the sample_sheet.csv file in parent_path with new information
# TO_DO - add idats to parent_path/idats/by_array/
# TO_DO - add the new array numbers to the queue file found in the parent_path

#2. process new array data with the command
./process_queue.csh GEO
#check the output in process_array_out/*.out files - it should say "process_array script - FINISHED"

#3. update tissue queue
Rscript create_tissue_queue.R /mnt/lustre/hms-01/fs01/joshua.moss/dor/meth_array_data/GEO

#4. Group by tissue type
./group_by_tissue_queue.csh GEO
#check the output in betas/by_tissue/group_by_tissue_out/*.out files - it should say "group_by_tissue_type script - FINISHED "

#5.Sum data by tissue type
./summary_scripts/sum_by_tissue_queue.csh GEO
# check the output in summary/by_tissue/group_by_tissue_out/*.out files - it should say "sum_by_tissue script - FINISHED"

#6.Sum data by statistics GEO
./summary_scripts/sum_by_stat.csh GEO
# check the output in summary/by_statistics/sum_by_stat.out file - it should say "sum_by_statistics script - FINISHED"



