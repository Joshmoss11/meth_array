#!/bin/csh
#SBATCH --mem=5000
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2

#adjust memory

if ($#argv != 3) then
    echo "Usage: $0 < sample_type, summary_files_dir, stat_files_dir #>"
        exit 0
endif
srun Rscript  /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_statistics_TCGA.R $1 $2 $3
