#!/bin/csh
#SBATCH --mem=15000
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2

#adjust memory

if ($#argv != 2) then
    echo "Usage: $0 < input_betas_path, output_summary_path #>"
        exit 0
endif
srun Rscript  /cs/icore/joshua.moss/scripts/meth_array_scripts/summary_scripts/sum_by_tissue.R $1 $2
