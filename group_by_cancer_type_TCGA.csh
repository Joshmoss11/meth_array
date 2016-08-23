#!/bin/csh
#SBATCH --mem=10000
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2

#adjust memory

if ($#argv != 2) then
    echo "Usage: $0 < parent_path, cancer_name #>"
        exit 0
endif
srun Rscript  /cs/icore/joshua.moss/scripts/meth_array_scripts/group_by_cancer_type_TCGA.R $1 $2 
