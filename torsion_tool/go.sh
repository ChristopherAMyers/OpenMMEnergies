#!/bin/bash

#SBATCH -n 16
#SBATCH -p minerva
#SBATCH --exclusive

# argument 1: Q-Chem-input-file
# argument 2: Q-Chem-out-file
# argument 3: must be "-save" if requesting to 
#             save job scratch files (optional)

source /network/rit/lab/ChenRNALab/bin/Q-Chem5.3/qcenv.sh
#scratch=/network/rit/home/cm114929/ChenRNALab/awesomeSauce/scratch
scratch=/tmp/${USER}
export QCSCRATCH=$scratch
mkdir -p $scratch

/network/rit/lab/ChenRNALab/bin/Q-Chem5.3/bin/qchem $3 -nt ${SLURM_NTASKS} $1 $2
