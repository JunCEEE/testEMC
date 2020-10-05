#!/bin/bash
#SBATCH --job-name=s11
#SBATCH --partition=exfel-theory
##SBATCH --partition=exfel
##SBATCH --partition=exfel-spb
#SBATCH --nodes=6
#SBATCH --time=100:00:00 
#SBATCH --output=log.diffr-%j

export LD_PRELOAD=""                 # useful on max-display nodes, harmless on others
source /etc/profile.d/modules.sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/gpfs/exfel/data/user/juncheng/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/gpfs/exfel/data/user/juncheng/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/gpfs/exfel/data/user/juncheng/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/gpfs/exfel/data/user/juncheng/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

#source $HOME/.zshrc
#conda activate simex
source $HOME/simex-0.5.profile

which conda
which python
module list
export OMP_NUM_THREADS=1

DATADIR=./
SCRATCH=/gpfs/exfel/data/scratch/juncheng
SPATH=$SCRATCH/$SLURM_JOB_NAME.$SLURM_JOB_ID
echo SCRATCH_PATH=$SPATH

# Show nodes hostnames
mpirun --map-by node -np $SLURM_JOB_NUM_NODES  hostname

# output file name without extension suffix
OUTFILE=diffr

### Clean up previous run
rm -r ${OUTFILE}.h5 ${OUTFILE}

# export SIMEX_VERBOSE=ON
mkdir -p $SPATH
python diffr.py -s $SPATH -o $OUTFILE
cp -r $SPATH/${OUTFILE}.h5 $SPATH/${OUTFILE} $DATADIR
rm -rf $SPATH
