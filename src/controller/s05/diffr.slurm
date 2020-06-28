#!/bin/bash
#SBATCH --job-name=s05
#SBATCH --partition=exfel
##SBATCH --partition=exfel-spb
#SBATCH --nodes=1
#SBATCH --mem=2kgb
#SBATCH --time=70:00:00 
#SBATCH --output=log.diffr   # Standard output and error log

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
conda activate simex

which conda
which python
module list

export OMP_NUM_THREADS=1

python diffr.py
