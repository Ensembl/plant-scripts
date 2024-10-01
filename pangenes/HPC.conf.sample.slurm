# cluster/farm configuration file, edit as needed (use spaces or tabs)
# PATH might be empty or set to a path/ ending with '/'
TYPE    slurm
SUBEXE  sbatch
CHKEXE  squeue
DELEXE  scancel
ERROR   F
# 70GB was enough for chr-split wheat analysis with minimap2
QARGS   -p production --time=24:00:00 --mem 70G
