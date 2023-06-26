#!/bin/bash
#SBATCH -p medium40
#SBATCH -t 48:00:00
#SBATCH --signal=B:12@1800

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
#dmtcp needs to be included in $PATH (https://dmtcp.sourceforge.io/)
#G-PhoCS needs to be included in $PATH (http://compgen.cshl.edu/GPhoCS/)

## Command-line args:
NT=$1
CTRL_FILE=$2
CHECK_OUT=$3

export OMP_NUM_THREADS=$NT # This line is necessary because otherwise G-PhoCS will only run with one core despite specifying -n

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gphocs_checkpoint.sh: Starting script."
echo -e "#### gphocs_checkpoint.sh: Number of threads: $NT"
echo -e "#### gphocs_checkpoint.sh: Control file: $CTRL_FILE"
echo -e "#### gphocs_checkpoint.sh: Checkpointing directory: $CHECK_OUT \n\n"

################################################################################
#### Coalescent modelling in G-PhoCS ####
################################################################################
trap 'echo -e "#### gphocs_checkpoint.sh: Checkpointing ..."; date; dmtcp_command --bcheckpoint; echo -e "#### gphocs_checkpoint.sh: Checkpointing done."; date; exit 12' 12

echo -e "#### gphocs_checkpoint.sh: Coalescent modelling in G-PhoCS ..."
dmtcp_launch --ckptdir $CHECK_OUT G-PhoCS $CTRL_FILE -n $NT -v &
wait




