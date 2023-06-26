#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
ALIGNMENT=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### concatenated_statistics.sh: Starting script."
echo -e "#### concatenated_statistics.sh: Alignment file: $ALIGNMENT \n\n"

################################################################################
#### CALCULATE STATISTICS FOR CONCATENATED ALIGNMENT ####
################################################################################
echo -e "#### concatenated_statistics.sh: Calculating statistics for concatenated alignment ..."
python $AMAS summary -f nexus -d dna -i $ALIGNMENT -o $ALIGNMENT.stats -s
mv $ALIGNMENT-seq-summary.txt $ALIGNMENT.taxon.stats

## Report:
echo -e "#### concatenated_statistics.sh: Done with script."
date





