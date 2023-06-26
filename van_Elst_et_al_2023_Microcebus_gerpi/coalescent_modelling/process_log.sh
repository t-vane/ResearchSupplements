#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
MODEL=$2
RUN_ID=$3
MCMC=$4
BURNIN=$5
LAST_SAMPLE=$6
MUTRATE=$7
MUTRATE_VAR=$8
GENTIME=$9
GENTIME_SD=$10
M_SCALE=$11
T_SCALE=$12
POPLIST_FILE=$13

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### process_logs.sh: Starting script."
echo -e "#### process_logs.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### process_logs.sh: Model: $MODEL"
echo -e "#### process_logs.sh: Run ID: $RUN_ID"
echo -e "#### process_logs.sh: MCMC file: $MCMC"
echo -e "#### process_logs.sh: Burn-in: $BURNIN"
echo -e "#### process_logs.sh: Last sample to keep: $LAST_SAMPLE"
echo -e "#### process_logs.sh: Gamma distribution of mutation rate will have mean $MUTRATE * 10e-8"
echo -e "#### process_logs.sh: Gamma distribution of mutation rate will have variance $MUTRATE_VAR * 10e-8"
echo -e "#### process_logs.sh: Lognormal distribution of generation time will have mean ln($GENTIME)"
echo -e "#### process_logs.sh: Lognormal distribution of generation time will have mean ln($GENTIME_SD)"
echo -e "#### process_logs.sh: Inverse scaling factor used in the G-PhoCS configuration file for migration parameter: $M_SCALE"
echo -e "#### process_logs.sh: Inverse scaling factor used in the G-PhoCS configuration file for tau and theta: $T_SCALE"
echo -e "#### process_logs.sh: File with information on parent and child populations: $POPLIST_FILE \n\n"

################################################################################
#### PROCESS LOGS ####
################################################################################
echo -e "#### process_logs.sh: Reformatting G-PhoCS output to remove empty column inserted by bug ..."
awk -v OFS="\t" '$1=$1' $MCMC > $MCMC.reform

echo -e "#### process_logs.sh: Processing logs and converting to demographic values ..."
Rscript $SCRIPTS_DIR/process_logs.R $MODEL $RUN_ID $MCMC.reform $BURNIN $LAST_SAMPLE $MUTRATE $MUTRATE_VAR $GENTIME $GENTIME_SD $M_SCALE $T_SCALE $POPLIST

## Report:
echo -e "\n#### process_logs.sh: Done with script."
date


