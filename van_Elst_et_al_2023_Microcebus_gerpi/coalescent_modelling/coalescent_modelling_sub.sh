################################################################################
#### COALESCENT MODELLING ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/coalescentModelling

SET_ID=gerpi
LOCUS_DIR=$PWORK/$SET_ID/locusExtraction/fasta/${SET_ID}_bylocus_final
COAL_DIR=$PWORK/$SET_ID/coalescentModelling

mkdir -p $COAL_DIR/logFiles

#################################################################
#### 0 PREPARATION ####
#################################################################
## Reformat per-locus files
mkdir -p $COAL_DIR/loci
sbatch --wait --account=nib00015 --output=$COAL_DIR/logFiles/prepareLoci.oe $SCRIPTS_DIR/prepareLoci.sh $LOCUS_DIR $COAL_DIR/loci

## Create sequence file containing all loci
cat $COAL_DIR/loci/loci.count > $COAL_DIR/$SET_ID.seq_file.txt
find $COAL_DIR/loci -maxdepth 1 -name "*locus" -type f -exec cat {} + >> $COAL_DIR/$SET_ID.seq_file.txt

## Remove undesired samples from sequence file
REMOVE_FILE=$COAL_DIR/remove_samples.txt # List with samples to remove
REMOVE_STRING=$(awk '$1=$1' RS= OFS='\\|' $REMOVE_FILE) # Creates string for subsequent sed command
sed "/$REMOVE_STRING/d" $COAL_DIR/$SET_ID.seq_file.txt > $COAL_DIR/$SET_ID.seq_file.reduced.txt

## Replace sample number in sequence file
NO_SAMPLES_TOTAL=$(sed -n 3p $COAL_DIR/$SET_ID.seq_file.txt | awk '{print $2}'); echo $NO_SAMPLES_TOTAL
NO_SAMPLES_REMOVED=$(cat $REMOVE_FILE | wc -l); echo $NO_SAMPLES_REMOVED
NO_SAMPLES_LEFT=$(( $NO_SAMPLES_TOTAL - $NO_SAMPLES_REMOVED )); echo $NO_SAMPLES_LEFT
sed -i "s/fa $NO_SAMPLES_TOTAL /fa $NO_SAMPLES_LEFT /g" $COAL_DIR/$SET_ID.seq_file.reduced.txt

#################################################################
#### 1 PRELIMINARY MODELS ####
#################################################################
MODELS="anc_sahaDobo_sakoVohiFina anc_sahaDobo_baoSakoVohiFina anc_bao_sakoVohiFina anc_maroJoll_gerpiRoot sako_mmaro vohiFina_sako vohiFina_dobo vohiFina_saha vohiFina_bao saha_dobo sako_bao mjoll_mmaro bao_mmaro" # Contains all model IDs

## After manually creating a control file for each preliminary model (ctrl_file.txt), submit jobs
NT=20
NJOBS=8 # Since there was walltime limit of 48 h on the server, we used DMTCP (https://dmtcp.sourceforge.io/) for checkpointing; $NJOBS jobs were then submitted, each of which would only run a maximum of 48 h 

for MODEL in $MODELS
do
	mkdir -p $COAL_DIR/models/$MODEL
	
	for i in $(seq 1 $NJOBS)
	do
		echo -e "#### Submitting job $i for model $MODEL"
		# If first script:
		if [ $i == 1 ]
		then
			# Declare directory to save checkpoint files
			CHECK_OUT=$COAL_DIR/models/$MODEL/checkpoint_$i
			mkdir -p $CHECK_OUT
			# Submit job and save submission ID
			JID=$(sbatch --account=nib00015 --output=$COAL_DIR/logFiles/$MODEL.$i.gphocs.oe $SCRIPTS_DIR/gphocs_checkpoint.sh $NT $COAL_DIR/models/$MODEL/ctrl_file.txt $CHECK_OUT)
			declare RUNID_$i=${JID##* }

		# If not first script:
		else
			# Assign input directory (which is output directory of previous iteration)
			CHECK_IN=$CHECK_OUT
			# Declare directory to save checkpoint files
			CHECK_OUT=$COAL_DIR/models/$MODEL/checkpoint_$i
			# Get submission ID of previous iteration
			VARNAME=RUNID_$(( $i - 1 ))
			# Submit next job and save submission ID
			JID=$(sbatch --account=nib00015 --output=$COAL_DIR/logFiles/$MODEL.$i.gphocs.oe --dependency=afterany:${!VARNAME} $SCRIPTS_DIR/gphocs_checkpoint_cont.sh $NT $CHECK_IN $CHECK_OUT)
			declare RUNID_$i=${JID##* }
		fi
	done
done

#################################################################
#### 2 FINAL MODELS ####
#################################################################
MODELS="noMig Mig" # Contains final model IDs

##################################################
#### 2.1 RUN FINAL MODELS ####
##################################################

## After manually creating a control file for each final model (ctrl_file.txt), submit jobs
NT=20
REPLICATES=4 # Number of replicate runs to be submitted for each model
NJOBS=15 # Since there was walltime limit of 48 h on the server, we used DMTCP (https://dmtcp.sourceforge.io/) for checkpointing; $NJOBS jobs were then submitted, each of which would only run a maximum of 48 h 

for j in $(seq 1 $REPLICATES)
	for MODEL in $MODELS
	do
		mkdir -p $COAL_DIR/models/$MODEL

		for i in $(seq 1 $NJOBS)
		do
			echo -e "#### Submitting job $i of run $j for model $MODEL"
			# If first script:
			if [ $i == 1 ]
			then
				# Declare directory to save checkpoint files
				CHECK_OUT=$COAL_DIR/models/$MODEL/checkpoint_run${j}_$i
				mkdir -p $CHECK_OUT
				# Submit job and save submission ID
				JID=$(sbatch --account=nib00015 --output=$COAL_DIR/logFiles/$MODEL.run$j.$i.gphocs.oe $SCRIPTS_DIR/gphocs_checkpoint.sh $NT $COAL_DIR/models/$MODEL/ctrl_file.run$j.txt $CHECK_OUT)
				declare RUNID_$i=${JID##* }

			# If not first script:
			else
				# Assign input directory (which is output directory of previous iteration)
				CHECK_IN=$CHECK_OUT
				# Declare directory to save checkpoint files
				CHECK_OUT=$COAL_DIR/models/$MODEL/checkpoint_run${j}_$i
				# Get submission ID of previous iteration
				VARNAME=RUNID_$(( $i - 1 ))
				# Submit next job and save submission ID
				JID=$(sbatch --account=nib00015 --output=$COAL_DIR/logFiles/$MODEL.run$j.$i.gphocs.oe --dependency=afterany:${!VARNAME} $SCRIPTS_DIR/gphocs_checkpoint_cont.sh $NT $CHECK_IN $CHECK_OUT)
				declare RUNID_$i=${JID##* }
			fi
		done
	done
done

##################################################
#### 2.2 PROCESS OUTPUT ####
##################################################
REPLICATES=4
BURNIN=200000 # Burn-in
LAST_SAMPLE=2000000 # Last MCMC sample to be considered
MUTRATE=1.236 # Gamma distribution of mutation rate will have mean $MUTRATE * 10e-8
MUTRATE_VAR=0.107 # Gamma distribution of mutation rate will have variance $MUTRATE_VAR * 10e-8
GENTIME=3.5 # Lognormal distribution of generation time will have mean ln($GENTIME)
GENTIME_SD=1.16 # Lognormal distriubtion of generation time will have standard deviation ln($GENTIME_SD)
M_SCALE=0.001 # Inverse scaling factor used in the G-PhoCS configuration file for migration parameter
T_SCALE=10000 # Inverse scaling factor used in the G-PhoCS configuration file for tau and theta
POPLIST=$COAL_DIR/parent_child_pops.txt # File with two columns containing information on parent and child populations (headers: "parent" and "child")

## Process output and convert to demographic values
for j in $(seq 1 $REPLICATES)
do
	for MODEL in $MODELS
	do
		sbatch --account=nib00015 --output=$COAL_DIR/logFiles/$MODEL.run$j.process_logs.oe $SCRIPTS_DIR/process_logs.sh $SCRIPTS_DIR $MODEL $j $COAL_DIR/models/$MODEL/$MODEL.run$j.mcmc.out $BURNIN $LAST_SAMPLE $MUTRATE $MUTRATE_VAR $GENTIME $GENTIME_SD $M_SCALE $T_SCALE $POPLIST_FILE
	done 
done

## Average log files across runs
for MODEL in $MODELS
do
	# Wait until processed logs are ready
	until [[ $(ls $COAL_DIR/models/$MODEL/cut/prep*mcmc.out | wc -l) == 4 ]]
	do
		sleep 5m
	done
	
	# Average
	sbatch --wait --acount=nib00015 --output=$COAL_DIR/logFiles/$MODEL.average_mcmcs.oe $SCRIPTS_DIR/average_mcmcs.sh $SCRIPTS_DIR \
		$COAL_DIR/models/$MODEL/cut/prep.$MODEL.run1.mcmc.out $COAL_DIR/models/$MODEL/cut/prep.$MODEL.run2.mcmc.out $COAL_DIR/models/$MODEL/cut/prep.$MODEL.run3.mcmc.out $COAL_DIR/models/$MODEL/cut/prep.$MODEL.run4.mcmc.out $COAL_DIR/models/$MODEL/cut/prep.$MODEL.average.mcmc.out
done

##################################################
#### 2.3 PLOT RESULTS ####
##################################################
## Plot main figures
SUMMARY=$COAL_DIR/models/$MODEL/cut/summary.xlsx # Excel file containing three sheets (div, mig, gdi) with columns with the following information (column headers are given in brackets): mean estimate (Mean), lower and upper 95 highest posterior density distribution limits (lower95 and upper95), underlying model (Model), and increasing integers (x) 
sbatch --acount=nib00015 --output=$COAL_DIR/logFiles/plot_coal_main.oe $SCRIPTS_DIR/plot_coal_main.sh $SCRIPTS_DIR $COAL_DIR/models/$MODEL/cut/prep.noMig.average.mcmc.out $COAL_DIR/models/$MODEL/cut/prep.Mig.average.mcmc.out $SUMMARY $COAL_DIR/models/$MODEL/cut/

## Plot supplementary figures (i.e., posterior distributions)
M_SCALE=0.001 # Inverse scaling factor used in the G-PhoCS configuration file for migration parameter
T_SCALE=10000 # Inverse scaling factor used in the G-PhoCS configuration file for tau and theta

sbatch --acount=nib00015 --output=$COAL_DIR/logFiles/plot_coal_posteriors.oe $SCRIPTS_DIR/plot_coal_posteriors.sh $SCRIPTS_DIR \
	$COAL_DIR/models/noMig/cut/cut.noMig.run1.mcmc.out $COAL_DIR/models/noMig/cut/cut.noMig.run2.mcmc.out $COAL_DIR/models/noMig/cut/cut.noMig.run3.mcmc.out $COAL_DIR/models/noMig/cut/cut.noMig.run4.mcmc.out \
	$COAL_DIR/models/noMig/cut/cut.Mig.run1.mcmc.out $COAL_DIR/models/noMig/cut/cut.Mig.run2.mcmc.out $COAL_DIR/models/noMig/cut/cut.Mig.run3.mcmc.out $COAL_DIR/models/noMig/cut/cut.Mig.run4.mcmc.out \
	$M_SCALE $T_SCALE $COAL_DIR/models/noMig/cut/
