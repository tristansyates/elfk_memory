#!/bin/sh

#SBATCH --account=psych
#SBATCH --job-name=fmriprep
#SBATCH --output=logs/fmriprep-%j
#SBATCH -c 12
#SBATCH --time=11:59:00
#SBATCH --mem-per-cpu=15gb

# This script runs one participant through fMRIPrep

# 1) define the participant object name to be input into the fMRIPrep function
# ID number of subject you want to run; $1 is used when running a batch of participants
# To run one subject, specify the specific ID number after the --participant flag (e.g. sub-PA001)
subj=$1

# 2) load the singularity software
module load singularity

# 3) print which participant fMRIPrep is running on for the purposes of tracking progress / errors
echo "Running fMRIPrep on participant $subj"

base_dir='' # ommitted for privacy 

# 4) Run fmriprep from the singularity container
singularity exec -e \
${base_dir}/fmriprep_23.1.3.sif fmriprep \
${base_dir}/bids/ \
${base_dir}/fmriprep_memory/ \
participant \
--participant_label $subj \
--skip-bids-validation \
--ignore fieldmaps \
--fd-spike-threshold 0.9 \
--bids-filter-file ${base_dir}/scripts_ginsburg/fmriprep_bids_filter.json \
--nthreads 12 \
--random-seed 0 \
--fs-license-file ${base_dir}/fs_license.txt \
--output-spaces MNI152NLin6Asym:res-2  MNI152NLin6Asym \
--skull-strip-t1w force \
--stop-on-first-crash \
--notrack \
--fs-no-reconall \
--work-dir ${base_dir}/work/ \


#################################
#### fMRIPrep basic function ####
#################################

# singularity function structure:
# singularity exec -e fmriprep_latest.sif
# exec = tells singularity software to run the fmriprep singularity image
# -e = tells singularity software to clean the environment before running fmriprep; necessary

# fmriprep function structure: 
# fmriprep bids-root/ output-folder/ participant --participant-label [additional flags]
# participant = keyword that tells fmriprep to run the first-level analysis 

# note: the \ at the end of each line allows the function to continue on different lines;
# 		in reality, it's really just one long line of code

#################################################
#### Basic flags used in fMRIPrep processing ####
#################################################

# --participant_label = specifies one particular subject

# --fd_thres = motion threshold for FD computation
#	set at 0.9mm for task-based ROI / whole brain analysis
#	set at 0.2mm (most stringent threshold, but can try 0.3mm and 0.5mm) for functional connectivity analyses

# --ants-nthreads = number of threads that will be set in ANTs processes; 2 is recommended

# --work-dir = stores intermediate results; need one for each participant, hence the work_$subj

# --no-sub = disables the submission of your study image quality metrics to MRIQC metrics repository

##########################################
#### Potential flags that can be used ####
##########################################

# --version = prints the MRIQC's version number and exits 
# 	sample code: singularity exec -e /moto/edu/psych/files/mriqc_latest.sif mriqc --version

# --session-id = use this flag when you only want mriqc to run on a specific session (e.g. ses-V1W1, ses-V1W2, ses-V2W2)
# 	without this flag, mriqc will process all sessions in the bids dataset as a default
# 	sample flag: --session-d ses-V1W1

# --task-id = use this flag when you only want mriqc to run on a specific task (e.g. matching, poke, REST)
# 	without this flag, mriqc will process all functional tasks in the bids dataset as a default
# 	sample flag: --task-id matching

# --run-id = use this flag when you only want mriqc to run on a specific task run (e.g. run-2, run-3)
# 	without this flag, mriqc will process all runs in the bids dataset as a default
# 	sample flag: --run-id run-1

# --deoblique = deoblique the functional scans during head motion correction preprocessing
#	necessary if fMRI acquisition was oblique (e.g. SB study)

# --hmc-fsl = use FSL MCFLIRT instead of AFNI 3dvolreg for head motion correction (HMC)
#	note: tried to apply this with PACCT and was not recognized as a command

# --ignore
