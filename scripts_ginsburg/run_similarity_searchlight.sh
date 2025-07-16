#!/bin/sh
# Run a searchlight analysis 
#SBATCH --account=psych
#SBATCH --job-name=sim_searchlight
#SBATCH --output=logs/sim_searchlight-%j
#SBATCH -c 5
#SBATCH --time=3:59:00
#SBATCH --mem=10gb

# get the inputs
sub=$1 # name of the subject
analysis_type=$2 # what type of analysis will you run? trialwise_detailed, trialwise_recognition, etc.  

# activate the conda environment we need (contains the requirements outlined in the README) 
conda activate brainiak

# run the actual script 
python scripts_ginsburg/Similarity_Searchlight.py $sub $analysis_type