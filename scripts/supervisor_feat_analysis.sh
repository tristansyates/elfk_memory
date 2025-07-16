#!/bin/bash
#
# Set up the FSF file needed for a given participant and analysis and then try to run feat without opening the GUI
#
# NOTE: this requires first creating an fsf template (inside the folder fsf_templates) with an example subject, saving an FSF file, and then replacing the parameters with placeholder names within that file.  
#
# example command: ./scripts/supervisor_feat_analysis.sh 002 Imm_Detail task-memory_run-1
# 
# V1 TY 10142023
# Updated 10182023

# source bashrc for module loading 
source ~/.bashrc

############################
####### Step 0: Hard code some parameters to use in the analysis IF they are an option 
# (i.e., high pass filtering is automatically turned off for MVPA preproc FSF templates) 
burn_in=4 
high_pass_cutoff=100 

preproc_suffix=_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz # use the version in native space
mask_suffix=space-MNI152NLin6Asym_desc-brain_mask.nii.gz # use the anatomical mask because fmriprep func mask is too strict

############################
####### Step 1: Get information and mask the 4D image we are using 

# take in the inputs 
sub=$1 # who is the subject? (note this should be JUST the number, e.g., 002)
analysis_type=$2 # what is the analysis type? FaceObject, Imm_Detail, Imm_Recognition

# check if they also provided a run type
if [ $# -ge 3 ]
then
	run=$3 # what run do you want to do? task-memory_run-1
else
	run=task-memory_run-1
fi

echo running ${analysis_type} analysis on sub-${sub} for run ${run}

# set these important things 
base_dir='' # ommitted for privacy
timing_dir='' # ommitted for privacy

# make the subject output folder if it doesn't exist
output_folder=${base_dir}/data/memory_feat_folders/sub-${sub}
mkdir -p $output_folder

# add a suffix if needed
if [ $# -ge 3 ]
then
	output_feat=${output_folder}/Memory_${analysis_type}_${run}
else
	output_feat=${output_folder}/Memory_${analysis_type}
fi

############################
####### Step 2: Mask the preproc data so the zstats look pretty in the end 
# Also don't forget to remove the burn-in so we don't have to think about it later 

# first get the subject specific files 
preproc_file=${base_dir}/data/fmriprep_memory/sub-${sub}/func/sub-${sub}_${run}${preproc_suffix}
preproc_file_masked=${preproc_file%%.nii.gz}_masked.nii.gz
mask_file=${base_dir}/data/fmriprep_memory/sub-${sub}/anat/sub-${sub}_${mask_suffix}

# what is the TR number? 
TR_Number=`fslnvols ${preproc_file}`

# Figure out the info for removing the burn in TRs
start_TR="$((burn_in-1))"
TR_Length="$((TR_Number-burn_in))"

# Create a masked (and trimmed) version of the preproc data if it doesn't exist
if [ ! -e ${preproc_file_masked} ]
then 
    echo creating ${preproc_file_masked} from ${mask_file}

    # put the mask in the example func space and bin it 
    flirt -in ${mask_file} -ref ${preproc_file} -applyxfm -usesqform -out ${mask_file%%.nii.gz}_reg_func.nii.gz
    fslmaths ${mask_file%%.nii.gz}_reg_func.nii.gz -bin ${mask_file%%.nii.gz}_reg_func.nii.gz

    # mask the preprocessed data 
    fslmaths ${preproc_file} -mas ${mask_file%%.nii.gz}_reg_func.nii.gz ${preproc_file_masked}
    
    # don't forget to remove the burn in though !!!! 
    echo removing first ${burn_in} TRs from ${preproc_file_masked}
    fslroi ${preproc_file_masked} ${preproc_file_masked} $start_TR $TR_Length
    
else
    echo masked data already created
fi

# what is the TR number? 
Trimmed_TR_Number=`fslnvols ${preproc_file_masked}`

# finally, include the confounds file created separately (e.g., checking if this is for mvpa or not)
if [[ "$analysis_type" == *"MVPA"* ]]
then
	confounds_file=${base_dir}/data/fmriprep_memory/sub-${sub}/func/sub-${sub}_${run}_mvpa_confounds.txt
else
	confounds_file=${base_dir}/data/fmriprep_memory/sub-${sub}/func/sub-${sub}_${run}_glm_confounds.txt
fi

############################
####### Step 3: Prepare the FSF template 
# check whether this has already been run 

if [ ! -e ${output_feat}.feat/filtered_func_data.nii.gz ]
then
    # remove if not 
	rm -rf ${output_feat}*.feat/
        
    # preset 
	fsf_template=${base_dir}/fsf_templates/Memory_${analysis_type}.fsf.template
	fsf_output=${output_folder}/Memory_${analysis_type}.fsf # what is the final output name

    # Put this info in the actual template we will use 
	cat $fsf_template \
	| sed "s:<?= \$FEAT_FOLDER ?>:$output_feat:g" \
	| sed "s:<?= \$TR_NUMBER ?>:$Trimmed_TR_Number:g" \
    | sed "s:<?= \$4D_DATA ?>:$preproc_file_masked:g" \
    | sed "s:<?= \$CONFOUNDS ?>:$confounds_file:g" \
   	| sed "s:<?= \$TIMING_DIR ?>:$timing_dir:g" \
	| sed "s:<?= \$HIGH_PASS_CUTOFF ?>:$high_pass_cutoff:g" \
        > ${fsf_output} #Output to this file
		
	echo Running $fsf_output
	feat $fsf_output
    
else
    echo sub-${sub} has already been analysed for $analysis_type
fi


