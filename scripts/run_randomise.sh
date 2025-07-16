#!/usr/bin/env bash
# Run FSL randomise to get group statistics across the whole brain 
# V1 10182023
# Add ability to have other regressors 
# Rework completely to be more foolproof
# 07222024

# source the bashrc for access to things
source ~/.bashrc

# set the project location
base_dir='' # ommitted for privacy

# get the inputs 
suffix=$1 # comp, pi, all 
analysis_type=$2 # For instance "FaceObject" or "trialwise_detailed"
contrast_num=$3 # For instance "zstat3" or "difference" for mvpa 
randomise_vers=$4 # default, age, group, or group_mem
mem_type=$5 # what mem type do you want to use if using group_mem ? detailed or recog

# use the info above to choose the subjects 
participants=`cat ${base_dir}/randomise_group_files/${suffix}_${analysis_type}.txt`

# what analysis are you running? 
echo running $randomise_vers randomise for ${analysis_type} on $suffix
echo $participants

# What is the intersect mask being used (assume all)
mask_file=$base_dir/data/intersect_mask.nii.gz # this is in 2mm iso standard space; 91 x 109 x 91

# Make the output folder if it doesn't already exist
mkdir -p $base_dir/data/randomise
mkdir -p $base_dir/data/randomise/${analysis_type}

# Specify the output root names
merged_file=$base_dir/data/randomise/${analysis_type}/${suffix}_merged_data_${contrast_num}.nii.gz 

if [ $randomise_vers == 'group' ]
then
    randomise_file=$base_dir/data/randomise/${analysis_type}/${suffix}_group_${contrast_num}
elif [ $randomise_vers == 'group_mem' ]
then
    randomise_file=$base_dir/data/randomise/${analysis_type}/${suffix}_group_mem_${contrast_num}
elif [ $randomise_vers == 'age' ]
then
    randomise_file=$base_dir/data/randomise/${analysis_type}/${suffix}_age_${contrast_num}
else
    randomise_file=$base_dir/data/randomise/${analysis_type}/${suffix}_${contrast_num}
fi

# Remove in case it exists
#rm -f $merged_file

echo Making $merged_file
echo Outputting $randomise_file
echo Using $mask_file

if [ -e $merged_file ] 
then
	echo merge file exists, continuing to randomise
else

# Cycle through the ppts
for ppt in ${participants}
do

    # if this is a glm output vs. mvpa output, we will have the outputs stored in different places 
    if [ $analysis_type == 'trialwise_detailed' ] || [ $analysis_type == 'trialwise_recognition' ]
    then
        # this is already aligned to standard, so just supply both as the same file 
        file_name=${base_dir}/data/similarity_searchlight/sub-${ppt}_${analysis_type}_${contrast_num}_summed_zscore_v2.nii.gz
        file_name_standard=${file_name}
    else
        secondlevel_path=$base_dir/data/memory_feat_folders/sub-${ppt}/
        file_name=${secondlevel_path}/Memory_${analysis_type}.feat/stats/${contrast_num}.nii.gz
        file_name_standard=${secondlevel_path}/Memory_${analysis_type}.feat/stats/${contrast_num}_registered_standard.nii.gz
    fi
    
   
    if [ -e $file_name ] 
    then      
    
        if [ ! -e $file_name_standard ]
        then
            echo registering zstat for sub-${ppt}
            # create the z-stat file in standard space -- not re-registering, just changing the voxel size 
            flirt -in ${file_name} -ref ${mask_file} -applyxfm -usesqform -out ${file_name_standard}
        fi

        # Either create or merge this file to the list
        if [ ! -e $merged_file ]
        then
            echo Initializing with $file_name_standard
            scp $file_name_standard $merged_file

        else
            fslmerge -t $merged_file $merged_file $file_name_standard
            echo Appending $file_name_standard
        fi
    
    else
        echo zstats not run for sub-${ppt} maybe check why?
    fi

done
fi

# Mask the merged file again so that all of the background values are zero
fslmaths $merged_file -mas $mask_file $merged_file

# Run randomise 
if [ $randomise_vers == 'default' ]
then
    # use the sign test contrast 
    randomise -i $merged_file -o $randomise_file -1 -n 1000 -x -T -t $base_dir/randomise_group_files/signtest.con -C 2.09 

elif [ $randomise_vers == 'group' ]
then
        design_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_group.mat
        contrast_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_group.con
        
        # TFCE analysis does not work with covariates, so save for later correction
        randomise -i $merged_file -o $randomise_file -d $design_file -t $contrast_file -m $mask_file -x -D -C 2.09 -n 1000

elif [ $randomise_vers == 'group_mem' ]
then
        design_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_group_${mem_type}_mem.mat
        contrast_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_group_${mem_type}_mem.con
        
        # TFCE analysis does not work with covariates, so save for later correction
        randomise -i $merged_file -o $randomise_file -d $design_file -t $contrast_file -m $mask_file -x -D -C 2.09 -n 1000
elif [ $randomise_vers == 'age' ]
then
        design_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_age.mat
        contrast_file=$base_dir/randomise_group_files/${suffix}_${analysis_type}_age.con
        
        randomise -i $merged_file -o $randomise_file -d $design_file -t $contrast_file -m $mask_file -x -T -C 2.09 -n 1000        
else
	echo not a recognized randomise version
fi

echo Finished
