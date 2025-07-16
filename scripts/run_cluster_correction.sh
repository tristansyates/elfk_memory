#!/usr/bin/env bash
# Code for manually running cluster correction with FSL

# take the inputs
suffix=$1 # for instance "all" or "comp" or "pi"
analysis_type=$2 # For instance "FaceObject" or "trialwise_detailed"
contrast=$3 # For instance "1" for zstat 1 or "difference" for mvpa 
thresh_type=$4 # For instance, "tstat1" 

# base directory
base_dir='' # ommitted for privacy

# where do randomise files live?
randomise_dir=${base_dir}/data/randomise/

# what is the intersect mask?
mask_file=${base_dir}/data/intersect_mask.nii.gz

# set some criteria for the analysis
pval='0.05' # p < .05
tval='2.41' # equal to p < .01 for an N of 46 

# find the tstat 
tstat_nii=${randomise_dir}/${analysis_type}/${suffix}_${contrast}_${thresh_type}.nii.gz 

echo $tstat_nii

# set the output names 
output_index=${randomise_dir}/${analysis_type}/${suffix}_${contrast}_${thresh_type}_cluster_index
output_thresh=${randomise_dir}/${analysis_type}/${suffix}_${contrast}_cluster_manual_${thresh_type}
output_text=${randomise_dir}/${analysis_type}/${suffix}_${contrast}_cluster_manual_${thresh_type}.txt

# run smoothest
smoothest_output=`smoothest -z ${tstat_nii} -m ${mask_file}`

echo $smoothest_output

keyword='DLH'
dval=$(awk -v key="$keyword" '{while(match($0,key" (.*)",a)){print a[1]; $0=substr($0,RSTART+RLENGTH)}}' <<< "$smoothest_output")

keyword='VOLUME'
volume=$(awk -v key="$keyword" '{while(match($0,key" (.*)",a)){print a[1]; $0=substr($0,RSTART+RLENGTH)}}' <<< "$smoothest_output")

# Finally run cluster correction 
fsl-cluster -i ${tstat_nii} -t ${tval} -p ${pval} -d ${dval} --volume=${volume} --mm -o ${output_index} --othresh=${output_thresh} > ${output_text}

cat ${output_text}
