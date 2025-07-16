# Script to convert large fmriprep confounds file into the relevant columns to be used for analyses
#
# example command: python scripts/convert_fmriprep_confounds.py 002 memory_run-1
# 
# Re-organized to be its own script 02/15/2025 TY  

import numpy as np
import pandas as pd
import os
import sys

# Set some paths 
base_dir='' # omitted for privacy 
fmriprep_path = '%s/data/fmriprep_memory/' % base_dir

# set some paramenters 
burn_in_TRs=4 # remove this number of TRs before analysis 
motion_threshold=0.9 # motion threshold for main analyses 

# take some inputs 
# first, which participant?
ppt=sys.argv[1] # numbers only, e.g., 002 
elfk_id='EL%s' % ppt # name with EL prefix (used in raw data)
sub_id='sub-%s' % ppt # name with sub- prefix (used in fmriprep) 

# then, which task? 
task=sys.argv[2] # memory_run-1, rest_run-1, rest_run-2

# what is the confound file suffix? 
confounds_suffix = 'task-%s_desc-confounds_timeseries.tsv' % task

# get the full path to the confounds 
# NOTE: for one subject, the rest run after encoding was actually rest run 3, so clarify that
if elfk_id =='EL141' and task=='rest_run-2':
    fmriprep_confs ='%s/%s/func/%s_%s' % (fmriprep_path,sub_id,sub_id,'task-rest_run-3_desc-confounds_timeseries.tsv')    
else:
    fmriprep_confs ='%s/%s/func/%s_%s' % (fmriprep_path,sub_id,sub_id,confounds_suffix)

# load as pandas data frame
all_confs_df = pd.read_csv(fmriprep_confs,delimiter='\t')

# remove 4 volumes of burn-in (as we did for the functionals)
all_confs_df = all_confs_df[burn_in_TRs:]

# get only the FD motion outliers (and not the dvars one) that exceed the threshold 
motion_outlier_columns = [val for val in all_confs_df.columns if 'motion_outlier' in val]
subset_motion_outliers = []
for motion_out in motion_outlier_columns:
    # Get the timepoints where this motion outlier is flagged
    flagged_timepoints = all_confs_df[motion_out] == 1
    if flagged_timepoints.any():
        # Get the FD values for those timepoints
        fd_values = all_confs_df.loc[flagged_timepoints, 'framewise_displacement']
        # Check if any of the FD values exceed the threshold
        if (fd_values > motion_threshold).any():
            subset_motion_outliers.append(motion_out)


#### GLM analysis confound file 
# append to the translation and rotation
motion_columns = np.hstack((['trans_x','trans_y','trans_z',
                 'rot_x','rot_y','rot_z',
                'trans_x_derivative1','trans_y_derivative1','trans_z_derivative1',
                'rot_x_derivative1','rot_y_derivative1','rot_z_derivative1',
                             'white_matter','csf'],subset_motion_outliers))
# make the df and save it
motion_confounds_df = all_confs_df[motion_columns]
np.savetxt('%s/%s/func/%s_task-%s_glm_confounds.txt' % (fmriprep_path,sub_id,sub_id,task),
                np.array(motion_confounds_df),)


#### Similarity analysis (MVPA) confound file 
cosines = [col for col in all_confs_df.columns if 'cosine' in col] ## cosine numbers could vary! 
# append to the translation and rotation
motion_columns = np.hstack((['trans_x','trans_y','trans_z',
                 'rot_x','rot_y','rot_z',
                'trans_x_derivative1','trans_y_derivative1','trans_z_derivative1',
                'rot_x_derivative1','rot_y_derivative1','rot_z_derivative1',
                             'white_matter','csf'],cosines,subset_motion_outliers))

# make the df and save it
motion_confounds_df = all_confs_df[motion_columns]
np.savetxt('%s/%s/func/%s_task-%s_mvpa_confounds.txt' % (fmriprep_path,sub_id,sub_id,task),
                np.array(motion_confounds_df),)
