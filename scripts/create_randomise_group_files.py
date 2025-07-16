# Create group randomise files without having to open the FSL GUI
# This creates a .txt file of participant names, and .mat and .con files to then provide to FSL randomise 
# These allow the possibility of looking at covariates for age and overall memory performance (z-scored within group) as well
#
# V1 02/17/2025 

# Import a few things 
import numpy as np
import pandas as pd
import os
import sys
import scipy.stats as stats

# set the base dir
base_dir='' # ommitted for privacy

# first read in the participant information (created separately) 
data_df=pd.read_csv('%s/data/memory_fmri_data.csv' %base_dir,index_col=0)

# Take the inputs from the command line 
suffix=sys.argv[1] # 'all' 'comp' 'pi' 
analysis_type=sys.argv[2] # 'FaceObject' 'Imm_Detail' 'Imm_Recognition' 'trialwise_detailed' 'trialwise_recognition'
mem_type=sys.argv[3] #'hit_rates_imm' 'detailed_assoc_imm'

# set a shorter suffix for save files 
if mem_type=='detailed_assoc_imm':
    mem_suffix='detailed' 
else:
    mem_suffix='recog'
    
# what is the value we are looking for in the pi_status column?
if suffix=='comp':
    group_val=0
elif suffix=='pi':
    group_val=1
    
# preset
memory=[]
ages=[]
usable_subs=[]

# cycle through the participants in the data frame
for s, sub_id in enumerate(data_df.sub_ids):

    # if you are getting all subjects, or this subject is in the group
    if suffix=='all' or float(data_df.pi_status[data_df.sub_ids==sub_id])==group_val:
        
        # we know we cannot include this subject for this analysis 
        if sub_id=='sub-134' and mem_suffix=='recog':
            print('skipping %s for having perfect memory' % sub_id)
            continue

        # check whether this is an encoding or reinstatement anlaysis 
        if 'trialwise' in analysis_type:
            inclusion = data_df.included_reinstatement[s]
        else:
            inclusion = data_df.included_encoding[s]
            
        # don't include this subject if they aren't supposed to be
        if inclusion ==0:
            print('not including %s ' % sub_id)
            continue

        # append memory, age, and id 
        memory.append(data_df[mem_type][s])
        ages.append(data_df.ages[s])
        usable_subs.append(sub_id.split('sub-')[-1])

# z-score the values for just this group of subjects 
zscored_mem=stats.zscore(memory)
zscored_age=stats.zscore(ages)

# print how many subjects are included 
print(len(zscored_age))

# First save a list of the participant IDs 
np.savetxt('%s/randomise_group_files/%s_%s.txt' %(base_dir,suffix,analysis_type),usable_subs,fmt='%s')

# Then we will fill out the randomise . mat files for having age as a covariate 
all_text = [] 

# The header includes information about the contrasts --- if you are looking at all subjects, also include the two groups
if suffix=='all':
    file_name='%s/randomise_group_files/%s_%s_group.mat' %(base_dir,suffix,analysis_type)
    # NOTE the ppheights don't matter for running randomise analysis 
    all_text.append('/NumWaves\t3\n/NumPoints\t%d\n/PPheights\n\n/Matrix' % (len(usable_subs)))

else:
    file_name='%s/randomise_group_files/%s_%s_age.mat' %(base_dir,suffix,analysis_type)
    # NOTE the ppheights don't matter for running randomise analysis 
    all_text.append('/NumWaves\t2\n/NumPoints\t%d\n/PPheights\n\n/Matrix' % (len(usable_subs)))

# Then, fill out the values for each subject included 
for s, sub_id in enumerate(usable_subs):

    # again, two columns for each group if you are looking at all subjecst
    if suffix=='all':
        if float(data_df.pi_status[data_df.sub_ids=='sub-%s' % sub_id])==0:
             all_text.append('%0.6f\t %0.6f\t%0.6f' %(1, 0, zscored_age[s])) # columns are comp, pi, age 
        else:
             all_text.append('%0.6f\t %0.6f\t%0.6f' %(0, 1, zscored_age[s])) # columns are comp, pi, age 
    else:
        
         all_text.append('%0.6f\t%0.6f' %(1, zscored_age[s])) # columns are group, age 

# Then save this out 
np.savetxt(file_name,all_text,fmt='%s')

# Also save the contrast files
all_text=[]
if suffix=='all':
    file_name='%s/randomise_group_files/%s_%s_group.con' %(base_dir,suffix,analysis_type)
    # These are the contrasts we care about 
    all_text = ['/ContrastName1\tcomp\n/ContrastName2\tpi',
            '/ContrastName3\tcomp > pi',
           '/ContrastName4\tpi > comp',
           '/ContrastName5\tpositive age',
           '/ContrastName6\tnegative age',
           '/NumWaves\t3',
           '/NumContrasts\t6',
           '/PPheights',
           '/RequiredEffect\n',
           '\n/Matrix',
           '1.000000e+00\t0.000000e+00\t0.000000e+00',
           '0.000000e+00\t1.000000e+00\t0.000000e+00',
           '1.000000e+00\t-1.000000e+00\t0.000000e+00',
           '-1.000000e+00\t1.000000e+00\t0.000000e+00',
           '0.000000e+00\t0.000000e+00\t1.000000e+00',
           '0.000000e+00\t0.000000e+00\t-1.000000e+00']
else:
    file_name='%s/randomise_group_files/%s_%s_age.con' %(base_dir,suffix,analysis_type)
    # These are the contrasts we care about 
    all_text = ['/ContrastName1\t group',
           '/ContrastName2\tpositive age',
           '/ContrastName3\tnegative age',
           '/NumWaves\t2',
           '/NumContrasts\t3',
           '/PPheights',
           '/RequiredEffect\n',
           '\n/Matrix',
           '1.000000e+00\t0.000000e+00',
           '0.000000e+00\t1.000000e+00',
           '0.000000e+00\t-1.000000e+00']
    
# Save it
np.savetxt(file_name,all_text,fmt='%s')

# Next we will fill out the mat files with memory as an additional column, only when we are looking at the contrast of all groups 
all_text=[]
if suffix =='all':

    all_text.append('/NumWaves\t4\n/NumPoints\t%d\n/PPheights\t\n\n/Matrix' % (len(usable_subs)))
    file_name='%s/randomise_group_files/%s_%s_group_%s_mem.mat' %(base_dir,suffix,analysis_type,mem_suffix)
    
    # First 
    for s, sub_id in enumerate(usable_subs):

        if float(data_df.pi_status[data_df.sub_ids=='sub-%s' % sub_id])==0:
            all_text.append('%0.6f\t %0.6f\t%0.6f\t %0.6f' %(1, 0, zscored_age[s],zscored_mem[s]))
        else:
            all_text.append('%0.6f\t %0.6f\t%0.6f\t %0.6f' %(0, 1, zscored_age[s],zscored_mem[s]))

    # Save it
    np.savetxt(file_name,all_text,fmt='%s')           

# also the contrast file 
if suffix=='all':
    file_name='%s/randomise_group_files/%s_%s_group_%s_mem.con' %(base_dir,suffix,analysis_type,mem_suffix)
    # These are the contrasts we care about 
    all_text = ['/ContrastName1\tcomp\n/ContrastName2\tpi',
            '/ContrastName3\tcomp > pi',
           '/ContrastName4\tpi > comp',
           '/ContrastName5\tpositive age',
           '/ContrastName6\tnegative age',
           '/ContrastName7\tmemory',
           '/NumWaves\t4',
           '/NumContrasts\t7',
           '/PPheights',
           '/RequiredEffect',
           '\n/Matrix',
           '1.000000e+00\t0.000000e+00\t0.000000e+00\t0.000000e+00',
           '0.000000e+00\t1.000000e+00\t0.000000e+00\t0.000000e+00',
           '1.000000e+00\t-1.000000e+00\t0.000000e+00\t0.000000e+00',
           '-1.000000e+00\t1.000000e+00\t0.000000e+00\t0.000000e+00',
           '0.000000e+00\t0.000000e+00\t1.000000e+00\t0.000000e+00',
           '0.000000e+00\t0.000000e+00\t-1.000000e+00\t0.000000e+00',
           '0.000000e+00\t0.000000e+00\t0.000000e+00\t1.000000e+00']

    # Save it
    np.savetxt(file_name,all_text,fmt='%s')
