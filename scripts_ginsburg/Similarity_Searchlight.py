# Calculate how much encoding items were reinstated at rest 
# Right now, this is done with trialwise encoded items that were either remembered or forgotten 
# 04082024
# bootstrap resample for z-scored outputs 
# 04182024
# sum instead of average; reverse the thresholding order 
# 04262024
# allow the possibility of also looking at pre- and post-encoding separately
# 07112024

######################################################################################
########## Step 1: Import stuff

# Import a few things 
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import sys 
import scipy.stats as stats
import nibabel as nib
from nibabel.processing import conform
from nilearn import image 
from brainiak.fcma.util import compute_correlation
from brainiak.searchlight.searchlight import Searchlight
from mpi4py import MPI

# Pull out the MPI information, make sure the rank is called rank
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

# define the path and some file names 
base_dir='' # ommitted for privacy 
mvpa_preproc_path='%s/mvpa_preproc_files/' % base_dir
output_path='%s/similarity_searchlight/' % base_dir
timing_path='%s/Ginsburg_Timing_Info/' % base_dir

# add a random seed for the bootstrapping stuff !! 
np.random.seed(0)

######################################################################################
####### Step 2 - load in the data

# which subject are you using?
sub_id=sys.argv[1]
elfk_id='EL%s'% sub_id #elfk ID form
sub_id='sub-%s'% sub_id #subject ID form

# which similarity analysis are you running? 
# persistence (vox x vox similarity post > pre) or trialwise (reinstatement of trials) by memory type 
analysis_type=sys.argv[2] 

# preset
bcvar=[]
data=[]

# first get the encoding trial information if you are running a trialwise analysis 
if 'trialwise' in analysis_type:

    # first get the timing that we care about 
    timing_all=np.loadtxt('%s/%s_trial_timing.txt' %(timing_path,elfk_id))
    bcvar.append(timing_all)
    
    memory_regs=pd.read_csv('%s/%s_memory_regressors.csv' %(timing_path,elfk_id))
    
    # then get the values for the specific type we are looking at (recognition, coarse, detailed)
    memory_reg = memory_regs[analysis_type.split('_')[1]]
    bcvar.append(np.array(memory_reg))
    
else:
    bcvar.append([None]) # Nothing in the first part of the bcvar 
    bcvar.append([None]) # Nothing in the second part of the bcvar 
    
# which tasks are we getting  -- note that one subject had rest run 3 instead of rest run 2
if sub_id == 'sub-141':
	tasks=['memory_run-1','rest_run-1','rest_run-3']
else:
	tasks=['memory_run-1','rest_run-1','rest_run-2']

# only load in the data on rank 0 though 
if rank ==0:
    for task in tasks:
        preproc_file='%s/%s_task-%s_filtered_func_data.nii.gz' %(mvpa_preproc_path, sub_id, task)
        data_4d=nib.load(preproc_file) # load 
        data_4d_zscore = stats.zscore(data_4d.get_fdata(),axis=3) # zscore over time
        data.append(data_4d_zscore) # append

else:
    data += [None]

print("Loaded participant %s \n" % (elfk_id))

# get brain mask (intersect of all participants)
brain_nii=nib.load('/burg/psych/users/elfk/tsy2105/intersect_mask.nii.gz')
affine_mat_standard = brain_nii.affine
dimsize = brain_nii.header.get_zooms()

# Now put it in native space for running the searchlight
brain_nii_native = conform(brain_nii, out_shape =data_4d.shape[:3],
                                    voxel_size = (data_4d.affine[0,0],
                                                  data_4d.affine[1,1],
                                                  data_4d.affine[2,2])) # put in the same shape as native brain 
brain_nii_native = image.binarize_img(brain_nii_native,threshold=0.1) # binarize the image 
affine_mat_native = brain_nii_native.affine

######################################################################################
####### Step 3 - set up the searchlight

mask=brain_nii_native.get_fdata()
sl_rad = 3
max_blk_edge = 5
pool_size = 1

# Create the searchlight object
sl = Searchlight(sl_rad=sl_rad,max_blk_edge=max_blk_edge)

# Distribute the information to the searchlights (note that data is already a list)
sl.distribute(data, mask)

sl.broadcast(bcvar)

######################################################################################
####### Step 4 - define kernel and its internal functions
# Note the most important function is reinstatement_kernel which is passed to the searchlight. It then calls the other functions

def bootstrap_summed_diff_zscore(data_1,data_2,nPerm=1000):
   # bootstrap resampling for the difference in two different data distributions
   # instead of averaging, though, sum over values 
    perm_dist_1=np.zeros(nPerm)
    perm_dist_2=np.zeros(nPerm)
    for perm in range(nPerm):
        # get it for data 1 
        sampidx=np.random.choice(np.arange(len(data_1)),len(data_1),replace=True)
        perm_dist_1[perm] = np.nansum(np.array(data_1)[sampidx])
    
        # get it for data 2
        sampidx=np.random.choice(np.arange(len(data_2)),len(data_2),replace=True)
        perm_dist_2[perm] = np.nansum(np.array(data_2)[sampidx])
        
    # turn it into a zscore = mean of diff minus the null (which is 0) divided by standard dev of the diff
    diff_dist = perm_dist_1 - perm_dist_2
    zscore=(np.mean(diff_dist))/np.std(diff_dist)
   
    return zscore, diff_dist

def convert_encode_data_to_trials(preproc_data,timing_all):
    # turn a vox by TR array into a trial by vox array based on timing information 
    
    all_trials = [] 
    
    for i in range(timing_all.shape[0]):
        # What is the time stamp
        time = timing_all[i,0]

        # What TR does this timepoint refer to?
        TR_start = int(time / 2) # TR length is 2 

        # What TR does this timepoint refer to?
        TR_end = int((time+timing_all[i,1])/2) # TR length is 2 

        # Add the condition label to this timepoint
        trial_pattern = np.nanmean(preproc_data[:,TR_start:TR_end],axis=1)

        # add the trials
        all_trials.append(trial_pattern)

    # then stack 
    all_trials = np.stack(all_trials)

    return all_trials.T


def calculate_trialwise_reinstatement(enc_trialwise_data,rest1_data,rest2_data,memory_reg):
    # correlate the activity pattern of each encoding trial with every possible TR in rest runs
    # calculate how much these patterns were activated on average (reinstated) and then subtract
    # post vs. pre rest run for remembered vs. forgotten encoding trials
    
    # correlate every event trial with every rest TR, using brainiak tools to speed up the process 
    tr_corr_pre = compute_correlation(enc_trialwise_data.T.copy(order='C'),rest1_data.T.copy(order='C'))
    tr_corr_post = compute_correlation(enc_trialwise_data.T.copy(order='C'),rest2_data.T.copy(order='C'))

    # preset 
    zscore_vals=dict()
    distribution_vals =dict()
    for mem_idx, mem_type in enumerate(['forgotten','remembered']):

        # Step 1: threshold the pre-rest by memory behavior for both pre and post rest 
        memory_idxs = np.array(memory_reg)==mem_idx # find the indices 
        
        # threshold the pre rest by memory behavior 
        tr_corr_pre_thresh=tr_corr_pre.copy()
        tr_corr_pre_thresh[~memory_idxs]=np.nan

        # and threshold the post rest by memory behavior 
        tr_corr_post_thresh=tr_corr_post.copy()
        tr_corr_post_thresh[~memory_idxs]=np.nan

        # Step 2: threshold the correlations that are greater than 1.5 sds from the mean
        pre_idxs = tr_corr_pre_thresh > np.nanmean(tr_corr_pre_thresh)+np.nanstd(tr_corr_pre_thresh)*1.5
        tr_corr_pre_thresh[~pre_idxs]=np.nan

        post_idxs = tr_corr_post_thresh > np.nanmean(tr_corr_post_thresh)+np.nanstd(tr_corr_post_thresh)*1.5
        tr_corr_post_thresh[~post_idxs]=np.nan

        # Step 3: turn the post - pre values into z-scores for later analysis 
        zscore, diff_dist = bootstrap_summed_diff_zscore(tr_corr_post_thresh[~np.isnan(tr_corr_post_thresh)],
                                             tr_corr_pre_thresh[~np.isnan(tr_corr_pre_thresh)])
        # save them out 
        zscore_vals[mem_type]=zscore
        distribution_vals[mem_type] = diff_dist
    
    # Step 4: get the difference of remem and forgotten 
    zscore_diff=(np.mean(distribution_vals['remembered'] -
                         distribution_vals['forgotten']))/np.std(distribution_vals['remembered'] - 
                                                                 distribution_vals['forgotten'])

    return [zscore_vals['remembered'], zscore_vals['forgotten'], zscore_diff]

def calculate_trialwise_reinstatement_uncorrected(enc_trialwise_data,rest1_data,rest2_data,memory_reg):
    # correlate the activity pattern of each encoding trial with every possible TR in rest runs
    # calculate how much these patterns were activated on average (reinstated) and then subtract
    # remembered vs. forgotten encoding trials (i.e., don't subtract post vs. pre here) 
    
    # correlate every event trial with every rest TR, using brainiak tools to speed up the process 
    tr_corr_pre = compute_correlation(enc_trialwise_data.T.copy(order='C'),rest1_data.T.copy(order='C'))
    tr_corr_post = compute_correlation(enc_trialwise_data.T.copy(order='C'),rest2_data.T.copy(order='C'))

    # preset 
    pre_vals=dict()
    post_vals=dict()
    for mem_idx, mem_type in enumerate(['forgotten','remembered']):

        # Step 1: threshold the pre-rest by memory behavior for both pre and post rest 
        memory_idxs = np.array(memory_reg)==mem_idx # find the indices 
        
        # threshold the pre rest by memory behavior 
        tr_corr_pre_thresh=tr_corr_pre.copy()
        tr_corr_pre_thresh[~memory_idxs]=np.nan

        # and threshold the post rest by memory behavior 
        tr_corr_post_thresh=tr_corr_post.copy()
        tr_corr_post_thresh[~memory_idxs]=np.nan

        # Step 2: threshold the correlations that are greater than 1.5 sds from the mean
        pre_idxs = tr_corr_pre_thresh > np.nanmean(tr_corr_pre_thresh)+np.nanstd(tr_corr_pre_thresh)*1.5
        tr_corr_pre_thresh[~pre_idxs]=np.nan

        post_idxs = tr_corr_post_thresh > np.nanmean(tr_corr_post_thresh)+np.nanstd(tr_corr_post_thresh)*1.5
        tr_corr_post_thresh[~post_idxs]=np.nan

        # !!! This part differs !!! 
        # Step 3: save out the values 
        pre_vals[mem_type]=tr_corr_pre_thresh[~np.isnan(tr_corr_pre_thresh)]
        post_vals[mem_type]=tr_corr_post_thresh[~np.isnan(tr_corr_post_thresh)]
    
    # Step 4: get the difference of remem and forgotten for post and pre encoding separately
    zscore_diff_pre, _ = bootstrap_summed_diff_zscore(pre_vals['remembered'],pre_vals['forgotten'])
    zscore_diff_post, _ = bootstrap_summed_diff_zscore(post_vals['remembered'],post_vals['forgotten'])

    return [zscore_diff_pre, zscore_diff_post]
    

def reinstatement_kernel(data,sl_mask,myrad,bcvar):
    #'''Searchlight kernel that reshapes the data and decides whether there are enough voxels to run the algorithm'''
    
    timing_all=bcvar[0] # timing of encoding trials during the memory task
    memory_reg=bcvar[1] # memory regressor
    
    # Make sure that we mask the data
    sl_mask_1d=sl_mask.reshape(sl_mask.shape[0] * sl_mask.shape[1] * sl_mask.shape[2]) # 1 dimensional sl mask
    
    #only run this operation if the number of brain voxels is greater than a certain amount
    if np.sum(sl_mask) >= 50:
        
        ## first reshape the data and preprocess it 
        # cycle through the different runs
        preproc_dict=dict()
        for r, run in enumerate(['memory','rest1','rest2']):
             
            # reshape 
            reshaped = data[r].reshape(sl_mask.shape[0] * sl_mask.shape[1] * sl_mask.shape[2], 
                                             data[r].shape[3]).T

            # Mask the data so that we only include data that is inside of the brain 
            # (otherwise we may get weird input from non-brain)
            masked_data=reshaped[:,sl_mask_1d==1]
            
            # change the orientation
            preproc_data = masked_data.T
            
            # add to the dictionary
            preproc_dict[run]=preproc_data
            
            # then if the run is the memory run, make the encode trial data
            if run =='memory':
                enc_trialwise_data = convert_encode_data_to_trials(preproc_data,timing_all)
        
#         # Get the output (post_remem - pre_remem) - (post_forg - post_remem)
        output = calculate_trialwise_reinstatement(enc_trialwise_data,preproc_dict['rest1'],
                                                  preproc_dict['rest2'],memory_reg)
        
        # Get the output (post_remem - post_forg) AND (pre_remem - pre_forg)
#         output = calculate_trialwise_reinstatement_uncorrected(enc_trialwise_data,preproc_dict['rest1'],
#                                                   preproc_dict['rest2'],memory_reg)
        
    else:
        output=np.nan
        
    return output

#############################################################
#############################################################
#####Step 4 - go go go searchlight !!

print("Begin SearchLight in rank %s\n" % rank)

all_sl_result = sl.run_searchlight(reinstatement_kernel,pool_size=pool_size)

print("End SearchLight in rank %s\n" % rank)

# save the data if on rank 0
if rank == 0:
    
    coords = np.where(mask)
    
    all_sl_result = all_sl_result[mask==1]
    all_sl_result = [3*[0] if not n else n for n in all_sl_result] # replace all None
    
    # cycle through remembered vs. forgotten 
    for i, label in enumerate(['remembered','forgotten','difference']):
    
    # cycle through post vs. pre 
    #for i, label in enumerate(['difference_pre','difference_post']):
    
        # get the values and turn them into a double, removing the nans 
        sl_result = [r[i] for r in all_sl_result]
        #print(label,sl_result)
        
        result_vol = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2]))  
        result_vol[coords[0], coords[1], coords[2]] = sl_result   
        result_vol = result_vol.astype('double')
        result_vol[np.isnan(result_vol)] = 0

        # Save the results!
        output_name = '%s/%s_%s_%s_summed_zscore_v2.nii.gz' %(output_path,sub_id,analysis_type,label)
        sl_nii_native = nib.Nifti1Image(result_vol, affine_mat_native)
        sl_nii_standard = conform(sl_nii_native, out_shape =brain_nii.shape,
                                        voxel_size = (affine_mat_standard[0,0],
                                                      affine_mat_standard[1,1],
                                                      affine_mat_standard[2,2]),order=0)

        hdr = sl_nii_standard.header
        hdr.set_zooms((dimsize[0], dimsize[1], dimsize[2]))
        nib.save(sl_nii_standard, output_name)  # Save

    print('Finished searchlight')
    
    
    