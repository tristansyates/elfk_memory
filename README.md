# Memory encoding and reinstatement in the developing brain  
#### TSY 2105

This folder contains scripts for examining memory behavioral and neural activity in children and adolescents from the Emotional Learning for Kids (ELFK) dataset in the DAN Lab. 

The notebook labeled `scripts/ELFK_Memory.ipynb` contains code that recreates figures and values from the manuscript. 

There were two main neural analyses that were run: a whole-brain univariate GLM analysis to examine activation during memory encoding, and a whole-brain multivariate pattern analysis to look at memory reinstatement during rest periods. 

### Requirements

The following packages are required to run the notebook / Python code:
1. [Numpy](https://numpy.org/) v.1.21.5
2. [Pandas](https://pandas.pydata.org/) v.1.3.5
3. [Matplotlib](https://matplotlib.org/) v.3.5.3
4. [Seaborn](https://seaborn.pydata.org/) v.0.12.2
5. [Scipy](https://scipy.org/) v.1.7.3
6. [Statsmodels](https://www.statsmodels.org/stable/index.html) v.0.13.5
7. [Nibabel](https://nipy.org/nibabel/) v.4.0.2
8. [Nilearn](https://nilearn.github.io/stable/index.html) v.0.10.1
9. [BrainIAK](https://brainiak.org/docs/) v.0.11

For general linear modeling (GLM) analyses, we used [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/) v.6.0.7.9.

For job scheduling, we used [slurm](https://slurm.schedmd.com/documentation.html). Analyses with a high processing load (i.e., memory reinstatement analyses) were run on a high-performance computing cluster for optimal performance, and are located in the `scripts_ginsburg/` folder (rather than the `scripts/` folder, which contains scripts that do not require the benefits of an HPC). 

### Data privacy

The folders labeled `data` and `randomise_group_files` are intentionally empty to preserve data anonymity, as files in these folders are de-identified but link participant IDs to potentially sensitive information (i.e., ages and adversity status). De-identified data files are available upon request.  

### Preprocessing

Behavioral data were transformed from the original E-Prime or Psychopy output files to timing files and summary behavioral memory metrics using Python code. Because this code includes some identifying information (i.e., dates of visits), it has been omitted from the code release, but is available upon request.

Neuroimaging data were preprocessed using [fMRIPrep](https://fmriprep.org/en/stable/) v.23.1.3 as a [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) container with the standard pipeline. The code was launched using `scripts_ginsburg/supervisor_fmriprep.sh` which runs each participant through `scripts_ginsburg/run_fmriprep.sh` and references a BIDs filter file `scripts_ginsburg/fmriprep_bids_filter.json`.

The templates for running the first level analyses using FEAT with FSL are located in the `fsf_templates` folder. These were created manually using the FSL GUI on an example subject, with placeholders inserted for lines in which data would change by the subject ID. These templates run some additional preprocessing steps that are not included within fMRIPrep (e.g., high-pass filtering). The templates used in the include: (1) Memory_FaceObject.fsf.template, which runs the control analysis looking at activation during face-scene and object-scene encoding trials, (2) Memory_Imm_Detail.fsf.template, which runs the main analysis looking at encoding activation for trials that will later be holistically remembered or forgotten, and (3) Memory_Preproc_MVPA.fsf.template, which runs minimal preprocessing on the memory and rest runs to later be used in the multivariate analyses.

### Memory encoding (general linear model analysis)
##### This analysis aims to look at subsequent memory encoding in the brain. In other words, which brain regions are more active during the encoding of trials that will later be remembered vs. forgotten? 

To run the first level FEAT analysis for a given participant, we ran `scripts/supervisor_feat_analysis.sh` with at least two inputs (participant ID number, e.g., '002' and the analysis type, e.g., 'Imm_Detail'), and an optional third input (the task type, e.g., 'task-memory_run-1') to specify which runs to analyse (defaults to the memory run). 

### Memory reinstatement (multivariate pattern analysis) 
##### This analysis aims to look at subsequent memory pattern reinstatement in the brain. In other words, which brain regions show greater pattern similarity between encoding and post-encoding rest (correcting for pre-encoding rest) for trials that will later be remembered vs. forgotten? 

To run the multivariate analysis for a given participant, we ran `scripts_ginsburg/run_similarity_searchlight.sh` with two inputs (participant ID number, e.g., '002' and the analysis type, e.g., 'trailwise_detailed') on a high-performance computing cluster. This script is a wrapper for the Python script (`scripts_ginsburg/Similarity_Searchlight.py`) that runs a whole-brain searchlight analysis.

### Group analyses 

The templates for running higher level analyses using FSL were created by running `scripts/create_randomise_group_files.py` with 3 inputs (suffix for the group type, e.g., 'all'; analysis type, e.g., 'Imm\_Detail'; memory type for what should be used as a covariate, e.g., 'detailed_assoc_imm').

To run the group analyses, we ran `scripts/run_randomise.sh` with at least 4 inputs (suffix for the group type, e.g., 'all'; analysis type, e.g., 'Imm\_Detail'; contrast number from the lower level FEAT, e.g., '3' or contrast type from similarity analysis, e.g., 'difference'; randomise version, e.g., 'default'), and an optional fifth input for the memory type to be used as a regressor if the randomise version is 'group\_mem' (e.g., 'detailed').

These whole brain statistical maps can be corrected for multiple comparisons corrections (using FSL's cluster tool) by running `scripts/run_cluster_correction.sh` with 4 inputs (suffix for the group type, e.g., 'all'; analysis type, e.g., 'Imm\_Detail'; contrast number from the lower level FEAT, e.g., '3' or contrast type from similarity analysis, e.g., 'difference'; contrast from randomise, e.g., 'tstat1').

### Contact

Questions about this repository can be directed to Tristan Yates at tsy2105@columbia.edu (long-term email: tristansyates@gmail.com). 

