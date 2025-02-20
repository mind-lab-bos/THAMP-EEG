# THAMP-EEG

Scripts for preprocessing and analysis of EEG Data for the THAMP project. More project-general information can be found on the [MIND Lab WIKI](https://github.com/mind-lab-bos/labwiki/wiki). 

Arun's preprocessing video tutorials are on the MINDLab dropbox

<details><summary>

## File Structure:

</summary>


- README.md

- 'THAMP Overview Document.pdf' -> original project procedure, guidelines, pipelines and data processing.
- **preprocessing**:
	1. THAMP_preprocess_updated.m
	2. THAMP_prune_trigs.m
	3. THAMP_compile_analysis_dir.m
	4. THAMP_standardize_EEG.m
- **analysis**: 
	1. THAMP_PLV_analyses.m
	2. THAMP_PLV_over_time.m
	- **PLV_R**:
		1. THAMP_PLV.m
		2. THAMP_R_EEG.R
- **example_EEG_data**: fully preprocessed data for 10 participants
	- `subID`
		- song_order.csv: 
		-	`finalEEGs`
			- EEG\[1-12\].set
			- EEG\[1-12\].fdt

- **mat_files**:

- **metadata**:
	1. THAMP_eeg_scored_qualtrics.csv
	2. qualtrics.csv
	3. THAMP Song Library.xlsx
- **old_scripts**:
	1. Jakob_THAMP_Preprocessing.m
	2. chanlabels64.m
	3. calcPSD.m
	4. THAMPcalcPSD_1stlvl.m
	5. THAMPcalcPSD_2ndlvl.m
	6. THAMPcalcPSD_2ndlvl_bytask.m

</details>

## Dependencies, Toolboxes and Versions:
* MATLAB -- analyses run on version 2024a
* RStudio
* MIRtoolbox (includes Auditory Toolbox)
* EEGlab
