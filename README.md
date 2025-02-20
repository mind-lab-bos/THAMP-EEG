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

<details><summary>

## Requirements and Dependencies:

* MATLAB -- analyses run on version 2024a
* RStudio
* MIRtoolbox (includes Auditory Toolbox)
* EEGlab

</summary>
</details>

<details><summary>

## Instructions for running examples:

1. Make sure toolboxes are installed (see above).
2. Run `THAMP_PLV.m` section by section in MATLAB. 
	* add MIRToolbox and EEGlab to your MATLAB path: `addpath(path/to/toolbox)`
	* This script uses the example preprocessed EEG data found in the `example_EEG_data` folder, iterating through 10 participants. You should be able to reproduce a phase-locking value over normalized figure and export data for further visualization and modeling in R `THAMP_R_EEG.R`. 
3. `THAMP_R_EEG.R`
	* You may run into packages that need to be installed through R Studio. Use the command: `install.packages("{package_name}")`
	* Reproduce plots that break down mod-PLV by ASRS and eBMRQ scores


</summary>
</details>
