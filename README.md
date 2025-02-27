# THAMP-EEG

Scripts for preprocessing and analysis of EEG Data for the THAMP project. More general information can be found on the [MIND Lab WIKI](https://github.com/mind-lab-bos/labwiki/wiki). 

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

- **Toolboxes**
	1. EEGlab
	2. MIRtoolbox1.8.2
</details>

<details><summary>

## Requirements and Dependencies:

</summary>
* MATLAB -- runs on version 2024a
* RStudio -- R Version 4.4.1
* [MIRtoolbox 1.8.2](https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mirtoolbox) (Lartillot & Toiviainen) (installation includes Auditory Toolbox)
* [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) (Delorme & Makeig, 2004)
* MATLAB add-on: raacampbell/shadedErrorBar -- for plotting PLVs.
	* install within MATLAB by searching for the package via "Add-Ons"


</details>

<details><summary>

## Instructions for running examples:

</summary>

1. Verify that toolboxes are installed (see above).
2. Download the `example_EEG_data` and `metadata` folders from the Dropbox link. 	Move both folders into the `THAMP-EEG` folder (the present Github project). 
3. Run `THAMP_PLV.m` section by section in MATLAB. 
	* This script uses the preprocessed example EEG data found in the `example_EEG_data` folder, iterating through 10 participants. You should be able to reproduce a "phase-locking value over normalized frequency" figure and export data to R before running `THAMP_R_EEG.R`.
	* navigate to the THAMP-EEG directory within MATLAB before running the script. 
	* remember to add MIRToolbox and EEGlab to your MATLAB path: `addpath(path/to/toolbox)`
	* if you want to visualize the results and skip running the analysis yourself, you can load in the `all_SART_PLVs.mat` file and move directly to the plotting section. 
	
4. `THAMP_R_EEG.R`
	* This script loads in PLV and RTCV data saved out from the MATLAB analysis for statistical models and plotting with ASRS and eBMRQ. Current data are already saved out so there is no need to run the MATLAB script before this one.
	* You may run into packages that need to be installed through R Studio. Use the command: `install.packages("{package_name}")`
	* Reproduce plots that model PLV+RTCV by ASRS and eBMRQ scores

</details>
