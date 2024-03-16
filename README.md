# THAMP-EEG
Instructions for recruiting/running participants can be found in the THAMP sub-page of the lab wiki. This repo contains instructions and code involved in the preprocessing and analysis of the data.


##Pre-processing
This is a guide for preprocessing using the script: Jakob_THAMP_Preprocessing.m. It is also very simple to do without the script, just click through the eeglab gui for each step.

Move the eeg files from the recording machine to discovery. The files should be stored in /work/mindlab/Projects/THAMP/EEG Data/raw
Download the desired participant folder to your local machine that you will be doing the preprocessing on. (Or on discovery)
If using a local machine ensure that eeglab and all the necessary plugins are installed and added to path. Either download plugins from the internet or copy from this folder /work/mindlab/Programs/eeglab2021.1/plugins. On discovery simply add eeglab to path from /work/mindlab/Programs/eeglab2021.1
Download all of the .m files found HERE or in /work/mindlab/Projects/THAMP/EEG Data. (Or start a virtual matlab session on discovery and open the .m files)
In the Jakob_THAMP_preprocessing change the path, ID, and ID_l variables (lines 19, 23, 24) to the correct folder and participant id that you would like to process.
Also change lines 53, and 62 to point the appropriate local location file, and output directory
Also line 99
Run Jakob_THAMP_Preprocessing. While running there are several manual stops that require input for it to keep running.
The full list of preprocessing steps are as follows
Load channel locations
Change song trigger labels, remove extra triggers.
Filter from .5 Hz to 50 Hz
Re-reference to TP9 and TP10
Resample from 5000 to 500
Automatic Rejection
Manual Rejection
Independent Component Analysis
Reject Components
Analysis Script
Epoching
PSD
Topos Plots
How to update Events
Look at already processed subject txt files for examples.
Export triggers in eeglab as txt. (the code automatically exports this file)
Copy and paste into excel or sheets.
Create a new column on the write.
Use this formula
	=if(AND((H2=H1),NOT(H2=H3)), TRUE, if(AND(H2=H4,not(H2=H5)),TRUE,FALSE))
The 'H' column should be referring to the trigger value column.
The new column header name should be called 'remove' and the values should all be TRUE/FALSE
In the sheet change the onset and offset triggers (invert for unmoved first)
	SART onset (modded first) 	(unmodded first)
	s111 -> s211			s_11 -> s101			
	s112 -> s212			s_12 -> s102	
	s113 -> s213			s_13 -> s103
	s_14 -> s104			s114 -> s214
	s_15 -> s105			s115 -> s215
	s_16 -> s106			s116 -> s216
	NBACK onset (modded first) 	(unmodded first)
	s111				s_21 -> s201
	s112				s_22 -> s202
	s113				s_23 -> s203
	s_24 -> S204			s114 
	s_25 -> S205			s115
	s_26 -> S206			s116
Make sure any extra triggers in the 100's and 200's are changed. (Change them to s240 or s255)
Save the sheet with the replaced triggers and put it in the folder. Save as FLLL_triggers_repl.tsv
Copy the values column into a new txt file and save it as FLLL_type.txt
Copy the remove column into a new txt file and save it as FLLL_remove.txt
Optional: remove the rows containing TRUE if response duplicate triggers need to be deleted
	go to edit/select_epochs_or_events
	in the keep field enter in the Selection box false
	check keep only selected events and remove all others

From this point on the code will be automated again and you can continue
Rejection
After automated rejection, the script will open the scroll plot. Visually inspect all 64 channels and look at the runlog to see if any channels were marked as noisy. Remove any channels that the automated rejection may have missed that seem obviously noisy.
Component Rejection
After ica runs, the script will plot components 1-30 with labels using ICLabel. Manually inspect the components by clicking on them to view more information. Enter the components by number to be rejected in the command window.
Manual Epoching
After pruning ICA components and before running the psd script
tools /extract epochs
Each epoch is 60800 ms long.
So the start value is 0 and end is 60.8 for the epoch limits
Select the correct 12 triggers based on being either modded first or unmodded first. Triggers can be found above in the events section.
Save the new dataset and add ‘_epochs’ to the end
EEG Analysis
Now open the THAMPcalcPSD_1stlevel.m script.
Edit the id variable, songorder variable, task order variable and mod order variable. Look at the thamp runlog for this information. Also change the file path to where you have the processed
Run the script to view the PSD and topos plots for each song. 
Plot additional plots as needed.
Next open THAMPcalcPSD_2ndlvl_bytask.m
Edit the id variable, songorder variable, task order variable and mod order variable for whichever participants are to be included
The code should be commented about what does what
The output are the psd plot with Modded vs unmodded plotted together and the two topos plots
There should be 1 plot for each task for each song so 60 plots total.
The script organizes and takes data from each participant that heard each specific song and averages each psd by song
Each subject has the 12 psd files for each song after the lvl 1 script

