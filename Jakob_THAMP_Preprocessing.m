clear all
clc
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%access the vhdr file in the path .../newbridge/ID/ID.vhdr
%path = '/Users/Kob/Documents/MINDLab/32_EEG_Preprocessing/newbridge/';
%THAMP path
path = '/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/';
disp(['Current path is set to: ',path]);
%ID = input('Input the Participant ID (YYMMDDFLLL#): ','s');


ID = '231024KKIM';
ID_l = 'kkim';


%mkdir
vhdr = [ID '.vhdr'];
fullpath = [path ID '/']
%vhdr = [insertAfter(ID,6,'_') '.vhdr']


%open vhdr file, use channels 1-64, name file
% CHANGE THE frames per epoch number for each participant??
EEG = pop_loadbv(fullpath,vhdr,[1 4592000],[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',ID,'gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[ID '.set'],'filepath',fullpath);

%%

chanlabels64


%open loc file and get channel locations
locpath = '/Users/Kob/Documents/MINDLab/HAPPE/HAPPE/Packages/eeglab2022.0/plugins/dipfit/standard_BEM/elec/';
loc = 'standard_1005.elc';
disp(['Current loc file path is set to: ',[locpath loc]]);
EEG = pop_chanedit(EEG, 'lookup', [locpath loc]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%EVENTS
pop_expevents(EEG, ['/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/',ID,'/',ID_l,'_events.txt'], 'samples');

%Copy and paste into excel or sheets.
% Create a new column on the write.
% Use this formula
% 	=if(AND((H2=H1),NOT(H2=H3)), TRUE, if(AND(H2=H4,not(H2=H5)),TRUE,FALSE))
% The 'H' column should be referring to the trigger value column.
% The new column header name should be called 'remove' and the values should all be TRUE/FALSE
% In the sheet change the onset and offset triggers (invert for unmoved first)
% 	SART onset (modded first) 	(unmodded first)
% 	s111 -> s211			s_11 -> s101			
% 	s112 -> s212			s_12 -> s102	
% 	s113 -> s213			s_13 -> s103
% 	s_14 -> s104			s114 -> s214
% 	s_15 -> s105			s115 -> s215
% 	s_16 -> s106			s116 -> s216
% 	NBACK onset (modded first) 	(unmodded first)
% 	s111				s_21 -> s201
% 	s112				s_22 -> s202
% 	s113				s_23 -> s203
% 	s_24 -> S204			s114 
% 	s_25 -> S205			s115
% 	s_26 -> S206			s116
% Make sure any extra triggers in the 100's and 200's are changed. (Change them to s240 or s255)
% Save the sheet with the replaced triggers and put it in the folder. Save as FLLL_triggers_repl.tsv
% 
% Copy the values column into a new txt file and save it as FLLL_type.txt
% Copy the remove column into a new txt file and save it as FLLL_remove.txt

check1 = input('Hit enter to unpause after manually removing and renaming triggers')

EEG = eeg_checkset( EEG );
EEG = pop_editeventfield( EEG, 'type',['/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/',ID,'/',ID_l,'_values.txt'],'remove',['/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/',ID,'/',ID_l,'_remove.txt']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);



%save
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[ID '.set'],'filepath',fullpath);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%filter and save
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',60,'plotfreqz',1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[ID '_filt'],'savenew',[fullpath ID '_filt.set'],'gui','off'); 
EEG = eeg_checkset( EEG );

%rereference to TP9 TP10 and save
EEG = pop_reref( EEG, [10 21] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[ID '_filt_reref'],'savenew',[fullpath ID '_filt_reref.set'],'gui','off'); 
EEG = eeg_checkset( EEG );

%resample and save
EEG = pop_resample( EEG, 500);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',[ID '_filt_reref_resamp'],'savenew',[fullpath ID '_filt_reref_resamp.set'],'gui','off'); 
EEG = eeg_checkset( EEG );

%reject automatic
%EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');

%channel rejection Manual
pop_eegplot( EEG, 1, 1, 1);
rej = input('Inspect the scroll plot for channels to manually reject and enter them seperated by a comma ex: "Fp1,Fp2,Cp2"\n','s');
rej2 = split(strip(rej),',');
EEG = pop_select( EEG, 'nochannel',rej2');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',[ID '_filt_reref_resamp_rej'],'savenew',[fullpath ID '_filt_reref_resamp_rej.set'],'gui','off'); 
EEG = eeg_checkset( EEG );


% SAVE DATA WITH DROPPED CHANNELS
%EEG = pop_saveset(EEG, 'filename',strcat(fileN,'_',condition,'_filt_badchan.set'),'filepath', savepath);

%channel interp do this after ICA maybe
%EEG = pop_interp(EEG, ALLEEG(4).chanlocs, 'spherical');
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'setname','230620SSCH3_filt_reref_resamp_rej_interp','savenew','/Users/Kob/Documents/MINDLab/32 EEG Preprocessing/newbridge/230620SSCH3/preprocessed/230620SSCH3_filt_reref_resamp_rej_interp.set','gui','off'); 
% EEG data structure that has all bad channels. Used for interpolating.
%EEG_filt = pop_loadset('filename',strcat(ID,'_',Condition,'_filt.set'),'filepath',EEG.filepath);
EEG_filt = pop_loadset('filename',strcat(ID,'_filt_reref_resamp.set'),'filepath',EEG.filepath);
EEG = pop_interp(EEG, EEG_filt.chanlocs, 'spherical');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname',[ID '_filt_reref_resamp_rej_interp'],'savenew',[fullpath ID '_filt_reref_resamp_rej_interp.set'],'gui','off'); EEG = eeg_checkset( EEG );

%ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

EEG = pop_iclabel(EEG, 'default')
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%pop_topoplot(EEG, 0, [1:29] ,EEG.setname,[5 6] ,0,'electrodes','on');
pop_selectcomps(EEG,1:30)
% gotta enter list of problematic subcomp
%EEG = pop_subcomp( EEG, [1  2  3  4], 0);  example: eyes = 1, muslces 2
%etc
comprej = input('Input the list of components to remove after reviewing (1,4,5,6,7,9,11 ...)/n','s');
comprej = split(comprej,',');
EEG = pop_subcomp( EEG, comprej, 0);
figure; topoplot([],EEG.chanlocs, 'style', 'blanEEk',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);



[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname',[ID '_filt_reref_resamp_rej_interp_prunedICA'],'savenew',[fullpath ID '_filt_reref_resamp_rej_interp_prunedICA.set'],'gui','off'); 
EEG = eeg_checkset( EEG );





