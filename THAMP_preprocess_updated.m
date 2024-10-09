%% EEG Preprocessing for THAMP
% Jakob pre-Jul 2024
% Arun Jul 2024
%% Hyperparameters
raw_path = "/work/mindlab/Projects/THAMP/EEG Data/raw/"; % path to raw data participant folders /THAMP/EEG_Data/Raw/

save_path = "/work/mindlab/Projects/THAMP/EEG Data/preprocessed/new_preprocess_batch"; % path to preprocessed data participant folders /THAMP/EEG_Data/Preprocessed/

eeglab_path = '/work/mindlab/Programs/eeglab2021.1/';
chanlocs_filepath = [eeglab_path,'plugins/dipfit4.3/standard_BEM/elec/standard_1005.elc']; %chanlocs path

chans_to_load = 1:64;
resample_rate = 500; % SHOULD THIS BE 1000Hz? Raw data seems to be at 5000Hz sampling rate
filter_cutoffs = [0.5 60];
rereference_chans = [10 21];

addpath(eeglab_path, raw_path);

%% determines subjects not yet preprocessed
to_process= {dir(raw_path).name}; to_process = string(to_process(~contains(to_process, ".")));
already_processed = [string({dir(save_path).name}), string({dir(fullfile(save_path, "..")).name})]; already_processed = already_processed(~contains(already_processed, "."));
to_process = setdiff(to_process, already_processed);
disp(to_process);
%% Change participant ID
participantID = "240424YDEN";
mkdir(fullfile(save_path, participantID));

% add eeglab to your matlab path :)

%% 1. Load EEG file
raw_data_path = fullfile(raw_path, participantID);
raw_data_filename = strcat(participantID, '.vhdr');

eeglab redraw; close; % make sure eeglab has loaded % if this doesn't work, make sure that EEGlab is added to your matlab path

EEG = pop_loadbv(raw_data_path, raw_data_filename);
%% 2. Load Channel Locations

EEG = pop_chanedit(EEG, 'lookup', chanlocs_filepath); 

%% 3. Select Channels
EEG = pop_select( EEG, 'channel',chans_to_load);

%% 4. Resample (and save)

EEG = pop_resample( EEG, resample_rate);
preprocess_steps_completed = 'resampled';

EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));

%% 5. Filter (and save)
EEG = pop_eegfiltnew(EEG, 'locutoff',filter_cutoffs(1),'hicutoff',filter_cutoffs(2));
preprocess_steps_completed = 'resampled_filtered';
EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));
%% 6. Rereference to TP9 and TP10 (and save)

EEG = pop_reref(EEG, rereference_chans); % DOUBLE CHECK
preprocess_steps_completed = 'resampled_filtered_reref';
EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));

%% 7. Reject Data w/ channel rejection + ASR
all_chanlocs = EEG.chanlocs; % save chanlocs before cleaning for interpolation
EEG = clean_artifacts(EEG,'ChannelCriterion',0.5,'LineNoiseCriterion',4,'BurstCriterion',10,'WindowCriterion','off','BurstRejection','off', 'ChannelCriterionMaxBadTime', .5);


%% 8. Visualize data
pop_eegplot(EEG, 1, 1, 1); figure; topoplot([],EEG.chanlocs, 'style', 'blanEEk',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);


%% 9. Manually reject bad channels (and save)
chans_to_remove = {'AF3', 'F7'}; % CHANGE THIS
EEG = pop_select( EEG, 'nochannel',chans_to_remove);

preprocess_steps_completed = 'resampled_filtered_reref_cleaned';
EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));

%% 10. Interpolate missing channels (and save)

EEG = pop_interp(EEG, all_chanlocs, 'spherical');
preprocess_steps_completed = 'resampled_filtered_reref_cleaned_interp';
EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));

%% repeat steps 8 through 10 until all bad channels are interpolated

%% 11. ICA 

EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
EEG = pop_iclabel(EEG, 'default');
% pop_viewprops(EEG, 0, 1:28); %% you may need to use the eeglab gui -- tools/Classify Components using ICLabel/Label Components
%%
classes = EEG.etc.ic_classification.ICLabel.classes;

classifications = EEG.etc.ic_classification.ICLabel.classifications;
% find components that are better explained by artifacts rather than brain data
components_to_remove = find(any(classifications(:, ~strcmp(string(classes), 'Brain')) > classifications(:, strcmp(string(classes), 'Brain')), 2));
%% Delete unwanted (any non-brain) components (and save)
EEG = pop_subcomp( EEG, components_to_remove, 0);

preprocess_steps_completed = 'resampled_filtered_reref_cleaned_interp_ica';
EEG = pop_saveset( EEG, 'filename',char(strcat(participantID, '_', preprocess_steps_completed, '.set')),'filepath',char(fullfile(save_path, participantID)));
