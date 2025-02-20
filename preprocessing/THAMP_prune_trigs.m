%% Thamp Preprocess Trigs (run after other preprocessing)
% Arun Jul 2024

%% Filepaths (Update only this section)
participantID = "240215KSOK"; 
preprocessed_data_path = "<path/to/thamp>/EEG_Data/preprocessed/";
thamp_runlog_path = "<path/to/thamp>/THAMP_runlog.xlsx"; % make sure a recent runlog is downloaded. 

% make sure eeglab is in your matlab path :)

%% Load post-ica EEG file
eeg_data_path = char(fullfile(preprocessed_data_path, participantID));

% get ica filename
all_files = dir(eeg_data_path);
filename = char(all_files(endsWith(string({all_files.name}), "ica.set", "IgnoreCase", true)).name);

eeglab redraw; close; % make sure eeglab has loaded % if this doesn't work, make sure that EEGlab is added to your matlab path

EEG = pop_loadset('filename',filename,'filepath',eeg_data_path);

%% Find participant song order

thamp_runlog = readtable(thamp_runlog_path, "Sheet",3, "VariableNamingRule","preserve");
first_task = thamp_runlog.("Block Order (first)")(thamp_runlog.("Participant ID") == participantID);
first_condition = thamp_runlog.("MOD Order (first)")(thamp_runlog.("Participant ID") == participantID);

% trigs by condition groups: each participant listens to 12 songs. there
% are 4 groups given the Task x Mod interaction -> 3 songs per group

trig_order_1 = [111, 112, 113, 14, 15, 16, 111, 112, 113, 24, 25, 26]; % SART+MOD First -> [(SARTMOD, SARTUNMOD), (2BACKMOD, 2BACKUNMOD)]
trig_order_2 = [11, 12, 13, 114, 115, 116, 21, 22, 23, 114, 115, 116]; % SART+UNMOD First -> [(SARTUNMOD, SARTMOD), (2BACKUNMOD, 2BACKMOD)]
trig_order_3 = [111, 112, 113, 24, 25, 26, 111, 112, 113, 14, 15, 16]; % 2BACK+MOD First -> [(2BACKMOD, 2BACKUNMOD), (SARTMOD, SARTUNMOD)]
trig_order_4 = [21, 22, 23, 114, 115, 116, 11, 12, 13, 114, 115, 116]; % 2BACK+UNMOD First -> [(2BACKUNMOD, 2BACKMOD), (SARTUNMOD, SARTMOD)]


switch strcat(string(first_task), string(first_condition))
    case "SARTMOD"
        trig_order = trig_order_1;
    case "SARTUNMOD"
        trig_order = trig_order_2;
    case "2BACKMOD"
        trig_order = trig_order_3;
    case "2BACKUNMOD"
        trig_order = trig_order_4;
    otherwise 
        error("No match found for song conditions")
end

%% Prune Trigs
eeg_trigs = string(regexprep({EEG.event(:).type}, '[A-Za-z ]', '')); % remove alphabet

pruned_trigs = [];

% get index and latencies
for i = 1:length(eeg_trigs)
    if eeg_trigs(i) == ""
        continue
    end
    pruned_trigs = [pruned_trigs; EEG.event(i).bvmknum, EEG.event(i).latency str2double(eeg_trigs(i))];
end

% remove consecutive duplicates
compiled_trigs = [];
index = 1;
while index <= length(pruned_trigs)
    compiled_trigs = [compiled_trigs; pruned_trigs(index,:)];
    if (index +1) > length(pruned_trigs) 
        break 
    end
    if pruned_trigs(index+1,3) == pruned_trigs(index,3) % check for consecutive duplicate
        index = index + 2;
    else
        index = index + 1;
    end
end
%% Find songs (first group)

song_trigs_1 = trig_order(1:6);
song_events_1 = compiled_trigs(find(compiled_trigs(:,3)==song_trigs_1(1),1, 'first'),:) ;
for i = 2:length(song_trigs_1)
    prev_trig_latency = song_events_1(end, 2);
    matched_trigs = compiled_trigs(compiled_trigs(:,3) == song_trigs_1(i), :);
    if isempty(matched_trigs)
        error(strcat("trig not found ", num2str(song_trigs_1(i))));
    end 
    [~, closest_trig] = min(abs(60.8 - (matched_trigs(:,2) - prev_trig_latency)/EEG.srate));
    song_events_1 = [song_events_1; matched_trigs(closest_trig,:)];
end

% check that all song lengths are around 60.8
disp(diff(song_events_1(:,2))/EEG.srate);

%% Find songs (second group)

song_trigs_2 = trig_order(7:12);
song_events_2 = compiled_trigs(find(compiled_trigs(:,3)==song_trigs_2(end),1, 'last'),:) ;
for i = flip(1:length(song_trigs_2)-1)
    prev_trig_latency = song_events_2(end, 2); 
    matched_trigs = compiled_trigs(compiled_trigs(:,3) ==  song_trigs_2(i), :); 
    if isempty(matched_trigs)
        error(strcat("trig not found ", num2str(song_trigs_2(i))));
    end 
    [~, closest_trig] = min(abs(60.8 - (matched_trigs(:,2) - prev_trig_latency)/EEG.srate));
    song_events_2 = [song_events_2; matched_trigs(closest_trig,:)];
end

song_events_2 = flip(song_events_2);

% check that all song lengths are around 60.8
disp(diff(song_events_2(:,2))/EEG.srate);

%% Epoch data

all_song_events = [song_events_1(:,1); song_events_2(:,1)];


EEG=pop_epoch(EEG, {},[0 60.8], 'eventindices', all_song_events');

pop_saveset(EEG, "filename", [erase(filename, ".set") '_epoched.set'], "filepath", eeg_data_path)

