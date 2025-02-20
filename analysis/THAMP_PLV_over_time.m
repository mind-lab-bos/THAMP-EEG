
%% PLV Over Time Analysis
%{
Saved with 61 channels (disgarded Aux1 earlier)
{"240418JDOR", "240621AWOO", "240709JMAR", "240730MANG", "240808CROC",
"240813IITS", "240815BZHE", "240816AJOL", "240820LWU"}
%}

% Saved with 62 channels (from Jakob's script renamed manually ---
% incorrect, should be rewritted with standard BEM chart)
% participantIDs = {"231107XHE","231114CWAY","231117SSAY","231205KKAT","231205NSAZ","231211TABO",...
%     "231212TNAR","231213DROD","240119ACHE","240119AWIL","240207FSUT","240215KSOK","240229LMAR"};


participantIDs = {"231107XHE","231114CWAY","231117SSAY","231205KKAT","231205NSAZ","231211TABO",...
    "231212TNAR","231213DROD","240119ACHE","240119AWIL","240207FSUT","240215KSOK","240229LMAR", "240418JDOR", ...
    "240621AWOO", "240709JMAR", "240730MANG", "240808CROC",...
    "240813IITS", "240815BZHE", "240816AJOL", "240820LWU"};

numParticipants = length(participantIDs);
path_to_data = '/Users/arun/Documents/MINDLab/THAMP/EEG_Data/combined_analyzed/';
song_folder = '/Users/arun/Documents/MINDLab/THAMP/Song Library/normalized_3_21_mp3_4_28';
thamp_song_library = readtable("/Users/arun/Documents/MINDLab/THAMP/THAMP Song Library.xlsx", "NumHeaderLines",0, "VariableNamingRule","preserve");
addpath(song_folder)

numBins = 101; % number of frequencies (frequency resolution)
lowFreq = 0.2;
low_norm = .1;
high_norm = 10.1; % normalization bounds
highFreq = 20.2;
% fois = linspace(lowFreq, highFreq, numBins); % frequency values at which we calculate PLV
norm_fois = linspace(low_norm, high_norm, numBins); % normalized frequency ratios at which we calculate PLV
cycles = 5; % determines resolution of wavelet

num_electrodes = 61;
% aggregateChannelPLVs = zeros(numParticipants,12, numBins, num_electrodes); % participants x songs x fois x electrodes
% mod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, num_electrodes); % participants x songs x fois x electrodes
% unmod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, num_electrodes); % participants x songs x fois x electrodes
% norm_aggregateChannelPLVs = zeros(numParticipants,12, numBins, num_electrodes);
% norm_mod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, num_electrodes);
% norm_unmod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, num_electrodes);

% window length (in secs)
window_length = 10;
hop_length = .5;

% take middle n seconds from song
excerpt_length = 50;

%%
all_norm_PLV_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 12, numBins, num_electrodes);
mod_norm_PLV_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 6, numBins, num_electrodes);
unmod_norm_PLV_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 6, numBins, num_electrodes);
all_alpha_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 12, num_electrodes);
%% Fix electrode locations for (write back to standard BEM)
participants_to_fix = {"231107XHE","231114CWAY","231117SSAY","231205KKAT","231205NSAZ","231211TABO",...
    "231212TNAR","231213DROD","240119ACHE","240119AWIL","240207FSUT","240215KSOK","240229LMAR"};
for participant_idx = 1:length(participants_to_fix)
    participantID = participants_to_fix{participant_idx};
    for song_idx = 1:12
        EEG = pop_loadset('filename',['EEG' num2str(song_idx) '.set'],'filepath', char(fullfile(path_to_data, participantID, "finalEEGs")));
        EEG = pop_saveset( EEG, 'filename',['EEG_badchan' num2str(song_idx) '.set'],'filepath',char(fullfile(path_to_data, participantID, "finalEEGs")));
        EEG.chanlocs = gt_chanlocs;
        EEG = pop_select( EEG, 'nochannel',{'Aux1'});
        EEG = pop_saveset( EEG, 'filename',['EEG' num2str(song_idx) '.set'],'filepath',char(fullfile(path_to_data, participantID, "finalEEGs")));
    end
end



%% Iterate through participants
for participant_idx = 1:numParticipants
    tic
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    mod_idxs = find(strcmp(string([songs.condition]), "Mod"));
    unmod_idxs = find(strcmp(string([songs.condition]), "Unmod"));
    norm_PLV_over_time = zeros(length(1:hop_length:(excerpt_length-window_length+1)), 12, numBins,num_electrodes);
    alpha_over_time = zeros(length(1:hop_length:(excerpt_length-window_length+1)), 12,num_electrodes);
    for song_idx = 1:12
        song_bpm_hz = thamp_song_library.BPM_Hz(strcmp(string(thamp_song_library.Song), regexprep(string(songs(song_idx,:).song), "_", " "))); % get BPM in Hz of song by searching for song name in THAMP_song_library
        norm_freqs = song_bpm_hz * norm_fois;
        switch string(songs(song_idx,:).condition)
            case "Mod"
                song_name = strcat(regexprep(string(songs(song_idx,:).song), "_", " "),".3.mp3");
            case "Unmod"
                song_name = strcat(regexprep(string(songs(song_idx,:).song), "_", " "),".05.mp3");
            otherwise
                error("unrecognized condition")
        end
        disp(strcat("   processing ", song_name));
        %%
        audio = miraudio(char(fullfile(song_folder, song_name)), 'Extract', 0, 60.8, 's'); % trim to 60.8 seconds
        filt = mirfilterbank(audio, 'Gammatone');
        env = mirenvelope(filt, 'Sampling', 500);
        sum_bands = mirsum(env);
        filterbank_audio_data = get(sum_bands,'Data'); filterbank_audio_data = filterbank_audio_data{1}{1};
        filterbank_sr = get(sum_bands, 'Sampling'); filterbank_sr = filterbank_sr{1};
        %%
        EEG = pop_loadset('filename',['EEG' num2str(song_idx) '.set'],'filepath', char(fullfile(path_to_data, participantID, "finalEEGs")));
        audio_signal = filterbank_audio_data;
        audio_fs = 500;
        if EEG.srate ~= audio_fs
            error("sampling rate mismatch")
        end
        %% laplacian
        [surf_lap,G,H] = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
        eegdata = surf_lap';
        %%
                
        % Trim
        if abs(length(eegdata) - length(audio_signal))/audio_fs > 1, error("length mismatch"), end % if the files are more than 1 second different in length
        if length(eegdata) > length(audio_signal)
            eegdata = eegdata(1:length(audio_signal),:);
        else
            audio_signal = audio_signal(1:length(eegdata));
        end
        
        % identify indeces for middle n seconds
        idx_middle = floor(size(eegdata,1)/2); % after trimming, indeces are same for eegdata and audio_signal
        idx_start = idx_middle - excerpt_length/2*audio_fs + 1;
        idx_end = idx_middle + excerpt_length/2*audio_fs;


        %%
        % calc PLVs
        % for freq_idx = 1:numBins
        %     [EEG_phases, ~] = eegConvolution(eegdata,cycles,fois(freq_idx),EEG.srate); % convolve with EEG
        %     [Audio_phases, ~] = eegConvolution(audio_signal,cycles,fois(freq_idx),audio_fs); % convolve with Audio
        %     aggregateChannelPLVs(participant_idx,song_idx, freq_idx,:) = abs(mean(exp(1i*(EEG_phases-Audio_phases)),1));
        % end

        % calc norm_PLVs
        window_idxs = idx_start:hop_length*audio_fs:idx_end-window_length*audio_fs+1;
        alpha_amplitude = abs(hilbert(bandpass(eegdata, [8 12], EEG.srate)));
        for window_idx = 1:length(window_idxs)
            alpha_over_time(window_idx,song_idx,:) = mean(alpha_amplitude(window_idxs(window_idx):window_idxs(window_idx)+window_length*EEG.srate-1,:), 1);
        end

        for freq_idx = 1:numBins
            %%
            [EEG_phases, EEG_amplitudes, EEG_filtered] = eegConvolution(eegdata,cycles,norm_freqs(freq_idx),EEG.srate); % convolve with EEG
            [Audio_phases, ~, Audio_filtered] = eegConvolution(audio_signal,cycles,norm_freqs(freq_idx),audio_fs); % convolve with Audio
            %%
            % plv = zeros(81, 62);
            for window_idx = 1:length(window_idxs)
            % EEG_phases = EEG_phases(idx_start:idx_end,:);
            % Audio_phases = Audio_phases(idx_start:idx_end);
                windowed_EEG_phases = EEG_phases(window_idxs(window_idx):window_idxs(window_idx)+window_length*EEG.srate-1,:);
                windowed_Audio_phases = Audio_phases(window_idxs(window_idx):window_idxs(window_idx)+window_length*audio_fs-1);
                norm_PLV_over_time(window_idx, song_idx, freq_idx, :) = abs(mean(exp(1i* (windowed_EEG_phases-windowed_Audio_phases)), 1));
                % plvs = abs(mean(exp(1i* (windowed_EEG_phases-windowed_Audio_phases)), 1));
                % [surr_m, surr_sd] = surrogate_PLV(windowed_EEG_phases, windowed_Audio_phases, 200);
                % plv(window_idx,:) = (plvs-surr_m)./surr_sd;
                % 
            end
            %%
        end
    
    end

    all_norm_PLV_over_time(participant_idx,:,:, :,:) = norm_PLV_over_time;
    mod_norm_PLV_over_time(participant_idx,:,:, :,:) = norm_PLV_over_time(:,mod_idxs, :,:);
    unmod_norm_PLV_over_time(participant_idx,:,:, :,:) = norm_PLV_over_time(:,unmod_idxs, :,:);
    all_alpha_over_time(participant_idx,:, :,:) = alpha_over_time;
    toc
end

% save("THAMP_PLV_over_time", "all_norm_PLV_over_time", "norm_fois","mod_norm_PLV_over_time", "unmod_norm_PLV_over_time", "all_alpha_over_time", "window_length","hop_length","excerpt_length")
%% surrogate PLV
function [surr_m, surr_sd] = surrogate_PLV(series1, series2, n_surrogate)
    surrogate = zeros(size(series1));
    for i = 1:n_surrogate
        shift = randi([round(size(series1, 1)/4) round(3*size(series1, 1)/4)]);
         
        surrogate(i,:) = abs(mean(exp(1i* (circshift(series1, shift, 1)-series2)), 1));
    end
    surr_m = mean(surrogate, 1);
    surr_sd = std(surrogate, [], 1);
end

%% testing plot
figure;
hold on;
plot(normalize(Audio_filtered(:)))
plot(normalize(EEG_filtered(:, 2)))
plot(window_idxs+250, normalize(norm_PLV_over_time(:,song_idx,freq_idx,10)), "LineWidth", 2)
% plot(window_idxs+250, normalize(plv(:,10)), "LineWidth", 2)
% plot(normalize(Audio_filtered - EEG_filtered(:, 2)))
hold off;
%%
participants = 1:22;

figure;
hold on;
plot(mean(mod_norm_PLV_over_time(participants, :, :, norm_fois==4, :), [3 5])', 'Color','r')
plot(mean(mod_norm_PLV_over_time(participants, :, :, norm_fois==4, :), [1, 3, 5]), 'Color','r', 'LineWidth',2)
plot(mean(unmod_norm_PLV_over_time(participants, :, :, norm_fois==4, :), [3 5])', 'Color','b')
plot(mean(unmod_norm_PLV_over_time(participants, :, :, norm_fois==4, :), [1, 3, 5]), 'Color','b', 'LineWidth',2)
hold off;

figure;
hold on;
plot(norm_fois, squeeze(mean(mod_norm_PLV_over_time(participants, :, :, :, :), [1,2, 3, 5]))', 'Color','r')
plot(norm_fois, squeeze(mean(unmod_norm_PLV_over_time(participants, :, :, :, :), [1, 2, 3, 5]))', 'Color','b')
hold off
%% RTCV Over Time Analysis


% window length (in secs)
window_length = 10;
hop_length = .5;

% take middle n seconds from song
excerpt_length = 50;

mod_SART_RTCV_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
mod_SART_RTCV = zeros(numParticipants, 3);
unmod_SART_RTCV_over_time  = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
unmod_SART_RTCV = zeros(numParticipants, 3);

mod_SART_RT_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
mod_SART_RT = zeros(numParticipants, 3);
unmod_SART_RT_over_time  = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
unmod_SART_RT = zeros(numParticipants, 3);

% Iterate through participants
for participant_idx = 1:numParticipants
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    for song_idx = 1:12
        RTCV_over_time = zeros(length(1:hop_length:(excerpt_length-window_length+1)), 1);
        RT_over_time = zeros(length(1:hop_length:(excerpt_length-window_length+1)), 1);
        if string(songs(song_idx,:).task) ~= "SART"
            continue
        end
        switch string(songs(song_idx,:).condition)
            case "Mod"
                song_name = strcat(regexprep(string(songs(song_idx,:).song), "_", " "),".3.mp3");
            case "Unmod"
                song_name = strcat(regexprep(string(songs(song_idx,:).song), "_", " "),".05.mp3");
            otherwise
                error("unrecognized condition")
        end
        disp(strcat("   processing ", song_name));
        
        EEG = pop_loadset('filename',['EEG' num2str(song_idx) '.set'],'filepath', char(fullfile(path_to_data, participantID, "finalEEGs")));
        
        if EEG.srate ~= 500
            error("sampling rate mismatch")
        end
        
        % prune EEG events
        event_numbers = str2double(string(regexprep({EEG.event(:).type}, '[A-Za-z ]', ''))); % remove letters
        event_latencies = [EEG.event.latency]/500; % in seconds % makes the assumption that units of EEG.event.latency is samples
        
        if abs(60.8 - max(event_latencies)) > 5
            error("event latencies weren't calculated correctly")
        end

        % remove consecutive duplicates
        compiled_trigs = [];
        index = 1;
        while index <= length(EEG.event)
            compiled_trigs = [compiled_trigs; event_numbers(index) event_latencies(index)];
            if (index +1) > length(EEG.event)
                break 
            end
            if event_numbers(index+1) == event_numbers(index) % check for consecutive duplicate
                index = index + 2;
            else
                index = index + 1;
            end
        end
        %% calc RTCV for whole song
        RTs = [];
        data = compiled_trigs;
        for trial_idx = 1:size(data,1)-1
            if data(trial_idx, 1) ~= 3 && data(trial_idx, 1) ~= 10 && data(trial_idx+1, 1) == 10
                RTs = [RTs; data(trial_idx+1, 2) - data(trial_idx, 2)];
            end
        end
        RT = mean(RTs);
        RTCV = std(RTs)/mean(RTs);

        %%


        % identify indeces for middle n seconds
        idx_middle = floor(size(EEG.data,2)/2)/EEG.srate; % in seconds for this analysis
        idx_start = idx_middle - excerpt_length/2;
        idx_end = idx_middle + excerpt_length/2;


        % calc RTCVs
        window_idxs = idx_start:hop_length:idx_end-window_length;

        for window_idx = 1:length(window_idxs)
            data_range = (compiled_trigs(:, 2) >= window_idxs(window_idx)) & (compiled_trigs(:, 2) < window_idxs(window_idx)+window_length);
            data = compiled_trigs(data_range, :);
            RTs = [];
            for trial_idx = 1:size(data,1)-1
                if data(trial_idx, 1) ~= 3 && data(trial_idx, 1) ~= 10 && data(trial_idx+1, 1) == 10
                    RTs = [RTs; data(trial_idx+1, 2) - data(trial_idx, 2)];
                end
            end
            RTCV_over_time(window_idx) = std(RTs)/mean(RTs);
            RT_over_time(window_idx) = mean(RTs);
        end

        switch string(songs(song_idx,:).condition)
            case "Mod"
                mod_SART_RTCV_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RTCV_over_time;
                mod_SART_RT_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RT_over_time;
                mod_SART_RTCV(participant_idx, mod(song_idx-1, 3)+1) = RTCV;
                mod_SART_RT(participant_idx, mod(song_idx-1, 3)+1) = RT;
            case "Unmod"
                unmod_SART_RTCV_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RTCV_over_time;
                unmod_SART_RT_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RT_over_time;
                unmod_SART_RTCV(participant_idx, mod(song_idx-1, 3)+1) = RTCV;
                unmod_SART_RT(participant_idx, mod(song_idx-1, 3)+1) = RT;
            otherwise
                error("unrecognized condition")
        end
    
    end
end
%%
% save("THAMP_RTCV_over_time_EEG", "mod_SART_RTCV_over_time", "unmod_SART_RTCV_over_time", "window_length","hop_length","excerpt_length")
%% get song order and group EEG files
mod_SART_songs = zeros(numParticipants, 3);
unmod_SART_songs = zeros(numParticipants, 3);
mod_SART_PLV_over_time = zeros(size(all_norm_PLV_over_time,1), size(all_norm_PLV_over_time, 2), 3, size(all_norm_PLV_over_time, 4), size(all_norm_PLV_over_time, 5));
unmod_SART_PLV_over_time = zeros(size(all_norm_PLV_over_time,1), size(all_norm_PLV_over_time, 2), 3, size(all_norm_PLV_over_time, 4), size(all_norm_PLV_over_time, 5));
mod_alpha_over_time = zeros(size(all_alpha_over_time,1), size(all_alpha_over_time, 2), 3, size(all_alpha_over_time, 4));
unmod_alpha_over_time = zeros(size(all_alpha_over_time,1), size(all_alpha_over_time, 2), 3, size(all_alpha_over_time, 4));
for participant_idx = 1:numParticipants
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    mod_SART_songs(participant_idx,:) = find(strcmp(songs.condition, "Mod") & strcmp(songs.task, "SART"));
    unmod_SART_songs(participant_idx,:) = find(strcmp(songs.condition, "Unmod") & strcmp(songs.task, "SART"));
    mod_SART_PLV_over_time(participant_idx, :,:,:,:) = squeeze(all_norm_PLV_over_time(participant_idx,:,mod_SART_songs(participant_idx,:),:,:));
    unmod_SART_PLV_over_time(participant_idx, :,:,:,:) = squeeze(all_norm_PLV_over_time(participant_idx,:,unmod_SART_songs(participant_idx,:),:,:));
    mod_alpha_over_time(participant_idx, :,:,:) = squeeze(all_alpha_over_time(participant_idx,:,mod_SART_songs(participant_idx,:),:));
    unmod_alpha_over_time(participant_idx, :,:,:) = squeeze(all_alpha_over_time(participant_idx,:,unmod_SART_songs(participant_idx,:),:));
end
%% RTCV PLV over time analysis

% save to R:
% mod_PLV_over_time = squeeze(mod_SART_PLV_over_time(:,:,:,norm_fois==4,:));
% unmod_PLV_over_time = squeeze(unmod_SART_PLV_over_time(:,:,:,norm_fois==4,:));
% save("PLV_over_time_to_R", "mod_PLV_over_time", "unmod_PLV_over_time", "mod_SART_RT_over_time", "unmod_SART_RT_over_time", "mod_SART_RTCV_over_time", "unmod_SART_RTCV_over_time")

%%
% participants = 1:13
figure; hold on;
plot(squeeze(mean(mod_SART_PLV_over_time(participants, :,song_number,  :, electrodes), [1, 2, 3, 5])))
plot(squeeze(mean(unmod_SART_PLV_over_time(participants, :,song_number,  :, electrodes), [1, 2, 3, 5])))
hold off;
%% one person (example)
participants = 1
song_number = 1:3
foi = norm_fois == 4;
electrodes = 1:30

figure;
tiledlayout('flow')
nexttile;
hold on;
plot(normalize([mean(mod_SART_PLV_over_time(participants, :,1,  foi, electrodes), [1, 3, 4, 5]), ...
    mean(mod_SART_PLV_over_time(participants, :,2,  foi, electrodes), [1, 3, 4, 5]),...
    mean(mod_SART_PLV_over_time(participants, :,3,  foi, electrodes), [1, 3, 4, 5])
    ]),"LineWidth",2, "DisplayName", "Phase-locking value (at mod rate)");

% plot(normalize([mean(mod_alpha_over_time(participants, :,1, electrodes), [1, 3, 4]), ...
%     mean(mod_alpha_over_time(participants, :,2, electrodes), [1, 3, 4]),...
%     mean(mod_alpha_over_time(participants, :,3, electrodes), [1, 3, 4])
%     ]));
% 
plot(normalize([mean(mod_SART_RT_over_time(participants, :,1), [1, 3]), ...
    mean(mod_SART_RT_over_time(participants, :,2), [1, 3]), ...
    mean(mod_SART_RT_over_time(participants, :,3), [1, 3])...
    ]),"LineWidth",2, "DisplayName", "Reaction Time");
plot(normalize([mean(mod_SART_RTCV_over_time(participants, :,1), [1, 3]), ...
    mean(mod_SART_RTCV_over_time(participants, :,2), [1, 3]), ...
    mean(mod_SART_RTCV_over_time(participants, :,3), [1, 3])...
    ]),"LineWidth",2, "DisplayName", "Reaction Time Variability");
yline(0, "HandleVisibility","off")
xline(82, "HandleVisibility","off")
xline(81+82, "HandleVisibility","off")
title("Example Timeseries of Variables of Interest")
% plot(normalize(mean(mod_SART_RTCV_over_time(participants, :,song_number), [1, 3])));
% 
% plot(normalize(mean(mod_SART_RT_over_time(participants, :,song_number), [1, 3])));
% plot(mean(unmod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [1, 3, 4, 5]));
hold off;
xlabel("Time [windows]")
ylabel("Normalized value [a.u.]")
legend
fontsize(16, "points")
set(gcf, "color", "w")


% ylim([0 .4])
% model=fitlm(normalize(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [1, 3, 4, 5])), normalize(mean(mod_SART_RT_over_time(participants, :,song_number), [1, 3])))

% nexttile;
% scatter(...
%     reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[]),...
%         (reshape(mod_SART_RTCV_over_time(participants,:,song_number), 1, []))...
%     )
% xlim([0 .3])
% plot(mean(unmod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [1, 3, 4, 5]));
% lsline

%% scatter


participants = 1
song_number = 1
foi = norm_fois == 4;
electrodes = 1
figure;
hold on;
plot(normalize(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [1, 3, 4, 5])));
plot(normalize(mean(mod_SART_RT_over_time(participants, :,song_number), [1, 3])));
% scatter(...
%     log10(reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[])),...
%         log10(reshape(mod_SART_RTCV_over_time(participants,:,song_number), 1, []))...
%     )
% [r, p] = corrcoef(...
%     log10(reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[])),...
%         log10(reshape(mod_SART_RT_over_time(participants,:,song_number), 1, []))...
%         )
% plot(mean(unmod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [1, 3, 4, 5]));
lsline
hold off;


%%
all_Rs = zeros(22, 3, 61);
all_Ps = zeros(22, 3, 61);
for participants = 1:22
    for song_number = 1:3
        for electrodes = 1:61
            % series1 = mean(mod_alpha_over_time(participants, :,song_number, electrodes), 4)';
            series1 = detrend(mean(unmod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5])');
            series2 = detrend(unmod_SART_RTCV_over_time(participants,:,song_number)');
            r = corr(series1, series2, "Type","Pearson");
            % [p, S] = polyfit(series1,series2,3);
            % r = S.rsquared;

            % 
            % surrogate_rs = surrogate_shifting(series1, series2, 200);
            % % r = sum(r > surrogate_rs)/length(surrogate_rs);
            % % 
            % [surr_m, surr_sd] = deal(mean(surrogate_rs), std(surrogate_rs));
            % r = (r-surr_m)/surr_sd;

        % [r, p] = corr(log10(reshape(mod_SART_RTCV_over_time(participants,:,song_number), 1, []))', ...
        %     log10(reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[]))',...
        %     "Type","Pearson"...
        %     );
        % [r, p] = corr(log10(reshape(mod_SART_RT_over_time(participants,:,song_number), 1, []))', ...
        %     log10(reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[]))',...
        %     "Type","Pearson"...
        %     )

        % [r, p] = corr(log10(reshape(mean(mod_alpha_over_time(participants, :,song_number, electrodes), [3, 4]), 1,[]))', ...
        %     log10(reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[]))',...
        %     "Type","Pearson"...
        %     );

        % [r, p] = corr(log10(reshape(mean(mod_alpha_over_time(participants, :,song_number, electrodes), [3, 4]), 1,[]))', ...
        %     log10(reshape(mod_SART_RT_over_time(participants,:,song_number), 1, []))',...
        %     "Type","Pearson"...
        %     )

        % [p, R] = polyfit(log10(reshape(mod_SART_RTCV_over_time(participants,:,song_number), 1, [])), ...
        %         reshape(mean(mod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[]), 2 ...
        %     )
        
        % r = R.rsquared;
        % [r, p] = corrcoef(log10(reshape(unmod_SART_RTCV_over_time(participants,:,song_number), 1, [])), ...
        %     reshape(mean(unmod_SART_PLV_over_time(participants, :,song_number,  foi, electrodes), [4, 5]), 1,[])...
        %     )

            % all_Ps(participants, song_number, electrodes) = p;
            all_Rs(participants, song_number, electrodes) = r;
        end
    end
end
%% surrogate
function surrogate_rs = surrogate_shifting(series1, series2, n_surrogate)
    surrogate_rs = zeros(n_surrogate, 1);
    
    for i = 1:n_surrogate
        shift = randi([round(length(series1)/4) round(3*length(series1)/4)]);
        surrogate_rs(i) = corr(circshift(series1, shift,1), series2, "Type","Pearson");
    end
    % [m, sd] = deal(mean(surrogate_rs), std(surrogate_rs));
end
%% Topo
% figure;
% topoplot(1- squeeze(median(all_Ps, [1, 2])), EEG.chanlocs);
% clim([.9 1])
% colormap("jet")
topos = load("/Users/arun/Documents/MINDLab/THAMP/EEG_Data/PLV_R/topos_out.mat").topos
figure;
topoplot(topos, EEG.chanlocs);
% clim([-1 1])
clim([-2 2])
% clim([0 .2])
c = colorbar
ylabel(c, "T-stat")
set(gcf, "color", "w")
title(["RTCV ~ Phase-locking over time, Unmod","T-Statistic"]) , %, "Permutation Shuffled Z-Score"]) % 
fontsize(20, "points")
%%
figure;
topoplot(squeeze(mean(all_Rs, [1, 2])), EEG.chanlocs);
clim([-1 1])
% clim([-.2 .2])
% clim([0 .2])
colorbar
set(gcf, "color", "w")
title("Phase-locking Reaction Time Correlations")
fontsize(16, "points")


%% Fischer's Z

z = mean(1/2 * log((1+all_Rs) ./ (1-all_Rs)));
R = (exp(2*z) -1)/(exp(2*z)+1)
%%
 % participants x songs x fois x electrodes
plot(fois, squeeze(mean(mod_aggregateChannelPLVs(:, :, :, :), [1,2, 4])), 'Color', 'b')
hold on;
plot(fois, squeeze(mean(unmod_aggregateChannelPLVs(:, :, :, :), [1,2 4])), 'Color','r')
hold off;
%% norm PLV
plot(norm_fois, squeeze(mean(norm_mod_aggregateChannelPLVs(:, :, :, :), [1,2, 4])), 'Color', 'b')
hold on;
plot(norm_fois, squeeze(mean(norm_unmod_aggregateChannelPLVs(:, :, :, :), [1,2 4])), 'Color','r')
hold off;
xlabel("multiple of beat rate")
ylabel("PLV")
%% shaded error plot


figure;
hold on;
% 
% hp = patch([8 8 9.6 9.6],[0 .14 .14 0],'k',...
%     'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs, [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs, [2 4]))/sqrt(size(norm_mod_aggregateChannelPLVs, 1)),'lineprops',{'Color',[.8 0 0],'LineWidth',4,'DisplayName',"Mod"})
shadedErrorBar(norm_fois,mean(norm_unmod_aggregateChannelPLVs, [1 2 4]),std(mean(norm_unmod_aggregateChannelPLVs, [2 4]))/sqrt(size(norm_mod_aggregateChannelPLVs, 1)),'lineprops',{'Color',[.7 .4 0],'LineWidth',4,'DisplayName',"Unmod"})
xline(4, "-", "mod rate", "LineWidth",1.5, "Color", "k", "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment","left", "HandleVisibility","off");
xline([.5 1 2], "--", "LineWidth",1, "HandleVisibility","off");
hold off
legend
xlim([0 10.2])
xticks([.5 1:10])
xlabel("multiple of beat rate")
ylabel("Phase-locking value")
title("Phase-locking Values at Multiples of the Beat Rate")
set(gcf, "color", "w")

fontsize(16, "points")

% print("THAMP_PLV.png", "-dpng", "-r500")

%%
freq_range = [4 4];
topoSlice = squeeze(mean(norm_mod_aggregateChannelPLVs(:,:,norm_fois>=freq_range(1) & norm_fois<=freq_range(2),:),[1 2 3]));
figure;
topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([.02 .17])
set(gcf, "color", "w")
colorbar

%% Plot Corrected



figure;
hold on;

shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs, [1 2 4]) - mean(noise_norm_mod_aggregateChannelPLVs, [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs, [2 4]))/sqrt(size(norm_mod_aggregateChannelPLVs, 1)),'lineprops',{'Color',[.8 0 0],'LineWidth',4,'DisplayName',"Mod"})
shadedErrorBar(norm_fois,mean(norm_unmod_aggregateChannelPLVs, [1 2 4]) - mean(noise_norm_unmod_aggregateChannelPLVs, [1 2 4]),std(mean(norm_unmod_aggregateChannelPLVs, [2 4]))/sqrt(size(norm_mod_aggregateChannelPLVs, 1)),'lineprops',{'Color',[.7 .4 0],'LineWidth',4,'DisplayName',"Unmod"})

xline(4, "-", "mod rate", "LineWidth",1.5, "Color", "k", "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment","left", "HandleVisibility","off");
xline([.5 1 2], "--", "LineWidth",1, "HandleVisibility","off");
yline(0, "k-", "LineWidth",1, "HandleVisibility","off");
hold off
legend
xlim([0 10.2])
xticks([.5 1:10])
xlabel("multiple of beat rate")
ylabel("Noise Subtracted Phase-locking value")
title("Noise Subtracted Phase-locking Values at Multiples of the Beat Rate")
set(gcf, "color", "w")

fontsize(16, "points")



 % participants x songs x fois x electrodes

%% Plot asrs vs not asrs
% [1, 1, x, 0, 0, 0, 0, 0, 0]
asrs_idxs = find([1, 1, 0, 0, 0, 0, 0, 0, 0]);
noasrs_idxs = find([0, 0, 0, 1, 1, 1, 1, 1, 1]);

% asrs_idxs = find([1, 0, 1, 0, 0, 0,1,0,0,1,1,0,1]);
% noasrs_idxs = find(~[1, 0, 1, 0, 0, 0,1,0,0,1,1,0,1]);


figure;
hold on;

shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs(asrs_idxs,:,:,:), [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs, [2 4]))/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .5 .1],'LineWidth',4,'DisplayName',"ASRS"})
shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs(noasrs_idxs,:,:,:), [1 2 4]),std(mean(norm_unmod_aggregateChannelPLVs, [2 4]))/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.1 .4 .6],'LineWidth',4,'DisplayName',"No ASRS"})

xline(4, "-", "mod rate", "LineWidth",1.5, "Color", "k", "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment","left", "HandleVisibility","off");
xline([.5 1 2], "--", "LineWidth",1, "HandleVisibility","off");
yline(0, "k-", "LineWidth",1, "HandleVisibility","off");
hold off
legend
xlim([0 10.2])
xticks([.5 1:10])
xlabel("multiple of beat rate")
ylabel("Phase-locking value")
title("Phase-locking Values at Multiples of the Beat Rate")
set(gcf, "color", "w")

fontsize(16, "points")



%% PLV over songs shaded error plot

% participant_idxs = find([1, 0, 1, 0, 0, 0,1,0,0,1,1,0,1]);
participant_idxs = find([1, 1, 0, 0, 0, 0, 0, 0, 0]);
figure;
hold on;
% 
% hp = patch([8 8 9.6 9.6],[0 .14 .14 0],'k',...
%     'facecolor',[.5 .5 .5],'facealpha', 0.2,'edgecolor','none', 'HandleVisibility', 'off') ;
shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs(participant_idxs,[3, 6], :, :), [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs(participant_idxs,[3, 6], :, :), [2 4]))/sqrt(length(participant_idxs)),'lineprops',{'Color',[.8 0 0],'LineWidth',4,'DisplayName',"Mod 3"})
shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs(participant_idxs, [2, 5], :,:), [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs(participant_idxs, [2, 5], :,:), [2 4]))/sqrt(length(participant_idxs)),'lineprops',{'Color',[.6 0 0],'LineWidth',4,'DisplayName',"Mod 2"})
shadedErrorBar(norm_fois,mean(norm_mod_aggregateChannelPLVs(participant_idxs, [1, 4], :,:), [1 2 4]),std(mean(norm_mod_aggregateChannelPLVs(participant_idxs, [1, 4], :,:), [2 4]))/sqrt(length(participant_idxs)),'lineprops',{'Color',[.4 0 0],'LineWidth',4,'DisplayName',"Mod 1"})
% shadedErrorBar(norm_fois,mean(norm_unmod_aggregateChannelPLVs, [1 2 4]),std(mean(norm_unmod_aggregateChannelPLVs, [2 4]))/sqrt(size(norm_mod_aggregateChannelPLVs, 1)),'lineprops',{'Color',[.7 .4 0],'LineWidth',4,'DisplayName',"Unmod"})
xline(4, "-", "mod rate", "LineWidth",1.5, "Color", "k", "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment","left", "HandleVisibility","off");
xline([.5 1 2], "--", "LineWidth",1, "HandleVisibility","off");
hold off
legend
xlim([0 10.2])
xticks([.5 1:10])
xlabel("multiple of beat rate")
ylabel("Phase-locking value")
title("Phase-locking Values at Multiples of the Beat Rate")
set(gcf, "color", "w")

fontsize(16, "points")
%% Aggregate 2 groups of PLV
% 
% mod_PLVs1 = mean(norm_mod_aggregateChannelPLVs, 4); unmod_PLVs1 = mean(norm_unmod_aggregateChannelPLVs, 4); 
% all_PLVs1 = mean(norm_aggregateChannelPLVs, 4);
% 
% mod_PLVs2 = mean(norm_mod_aggregateChannelPLVs, 4); unmod_PLVs2 = mean(norm_unmod_aggregateChannelPLVs, 4);
% all_PLVs2 = mean(norm_aggregateChannelPLVs, 4);
% 
mod_PLVs = [mod_PLVs1; mod_PLVs2];
unmod_PLVs = [unmod_PLVs1; unmod_PLVs2];
all_PLVs = [all_PLVs1; all_PLVs2];

%% Find sart conditions
participants = {"240418JDOR", "240621AWOO", "240709JMAR", "240730MANG", "240808CROC", "240813IITS", "240815BZHE", "240816AJOL", "240820LWU",...
    "231107XHE","231114CWAY","231117SSAY","231205KKAT","231205NSAZ","231211TABO",...
    "231212TNAR","231213DROD","240119ACHE","240119AWIL","240207FSUT","240215KSOK","240229LMAR"};
path_to_data = '/Users/arun/Documents/MINDLab/THAMP/EEG_Data/combined_analyzed/';
sart_mod_songs = zeros(3, length(participants));
sart_unmod_songs = zeros(3, length(participants));
for participant_idx = 1:length(participants)
    songs = readtable(fullfile(path_to_data, participants{participant_idx}, "song_order.csv"));
    sart_mod_songs(:, participant_idx) = find(songs.condition == "Mod" & songs.task == "SART");
    sart_unmod_songs(:, participant_idx) = find(songs.condition == "Unmod" & songs.task == "SART");
end

%% aggregate PLVs over different electrode numbers
sart_mod_PLV = zeros(size(all_PLVs,1), 3, size(all_PLVs,3));
sart_unmod_PLV = zeros(size(all_PLVs,1), 3, size(all_PLVs,3));

for idx = 1:size(all_PLVs,1)
    sart_mod_PLV(idx,:,:) = all_PLVs(idx,sart_mod_songs(:,idx), :);
    sart_unmod_PLV(idx,:,:) = all_PLVs(idx,sart_unmod_songs(:,idx), :);
end


%%
asrs_idxs = find([[1, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0,1,0,0,1,1,0,1]]);

noasrs_idxs = find([[0, 0, 0, 1, 1, 1, 1, 1, 1], ~[1, 0, 1, 0, 0, 0,1,0,0,1,1,0,1]]);
% 
% asrs_idxs = 1:22;



figure;
hold on;
% shadedErrorBar(norm_fois,mean(mod_PLVs(:,:,:), [1 2]),std(mean(mod_PLVs(:,:,:), 2))/sqrt(size(mod_PLVs,1)),'lineprops',{'Color',[.7 .1 .1],'LineWidth',4,'DisplayName',"Mod"})
% shadedErrorBar(norm_fois,mean(unmod_PLVs(:,:,:), [1 2]),std(mean(unmod_PLVs(:,:,:), 2))/sqrt(size(unmod_PLVs,1)),'lineprops',{'Color',[.6 .2 .1],'LineWidth',4,'DisplayName',"Unmod"})

% shadedErrorBar(norm_fois,mean(sart_mod_PLV(asrs_idxs,3,:), 1),std(sart_mod_PLV(asrs_idxs,3,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .5 .1],'LineWidth',4,'DisplayName',"Song3"})
% shadedErrorBar(norm_fois,mean(sart_mod_PLV(asrs_idxs,2,:), 1),std(sart_mod_PLV(asrs_idxs,2,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .2 .1],'LineWidth',4,'DisplayName',"Song2"})
% shadedErrorBar(norm_fois,mean(sart_mod_PLV(asrs_idxs,1,:), 1),std(sart_mod_PLV(asrs_idxs,1,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.4 0 .1],'LineWidth',4,'DisplayName',"Song1"})
% 
shadedErrorBar(norm_fois,mean(sart_mod_PLV(noasrs_idxs,3,:), 1),std(sart_mod_PLV(noasrs_idxs,3,:),1)/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.6 .5 .8],'LineWidth',4,'DisplayName',"Song3"})
shadedErrorBar(norm_fois,mean(sart_mod_PLV(noasrs_idxs,2,:), 1),std(sart_mod_PLV(noasrs_idxs,2,:),1)/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.6 .2 .6],'LineWidth',4,'DisplayName',"Song2"})
shadedErrorBar(norm_fois,mean(sart_mod_PLV(noasrs_idxs,1,:), 1),std(sart_mod_PLV(noasrs_idxs,1,:),1)/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.2 0 .6],'LineWidth',4,'DisplayName',"Song1"})


% shadedErrorBar(norm_fois,mean(sart_unmod_PLV(asrs_idxs,3,:), 1),std(sart_unmod_PLV(asrs_idxs,3,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .5 .1],'LineWidth',4,'DisplayName',"ASRS3"})
% shadedErrorBar(norm_fois,mean(sart_unmod_PLV(asrs_idxs,2,:), 1),std(sart_unmod_PLV(asrs_idxs,2,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .2 .1],'LineWidth',4,'DisplayName',"ASRS2"})
% shadedErrorBar(norm_fois,mean(sart_unmod_PLV(asrs_idxs,1,:), 1),std(sart_unmod_PLV(asrs_idxs,1,:),1)/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.4 0 .1],'LineWidth',4,'DisplayName',"ASRS1"})


% 
% shadedErrorBar(norm_fois,mean(mod_PLVs(asrs_idxs,:,:), [1 2]),std(mean(mod_PLVs(asrs_idxs,:,:), 2))/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.7 .1 .1],'LineWidth',4,'DisplayName',"Mod ASRS Positive"})
% shadedErrorBar(norm_fois,mean(unmod_PLVs(asrs_idxs,:,:), [1 2]),std(mean(unmod_PLVs(asrs_idxs,:,:), 2))/sqrt(length(asrs_idxs)),'lineprops',{'Color',[.6 .2 .1],'LineWidth',4,'DisplayName',"Unmod ASRS Positive"})
% 
% shadedErrorBar(norm_fois,mean(mod_PLVs(noasrs_idxs,:,:), [1 2]),std(mean(mod_PLVs(noasrs_idxs,:,:), 2))/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.1 .1 .6],'LineWidth',4,'DisplayName',"Mod ASRS Negative"})
% shadedErrorBar(norm_fois,mean(unmod_PLVs(noasrs_idxs,:,:), [1 2]),std(mean(unmod_PLVs(noasrs_idxs,:,:), 2))/sqrt(length(noasrs_idxs)),'lineprops',{'Color',[.1 .4 .6],'LineWidth',4,'DisplayName',"Unmod ASRS Negative"})

xline(4, "-", "mod rate", "LineWidth",1.5, "Color", "k", "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment","left", "HandleVisibility","off");
xline([.5 1 2], "--", "LineWidth",1, "HandleVisibility","off");
yline(0, "k-", "LineWidth",1, "HandleVisibility","off");
hold off
legend
ylim([0 .25])
xlim([0 10.2])
xticks([.5 1:10])
xlabel("multiple of beat rate")
ylabel("Phase-locking value")
title("Phase-locking Over SART Songs: ASRS-")
set(gcf, "color", "w")

fontsize(16, "points")
%%





%% test the effect of amplitude modulation on ASRS positive and negative

freq_idx = norm_fois == 4; % mod rate

% [h, p, ci, stats1] = ttest2(mean(mod_PLVs(asrs_idxs,:,freq_idx), 3)', mean(unmod_PLVs(asrs_idxs,:,freq_idx), 3)');
% [h, p, ci, stats2] = ttest2(mean(mod_PLVs(noasrs_idxs,:,freq_idx), 3)', mean(unmod_PLVs(noasrs_idxs,:,freq_idx), 3)');


[h, p, ci, stats1] = ttest2(mean(mod_PLVs(asrs_idxs,:,freq_idx), 2), mean(mod_PLVs(noasrs_idxs,:,freq_idx), 2));



%%
% expects each column to be a channel timeseries 
function [phase_timeSeries, amplitude_timeSeries, filtered_timeSeries] = eegConvolution(eegdata, cycles, frequency, srate)
    time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
    % create complex sine wave
    sine_wave = exp( 1i*2*pi*frequency.*time );
    % create Gaussian window
    s = cycles / (2*pi*frequency); % this is the standard deviation of the gaussian
    gaus_win  = exp( (-time.^2) ./ (2*s^2) );
    % now create Morlet wavelet
    cmw = gaus_win.*sine_wave;

    nData = length(eegdata);
    nKern = length(cmw);
    nConv = nData + nKern - 1;
    %FFTs:
    % note that the "N" parameter is the length of convolution, NOT the length
    % of the original signals! Super-important!
    % FFT of wavelet, and amplitude-normalize in the frequency domain
    cmwX = fft(cmw,nConv);
    cmwX = cmwX ./ max(cmwX);
    % FFT of data
    dataX = fft(eegdata,nConv);
    conv_res = dataX.*cmwX.'; %% --> we have just performed a convolution in the frequency domain (potentially the coolest thing ever) 
    % now back to the time domain
    % cut 1/2 of the length of the wavelet from the beginning and from the end
    half_wav = floor( length(cmw)/2 )+1;
    % take inverse Fourier transform
    conv_res_timedomain = ifft(conv_res);
    conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav,:);
    % conv_res_timedomain = conv_res_timedomain(100:end-100,:); %remove edge artifacts from convolution by trimming 100 samples from both ends
    phase_timeSeries = angle(conv_res_timedomain);
    amplitude_timeSeries = abs(conv_res_timedomain);
    filtered_timeSeries = real(conv_res_timedomain);
    % power_timeSeries = abs(conv_res_timedomain).^2;
end
%%

% CITE: Information Theory Toolbox by Mo Chen
function z = mutInfo(x, y)
% Compute mutual information I(x,y) of two discrete variables x and y.
% Input:
%   x, y: two integer vector of the same length 
% Output:
%   z: mutual information z=I(x,y)
% Written by Mo Chen (sth4nth@gmail.com).
assert(numel(x) == numel(y));
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);
l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));
idx = 1:n;
Mx = sparse(idx,x,1,n,k,n);
My = sparse(idx,y,1,n,k,n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log2(Pxy));
Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));
% entropy of Py and Px
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));
% mutual information
z = Hx+Hy-Hxy;
z = max(0,z);
end

% CITE: Information Theory Toolbox by Mo Chen
function z = nmi(x, y)
% Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
% Input:
%   x, y: two integer vector of the same length 
% Ouput:
%   z: normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))
% Written by Mo Chen (sth4nth@gmail.com).
assert(numel(x) == numel(y));
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);
l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));
idx = 1:n;
Mx = sparse(idx,x,1,n,k,n);
My = sparse(idx,y,1,n,k,n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log2(Pxy));
% hacking, to elimative the 0log0 issue
Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));
% entropy of Py and Px
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));
% mutual information
MI = Hx + Hy - Hxy;
% normalized mutual information
z = sqrt((MI/Hx)*(MI/Hy));
z = max(0,z);
end