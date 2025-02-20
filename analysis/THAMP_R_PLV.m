%% THAMP PLV analysis

% Arun Asthagiri


path_to_data = '/Users/arun/Documents/MINDLab/THAMP/EEG_Data/analyzed/final_analyzed';
song_folder = '/Users/arun/Documents/MINDLab/THAMP/Song Library/normalized_3_21_mp3_4_28';
qualtrics_table = readtable('/Users/arun/Documents/MINDLab/THAMP/THAMP_eeg_scored_qualtrics.csv');
thamp_song_library = readtable("/Users/arun/Documents/MINDLab/THAMP/THAMP Song Library.xlsx", "NumHeaderLines",0, "VariableNamingRule","preserve");
addpath(song_folder)

numBins = 101; % number of frequencies (frequency resolution)
lowFreq = 0.2;
low_norm = .1;
high_norm = 10.1; % normalization bounds
highFreq = 20.2;
fois = linspace(lowFreq, highFreq, numBins); % frequency values at which we calculate PLV
norm_fois = linspace(low_norm, high_norm, numBins); % normalized frequency ratios at which we calculate PLV
cycles = 5; % determines resolution of wavelet

num_electrodes = 61;

%%
participantIDs = string({dir(path_to_data).name});
participantIDs = participantIDs(~contains(participantIDs, "."));
numParticipants = length(participantIDs);
% SART ONLY
% norm_aggregateChannelPLVs = zeros(numParticipants,6, numBins, num_electrodes); % participants x songs x fois x electrodes
norm_mod_aggregateChannelPLVs = nan(numParticipants,3, numBins, num_electrodes); 
norm_unmod_aggregateChannelPLVs = nan(numParticipants,3, numBins, num_electrodes);
asrs_scores = nan(numParticipants,1);
ebmrq_scores = nan(numParticipants,1);

%% Iterate through participants
for participant_idx = 1:numParticipants
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
    try
        asrs_scores(participant_idx)=qualtrics_table(strcmp([qualtrics_table.SUBJECTID], participantID),:).ASRS_Positive;
        ebmrq_scores(participant_idx)=qualtrics_table(strcmp([qualtrics_table.SUBJECTID], participantID),:).eBMRQ_Total_Score;
    catch
        disp("qualtrics data not found")
    end
    
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    sart_conditions = strcmp(string([songs.task]), "SART");
    mod_idxs = strcmp(string([songs.condition]), "Mod"); mod_idxs = mod_idxs & sart_conditions;
    unmod_idxs = strcmp(string([songs.condition]), "Unmod"); unmod_idxs = unmod_idxs & sart_conditions;
    
    norm_aggregateChannelPLVs= nan(12, numBins,num_electrodes);
    
    for song_idx = 1:12
        try
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
        eegdata = EEG.data';
        % %% laplacian
        % [surf_lap,G,H] = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
        % eegdata = surf_lap';
        %%
        
        % Trim
        if abs(length(eegdata) - length(audio_signal))/audio_fs > 1, error("length mismatch"), end % if the files are more than 1 second different in length
        if length(eegdata) > length(audio_signal)
            eegdata = eegdata(1:length(audio_signal),:);
        else
            audio_signal = audio_signal(1:length(eegdata));
        end
        % calc norm_PLVs
        for freq_idx = 1:numBins
            [EEG_phases, ~] = eegConvolution(eegdata,cycles,norm_freqs(freq_idx),EEG.srate); % convolve with EEG
            [Audio_phases, ~] = eegConvolution(audio_signal,cycles,norm_freqs(freq_idx),audio_fs); % convolve with Audio
            norm_aggregateChannelPLVs(song_idx, freq_idx,:) = abs(mean(exp(1i*(EEG_phases-Audio_phases)),1));
        end
        catch
            disp("skipping song")
        end
    
    end
    norm_mod_aggregateChannelPLVs(participant_idx,:, :,:) = norm_aggregateChannelPLVs(mod_idxs, :,:);
    norm_unmod_aggregateChannelPLVs(participant_idx,:,:,:) = norm_aggregateChannelPLVs(unmod_idxs, :,:);
end

%%
% writematrix([participantIDs' asrs_scores ebmrq_scores],"bmrq.csv")
% save("all_SART_PLVs", "norm_mod_aggregateChannelPLVs", "norm_unmod_aggregateChannelPLVs")
%% RTCV

mod_SART_RT = nan(numParticipants, 3);
unmod_SART_RT = nan(numParticipants, 3);

mod_SART_RTCV = nan(numParticipants, 3);
unmod_SART_RTCV = nan(numParticipants, 3);


% Iterate through participants
for participant_idx = 1:numParticipants
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    for song_idx = 1:12
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
        % calc RTCV for whole song
        RTs = [];
        data = compiled_trigs;
        for trial_idx = 1:size(data,1)-1
            if data(trial_idx, 1) ~= 3 && data(trial_idx, 1) ~= 10 && data(trial_idx+1, 1) == 10
                RTs = [RTs; data(trial_idx+1, 2) - data(trial_idx, 2)];
            end
        end
        RT = mean(RTs);
        RTCV = std(RTs)/mean(RTs);


        switch string(songs(song_idx,:).condition)
            case "Mod"
                mod_SART_RTCV(participant_idx, mod(song_idx-1, 3)+1) = RTCV;
                mod_SART_RT(participant_idx, mod(song_idx-1, 3)+1) = RT;
            case "Unmod"
                unmod_SART_RTCV(participant_idx, mod(song_idx-1, 3)+1) = RTCV;
                unmod_SART_RT(participant_idx, mod(song_idx-1, 3)+1) = RT;
            otherwise
                error("unrecognized condition")
        end
    
    end
end
%% save out to R
% writematrix([participantIDs' asrs_scores ebmrq_scores],"qualtrics.csv")
sart_mod_PLV_modfreq = squeeze(norm_mod_aggregateChannelPLVs(:,:,norm_fois==4,:));
sart_unmod_PLV_modfreq = squeeze(norm_unmod_aggregateChannelPLVs(:,:,norm_fois==4,:));
% save("PLV_to_R", "sart_mod_PLV_modfreq", "sart_unmod_PLV_modfreq", "mod_SART_RT", "mod_SART_RTCV", "unmod_SART_RT", "unmod_SART_RTCV");

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


%% PLV over Time
% window length (in secs)
window_length = 10;
hop_length = .5;

% take middle n seconds from song
excerpt_length = 50;

% SART Only
% mod_norm_PLV_over_time = nan(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3, numBins, num_electrodes);
% unmod_norm_PLV_over_time = nan(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3, numBins, num_electrodes);
mod_norm_PLV_over_time = nan(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3, num_electrodes);
unmod_norm_PLV_over_time = nan(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3, num_electrodes);

norm_foi = 4;

%% Iterate through participants
for participant_idx = 1:numParticipants
    tic
    participantID = participantIDs{participant_idx};
    disp(strcat("processing participant: ", participantID))
   
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    sart_conditions = strcmp(string([songs.task]), "SART");
    mod_idxs = strcmp(string([songs.condition]), "Mod"); mod_idxs = mod_idxs & sart_conditions;
    unmod_idxs = strcmp(string([songs.condition]), "Unmod"); unmod_idxs = unmod_idxs & sart_conditions;
    
    norm_PLV_over_time = nan(length(1:hop_length:(excerpt_length-window_length+1)), 12,num_electrodes);
    for song_idx = 1:12
        try
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
            % %% laplacian
            % [surf_lap,G,H] = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]); 
            % eegdata = surf_lap';
            
                    
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
    
    
            % calc norm_PLVs
            window_idxs = idx_start:hop_length*audio_fs:idx_end-window_length*audio_fs+1;
            % alpha_amplitude = abs(hilbert(bandpass(eegdata, [8 12], EEG.srate)));
            % for window_idx = 1:length(window_idxs)
            %     alpha_over_time(window_idx,song_idx,:) = mean(alpha_amplitude(window_idxs(window_idx):window_idxs(window_idx)+window_length*EEG.srate-1,:), 1);
            % end
    
            
            %%
            [EEG_phases, EEG_amplitudes, EEG_filtered] = eegConvolution(eegdata,cycles,norm_freqs(norm_foi),EEG.srate); % convolve with EEG
            [Audio_phases, ~, Audio_filtered] = eegConvolution(audio_signal,cycles,norm_freqs(norm_foi),audio_fs); % convolve with Audio
            %%
            for window_idx = 1:length(window_idxs)
                windowed_EEG_phases = EEG_phases(window_idxs(window_idx):window_idxs(window_idx)+window_length*EEG.srate-1,:);
                windowed_Audio_phases = Audio_phases(window_idxs(window_idx):window_idxs(window_idx)+window_length*audio_fs-1);
                norm_PLV_over_time(window_idx, song_idx, :) = abs(mean(exp(1i* (windowed_EEG_phases-windowed_Audio_phases)), 1));
            end
            
        catch
            disp("skipped song")
        end
    end

    mod_norm_PLV_over_time(participant_idx,:,:,:) = norm_PLV_over_time(:,mod_idxs,:);
    unmod_norm_PLV_over_time(participant_idx,:,:,:) = norm_PLV_over_time(:,unmod_idxs,:);
    % all_alpha_over_time(participant_idx,:, :,:) = alpha_over_time;
    toc
end

%% RTCV Over Time


unmod_SART_RTCV_over_time  = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
mod_SART_RTCV_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
unmod_SART_RT_over_time  = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);
mod_SART_RT_over_time = zeros(numParticipants,length(1:hop_length:(excerpt_length-window_length+1)), 3);


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
            case "Unmod"
                unmod_SART_RTCV_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RTCV_over_time;
                unmod_SART_RT_over_time(participant_idx,:, mod(song_idx-1, 3)+1) = RT_over_time;
            otherwise
                error("unrecognized condition")
        end
    
    end
end
%%
%save("THAMP_PLV_RTCV_over_time", "unmod_SART_RTCV_over_time", "mod_SART_RTCV_over_time", "unmod_SART_RT_over_time", "mod_SART_RT_over_time", "mod_norm_PLV_over_time", "unmod_norm_PLV_over_time")
%% compute correlation
all_Rs = zeros(numParticipants, 3, 61);
for participants = 1:numParticipants
    for song_number = 1:3
        for electrodes = 1:61
            series1 = mean(mod_norm_PLV_over_time(participants, :,song_number, electrodes), 4)';
            series2 = mod_SART_RT_over_time(participants,:,song_number)';
            r = corr(series1, series2, "Type","Pearson");
            % surrogate_rs = surrogate_shifting(series1, series2, 200);
            % [surr_m, surr_sd] = deal(mean(surrogate_rs), std(surrogate_rs));
            % r = (r-surr_m)/surr_sd;

            all_Rs(participants, song_number, electrodes) = r;
        end
    end
end
%%
% [h, p, ci, stats] = ttest(reshape(mean(all_Rs, [ 3]), 1,[]));

%%
figure;
topoplot(squeeze(mean(all_Rs, [1, 2])), EEG.chanlocs);
% clim([-1 1])
clim([-.2 .2])
% clim([0 .2])
colorbar
set(gcf, "color", "w")
title("Mod Phase-locking vs RTCV Correlations over Time")
fontsize(16, "points")
%% surrogate
function surrogate_rs = surrogate_shifting(series1, series2, n_surrogate)
    surrogate_rs = zeros(n_surrogate, 1);
    
    for i = 1:n_surrogate
        shift = randi([round(length(series1)/4) round(3*length(series1)/4)]);
        surrogate_rs(i) = corr(circshift(series1, shift,1), series2, "Type","Pearson");
    end
    % [m, sd] = deal(mean(surrogate_rs), std(surrogate_rs));
end

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

