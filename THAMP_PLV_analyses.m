%% PLV Analysis
participantIDs = {"240418JDOR", "240621AWOO", "240709JMAR", "240730MANG", "240808CROC", "240813IITS", "240815BZHE", "240816AJOL", "240820LWU"};
numParticipants = length(participantIDs);
path_to_data = '/path/to/THAMP/EEG_Data/analyzed/'; % on discovery
song_folder = '/path/to/THAMP/Song Library/normalized_3_21_mp3_4_28'; % download from dropbox
thamp_song_library = readtable("/path/to/THAMP/THAMP Song Library.xlsx", "NumHeaderLines",0, "VariableNamingRule","preserve"); % download from dropbox
addpath(song_folder)

% parameters and initialization
numBins = 101; % number of frequencies (frequency resolution)
low_norm = .1; % normalization bounds
high_norm = 10.1; % normalization bounds
lowFreq = 0.2;
highFreq = 20.2;
fois = linspace(lowFreq, highFreq, numBins); % frequency values at which we calculate PLV
norm_fois = linspace(low_norm, high_norm, numBins); % normalized frequency ratios at which we calculate PLV
cycles = 5; % determines resolution of wavelet
aggregateChannelPLVs = zeros(numParticipants,12, numBins, 61); % participants x songs x fois x electrodes
mod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, 61); % participants x songs x fois x electrodes
unmod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, 61); % participants x songs x fois x electrodes
norm_aggregateChannelPLVs = zeros(numParticipants,12, numBins, 61); 
norm_mod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, 61); 
norm_unmod_aggregateChannelPLVs = zeros(numParticipants,6, numBins, 61);

%%
for participant_idx = 1:numParticipants
    participantID = participantIDs{participant_idx};
    
    disp(strcat("processing participant: ", participantID))
    songs = readtable(fullfile(path_to_data,participantID, "song_order.csv"));
    mod_idxs = find(strcmp(string([songs.condition]), "Mod"));
    unmod_idxs = find(strcmp(string([songs.condition]), "Unmod"));
    
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
        %% MIR Toolbox preprocessing audio (cochlear filterbank)
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
        eegdata = EEG.data.';
        
        % Trim
        if abs(length(eegdata) - length(audio_signal))/audio_fs > 1, error("length mismatch"), end % if the files are more than 1 second different in length
        if length(eegdata) > length(audio_signal)
            eegdata = eegdata(1:length(audio_signal),:);
        else
            audio_signal = audio_signal(1:length(eegdata));
        end
        %%
        % calc PLVs
        for freq_idx = 1:numBins
            [EEG_phases, ~] = eegConvolution(eegdata,cycles,fois(freq_idx),EEG.srate); % convolve with EEG
            [Audio_phases, ~] = eegConvolution(audio_signal,cycles,fois(freq_idx),audio_fs); % convolve with Audio
            aggregateChannelPLVs(participant_idx,song_idx, freq_idx,:) = abs(mean(exp(1i*(EEG_phases-Audio_phases)),1));
        end
        % calc norm_PLVs
        for freq_idx = 1:numBins
            [EEG_phases, ~] = eegConvolution(eegdata,cycles,norm_freqs(freq_idx),EEG.srate); % convolve with EEG
            [Audio_phases, ~] = eegConvolution(audio_signal,cycles,norm_freqs(freq_idx),audio_fs); % convolve with Audio
            norm_aggregateChannelPLVs(participant_idx,song_idx, freq_idx,:) = abs(mean(exp(1i*(EEG_phases-Audio_phases)),1));
        end
    end
    mod_aggregateChannelPLVs(participant_idx,:, :,:) = aggregateChannelPLVs(participant_idx,mod_idxs, :,:);
    unmod_aggregateChannelPLVs(participant_idx,:,:,:) = aggregateChannelPLVs(participant_idx,unmod_idxs, :,:);
    norm_mod_aggregateChannelPLVs(participant_idx,:, :,:) = norm_aggregateChannelPLVs(participant_idx,mod_idxs, :,:);
    norm_unmod_aggregateChannelPLVs(participant_idx,:,:,:) = norm_aggregateChannelPLVs(participant_idx,unmod_idxs, :,:);
end
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
hold off
legend
xlim([0 10.2])
xlabel("multiple of beat rate")
ylabel("Phase-locking value")
title("Phase-locking Values at Multiples of the Beat Rate")
set(gcf, "color", "w")

fontsize(16, "points")

% print("THAMP_PLV.png", "-dpng", "-r500")

%%
freq_range = [3.8 4.2];
topoSlice = squeeze(mean(norm_unmod_aggregateChannelPLVs(:,:,norm_fois>=freq_range(1) & norm_fois<=freq_range(2),:),[1 2 3]));
figure;
topoplot(topoSlice, EEG.chanlocs,  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
clim([.02 .15])
colorbar

%%
% expects each column to be a channel timeseries 
function [phase_timeSeries, amplitude_timeSeries] = eegConvolution(eegdata, cycles, frequency, srate)
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
    % filtered_timeSeries = real(conv_res_timedomain);
    % power_timeSeries = abs(conv_res_timedomain).^2;
end

