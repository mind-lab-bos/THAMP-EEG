%THAMPcalcPSD

% Open the dataset ending in
% '_resamp_filt_reref_rej_interp_prunedICA_epochs' for the desired
% participant

% could do: code the creation of epochs to make epoch 1,2,3,4 for the avgs
% of each condition. will have to redo the epoching as pop_select(EEG,'trial', [1:3]) 
% doesnt work in calcPSD, need to combine the 3 epochs into 1 epoch


%%
clear all
clc 

%SINGLE PARTICIPANT PROCESSING
%EDIT THESE EACH TIME, check the runlog for song order and taskorder and
%mod order
ID = '240117LAYR';
ID_l = 'layr';
id = ID;
songorder1 = 1;
songorder2 = 5;
taskorder = 2; %2 for 2back or 1 for sart depending on which is first in the runlog
modorder = 1; %1 for mod or 2 for unmod depending on which is first in the runlog

% also edit the appropriate file path to where the epoched/processed eeg
% file is
path = '/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/';
fullpath = [path ID '/'];
mkdir(fullpath,'analysis');
savepath = [fullpath 'analysis/'];

%opens eeglab and the epoched processed file
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[ID '_filt_reref_resamp_rej_interp_prunedICA_epochs.set'],'filepath',fullpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%looks at song order and finds the specific songs heard
if songorder1 == 1 || songorder1 == 2
    ref1 = [18,14,8,24,23,21];
elseif songorder1 == 3 || songorder1 == 4
    ref1 = [1,7,20,30,5,16];
elseif songorder1 == 5 || songorder1 == 6
    ref1 = [29,11,31,32,12,28];
elseif songorder1 == 7 || songorder1 == 8
    ref1 = [19,17,13,9,2,4];
elseif songorder1 == 9 || songorder1 == 10
    ref1 = [26,22,27,10,3,25];
end
if songorder2 == 1 || songorder2 == 2
    ref2 = [18,14,8,24,23,21];
elseif songorder2 == 3 || songorder2 == 4
    ref2 = [1,7,20,30,5,16];
elseif songorder2 == 5 || songorder2 == 6
    ref2 = [29,11,31,32,12,28];
elseif songorder2 == 7 || songorder2 == 8
    ref2 = [19,17,13,9,2,4];
elseif songorder2 == 9 || songorder2 == 10
    ref2 = [26,22,27,10,3,25];
end

%labelling the order of the blocks
block1 = {'sartmod1','sartmod2','sartmod3'};
block2 = {'sartunmod1','sartunmod2','sartunmod3'};
block3 = {'nbackmod1','nbackmod2','nbackmod3'};
block4 = {'nbackunmod1','nbackunmod2','nbackunmod3'};
if taskorder == 1 && modorder == 1
    titles = {block1,block2,block3,block4};
elseif taskorder == 1 && modorder == 2
    titles = {block2,block1,block4,block3};
elseif taskorder == 2 && modorder == 1
    titles = {block3,block4,block1,block2};
elseif taskorder == 2 && modorder == 2
    titles = {block4,block3,block2,block1};
end

% lookup variable for the theta freq relative to the bpm previously
% calculated for each song
thetasong = [5.266666667,4.533333333,6.2,5.2,7.266666667,8.2,4.8,5.8,4.866666667,5.333333333,5.066666667,5.333333333,6.466666667,5.266666667,7.733333333,8,7.6,4.666666667,8.733333333,4.933333333,7.2,6.4,8,6.333333333,4.866666667,5.6,8.2,7.266666667,6.933333333,6.133333333,4.533333333,5.133333333];

% labels whether each song for this participant was sart or nback or mod or
% unmod and the theta freq value
song1 = {titles{1}{1},num2str(thetasong(ref1(1)))};
song2 = {titles{1}{2},num2str(thetasong(ref1(2)))};
song3 = {titles{1}{3},num2str(thetasong(ref1(3)))};
song4 = {titles{2}{1},num2str(thetasong(ref1(4)))};
song5 = {titles{2}{2},num2str(thetasong(ref1(5)))};
song6 = {titles{2}{3},num2str(thetasong(ref1(6)))};
song7 = {titles{3}{1},num2str(thetasong(ref2(1)))};
song8 = {titles{3}{2},num2str(thetasong(ref2(2)))};
song9 = {titles{3}{3},num2str(thetasong(ref2(3)))};
song10 = {titles{4}{1},num2str(thetasong(ref2(4)))};
song11 = {titles{4}{2},num2str(thetasong(ref2(5)))};
song12 = {titles{4}{3},num2str(thetasong(ref2(6)))};

%selects each of the epochs
EEG1 = pop_select(EEG, 'trial', 1);
EEG2 = pop_select(EEG, 'trial', 2);
EEG3 = pop_select(EEG, 'trial', 3);
EEG4 = pop_select(EEG, 'trial', 4);
EEG5 = pop_select(EEG, 'trial', 5);
EEG6 = pop_select(EEG, 'trial', 6);
EEG7 = pop_select(EEG, 'trial', 7);
EEG8 = pop_select(EEG, 'trial', 8);
EEG9 = pop_select(EEG, 'trial', 9);
EEG10 = pop_select(EEG, 'trial', 10);
EEG11 = pop_select(EEG, 'trial', 11);
%EEG12 = pop_select(EEG, 'trial', 12);

% alternative way of grouping the epochs
% EEGa = pop_select(EEG, 'trial', [1:3]);
% EEGb= pop_select(EEG, 'trial', [4:6]);
% EEGc = pop_select(EEG, 'trial', [7:9]);
% EEGd = pop_select(EEG, 'trial', [10:12]);

% Constants
extendSecs = 0;                         % Pre- and post-window in secs
lfT   =  1;                              % lowest frequency bin
hfT   =  55;                              % highest frequency bin
windowT   = 16000; %8000/2 %15000*2?             % PSD analysis window size (samples)
overlapT  = 8000; %4000/2 % 7500*2?             % PSD overlap (samples) 
nperT = 120;                               % number of bins per octave
noctT = log2(hfT)-log2(lfT);               % number of octaves
NT = noctT*nperT+1;                        % Total number of bins
cfT = logspace(log10(lfT),log10(hfT),NT);  % Log spaced bin center frequencies

% DOCUMENT 
% make sure calcPSD.m function is added to path (is in the github)

% [EEGPSDa, cfT] = calcPSD(EEGa, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDb, cfT] = calcPSD(EEGb, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDc, cfT] = calcPSD(EEGc, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDd, cfT] = calcPSD(EEGd, lfT, hfT, windowT, overlapT); % write PSDs

%calculates the psd for each epoch
[EEGPSD1, cfT] = calcPSD(EEG1, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD2, cfT] = calcPSD(EEG2, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD3, cfT] = calcPSD(EEG3, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD4, cfT] = calcPSD(EEG4, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD5, cfT] = calcPSD(EEG5, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD6, cfT] = calcPSD(EEG6, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD7, cfT] = calcPSD(EEG7, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD8, cfT] = calcPSD(EEG8, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD9, cfT] = calcPSD(EEG9, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD10, cfT] = calcPSD(EEG10, lfT, hfT, windowT, overlapT); % write PSDs
[EEGPSD11, cfT] = calcPSD(EEG11, lfT, hfT, windowT, overlapT); % write PSDs
%[EEGPSD12, cfT] = calcPSD(EEG12, lfT, hfT, windowT, overlapT); % write PSDs

%saves each psd as a different file
save(fullfile(savepath,strcat(ID,'_EEGPSD1.mat')),'EEGPSD1')
save(fullfile(savepath,strcat(ID,'_EEGPSD2.mat')),'EEGPSD2')
save(fullfile(savepath,strcat(ID,'_EEGPSD3.mat')),'EEGPSD3')
save(fullfile(savepath,strcat(ID,'_EEGPSD4.mat')),'EEGPSD4')
save(fullfile(savepath,strcat(ID,'_EEGPSD5.mat')),'EEGPSD5')
save(fullfile(savepath,strcat(ID,'_EEGPSD6.mat')),'EEGPSD6')
save(fullfile(savepath,strcat(ID,'_EEGPSD7.mat')),'EEGPSD7')
save(fullfile(savepath,strcat(ID,'_EEGPSD8.mat')),'EEGPSD8')
save(fullfile(savepath,strcat(ID,'_EEGPSD9.mat')),'EEGPSD9')
save(fullfile(savepath,strcat(ID,'_EEGPSD10.mat')),'EEGPSD10')
save(fullfile(savepath,strcat(ID,'_EEGPSD11.mat')),'EEGPSD11')
%save(fullfile(savepath,strcat(ID,'_EEGPSD12.mat')),'EEGPSD12')



%%
% plotting each song
% Uses the appropriate hz value between 4-8 for each song 1-12. Then
% look at the cfT variable in calcPSD.m to see which bin that freq is in.
% Look at the songlists to see what the 12 songs are for each participant.

% need to add clim statements and colorbar for the topoplots

%SONG 1
figure; 
plot(cfT,mean(EEGPSD1,1))
legend(song1(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song1(1)));
[ d, ix ] = min( abs( cfT-str2double(song1(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD1(:,ix),EEG1.chanlocs,'electrodes','off');
title(id+" "+string(song1(1))+" "+string(song1(2))+'Hz (bin '+string(ix)+')');

%SONG 2
figure; 
plot(cfT,mean(EEGPSD2,1))
legend(song2(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song2(1)));
[ d, ix ] = min( abs( cfT-str2double(song2(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD2(:,ix),EEG2.chanlocs,'electrodes','off');
title(id+" "+string(song2(1))+" "+string(song2(2))+'Hz (bin '+string(ix)+')');

%SONG 3
figure; 
plot(cfT,mean(EEGPSD3,1))
legend(song3(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song3(1)));
[ d, ix ] = min( abs( cfT-str2double(song3(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD3(:,ix),EEG3.chanlocs,'electrodes','off');
title(id+" "+string(song3(1))+" "+string(song3(2))+'Hz (bin '+string(ix)+')');

%SONG 4
figure; 
plot(cfT,mean(EEGPSD4,1))
legend(song4(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song4(1)));
[ d, ix ] = min( abs( cfT-str2double(song4(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD4(:,ix),EEG4.chanlocs,'electrodes','off');
title(id+" "+string(song4(1))+" "+string(song4(2))+'Hz (bin '+string(ix)+')');

%SONG 5
figure; 
plot(cfT,mean(EEGPSD5,1))
legend(song5(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song5(1)));
[ d, ix ] = min( abs( cfT-str2double(song5(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD5(:,ix),EEG5.chanlocs,'electrodes','off');
title(id+" "+string(song5(1))+" "+string(song5(2))+'Hz (bin '+string(ix)+')');

%SONG 6
figure; 
plot(cfT,mean(EEGPSD6,1))
legend(song6(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song6(1)));
[ d, ix ] = min( abs( cfT-str2double(song6(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD6(:,ix),EEG6.chanlocs,'electrodes','off');
title(id+" "+string(song6(1))+" "+string(song6(2))+'Hz (bin '+string(ix)+')');

%all sart plot
figure;
plot(cfT,mean(EEGPSD1,1),'b',cfT,mean(EEGPSD2,1),'b',cfT,mean(EEGPSD3,1),'b',cfT,mean(EEGPSD4,1),'r',cfT,mean(EEGPSD5,1),'r',cfT,mean(EEGPSD6,1),'r')
xlabel('Hz'),ylabel('PSD');
title(id+" "+'Sart, all songs')
xlim([0 45])
legend(string(song1(1)),string(song2(1)),string(song3(1)),string(song4(1)),string(song5(1)),string(song6(1)))

%SONG 7
figure; 
plot(cfT,mean(EEGPSD7,1))
legend(song7(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song7(1)));
[ d, ix ] = min( abs( cfT-str2double(song7(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD7(:,ix),EEG7.chanlocs,'electrodes','off');
title(id+" "+string(song7(1))+" "+string(song7(2))+'Hz (bin '+string(ix)+')');

%SONG 8
figure; 
plot(cfT,mean(EEGPSD8,1))
legend(song8(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song8(1)));
[ d, ix ] = min( abs( cfT-str2double(song8(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD8(:,ix),EEG8.chanlocs,'electrodes','off');
title(id+" "+string(song8(1))+" "+string(song8(2))+'Hz (bin '+string(ix)+')');

%SONG 9
figure; 
plot(cfT,mean(EEGPSD9,1))
legend(song9(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song9(1)));
[ d, ix ] = min( abs( cfT-str2double(song9(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD9(:,ix),EEG9.chanlocs,'electrodes','off');
title(id+" "+string(song9(1))+" "+string(song9(2))+'Hz (bin '+string(ix)+')');

%SONG 10
figure; 
plot(cfT,mean(EEGPSD10,1))
legend(song10(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song10(1)));
[ d, ix ] = min( abs( cfT-str2double(song10(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD10(:,ix),EEG10.chanlocs,'electrodes','off');
title(id+" "+string(song10(1))+" "+string(song10(2))+'Hz (bin '+string(ix)+')');

%SONG 11
figure; 
plot(cfT,mean(EEGPSD11,1))
legend(song11(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song11(1)));
[ d, ix ] = min( abs( cfT-str2double(song11(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD11(:,ix),EEG11.chanlocs,'electrodes','off');
title(id+" "+string(song11(1))+" "+string(song11(2))+'Hz (bin '+string(ix)+')');

%SONG 12
figure; 
plot(cfT,mean(EEGPSD12,1))
legend(song12(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song12(1)));
[ d, ix ] = min( abs( cfT-str2double(song12(2)) ) ); %choose the closest bin to the songs theta Hz value
figure
topoplot(EEGPSD12(:,ix),EEG12.chanlocs,'electrodes','off');
title(id+" "+string(song12(1))+" "+string(song12(2))+'Hz (bin '+string(ix)+')');

%all nback plot
figure;
plot(cfT,mean(EEGPSD7,1),'b',cfT,mean(EEGPSD8,1),'b',cfT,mean(EEGPSD9,1),'b',cfT,mean(EEGPSD10,1),'r',cfT,mean(EEGPSD11,1),'r',cfT,mean(EEGPSD12,1),'r')
xlabel('Hz'),ylabel('PSD');
title(id+" "+'Nback, all songs')
xlim([0 45])
legend(string(song7(1)),string(song8(1)),string(song9(1)),string(song10(1)),string(song11(1)),string(song12(1)))


% plotting avgs
% be sure to edit legend based on sart/nback first and mod/unmod first
% figure; 
% plot(cfT,mean(EEGPSDa,1),cfT,mean(EEGPSDb,1))
% legend({'sartMod','sartUnmod'});
% xlabel('Hz'),ylabel('PSD');
% 
% figure; 
% plot(cfT,mean(EEGPSDc,1),cfT,mean(EEGPSDd,1))
% legend({'nbackMod','nbackUnmod'});
% xlabel('Hz'),ylabel('PSD');
% 
% DiffEEGPSDsart = EEGPSDa-EEGPSDb;
% figure; 
% topoplot(mean(DiffEEGPSDsart(:,638:642),2),EEG.chanlocs,'electrodes','off');
% title('40 Hz');
% 
% DiffEEGPSDnback = EEGPSDc-EEGPSDd;
% figure; 
% topoplot(mean(DiffEEGPSDnback(:,638:642),2),EEG.chanlocs,'electrodes','off');
% title('40 Hz');

%close all

