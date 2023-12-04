%THAMPcalcPSD

% Open the dataset ending in
% '_resamp_filt_reref_rej_interp_prunedICA_epochs' for the desired
% participant

% could do: code the creation of epochs to make epoch 1,2,3,4 for the avgs
% of each condition. will have to redo the epoching as pop_select(EEG,'trial', [1:3]) 
% doesnt work in calcPSD, need to combine the 3 epochs into 1 epoch


% automate the theta values, songnames, condition names from songlist input
id = '230926TELZ';
songlist = [1,5];
song1 = {'sartmod1','4.6666666'};
song2 = {'sartmod2','5.26666666'};
song3 = {'sartmod3','5.8'};
song4 = {'sartunmod1','6.333'};
song5 = {'sartunmod2','8'};
song6 = {'sartunmod3','7.2'};
song7 = {'nbackmod1','6.93333333'};
song8 = {'nbackmod2','5.066666'};
song9 = {'nbackmod3','4.533333'};
song10 = {'nbackunmod1','5.13333'};
song11 = {'nbackunmod2','5.33333'};
song12 = {'nbackunmod3','7.26666'};

% id = '230929SCOE';
% songlist = [2,6];
% song1 = {'nbackunmod1','4.6666666'};
% song2 = {'nbackunmod2','5.26666666'};
% song3 = {'nbackunmod3','5.8'};
% song4 = {'nbackmod1','6.333'};
% song5 = {'nbackmod2','8'};
% song6 = {'nbackmod3','7.2'};
% song7 = {'sartunmod1','6.93333333'};
% song8 = {'sartunmod2','5.066666'};
% song9 = {'sartunmod3','4.533333'};
% song10 = {'sartmod1','5.13333'};
% song11 = {'sartmod2','5.33333'};
% song12 = {'sartmod3','7.26666'};

% id = '231004LLE';
% songlist = [3,7];
% song1 = {'sartmod1','5.2666666'};
% song2 = {'sartmod2','4.8'};
% song3 = {'sartmod3','4.933333'};
% song4 = {'sartunmod1','6.1333333'};
% song5 = {'sartunmod2','7.2666666'};
% song6 = {'sartunmod3','8'};
% song7 = {'nbackmod1','8.733333'};
% song8 = {'nbackmod2','7.6'};
% song9 = {'nbackmod3','6.46666'};
% song10 = {'nbackunmod1','4.8666666'};
% song11 = {'nbackunmod2','4.53333'};
% song12 = {'nbackunmod3','5.2'};

%EEG = pop_select(EEG, 'point', epochPnts)
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
EEG12 = pop_select(EEG, 'trial', 12);
% 
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

% [EEGPSDa, cfT] = calcPSD(EEGa, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDb, cfT] = calcPSD(EEGb, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDc, cfT] = calcPSD(EEGc, lfT, hfT, windowT, overlapT); % write PSDs
% [EEGPSDd, cfT] = calcPSD(EEGd, lfT, hfT, windowT, overlapT); % write PSDs

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
[EEGPSD12, cfT] = calcPSD(EEG12, lfT, hfT, windowT, overlapT); % write PSDs


%%
% plotting each song
% look at https://www.dropbox.com/scl/fi/uwppvegwe1k2bal1r78py/THAMP-Song-Library.xlsx?dl=0&rlkey=r3foy9pl26jdn7plgnpuohi0h
% to find the appropriate hz value between 4-8 for each song 1-12. Then
% look at the cfT variable in calcPSD.m to see which bin that freq is in.
% Look at the songlists to see what the 12 songs are for each participant.

%TOPOS plot scale option?
%'maplimits' 'absmax' +/- the absolute-max 'maxmin' scale to data range [clim1,clim2] user-definined lo/hi

%SONG 1
figure; 
plot(cfT,mean(EEGPSD1,1))
legend(song1(1));
xlabel('Hz'),ylabel('PSD');
title(id + " " + string(song1(1)));

% figure;
% topoplot(mean(EEGPSD1(:,241:360),2),EEG1.chanlocs,'electrodes','off');
% title('4-8 Hz');

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

%%
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