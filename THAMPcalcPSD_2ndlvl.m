
% second level 
% opens all existing psd files from each participant
% AVERAGING ACROSS ALL PARTICIPANTS

clc
clear all

% edit this ID list to have all the fully analyzed participants
% edit the song order and mod order list to match the order of the idlist
ID_list = {'230926TELZ' '230929SCOE' '231004LLE' '231114CWAY' '231117SSAY' '231120ICON' '231128MYAN' '231205KKAT' '231205NSAZ' '231211TABO' '231212TNAR' '231213DROD'};
songorders = [1,	5; 2,	6;3,	7;1,	5;3,	7;2,	6;4,	8;5,	9;7,	1;6,	10;8,	2;9,	3];
modorder = {'MOD','UNMOD','MOD','MOD','MOD','UNMOD','UNMOD','MOD','MOD','UNMOD','UNMOD','MOD'};

% edit the songorder for each participant, copied in from the spreadsheet
songorder=[18	14	8	24	23	21	29	11	31	32	12	28;
18	14	8	24	23	21	29	11	31	32	12	28;
1	7	20	30	5	16	19	17	13	9	2	4;
18	14	8	24	23	21	29	11	31	32	12	28;
1	7	20	30	5	16	19	17	13	9	2	4;
18	14	8	24	23	21	29	11	31	32	12	28;
1	7	20	30	5	16	19	17	13	9	2	4;
29	11	31	32	12	28	26	22	27	10	3	25;
19	17	13	9	2	4	18	14	8	24	23	21;
29	11	31	32	12	28	26	22	27	10	3	25;
19	17	13	9	2	4	18	14	8	24	23	21
26	22	27	10	3	25	1	7	20	30	5	16];

% initializes all the variables set to zero
psd1m = zeros(62, 694);  psd1u = zeros(62, 694);
psd2m = zeros(62, 694);  psd2u = zeros(62, 694);
psd3m = zeros(62, 694);  psd3u = zeros(62, 694);
psd4m = zeros(62, 694);  psd4u = zeros(62, 694);
psd5m = zeros(62, 694);  psd5u = zeros(62, 694);
psd6m = zeros(62, 694);  psd6u = zeros(62, 694);
psd7m = zeros(62, 694);  psd7u = zeros(62, 694);
psd8m = zeros(62, 694);  psd8u = zeros(62, 694);
psd9m = zeros(62, 694);  psd9u = zeros(62, 694);
psd10m = zeros(62, 694); psd10u = zeros(62, 694);
psd11m = zeros(62, 694); psd11u = zeros(62, 694);
psd12m = zeros(62, 694); psd12u = zeros(62, 694);
psd13m = zeros(62, 694); psd13u = zeros(62, 694);
psd14m = zeros(62, 694); psd14u = zeros(62, 694);
psd15m = zeros(62, 694); psd15u = zeros(62, 694);
psd16m = zeros(62, 694); psd16u = zeros(62, 694);
psd17m = zeros(62, 694); psd17u = zeros(62, 694);
psd18m = zeros(62, 694); psd18u = zeros(62, 694);
psd19m = zeros(62, 694); psd19u = zeros(62, 694);
psd20m = zeros(62, 694); psd20u = zeros(62, 694);
psd21m = zeros(62, 694); psd21u = zeros(62, 694);
psd22m = zeros(62, 694); psd22u = zeros(62, 694);
psd23m = zeros(62, 694); psd23u = zeros(62, 694);
psd24m = zeros(62, 694); psd24u = zeros(62, 694);
psd25m = zeros(62, 694); psd25u = zeros(62, 694);
psd26m = zeros(62, 694); psd26u = zeros(62, 694);
psd27m = zeros(62, 694); psd27u = zeros(62, 694);
psd28m = zeros(62, 694); psd28u = zeros(62, 694);
psd29m = zeros(62, 694); psd29u = zeros(62, 694);
psd30m = zeros(62, 694); psd30u = zeros(62, 694);
psd31m = zeros(62, 694); psd31u = zeros(62, 694);
psd32m = zeros(62, 694); psd32u = zeros(62, 694);

% initializes all the coounter variables at 0
counter1u = 0; counter1m = 0;
counter2u = 0; counter2m = 0;
counter3u = 0; counter3m = 0;
counter4u = 0; counter4m = 0;
counter5u = 0; counter5m = 0;
counter6u = 0; counter6m = 0;
counter7u = 0; counter7m = 0;
counter8u = 0; counter8m = 0;
counter9u = 0; counter9m = 0;
counter10u = 0; counter10m = 0;
counter11u = 0; counter11m = 0;
counter12u = 0; counter12m = 0;
counter13u = 0; counter13m = 0;
counter14u = 0; counter14m = 0;
counter15u = 0; counter15m = 0;
counter16u = 0; counter16m = 0;
counter17u = 0; counter17m = 0;
counter18u = 0; counter18m = 0;
counter19u = 0; counter19m = 0;
counter20u = 0; counter20m = 0;
counter21u = 0; counter21m = 0;
counter22u = 0; counter22m = 0;
counter23u = 0; counter23m = 0;
counter24u = 0; counter24m = 0;
counter25u = 0; counter25m = 0;
counter26u = 0; counter26m = 0;
counter27u = 0; counter27m = 0;
counter28u = 0; counter28m = 0;
counter29u = 0; counter29m = 0;
counter30u = 0; counter30m = 0;
counter31u = 0; counter31m = 0;
counter32u = 0; counter32m = 0;

% there are 32 total songs, 30 of those were used
% loop through each song
% and then for each song loop through each participant
% and then for each participant loop through their 12 songs
% 3 nested for loops
for j = 1:32 
    for i = 1:length(ID_list)
        ID = ID_list{i}; % current participant in loop
        path = '/Users/Kob/Documents/MINDLab/eeg stuff/thamp_eeg-data/'; %change path for other machines
        fullpath = [path ID '/'];
        savepath = [fullpath 'analysis/'];
        % the files are found in /Users/Kob/Documents/MINDLab/eeg
        % stuff/thamp_eeg-data/(id)/analysis on jakobs machine, change this
        for k = 1:12 
            psdm = [];
            psdu = [];
            %checking to see if modded or unmodded
            if strcmp(modorder{i}, 'MOD') && (k == 1 || k == 2 || k == 3 || k == 7 || k == 8 || k == 9)
                modded = true;
            else, modded = false;
            end
            if strcmp(modorder{i}, 'UNMOD') && (k == 1 || k == 2 || k == 3 || k == 7 || k == 8 || k == 9)
                modded = false;
            elseif strcmp(modorder{i}, 'UNMOD') && ~(k == 1 || k == 2 || k == 3 || k == 7 || k == 8 || k == 9)
                modded = true;
            end
            % j is the current song (1-32), i is the current participant,
            % and k is the current song from that participant (1-12)
            % the two if statements check to see if this current
            % participant listened to the song and if it was modded or
            % unmodded
            if songorder(i,k) == j && modded
                %if the participant listed to that song then add that song
                %to temp variable psdm
                psdm = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                % add loaded variable to psdjm
                %eval([strcat('psd',string(j),'m')]) = eval([strcat('psd',string(j),'m')]) + psdm;
                
                if k == 1, psdm = psdm.EEGPSD1;end
                if k == 2, psdm = psdm.EEGPSD2; end
                if k == 3, psdm = psdm.EEGPSD3; end
                if k == 4, psdm = psdm.EEGPSD4; end
                if k == 5, psdm = psdm.EEGPSD5; end
                if k == 6, psdm = psdm.EEGPSD6; end
                if k == 7, psdm = psdm.EEGPSD7; end
                if k == 8, psdm = psdm.EEGPSD8; end
                if k == 9, psdm = psdm.EEGPSD9; end
                if k == 10, psdm = psdm.EEGPSD10; end
                if k == 11, psdm = psdm.EEGPSD11; end
                if k == 12, psdm = psdm.EEGPSD12; end

                %some of the eegpsd.mat files have the last row as all -inf
                %which breaks the plots. this checks to see if the last row
                %is -infs and replaces that row with the row above.
                if all(psdm(end, :) == -inf),psdm(end, :) = psdm(end - 1, :);end

                % then add the temp variable psdm to the psd variable for
                % each song
                if j == 1, psd1m = psd1m + psdm; counter1m = counter1m + 1;psd1m(62,1), end
                if j == 2, psd2m = psd2m + psdm; counter2m = counter2m + 1; end
                if j == 3, psd3m = psd3m + psdm; counter3m = counter3m + 1; end
                if j == 4, psd4m = psd4m + psdm; counter4m = counter4m + 1; end
                if j == 5, psd5m = psd5m + psdm; counter5m = counter5m + 1; end
                if j == 6, psd6m = psd6m + psdm; counter6m = counter6m + 1; end
                if j == 7, psd7m = psd7m + psdm; counter7m = counter7m + 1; end
                if j == 8, psd8m = psd8m + psdm; counter8m = counter8m + 1; end
                if j == 9, psd9m = psd9m + psdm; counter9m = counter9m + 1; end
                if j == 10, psd10m = psd10m + psdm; counter10m = counter10m + 1; end
                if j == 11, psd11m = psd11m + psdm; counter11m = counter11m + 1; end
                if j == 12, psd12m = psd12m + psdm; counter12m = counter12m + 1; end
                if j == 13, psd13m = psd13m + psdm; counter13m = counter13m + 1; end
                if j == 14, psd14m = psd14m + psdm; counter14m = counter14m + 1; end
                if j == 15, psd15m = psd15m + psdm; counter15m = counter15m + 1; end
                if j == 16, psd16m = psd16m + psdm; counter16m = counter16m + 1; end
                if j == 17, psd17m = psd17m + psdm; counter17m = counter17m + 1; end
                if j == 18, psd18m = psd18m + psdm; counter18m = counter18m + 1; end
                if j == 19, psd19m = psd19m + psdm; counter19m = counter19m + 1; end
                if j == 20, psd20m = psd20m + psdm; counter20m = counter20m + 1; end
                if j == 21, psd21m = psd21m + psdm; counter21m = counter21m + 1; end
                if j == 22, psd22m = psd22m + psdm; counter22m = counter22m + 1; end
                if j == 23, psd23m = psd23m + psdm; counter23m = counter23m + 1; end
                if j == 24, psd24m = psd24m + psdm; counter24m = counter24m + 1; end
                if j == 25, psd25m = psd25m + psdm; counter25m = counter25m + 1; end
                if j == 26, psd26m = psd26m + psdm; counter26m = counter26m + 1; end
                if j == 27, psd27m = psd27m + psdm; counter27m = counter27m + 1; end
                if j == 28, psd28m = psd28m + psdm; counter28m = counter28m + 1; end
                if j == 29, psd29m = psd29m + psdm; counter29m = counter29m + 1; end
                if j == 30, psd30m = psd30m + psdm; counter30m = counter30m + 1; end
                if j == 31, psd31m = psd31m + psdm; counter31m = counter31m + 1; end
                if j == 32, psd32m = psd32m + psdm; counter32m = counter32m + 1; end

            end

            % this whole if block is the same as above but for unmodded
            % songs
            if songorder(i,k) == j && ~modded
                 %psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                 psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
%                 temp = eval([strcat('psd',string(j),'u')]);
%                 temp = double(temp);
%                 eval([temp ' = ' temp ' + ' psdu]);
%                 eval([strcat('psd',string(j),'u')]) = temp;

                if k == 1, psdu = psdu.EEGPSD1; end
                if k == 2, psdu = psdu.EEGPSD2; end
                if k == 3, psdu = psdu.EEGPSD3; end
                if k == 4, psdu = psdu.EEGPSD4; end
                if k == 5, psdu = psdu.EEGPSD5; end
                if k == 6, psdu = psdu.EEGPSD6; end
                if k == 7, psdu = psdu.EEGPSD7; end
                if k == 8, psdu = psdu.EEGPSD8; end
                if k == 9, psdu = psdu.EEGPSD9; end
                if k == 10, psdu = psdu.EEGPSD10; end
                if k == 11, psdu = psdu.EEGPSD11; end
                if k == 12, psdu = psdu.EEGPSD12; end

                %some of the eegpsd.mat files have the last row as all -inf
                %which breaks the plots. this checks to see if the last row
                %is -infs and replaces that row with the row above.
                if all(psdu(end, :) == -inf),psdu(end, :) = psdu(end - 1, :);end

                if j == 1, psd1u = psd1u + psdu; counter1u = counter1u + 1; end
                if j == 2, psd2u = psd2u + psdu; counter2u = counter2u + 1; end
                if j == 3, psd3u = psd3u + psdu; counter3u = counter3u + 1; end
                if j == 4, psd4u = psd4u + psdu; counter4u = counter4u + 1; end
                if j == 5, psd5u = psd5u + psdu; counter5u = counter5u + 1; end
                if j == 6, psd6u = psd6u + psdu; counter6u = counter6u + 1; end
                if j == 7, psd7u = psd7u + psdu; counter7u = counter7u + 1; end
                if j == 8, psd8u = psd8u + psdu; counter8u = counter8u + 1; end
                if j == 9, psd9u = psd9u + psdu; counter9u = counter9u + 1; end
                if j == 10, psd10u = psd10u + psdu; counter10u = counter10u + 1; end
                if j == 11, psd11u = psd11u + psdu; counter11u = counter11u + 1; end
                if j == 12, psd12u = psd12u + psdu; counter12u = counter12u + 1; end
                if j == 13, psd13u = psd13u + psdu; counter13u = counter13u + 1; end
                if j == 14, psd14u = psd14u + psdu; counter14u = counter14u + 1; end
                if j == 15, psd15u = psd15u + psdu; counter15u = counter15u + 1; end
                if j == 16, psd16u = psd16u + psdu; counter16u = counter16u + 1; end
                if j == 17, psd17u = psd17u + psdu; counter17u = counter17u + 1; end
                if j == 18, psd18u = psd18u + psdu; counter18u = counter18u + 1; end
                if j == 19, psd19u = psd19u + psdu; counter19u = counter19u + 1; end
                if j == 20, psd20u = psd20u + psdu; counter20u = counter20u + 1; end
                if j == 21, psd21u = psd21u + psdu; counter21u = counter21u + 1; end
                if j == 22, psd22u = psd22u + psdu; counter22u = counter22u + 1; end
                if j == 23, psd23u = psd23u + psdu; counter23u = counter23u + 1; end
                if j == 24, psd24u = psd24u + psdu; counter24u = counter24u + 1; end
                if j == 25, psd25u = psd25u + psdu; counter25u = counter25u + 1; end
                if j == 26, psd26u = psd26u + psdu; counter26u = counter26u + 1; end
                if j == 27, psd27u = psd27u + psdu; counter27u = counter27u + 1; end
                if j == 28, psd28u = psd28u + psdu; counter28u = counter28u + 1; end
                if j == 29, psd29u = psd29u + psdu; counter29u = counter29u + 1; end
                if j == 30, psd30u = psd30u + psdu; counter30u = counter30u + 1; end
                if j == 31, psd31u = psd31u + psdu; counter31u = counter31u + 1; end
                if j == 32, psd32u = psd32u + psdu; counter32u = counter32u + 1; end

            end
        end
    end
end

% check to see if any of the song versions were never listened to, if they
% were listened to, then calculate the average by dividing the summed up
% psds for each song and dividing by the counter
if counter1m ~= 0, psd1m = psd1m/counter1m; end
if counter1u ~= 0, psd1u = psd1u/counter1u; end
if counter2m ~= 0, psd2m = psd2m/counter2m; end
if counter2u ~= 0, psd2u = psd2u/counter2u; end
if counter3m ~= 0, psd3m = psd3m/counter3m; end
if counter3u ~= 0, psd3u = psd3u/counter3u; end
if counter4m ~= 0, psd4m = psd4m/counter4m; end
if counter4u ~= 0, psd4u = psd4u/counter4u; end
if counter5m ~= 0, psd5m = psd5m/counter5m; end
if counter5u ~= 0, psd5u = psd5u/counter5u; end
if counter6m ~= 0, psd6m = psd6m/counter6m; end
if counter6u ~= 0, psd6u = psd6u/counter6u; end
if counter7m ~= 0, psd7m = psd7m/counter7m; end
if counter7u ~= 0, psd7u = psd7u/counter7u; end
if counter8m ~= 0, psd8m = psd8m/counter8m; end
if counter8u ~= 0, psd8u = psd8u/counter8u; end
if counter9m ~= 0, psd9m = psd9m/counter9m; end
if counter9u ~= 0, psd9u = psd9u/counter9u; end
if counter10m ~= 0, psd10m = psd10m/counter10m; end
if counter10u ~= 0, psd10u = psd10u/counter10u; end
if counter11m ~= 0, psd11m = psd11m/counter11m; end
if counter11u ~= 0, psd11u = psd11u/counter11u; end
if counter12m ~= 0, psd12m = psd12m/counter12m; end
if counter12u ~= 0, psd12u = psd12u/counter12u; end
if counter13m ~= 0, psd13m = psd13m/counter13m; end
if counter13u ~= 0, psd13u = psd13u/counter13u; end
if counter14m ~= 0, psd14m = psd14m/counter14m; end
if counter14u ~= 0, psd14u = psd14u/counter14u; end
if counter15m ~= 0, psd15m = psd15m/counter15m; end
if counter15u ~= 0, psd15u = psd15u/counter15u; end
if counter16m ~= 0, psd16m = psd16m/counter16m; end
if counter16u ~= 0, psd16u = psd16u/counter16u; end
if counter17m ~= 0, psd17m = psd17m/counter17m; end
if counter17u ~= 0, psd17u = psd17u/counter17u; end
if counter18m ~= 0, psd18m = psd18m/counter18m; end
if counter18u ~= 0, psd18u = psd18u/counter18u; end
if counter19m ~= 0, psd19m = psd19m/counter19m; end
if counter19u ~= 0, psd19u = psd19u/counter19u; end
if counter20m ~= 0, psd20m = psd20m/counter20m; end
if counter20u ~= 0, psd20u = psd20u/counter20u; end
if counter21m ~= 0, psd21m = psd21m/counter21m; end
if counter21u ~= 0, psd21u = psd21u/counter21u; end
if counter22m ~= 0, psd22m = psd22m/counter22m; end
if counter22u ~= 0, psd22u = psd22u/counter22u; end
if counter23m ~= 0, psd23m = psd23m/counter23m; end
if counter23u ~= 0, psd23u = psd23u/counter23u; end
if counter24m ~= 0, psd24m = psd24m/counter24m; end
if counter24u ~= 0, psd24u = psd24u/counter24u; end
if counter25m ~= 0, psd25m = psd25m/counter25m; end
if counter25u ~= 0, psd25u = psd25u/counter25u; end
if counter26m ~= 0, psd26m = psd26m/counter26m; end
if counter26u ~= 0, psd26u = psd26u/counter26u; end
if counter27m ~= 0, psd27m = psd27m/counter27m; end
if counter27u ~= 0, psd27u = psd27u/counter27u; end
if counter28m ~= 0, psd28m = psd28m/counter28m; end
if counter28u ~= 0, psd28u = psd28u/counter28u; end
if counter29m ~= 0, psd29m = psd29m/counter29m; end
if counter29u ~= 0, psd29u = psd29u/counter29u; end
if counter30m ~= 0, psd30m = psd30m/counter30m; end
if counter30u ~= 0, psd30u = psd30u/counter30u; end
if counter31m ~= 0, psd31m = psd31m/counter31m; end
if counter31u ~= 0, psd31u = psd31u/counter31u; end
if counter32m ~= 0, psd32m = psd32m/counter32m; end
if counter32u ~= 0, psd32u = psd32u/counter32u; end

%%
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

% relative theta freq values for each of the 32 songs in order
thetasong = [5.266666667,4.533333333,6.2,5.2,7.266666667,8.2,4.8,5.8,4.866666667,5.333333333,5.066666667,5.333333333,6.466666667,5.266666667,7.733333333,8,7.6,4.666666667,8.733333333,4.933333333,7.2,6.4,8,6.333333333,4.866666667,5.6,8.2,7.266666667,6.933333333,6.133333333,4.533333333,5.133333333];

songtitles = {'Hello'  %titles of all the songs
'Someone Like You'
'Girl on Fire'
'Complicated'
"(I've Had) The Time of My Life"
'Livin on Prayer'
"I'll Make Love To You"
'Versace On The Floor'
'When I Was Your Man'
'See You Again'
'Beautiful'
'I Will Follow You Into The Dark'
'Perfect'
'Thinking Aloud'
'All By Myself'
'All Of Me'
'Issues'
'Because Of You'
'Since You Been Gone'
'Praying'
'Need You Now'
'Shallow'
'Wrecking Ball'
'Just Give Me A Reason'
'Hallelujah'
'Stay With Me'
'Not Gonna Write You A Love Song'
'Eye of the Tiger'
'Hey There Delilah'
'Africa'
'I Will Always Love You'
'I Have Nothing'};

%%
%plotting

%opens eeglab and opens one of the participants epoched file. this is not
%for plotting just so that the channel locations are loaded in the EEG
%varible
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[ID '_filt_reref_resamp_rej_interp_prunedICA_epochs.set'],'filepath',fullpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% song plot loop
for k = 1:length(songtitles)
    % if current song matches then load the psd into a temp variable
    if k == 1, psdkm = psd1m; psdku = psd1u; end
    if k == 2, psdkm = psd2m; psdku = psd2u; end
    if k == 3, psdkm = psd3m; psdku = psd3u; end
    if k == 4, psdkm = psd4m; psdku = psd4u; end
    if k == 5, psdkm = psd5m; psdku = psd5u; end
    if k == 6, continue; end %this song is not included
    if k == 7, psdkm = psd7m; psdku = psd7u; end
    if k == 8, psdkm = psd8m; psdku = psd8u; end
    if k == 9, psdkm = psd9m; psdku = psd9u; end
    if k == 10, psdkm = psd10m; psdku = psd10u; end
    if k == 11, psdkm = psd11m; psdku = psd11u; end
    if k == 12, psdkm = psd12m; psdku = psd12u; end
    if k == 13, psdkm = psd13m; psdku = psd13u; end
    if k == 14, psdkm = psd14m; psdku = psd14u; end
    if k == 15, continue; end %this song is not included
    if k == 16, psdkm = psd16m; psdku = psd16u; end
    if k == 17, psdkm = psd17m; psdku = psd17u; end
    if k == 18, psdkm = psd18m; psdku = psd18u; end
    if k == 19, psdkm = psd19m; psdku = psd19u; end
    if k == 20, psdkm = psd20m; psdku = psd20u; end
    if k == 21, psdkm = psd21m; psdku = psd21u; end
    if k == 22, psdkm = psd22m; psdku = psd22u; end
    if k == 23, psdkm = psd23m; psdku = psd23u; end
    if k == 24, psdkm = psd24m; psdku = psd24u; end
    if k == 25, psdkm = psd25m; psdku = psd25u; end
    if k == 26, psdkm = psd26m; psdku = psd26u; end
    if k == 27, psdkm = psd27m; psdku = psd27u; end
    if k == 28, psdkm = psd28m; psdku = psd28u; end
    if k == 29, psdkm = psd29m; psdku = psd29u; end
    if k == 30, psdkm = psd30m; psdku = psd30u; end
    if k == 31, psdkm = psd31m; psdku = psd31u; end
    if k == 32, psdkm = psd32m; psdku = psd32u; end

    figure; 
    subplot(2,2,[1,2])           %plot all plots on one figure
    semilogx(cfT,mean(psdkm,1),'b')   % plot the modded psd                                              
    hold on
    semilogx(cfT,mean(psdku,1),'r')    %plot the unmodded psd
    legend('Modded','Unmodded');
    xlabel('Hz'),ylabel('PSD');
    title(strcat('Average PSD of Song ',string(k),' : ', songtitles{k}))
    [ d, ix ] = min(abs(cfT-thetasong(k))); %choose the closest bin to the songs theta Hz value
    subplot(2,2,3)
    if ~all(psdkm(:) == 0),topoplot(psdkm(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for modded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])]) %*
    colorbar
    title(strcat("Avg Modded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin ',string(ix),')'));
    subplot(2,2,4)
    if ~all(psdku(:) == 0),topoplot(psdku(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for unmodded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])])
    colorbar
    title(strcat("Avg Unmodded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin '+string(ix),')'));
   
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'.fig')) %save figures as matlab fig
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'.png')) %save figure as png
end


%%

% no loop Each song 1 through 32
% %SONG 1
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd1m,1),'b')
% 
% hold on
% plot(cfT,mean(psd1u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 1 :', songtitles{1}))
% [ d, ix ] = min(abs(cfT-thetasong(1))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd1m(:) == 0),topoplot(psd1m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd1u), mean(psd1m)]), max([mean(psd1u), mean(psd1m)])])
% colorbar
% title("Avg Modded Song 1 at "+string(thetasong(1))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd1u(:) == 0),topoplot(psd1u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd1u), mean(psd1m)]), max([mean(psd1u), mean(psd1m)])])
% colorbar
% title("Avg Unmodded Song 1 at "+string(thetasong(1))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song1.fig')
% saveas(gcf,'Song1.png')
% 
% %song 2
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd2m,1))
% hold on
% plot(cfT,mean(psd2u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 2 :', songtitles{2}))
% [ d, ix ] = min(abs(cfT-thetasong(2))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd2m(:) == 0),topoplot(psd2m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd2u), mean(psd2m)]), max([mean(psd2u), mean(psd2m)])])
% colorbar
% title("Avg Modded Song 2 at "+string(thetasong(2))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd2u(:) == 0),topoplot(psd2u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd2u), mean(psd2m)]), max([mean(psd2u), mean(psd2m)])])
% colorbar
% title("Avg Unmodded Song 2 at "+string(thetasong(2))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song2.fig')
% saveas(gcf,'Song2.png')
% 
% %song 3
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd3m,1))
% hold on
% plot(cfT,mean(psd3u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 3 :', songtitles{3}))
% [ d, ix ] = min(abs(cfT-thetasong(3))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd3m(:) == 0),topoplot(psd3m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd3u), mean(psd3m)]), max([mean(psd3u), mean(psd3m)])])
% colorbar
% title("Avg Modded Song 3 at "+string(thetasong(3))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd3u(:) == 0),topoplot(psd3u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd3u), mean(psd3m)]), max([mean(psd3u), mean(psd3m)])])
% colorbar
% title("Avg Unmodded Song 3 at "+string(thetasong(3))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song3.fig')
% saveas(gcf,'Song3.png')
% 
% %song 4
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd4m,1))
% hold on
% plot(cfT,mean(psd4u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 4 :', songtitles{4}))
% [ d, ix ] = min(abs(cfT-thetasong(4))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd4m(:) == 0),topoplot(psd4m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd4u), mean(psd4m)]), max([mean(psd4u), mean(psd4m)])])
% colorbar
% title("Avg Modded Song 4 at "+string(thetasong(4))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd4u(:) == 0),topoplot(psd4u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd4u), mean(psd4m)]), max([mean(psd4u), mean(psd4m)])])
% colorbar
% title("Avg Unmodded Song 4 at "+string(thetasong(4))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song4.fig')
% saveas(gcf,'Song4.png')
% 
% %
% %SONG 5
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd5m,1),'b')
% hold on
% plot(cfT,mean(psd5u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 5 :', songtitles{5}))
% [ d, ix ] = min(abs(cfT-thetasong(5))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd5m(:) == 0),topoplot(psd5m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd5u), mean(psd5m)]), max([mean(psd5u), mean(psd5m)])])
% colorbar
% title("Avg Modded Song 5 at "+string(thetasong(5))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd5u(:) == 0),topoplot(psd5u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd5u), mean(psd5m)]), max([mean(psd5u), mean(psd5m)])])
% colorbar
% title("Avg Unmodded Song 5 at "+string(thetasong(5))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song5.fig')
% saveas(gcf,'Song5.png')
% 
% %SONG 6
% % figure; 
% % subplot(2,2,[1,2])
% % plot(cfT,mean(psd6m,1),'b')
% % hold on
% % plot(cfT,mean(psd6u,1),'r')
% % legend('Modded','Unmodded');
% % xlabel('Hz'),ylabel('PSD');
% % title(strcat('Average PSD of Song 6 :', songtitles{6}))
% % [ d, ix ] = min(abs(cfT-thetasong(6))); %choose the closest bin to the songs theta Hz value
% % subplot(2,2,3)
% % if ~all(psd6m(:) == 0),topoplot(psd6m(:,ix),EEG.chanlocs,'electrodes','off');clim([min([mean(psd6u), mean(psd6m)]), max([mean(psd6u), mean(psd6m)])]);colorbar;end
% % title("Avg Modded Song 6 at "+string(thetasong(6))+'Hz (bin '+string(ix)+')');
% % subplot(2,2,4)
% % if ~all(psd6u(:) == 0),topoplot(psd6u(:,ix),EEG.chanlocs,'electrodes','off');clim([min([mean(psd6u), mean(psd6m)]), max([mean(psd6u), mean(psd6m)])]);colorbar;end
% % title("Avg Unmodded Song 6 at "+string(thetasong(6))+'Hz (bin '+string(ix)+')');
% % saveas(gcf,'Song6.fig')
% % saveas(gcf,'Song6.png')
% 
% %SONG 7
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd7m,1),'b')
% hold on
% plot(cfT,mean(psd7u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 7 :', songtitles{7}))
% [ d, ix ] = min(abs(cfT-thetasong(7))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd7m(:) == 0),topoplot(psd7m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd7u), mean(psd5m)]), max([mean(psd7u), mean(psd7m)])])
% colorbar
% title("Avg Modded Song 7 at "+string(thetasong(7))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd7u(:) == 0),topoplot(psd7u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd7u), mean(psd7m)]), max([mean(psd7u), mean(psd7m)])])
% colorbar
% title("Avg Unmodded Song 7 at "+string(thetasong(7))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song7.fig')
% saveas(gcf,'Song7.png')
% 
% %SONG 8
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd8m,1),'b')
% hold on
% plot(cfT,mean(psd8u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 8 :', songtitles{8}))
% [ d, ix ] = min(abs(cfT-thetasong(8))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd8m(:) == 0),topoplot(psd8m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd8u), mean(psd8m)]), max([mean(psd8u), mean(psd8m)])])
% colorbar
% title("Avg Modded Song 8 at "+string(thetasong(8))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd8u(:) == 0),topoplot(psd8u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd8u), mean(psd8m)]), max([mean(psd8u), mean(psd8m)])])
% colorbar
% title("Avg Unmodded Song 8 at "+string(thetasong(8))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song8.fig')
% saveas(gcf,'Song8.png')
% 
% %SONG 9
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd9m,1),'b')
% hold on
% plot(cfT,mean(psd9u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 9 :', songtitles{9}))
% [ d, ix ] = min(abs(cfT-thetasong(9))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd9m(:) == 0),topoplot(psd9m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd9u), mean(psd9m)]), max([mean(psd9u), mean(psd9m)])])
% colorbar
% title("Avg Modded Song 9 at "+string(thetasong(9))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd9u(:) == 0),topoplot(psd9u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd9u), mean(psd9m)]), max([mean(psd9u), mean(psd9m)])])
% colorbar
% title("Avg Unmodded Song 9 at "+string(thetasong(9))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song9.fig')
% saveas(gcf,'Song9.png')
% 
% %SONG 10
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd10m,1),'b')
% hold on
% plot(cfT,mean(psd10u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 10 :', songtitles{10}))
% [ d, ix ] = min(abs(cfT-thetasong(10))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd10m(:) == 0),topoplot(psd10m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd10u), mean(psd10m)]), max([mean(psd10u), mean(psd10m)])])
% colorbar
% title("Avg Modded Song 10 at "+string(thetasong(10))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd10u(:) == 0),topoplot(psd10u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd10u), mean(psd10m)]), max([mean(psd10u), mean(psd10m)])])
% colorbar
% title("Avg Unmodded Song 10 at "+string(thetasong(10))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song10.fig')
% saveas(gcf,'Song10.png')
% 
% %SONG 11
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd11m,1),'b')
% hold on
% plot(cfT,mean(psd11u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 11 :', songtitles{11}))
% [ d, ix ] = min(abs(cfT-thetasong(11))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd11m(:) == 0),topoplot(psd11m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd11u), mean(psd11m)]), max([mean(psd11u), mean(psd11m)])])
% colorbar
% title("Avg Modded Song 11 at "+string(thetasong(11))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd11u(:) == 0),topoplot(psd11u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd11u), mean(psd11m)]), max([mean(psd11u), mean(psd11m)])])
% colorbar
% title("Avg Unmodded Song 11 at "+string(thetasong(11))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song11.fig')
% saveas(gcf,'Song11.png')
% 
% %SONG 12
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd12m,1),'b')
% hold on
% plot(cfT,mean(psd12u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 12 :', songtitles{12}))
% [ d, ix ] = min(abs(cfT-thetasong(12))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd12m(:) == 0),topoplot(psd12m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd12u), mean(psd12m)]), max([mean(psd12u), mean(psd12m)])])
% colorbar
% title("Avg Modded Song 12 at "+string(thetasong(12))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd12u(:) == 0),topoplot(psd12u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd12u), mean(psd12m)]), max([mean(psd12u), mean(psd12m)])])
% colorbar
% title("Avg Unmodded Song 12 at "+string(thetasong(12))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song12.fig')
% saveas(gcf,'Song12.png')
% 
% %SONG 13
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd13m,1),'b')
% hold on
% plot(cfT,mean(psd13u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 13 :', songtitles{13}))
% [ d, ix ] = min(abs(cfT-thetasong(13))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd13m(:) == 0),topoplot(psd13m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd13u), mean(psd13m)]), max([mean(psd13u), mean(psd13m)])])
% colorbar
% title("Avg Modded Song 13 at "+string(thetasong(13))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd13u(:) == 0),topoplot(psd13u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd13u), mean(psd13m)]), max([mean(psd13u), mean(psd13m)])])
% colorbar
% title("Avg Unmodded Song 13 at "+string(thetasong(13))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song13.fig')
% saveas(gcf,'Song13.png')
% 
% %SONG 14
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd14m,1),'b')
% hold on
% plot(cfT,mean(psd14u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 14 :', songtitles{14}))
% [ d, ix ] = min(abs(cfT-thetasong(14))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd14m(:) == 0),topoplot(psd14m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd14u), mean(psd14m)]), max([mean(psd14u), mean(psd14m)])])
% colorbar
% title("Avg Modded Song 14 at "+string(thetasong(14))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd14u(:) == 0),topoplot(psd14u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd14u), mean(psd14m)]), max([mean(psd14u), mean(psd14m)])])
% colorbar
% title("Avg Unmodded Song 14 at "+string(thetasong(14))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song14.fig')
% saveas(gcf,'Song14.png')
% 
% %SONG 15
% % figure; 
% % subplot(2,2,[1,2])
% % plot(cfT,mean(psd15m,1),'b')
% % hold on
% % plot(cfT,mean(psd15u,1),'r')
% % legend('Modded','Unmodded');
% % xlabel('Hz'),ylabel('PSD');
% % title(strcat('Average PSD of Song 15 :', songtitles{15}))
% % [ d, ix ] = min(abs(cfT-thetasong(15))); %choose the closest bin to the songs theta Hz value
% % subplot(2,2,3)
% % if ~all(psd15m(:) == 0),topoplot(psd15m(:,ix),EEG.chanlocs,'electrodes','off');clim([min([mean(psd15u), mean(psd15m)]), max([mean(psd15u), mean(psd15m)])]);colorbar;end
% % title("Avg Modded Song 15 at "+string(thetasong(15))+'Hz (bin '+string(ix)+')');
% % subplot(2,2,4)
% % if ~all(psd15u(:) == 0),topoplot(psd15u(:,ix),EEG.chanlocs,'electrodes','off');clim([min([mean(psd15u), mean(psd15m)]), max([mean(psd15u), mean(psd15m)])]);colorbar;end
% % title("Avg Unmodded Song 15 at "+string(thetasong(15))+'Hz (bin '+string(ix)+')');
% % saveas(gcf,'Song15.fig')
% % saveas(gcf,'Song15.png')
% 
% %SONG 16
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd16m,1),'b')
% hold on
% plot(cfT,mean(psd16u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 16 :', songtitles{16}))
% [ d, ix ] = min(abs(cfT-thetasong(16))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd16m(:) == 0),topoplot(psd16m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd16u), mean(psd16m)]), max([mean(psd16u), mean(psd16m)])])
% colorbar
% title("Avg Modded Song 16 at "+string(thetasong(16))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd16u(:) == 0),topoplot(psd16u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd16u), mean(psd16m)]), max([mean(psd16u), mean(psd16m)])])
% colorbar 
% title("Avg Unmodded Song 16 at "+string(thetasong(16))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song16.fig')
% saveas(gcf,'Song16.png')
% 
% %SONG 17
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd17m,1),'b')
% hold on
% plot(cfT,mean(psd17u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 17 :', songtitles{17}))
% [ d, ix ] = min(abs(cfT-thetasong(17))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd17m(:) == 0),topoplot(psd17m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd17u), mean(psd17m)]), max([mean(psd17u), mean(psd17m)])])
% colorbar
% title("Avg Modded Song 17 at "+string(thetasong(17))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd17u(:) == 0),topoplot(psd17u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd17u), mean(psd17m)]), max([mean(psd17u), mean(psd17m)])])
% colorbar
% title("Avg Unmodded Song 17 at "+string(thetasong(17))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song17.fig')
% saveas(gcf,'Song17.png')
% 
% %SONG 18
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd18m,1),'b')
% hold on
% plot(cfT,mean(psd18u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 18 :', songtitles{18}))
% [ d, ix ] = min(abs(cfT-thetasong(18))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd18m(:) == 0),topoplot(psd18m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd18u), mean(psd18m)]), max([mean(psd18u), mean(psd18m)])])
% colorbar
% title("Avg Modded Song 18 at "+string(thetasong(18))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd18u(:) == 0),topoplot(psd18u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd18u), mean(psd18m)]), max([mean(psd18u), mean(psd18m)])])
% colorbar
% title("Avg Unmodded Song 18 at "+string(thetasong(18))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song18.fig')
% saveas(gcf,'Song18.png')
% 
% %SONG 19
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd19m,1),'b')
% hold on
% plot(cfT,mean(psd19u,1),'r')
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 19 :', songtitles{19}))
% [ d, ix ] = min(abs(cfT-thetasong(19))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd19m(:) == 0),topoplot(psd19m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd19u), mean(psd19m)]), max([mean(psd19u), mean(psd19m)])])
% colorbar
% title("Avg Modded Song 19 at "+string(thetasong(19))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd19u(:) == 0),topoplot(psd19u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd19u), mean(psd19m)]), max([mean(psd19u), mean(psd19m)])])
% colorbar
% title("Avg Unmodded Song 19 at "+string(thetasong(19))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song19.fig')
% saveas(gcf,'Song19.png')
% 
% %song 20
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd20m,1))
% hold on
% plot(cfT,mean(psd20u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 20 :', songtitles{20}))
% [ d, ix ] = min(abs(cfT-thetasong(20))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd20m(:) == 0),topoplot(psd20m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd20u), mean(psd20m)]), max([mean(psd20u), mean(psd20m)])])
% colorbar
% title("Avg Modded Song 20 at "+string(thetasong(20))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd20u(:) == 0),topoplot(psd20u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd20u), mean(psd20m)]), max([mean(psd20u), mean(psd20m)])])
% colorbar
% title("Avg Unmodded Song 20 at "+string(thetasong(20))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song20.fig')
% saveas(gcf,'Song20.png')
% 
% %song 21
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd21m,1))
% hold on
% plot(cfT,mean(psd21u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 21 :', songtitles{21}))
% [ d, ix ] = min(abs(cfT-thetasong(21))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd21m(:) == 0),topoplot(psd21m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd21u), mean(psd21m)]), max([mean(psd21u), mean(psd21m)])])
% colorbar
% title("Avg Modded Song 21 at "+string(thetasong(21))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd21u(:) == 0),topoplot(psd21u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd21u), mean(psd21m)]), max([mean(psd21u), mean(psd21m)])])
% colorbar
% title("Avg Unmodded Song 21 at "+string(thetasong(21))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song21.fig')
% saveas(gcf,'Song21.png')
% 
% %song 22
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd22m,1))
% hold on
% plot(cfT,mean(psd22u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 22 :', songtitles{22}))
% [ d, ix ] = min(abs(cfT-thetasong(22))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd22m(:) == 0),topoplot(psd22m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd22u), mean(psd22m)]), max([mean(psd22u), mean(psd22m)])])
% colorbar
% title("Avg Modded Song 22 at "+string(thetasong(22))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd22u(:) == 0),topoplot(psd22u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd22u), mean(psd22m)]), max([mean(psd22u), mean(psd22m)])])
% colorbar
% title("Avg Unmodded Song 22 at "+string(thetasong(22))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song22.fig')
% saveas(gcf,'Song22.png')
% 
% %song 23
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd23m,1))
% hold on
% plot(cfT,mean(psd23u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 23 :', songtitles{23}))
% [ d, ix ] = min(abs(cfT-thetasong(23))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd23m(:) == 0),topoplot(psd23m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd23u), mean(psd23m)]), max([mean(psd23u), mean(psd23m)])])
% colorbar
% title("Avg Modded Song 23 at "+string(thetasong(23))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd23u(:) == 0),topoplot(psd23u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd23u), mean(psd23m)]), max([mean(psd23u), mean(psd23m)])])
% colorbar
% title("Avg Unmodded Song 23 at "+string(thetasong(23))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song23.fig')
% saveas(gcf,'Song23.png')
% 
% %song 24
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd24m,1))
% hold on
% plot(cfT,mean(psd24u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 24 :', songtitles{24}))
% [ d, ix ] = min(abs(cfT-thetasong(24))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd24m(:) == 0),topoplot(psd24m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd24u), mean(psd24m)]), max([mean(psd24u), mean(psd24m)])])
% colorbar
% title("Avg Modded Song 24 at "+string(thetasong(24))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd24u(:) == 0),topoplot(psd24u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd24u), mean(psd24m)]), max([mean(psd24u), mean(psd24m)])])
% colorbar
% title("Avg Unmodded Song 24 at "+string(thetasong(24))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song24.fig')
% saveas(gcf,'Song24.png')
% 
% %song 25
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd25m,1))
% hold on
% plot(cfT,mean(psd25u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 25 :', songtitles{25}))
% [ d, ix ] = min(abs(cfT-thetasong(25))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd25m(:) == 0),topoplot(psd25m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd25u), mean(psd25m)]), max([mean(psd25u), mean(psd25m)])])
% colorbar
% title("Avg Modded Song 25 at "+string(thetasong(2))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd25u(:) == 0),topoplot(psd25u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd25u), mean(psd25m)]), max([mean(psd25u), mean(psd25m)])])
% colorbar
% title("Avg Unmodded Song 25 at "+string(thetasong(25))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song25.fig')
% saveas(gcf,'Song25.png')
% 
% %song 26
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd26m,1))
% hold on
% plot(cfT,mean(psd26u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 26 :', songtitles{26}))
% [ d, ix ] = min(abs(cfT-thetasong(26))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd26m(:) == 0),topoplot(psd26m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd26u), mean(psd26m)]), max([mean(psd26u), mean(psd26m)])])
% colorbar
% title("Avg Modded Song 26 at "+string(thetasong(26))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd26u(:) == 0),topoplot(psd26u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd26u), mean(psd26m)]), max([mean(psd26u), mean(psd26m)])])
% colorbar
% title("Avg Unmodded Song 26 at "+string(thetasong(26))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song26.fig')
% saveas(gcf,'Song26.png')
% 
% %song 27
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd27m,1))
% hold on
% plot(cfT,mean(psd27u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 27 :', songtitles{27}))
% [ d, ix ] = min(abs(cfT-thetasong(27))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd27m(:) == 0),topoplot(psd27m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd27u), mean(psd27m)]), max([mean(psd27u), mean(psd27m)])])
% colorbar
% title("Avg Modded Song 27 at "+string(thetasong(27))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd27u(:) == 0),topoplot(psd27u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd27u), mean(psd27m)]), max([mean(psd27u), mean(psd27m)])])
% colorbar
% title("Avg Unmodded Song 27 at "+string(thetasong(27))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song27.fig')
% saveas(gcf,'Song27.png')
% 
% %song 28
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd28m,1))
% hold on
% plot(cfT,mean(psd28u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 28 :', songtitles{28}))
% [ d, ix ] = min(abs(cfT-thetasong(28))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd28m(:) == 0),topoplot(psd28m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd28u), mean(psd28m)]), max([mean(psd28u), mean(psd28m)])])
% colorbar
% title("Avg Modded Song 28 at "+string(thetasong(28))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd28u(:) == 0),topoplot(psd28u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd28u), mean(psd28m)]), max([mean(psd28u), mean(psd28m)])])
% colorbar
% title("Avg Unmodded Song 28 at "+string(thetasong(28))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song28.fig')
% saveas(gcf,'Song28.png')
% 
% %song 29
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd29m,1))
% hold on
% plot(cfT,mean(psd29u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 29 :', songtitles{29}))
% [ d, ix ] = min(abs(cfT-thetasong(29))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd29m(:) == 0),topoplot(psd29m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd29u), mean(psd29m)]), max([mean(psd29u), mean(psd29m)])])
% colorbar
% title("Avg Modded Song 29 at "+string(thetasong(29))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd29u(:) == 0),topoplot(psd29u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd29u), mean(psd29m)]), max([mean(psd29u), mean(psd29m)])])
% colorbar
% title("Avg Unmodded Song 29 at "+string(thetasong(29))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song29.fig')
% saveas(gcf,'Song29.png')
% 
% %song 30
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd30m,1))
% hold on
% plot(cfT,mean(psd30u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 30 :', songtitles{30}))
% [ d, ix ] = min(abs(cfT-thetasong(30))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd30m(:) == 0),topoplot(psd30m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd30u), mean(psd30m)]), max([mean(psd30u), mean(psd30m)])])
% colorbar
% title("Avg Modded Song 30 at "+string(thetasong(30))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd30u(:) == 0),topoplot(psd30u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd30u), mean(psd30m)]), max([mean(psd30u), mean(psd30m)])])
% colorbar
% title("Avg Unmodded Song 30 at "+string(thetasong(30))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song30.fig')
% saveas(gcf,'Song30.png')
% 
% %song 31
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd31m,1))
% hold on
% plot(cfT,mean(psd31u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 31 :', songtitles{31}))
% [ d, ix ] = min(abs(cfT-thetasong(31))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd31m(:) == 0),topoplot(psd31m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd31u), mean(psd31m)]), max([mean(psd31u), mean(psd31m)])])
% colorbar
% title("Avg Modded Song 31 at "+string(thetasong(31))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd31u(:) == 0),topoplot(psd31u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd31u), mean(psd31m)]), max([mean(psd31u), mean(psd31m)])])
% colorbar
% title("Avg Unmodded Song 31 at "+string(thetasong(31))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song31.fig')
% saveas(gcf,'Song31.png')
% 
% %song 32
% figure; 
% subplot(2,2,[1,2])
% plot(cfT,mean(psd32m,1))
% hold on
% plot(cfT,mean(psd32u,1))
% legend('Modded','Unmodded');
% xlabel('Hz'),ylabel('PSD');
% title(strcat('Average PSD of Song 32 :', songtitles{32}))
% [ d, ix ] = min(abs(cfT-thetasong(32))); %choose the closest bin to the songs theta Hz value
% subplot(2,2,3)
% if ~all(psd32m(:) == 0),topoplot(psd32m(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd32u), mean(psd32m)]), max([mean(psd32u), mean(psd32m)])])
% colorbar
% title("Avg Modded Song 32 at "+string(thetasong(32))+'Hz (bin '+string(ix)+')');
% subplot(2,2,4)
% if ~all(psd3u(:) == 0),topoplot(psd32u(:,ix),EEG.chanlocs,'electrodes','off');end
% clim([min([mean(psd32u), mean(psd32m)]), max([mean(psd32u), mean(psd32m)])])
% colorbar
% title("Avg Unmodded Song 32 at "+string(thetasong(32))+'Hz (bin '+string(ix)+')');
% saveas(gcf,'Song32.fig')
% saveas(gcf,'Song32.png')
% 
% %/


