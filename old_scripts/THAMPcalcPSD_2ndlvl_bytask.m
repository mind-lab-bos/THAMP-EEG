
% second level 
% opens all existing psd files from each participant
% AVERAGING ACROSS ALL PARTICIPANTS for both tasks combined

clc
clear all

% edit this ID list to have all the fully analyzed participants
% edit the song order and mod order list to match the order of the idlist
ID_list = {'230926TELZ' '230929SCOE' '231004LLE' '231114CWAY' '231117SSAY' '231120ICON' '231128MYAN' '231205KKAT' '231205NSAZ' '231211TABO' '231212TNAR' '231213DROD'...
    '240117LAYR'...
'240119ACHE'...
'240119AWIL'...
'240124LYUA'...
'240125CSMI'...
'240126ISHI'...
'240206DLY'...
'240207FSUT'...
'240208AFAS'...
'240209DKIM'...
'240215KSOK'...
'240229LMAR'};
songorders = [1,	5; 2,	6;3,	7;1,	5;3,	7;2,	6;4,	8;5,	9;7,	1;6,	10;8,	2;9,	3;...
    1	5;...
10	4;...
2	6;...
3	7;...
5	9;...
4	8;...
6	10;...
7	1;...
1	3;...
8	2;...
10	4;...
1	5];
modorder = {'MOD','UNMOD','MOD','MOD','MOD','UNMOD','UNMOD','MOD','MOD','UNMOD','UNMOD','MOD'...
'MOD'...
'UNMOD'...
'UNMOD'...
'MOD'...
'MOD'...
'UNMOD'...
'UNMOD'...
'MOD'...
'MOD'...
'UNMOD'...
'UNMOD'...
'MOD'};
blockorder = {'SART','2BACK','SART','SART','2BACK','SART','2BACK','SART','2BACK','SART','2BACK','SART',...
    '2BACK',...
'SART',...
'2BACK',...
'SART',...
'2BACK',...
'SART',...
'2BACK',...
'SART',...
'2BACK',...
'SART',...
'2BACK'...
'SART'}; %sart first or nback first

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
19	17	13	9	2	4	18	14	8	24	23	21;
26	22	27	10	3	25	1	7	20	30	5	16;
18	14	8	24	23	21  29	11	31	32	12	28;
26	22	27	10	3	25  1	7	20	30	5	16;
18	14	8	24	23	21  29	11	31	32	12	28;
1	7	20	30	5	16  19	17	13	9	2	4;
29	11	31	32	12	28  26	22	27	10	3	25;
1	7	20	30	5	16  19	17	13	9	2	4;
29	11	31	32	12	28  26	22	27	10	3	25;
19	17	13	9	2	4  18	14	8	24	23	21;
18	14	8	24	23	21  1	7	20	30	5	16;
19	17	13	9	2	4  18	14	8	24	23	21;
26	22	27	10	3	25  1	7	20	30	5	16;
18	14	8	24	23	21  29	11	31	32	12	28];

thetasong = [5.266666667,4.533333333,6.2,5.2,7.266666667,8.2,4.8,5.8,4.866666667,5.333333333,5.066666667,5.333333333,6.466666667,5.266666667,7.733333333,8,7.6,4.666666667,8.733333333,4.933333333,7.2,6.4,8,6.333333333,4.866666667,5.6,8.2,7.266666667,6.933333333,6.133333333,4.533333333,5.133333333];

%CONSTANTS
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
%SART

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

            if strcmp(blockorder{i}, 'SART') && (k == 1 || k == 2 || k == 3 || k == 4 || k == 5 || k == 6)
                task = 'SART';
            elseif strcmp(blockorder{i}, 'SART') && (k == 7 || k == 8 || k == 9 || k == 10 || k == 11 || k == 12)
                task = '2BACK';
            end
            if strcmp(blockorder{i}, '2BACK') && (k == 1 || k == 2 || k == 3 || k == 4 || k == 5 || k == 6)
                task = '2BACK';
            elseif strcmp(blockorder{i}, '2BACK') && (k == 7 || k == 8 || k == 9 || k == 10 || k == 11 || k == 12)
                task = 'SART';
            end

            % j is the current song (1-32), i is the current participant,
            % and k is the current song from that participant (1-12)
            % the two if statements check to see if this current
            % participant listened to the song and if it was modded or
            % unmodded
            if songorder(i,k) == j && modded && strcmp(task,'SART')
                %if the participant listed to that song then add that song
                %to temp variable psdm
                if exist(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat'))) == 2
                    psdm = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                    exist12 = 1;
                else
                    exist12 = 0;
                end
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
                if k == 12 && exist12==1, psdm = psdm.EEGPSD12; end

                %some of the eegpsd.mat files have the last row as all -inf
                %which breaks the plots. this checks to see if the last row
                %is -infs and replaces that row with the row above.
                if exist12 == 1
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

            end

            % this whole if block is the same as above but for unmodded
            % songs
            if songorder(i,k) == j && ~modded && strcmp(task,'SART')
                 %psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                 if exist(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat'))) == 2 %if psd12 doesnt exist it will skip it
                    psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                    exist12 = 1;
                 else
                     exist12 = 0;
                 end

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
                if k == 12 && exist12==1, psdu = psdu.EEGPSD12; end

                %some of the eegpsd.mat files have the last row as all -inf
                %which breaks the plots. this checks to see if the last row
                %is -infs and replaces that row with the row above.
                if exist12 == 1
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
%plotting

%opens eeglab and opens one of the participants epoched file. this is not
%for plotting just so that the channel locations are loaded in the EEG
%varible
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[ID '_filt_reref_resamp_rej_interp_prunedICA_epochs.set'],'filepath',fullpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

close all

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
    hold on
    xline(thetasong(k),'--')%plot the vertical dotted line at the freq bin for this song
    
    legend('Modded','Unmodded');
    xlabel('Hz'),ylabel('PSD');
    title(strcat('SART Average PSD of Song ',string(k),' : ', songtitles{k}))
    [ d, ix ] = min(abs(cfT-thetasong(k))); %choose the closest bin to the songs theta Hz value
    subplot(2,2,3)
    if ~all(psdkm(:) == 0),topoplot(psdkm(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for modded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])]) %*
    colorbar
    title(strcat("SART Avg Modded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin ',string(ix),')'));
    subplot(2,2,4)
    if ~all(psdku(:) == 0),topoplot(psdku(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for unmodded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])])
    colorbar
    title(strcat("SART Avg Unmodded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin '+string(ix),')'));
   
    songtitles{k}

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'_SART.fig')) %save figures as matlab fig
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'_SART.png')) %save figure as png
end

%%
%2BACK

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

            if strcmp(blockorder{i}, 'SART') && (k == 1 || k == 2 || k == 3 || k == 4 || k == 5 || k == 6)
                task = 'SART';
            elseif strcmp(blockorder{i}, 'SART') && (k == 7 || k == 8 || k == 9 || k == 10 || k == 11 || k == 12)
                task = '2BACK';
            end
            if strcmp(blockorder{i}, '2BACK') && (k == 1 || k == 2 || k == 3 || k == 4 || k == 5 || k == 6)
                task = '2BACK';
            elseif strcmp(blockorder{i}, '2BACK') && (k == 7 || k == 8 || k == 9 || k == 10 || k == 11 || k == 12)
                task = 'SART';
            end

            % j is the current song (1-32), i is the current participant,
            % and k is the current song from that participant (1-12)
            % the two if statements check to see if this current
            % participant listened to the song and if it was modded or
            % unmodded
            if songorder(i,k) == j && modded && strcmp(task,'2BACK')
                %if the participant listed to that song then add that song
                %to temp variable psdm
                if exist(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat'))) == 2
                    psdm = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                    exist12 = 1
                else
                    exist12 = 0;
                end

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
                if k == 12 && exist12 ==1, psdm = psdm.EEGPSD12; end
                
                if exist12 == 1
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

            end

            % this whole if block is the same as above but for unmodded
            % songs
            if songorder(i,k) == j && ~modded && strcmp(task,'2BACK')
                 %psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                 if exist(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat'))) == 2
                    psdu = load(fullfile(savepath,strcat(ID,'_EEGPSD',string(k),'.mat')));
                    exist12 = 1;
                 else
                     exist12 = 0;
                 end
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
                if k == 12 && exist12 == 1, psdu = psdu.EEGPSD12; end

                %some of the eegpsd.mat files have the last row as all -inf
                %which breaks the plots. this checks to see if the last row
                %is -infs and replaces that row with the row above.
                if exist12 == 1
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
%plotting

%opens eeglab and opens one of the participants epoched file. this is not
%for plotting just so that the channel locations are loaded in the EEG
%varible
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[ID '_filt_reref_resamp_rej_interp_prunedICA_epochs.set'],'filepath',fullpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

close all

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
    hold on
    xline(thetasong(k),'--');%plot the vertical dotted line at the freq bin for this song%plot the vertical dotted line at the freq bin for this song
    legend('Modded','Unmodded');
    xlabel('Hz'),ylabel('PSD');
    title(strcat('NBACK Average PSD of Song ',string(k),' : ', songtitles{k}))
    [ d, ix ] = min(abs(cfT-thetasong(k))); %choose the closest bin to the songs theta Hz value
    subplot(2,2,3)
    if ~all(psdkm(:) == 0),topoplot(psdkm(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for modded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])]) %*
    colorbar
    title(strcat("NBACK Avg Modded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin ',string(ix),')'));
    subplot(2,2,4)
    if ~all(psdku(:) == 0),topoplot(psdku(:,ix),EEG.chanlocs,'electrodes','off');end %plot the topos for unmodded if not empty
    clim([min([mean(psdku), mean(psdkm)]), max([mean(psdku), mean(psdkm)])])
    colorbar
    title(strcat("NBACK Avg Unmodded Song ",string(k)," at ",string(thetasong(k)),'Hz (bin '+string(ix),')'));
   
    songtitles{k}

    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'_2BACK.fig')) %save figures as matlab fig
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,strcat(songtitles{k},'_2BACK.png')) %save figure as png
end


