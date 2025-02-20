function [eegPSD, cfT] = calcPSD(EEG, lfT, hfT, windowT, overlapT)
    nperT = 120;                               % number of bins per octave
    noctT = log2(hfT)-log2(lfT);               % number of octaves
    NT = noctT*nperT+1;                        % Total number of bins
    cfT = logspace(log10(lfT),log10(hfT),NT);  % Log spaced bin center frequencies
    
    eegPSD = pwelch(EEG.data', windowT, overlapT, cfT, EEG.srate)'; % PSD
    eegPSD = 10*log10(eegPSD); % convert to dB
end
