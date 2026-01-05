function [aSC] = f_aSC_SeedtoAll(abspower, freq, centerFreq, seedChanIdx, allChanLab, HalfHz)

    % inputs:
    % relpower - matrix with powerspectra channel & frequency 
    % freq - frequency vector for the power matrix
    % centerFreq - for that subject. E.g. IAF. If none - 10Hz default 
    % chans = a pair or matrix of chanLab for power matrix for which to
    % compute aSC - must be nx2 matrix/vector
    % chanLab = channel label for all channels
    
    if ~exist('centerFreq', 'var')
        centerFreq = 10;
    end

    fmin = centerFreq-HalfHz;
    fmax = centerFreq+HalfHz;
    idxf = find(freq >= fmin & freq <= fmax);
    %idxfmax = find(freq <= fmax, 1, 'last');
    
    indices_to_run = [seedChanIdx+1:size(abspower,2)];
    
    aSC = nan(1,size(abspower,2));
    
    spectWave1 = abspower(idxf,seedChanIdx);
    
    tmp = corr([spectWave1, abspower(idxf,indices_to_run)]);
    aSC(indices_to_run) = tmp(2:end,1)';
    
end
