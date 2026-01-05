function [detrended, end_sample, stds] = f_DetrendLinearStdMax(signal)
% f_DetrendLogisticStdMax - Detrends the signal with a linear model. End
% sample of detrending is determined by the maximum of std(signal(1:N)). We
% call this function on the signal that has already been quadratically
% detrended. This is the "last resort" if all the other models failed to
% fit. 
%
% Input:
%   signal:         1*N array representing the input signal. Usually
%                   already detrended in some way
%                   
%
% Output:
%   detrended:      The signal detrended up until end_sample.
%   end_sample:     The sample that the signal was detrended up until
%   stds:           Array of std(signal(1:N)). N starts at 50 and is
%                   incremented by steps of 10
%
% T. Valles 2024

    start = 50;
    curr = start;
    step = 10;
    stds = [];
    while curr + step <= length(signal)
        std_curr = std(signal(1:curr));
        stds = [stds std_curr];
        curr = curr + step;
    end
    
    [m, idx_max] = max(stds);
    end_sample = min((idx_max - 1) * step + start, 2000);
    detrended = [detrend(signal(1:end_sample)) signal(end_sample + 1:end)];
end