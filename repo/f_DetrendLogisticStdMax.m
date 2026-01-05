function [detrended, end_sample, stds] = f_DetrendLogisticStdMax(signal)
% f_DetrendLogisticStdMax - Detrends the signal with a logistic model. End
% sample of detrending is determined by the maximum of std(signal(1:N)). We
% call this function on the signal that has already been quadratically
% detrended
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
    X = 1:size(signal(1:end_sample),2);

    %Force the start of the logistic fit to be close to that actual
    %starting point of the signal
    f = fit(X', double(signal(1:end_sample))', 'logistic4',...
        'Lower', [signal(1), -Inf, -Inf, -Inf],...
        'Upper', [signal(1), Inf, Inf, Inf]);


    % Subtract the fitted values from the original data
    fitted_values = feval(f, X);
    residuals = signal(1:end_sample)' - fitted_values;

    detrended = [residuals' signal(end_sample + 1:end)];