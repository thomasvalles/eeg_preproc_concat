function detrended = f_DetrendLogFull(signal)
% f_DetrendLogFull - Detrends the entire signal with a logistic model. No
% prior detrending necessary
%
% Input:
%   signal:         1*N array representing the input signal.
%                   
%
% Output:
%   detrended:      The signal detrended up until end_sample.
%
% T. Valles 2024

    X = 1:size(signal,2);

    % Force the starting point of the fit to be close to the starting point
    % of the actual data
    f = fit(X', double(signal)', 'logistic4',...
        'Lower', [signal(1), -Inf, -Inf, -Inf],...
        'Upper', [signal(1), Inf, Inf, Inf]);

    fitted_values = feval(f, X);
    detrended = signal' - fitted_values;
       
end