function [detrended, end_sample] = f_DetrendQuadraticSlopeChange(signal, chunk_size, degree)
% f_DetrendLogisticSlopeChange - Detrends the signal with a logistic model. End
% sample of detrending is determined by when the signal starts to change. We
% call this function on the signal that has already been quadratically
% detrended. The end sample for this method is usually less than that of
% the StdMax method, so it may have a better chance at fitting if the first
% choice method fails.
%
% Input:
%   signal:         1*N array representing the input signal. Usually
%                   already detrended in some way
%                   
%
% Output:
%   detrended:      The signal detrended up until end_sample.
%   end_sample:     The sample that the signal was detrended up until
%
% T. Valles 2024
    condition = true;
    %chunk_size = 50;
    start = 1;
    slide_step = 10;
    while condition
        r_curr = corr([1:chunk_size; signal(start:start + chunk_size - 1)]');
        r_curr = r_curr(1, 2);

        r_next = corr([1:chunk_size; signal(start + chunk_size: start + 2*chunk_size -1)]');
        r_next = r_next(1, 2); 
      
        start = start + slide_step;

        condition = start + 2*chunk_size <= length(signal);

        % Require some moderate correlation with sample number in one
        % direction followed by a moderate correlation in the opposite
        % direction
        condition = condition && ( ...
            ~( (abs(r_curr) >= 0.3) && ...
            (abs(r_next) >= 0.3) && ...
            sign(r_curr) * sign(r_next) < 0 ...
            ));
    end
    end_sample = min(start - slide_step + chunk_size, 100); %should not got very far here

    first_part = signal(1:end_sample);
    residuals = detrendnonlin(first_part, degree);

    detrended = [residuals' signal(end_sample + 1:end)];
end

