function [YPre2, YPost21] = f_DemeanAndDetrend(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, order, device)
% f_DemeanAndDetrend - Demeans and detrends EEG data in a specific epoch
%
% Input:
%   EEGdata:        The signal (one channel, entire recording) to detrend.
%
%   pre_sec:        The amount of time (in seconds) to cut before
%                   stimulation
%
%   post_sec:       The amount of time (in seconds) to cut after
%                   stimulation (first one to try)
%
%   burst_start:    Sample at which the burst of interest starts
%
%   burst_end:      Sample at which the burst ends
%
%
%   srate:          The sampling rate of the EEG
%
%   order:          The order of detrending methods to try. Order 1 is
%                   quadratic+log (stdmax) -> quadratic+log (slope change) 
%                   -> full log -> quadratic+linear (stdmax). Order 2 is
%                   quadratic+log (stdmax) -> log+log (stdmax) ->
%                   quadratic+log (slope change) -> full log ->
%                   quadratic+linear(stdmax). Order 3 is pchip
%
%   device:         Device used. will use different pchip parameters for
%                   magstim vs. magventure
%
% Output:
%   YPre2:          The pre-stimulation signal, demeaned and quadratically 
%                   detrended
%
%   YPost21:        The post-stimulation signal, demeaned and detrended
%
% T. Valles 2024
    
    % 1) demean data for interpolation
    [YPre, YPost] = f_DemeanPreAndPost(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate);
    
    % 2) detrend for interpolation
    YPre2 = detrendnonlin(YPre,2)';
    if order == 1
        YPost21 = f_DetrendOrder1(YPost);
    elseif order == 2
        YPost21 = f_DetrendOrder2(YPost);
    elseif order == 3
        if strcmp(device, 'Magventure')
            step = round(srate / 5);%400;
            small_step = round(srate / 40);%50;
            small_end = round(srate / 20) + 1;%101;
        elseif strcmp(device, 'Magstim')
            step = round(srate / 5);%400;
            small_step = round(srate / (200/3));%30;
            small_end = round(srate / 20);%100
        end

        large_step = srate; % try one second
        large_start = srate * 2; % starting at the second second?
        YPost21 = f_DetrendHermite(YPost, step, small_step, small_end, large_step, large_start);
    end
    
    YPost21 = reshape(YPost21, 1, length(YPost21));
    
end