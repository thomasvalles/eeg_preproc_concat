function [YPre, YPost] = f_DemeanPreAndPost(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate)
% f_DemeanPreAndPost- Demeans pre and post burst EEG data in a specific 
%                     epoch separately
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
%   srate:          The sampling rate of the EEG
%
% Output:
%   YPre:          The pre-stimulation signal, demeaned
%
%   YPost2:        The post-stimulation signal, demeaned
%
% T. Valles 2024
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    if(srate == 2000 || srate == 16000 || srate == 500)
        % concatenate 1s before & 2s after pulse
        idxArt2Remove = [burst_start-srate*pre_sec burst_end+srate*post_sec]; % the pulse to be remove
        idxWholeEpoch = [burst_start-srate*1+1 burst_end+srate*3-1]; % whole new epoch period
    elseif(srate == 2048)
        % I think this branch exists bc of an off by one error that will
        % occur if srate isn't 2000.
        
        % concatenate 1s before & 2s after pulse
        idxArt2Remove = [burst_start-srate*pre_sec burst_end+srate*post_sec]; % the pulse to be remove
        idxWholeEpoch = [burst_start-srate*1 burst_end+srate*3-1]; % whole new epoch period
    end

    % demean data for interpolation
    YPre = EEGdata(idxWholeEpoch(1):idxArt2Remove(1)) - mean(EEGdata(idxWholeEpoch(1):idxArt2Remove(1)));
    YPost = EEGdata(idxArt2Remove(2):idxWholeEpoch(2)) - mean(EEGdata(idxArt2Remove(2):idxWholeEpoch(2)));
    

    YPre_diffs = diff(YPre);

    % if there's some big jumping (>30) in the last 10 samples
    bad_samples = find(YPre_diffs > 30);
    bad_sample = min(bad_samples(bad_samples > numel(YPre_diffs) - 10));
    if ~isempty(bad_sample)
        num_rewind = (numel(YPre) - bad_sample) + 5;
        YPre = EEGdata((idxWholeEpoch(1) - num_rewind):idxArt2Remove(1) - num_rewind) - mean(EEGdata((idxWholeEpoch(1) - num_rewind):idxArt2Remove(1) - num_rewind));
        %disp(['bad sample ' num2str(bad_sample) '. num rewind: ' num2str(num_rewind)]);
    end

end