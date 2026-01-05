function [interpolData] = f_FixArtifact(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, device)
    order = 3; % always use the pchip method.
    [YPre2, YPost21] = f_DemeanAndDetrend(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, order, device);
    % if strcmp(device, 'Magventure')
    %     order = 1;
    %     [YPre2, YPost21] = f_DemeanAndDetrend(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, order);
    % elseif strcmp(device, 'Magstim')
    %     [YPre2, YPost21] = f_FixArtifactMagstim(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate);
    % end

    nInterpol = pre_sec*srate+post_sec*srate;
    nPrePost = 2;

    % interpolate the region we cut out with a line (we cut 5ms pre pulse
    % and 20 or 35ms post pulse)
    [~, ~, ~, YY] = f_CutInterpolArtifact(YPre2, YPost21, nInterpol, nPrePost);
    
    interpolData = [YPre2, YY, YPost21];
end