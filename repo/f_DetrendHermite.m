function [ddetrended] = f_DetrendHermite(signal, step, small_step, small_end, large_step, large_start)
    % step = 100;
    % x = linspace(0, numel(signal) - 1, numel(signal))';
    % signal_mean = movmean(signal, step)';
    % x_sample = x(1:step:end);
    % y_sample = signal_mean(1:step:end);
    % p = pchip(x_sample, y_sample, x);
    % 
    % detrended = signal' - p;

    % step = 400; the node spacing later in the signal
    % small_step = 30; the node spacing early in the signal, should be
    % short bc the signal is so steep/variable early on
    % small_end = 100; when to switch to the bigger spacing
    xp = linspace(1, numel(signal), numel(signal))';
    
    % take a moving average over the 3 step sizes
    signal_mean_smallstep = movmean(signal, small_step)';
    signal_mean = movmean(signal, step)';
    signal_mean_largestep = movmean(signal, large_step)';
    
    % lay out the small-, medium-, large- spaced nodes
    x_sample = xp(1:step:end);
    x_sample_smallstep = xp(1:small_step:end);
    x_sample_largestep = xp(1:large_step:end);

    % concatenate nodes based on small_end/large_start
    x_sample = [x_sample_smallstep(1:ceil(small_end/small_step));
                x_sample( (x_sample > x_sample_smallstep(ceil(small_end/small_step))) & (x_sample < x_sample(ceil(large_start/step))));
                x_sample_largestep(x_sample_largestep > x_sample(ceil(large_start/step)))];
    
    % add a node at the end
    x_sample(end + 1) = xp(end);
    
    % get the y-values at the nodes
    y_sample = [signal_mean_smallstep(1:small_step:small_end);
                signal_mean(x_sample(ceil(small_end/small_step) + 1 : ceil(large_start/step)));
                signal_mean_largestep(x_sample(ceil(large_start/step) + 1 : end))];
    
    p = pchip(x_sample, y_sample, xp);


    detrended = signal' - p;
    slope_change_chunk_size = 10; % should be pretty small. after pchip the steep rise doesn't go
                                  % for more than 30is samples

    [ddetrended, end_sample] = f_DetrendQuadraticSlopeChange(detrended', slope_change_chunk_size, 5);
                           

    % try
    %     [detrended, end_sample] = f_DetrendLogisticSlopeChange(detrended', slope_change_chunk_size);
    % catch
    %     try
    %         [detrended, end_sample] = f_DetrendLogisticStdMax(detrended');
    %     catch
    %         [detrended, end_sample] = f_DetrendLinearStdMax(detrended');
    %     end
    % end

     % gap = detrended(end_sample+1) - detrended(end_sample);
     % detrended(end_sample + 1:end) = detrended(end_sample + 1:end) - gap;
end