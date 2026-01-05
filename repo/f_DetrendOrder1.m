function [YPost21] = f_DetrendOrder1(YPost)
% f_DetrendOrder1 - Detrends YPost, trying some different methods until one
%                   converges.
%
% Input:
%   YPost:          The signal (one channel, one epoch, post stim) to 
%                   detrend.
%
% Output:
%
%   YPost21:        The post-stimulation signal, demeaned and detrended
%
%
% T. Valles 2024

    YPost2 = detrendnonlin(YPost,2)';
    slope_change_chunk_size = 50;

    try % first try log detrend on quad detrend up to max std dev
        [YPost21, end_sample, ~] = f_DetrendLogisticStdMax(YPost2);

        % always stitch the gap between when we stop detrending and the
        % subsequent sample
        gap = YPost21(end_sample+1) - YPost21(end_sample);
        YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;

    catch 
        try % try log detrend on quad detrend up to slope change (usually shorter than prev.)
            [YPost21, end_sample] = f_DetrendLogisticSlopeChange(YPost2, slope_change_chunk_size);
            gap = YPost21(end_sample+1) - YPost21(end_sample);
            YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;
            
        catch
            try % try log detrend on all of post
                YPost21 = f_DetrendLogFull(YPost);
            catch
                [YPost21, end_sample, ~] = f_DetrendLinearStdMax(YPost2);
                gap = YPost21(end_sample+1) - YPost21(end_sample);
                YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;

            end
        end
    end
end