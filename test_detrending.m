sub = 'SCC24035';
ses = '01a01';


mat_files = dir(['/Volumes/Leuchter/BIDS_DATABASE/derivatives/int_rs_concat/sub-' sub filesep 'ses-' ses filesep 'Step3_CutFiles/TriggersInterrog.mat']);
eeg_files = dir(['/Volumes/Leuchter/BIDS_DATABASE/rawdata/sub-' sub filesep 'ses-' ses filesep 'eeg/*.set']);

EEG = pop_loadset(fullfile(eeg_files(1).folder, eeg_files(1).name));
load(fullfile(mat_files(1).folder, mat_files(1).name));

EEG_int = pop_select(EEG, 'point', [segments(2).start, segments(2).end]);
EEGdata = EEG_int.data(5, :);
srate = EEG_int.srate;

%%
iFreq = 50;
burst_start = burstStartEnd(iFreq, 1);
burst_end = burstStartEnd(iFreq, 2);

% duration of the artifact to be cut
pre_sec = 0.005;
post_sec = 0.020;

%%
[YPre, YPost] = f_DemeanPreAndPost(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate);

step = round(srate / 5);%400;
small_step = round(srate / (200/3));%30;
small_end = round(srate / 20);%100
large_step = srate; % try one second
large_start = srate * 2; % starting at the second second?

x = (1:numel(YPre)) / srate;
a1 = subplot(2, 4, 1);
hold on
plot(x, YPre, 'Color', 'b');
p = polyfit(((1:numel(YPre)) / srate)', YPre(:), 2);
plot(x, polyval(p, x), 'Color', 'r', 'LineWidth', 2);
title('4s pre demeaned, quadratic fit');

a2 = subplot(2, 4, 2);
hold on
pred = YPre - polyval(p, x);
plot(x, pred, 'Color', 'b');
title('4s pre demeaned, detrended')

x = (1:numel(YPost)) / srate;
a3 = subplot(2, 4, 5);
hold on
plot(x, YPost, 'Color', 'b')

xp = linspace(1, numel(YPost), numel(YPost))';

YPost_mean_smallstep = movmean(YPost, small_step)';
YPost_mean = movmean(YPost, step)';
YPost_mean_largestep = movmean(YPost, large_step)';

x_sample = xp(1:step:end);
x_sample_smallstep = xp(1:small_step:end);
x_sample_largestep = xp(1:large_step:end);
x_sample = [x_sample_smallstep(1:ceil(small_end/small_step));
            x_sample( (x_sample > x_sample_smallstep(ceil(small_end/small_step))) & (x_sample < x_sample(ceil(large_start/step))));
            x_sample_largestep(x_sample_largestep > x_sample(ceil(large_start/step)))];

x_sample(end + 1) = xp(end);

y_sample = [YPost_mean_smallstep(1:small_step:small_end);
            YPost_mean(x_sample(ceil(small_end/small_step) + 1 : ceil(large_start/step)));
            YPost_mean_largestep(x_sample(ceil(large_start/step) + 1 : end))];

scatter(x_sample / srate, y_sample, 'filled', 'red');
p = pchip(x_sample, y_sample, xp);

plot(x, p, 'Color', 'r');
title('YPost raw')


detrended = YPost' - p;
a4 = subplot(2, 4, 6);
plot(x, detrended, 'Color', 'b');
title('YPost pchip detrend')

slope_change_chunk_size = 10; % should be pretty small. after pchip the steep rise doesn't go
                              % for more than 30is samples

[ddetrended, end_sample] = f_DetrendQuadraticSlopeChange(detrended', slope_change_chunk_size, 5);
a5 = subplot(2, 4, 7);
plot(x, ddetrended, 'Color', 'b');

linkaxes([a1, a2, a3, a4, a5], 'x');
title('YPost detrended final')

% [interpolData] = f_FixArtifact(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, device);
 srate = EEG.srate;
winlength = srate; % at least 4 seconds of data to have a reasonable resolution.
fftlength = 2^nextpow2(winlength);
s_FreqRes = srate/fftlength;
nFreq = numel(0:s_FreqRes:srate/2);
overlap = 0;

% % pre power
% [v_PSD_pre,v_Freq_pre] = pwelch(YPre',...
%     winlength,overlap,fftlength,srate);
% 
% [v_PSD_pred,v_Freq_pred] = pwelch(pred',...
%     winlength,overlap,fftlength,srate);
% 
% v_PSD_pre = v_PSD_pre((v_Freq_pre >= 1) & (v_Freq_pre <= 20));
% v_PSD_pred = v_PSD_pred((v_Freq_pred >= 1) & (v_Freq_pred <= 20));
% v_PSD_pre = 100 * v_PSD_pre / sum(v_PSD_pre);
% v_PSD_pred = 100 * v_PSD_pred / sum(v_PSD_pred);
% subplot(2, 4, 4)
% hold on
% plot(v_Freq_pre((v_Freq_pre >= 1) & (v_Freq_pre <= 20)), v_PSD_pre);
% plot(v_Freq_pred((v_Freq_pred >= 1) & (v_Freq_pred < 20)), v_PSD_pred);
% legend('Raw', 'Detrended')
% 
% % post power 
% [v_PSD_post,v_Freq_post] = pwelch(YPost',...
%     winlength,overlap,fftlength,srate);
% 
% [v_PSD_postd,v_Freq_postd] = pwelch(ddetrended',...
%     winlength,overlap,fftlength,srate);
% 
% subplot(2, 4, 8)
% hold on
% v_PSD_post = v_PSD_post((v_Freq_post >= 1) & (v_Freq_post <= 20));
% v_PSD_postd = v_PSD_postd((v_Freq_postd >= 1) & (v_Freq_postd <= 20));
% v_PSD_post = 100 * v_PSD_post / sum(v_PSD_post);
% v_PSD_postd = 100 * v_PSD_postd / sum(v_PSD_postd);
% plot(v_Freq_post((v_Freq_post >= 1) & (v_Freq_post <= 20)), v_PSD_post);
% plot(v_Freq_postd((v_Freq_postd >= 1) & (v_Freq_postd < 20)), v_PSD_postd);
legend('Raw', 'Detrended')
sgtitle("sub-" + sub + ", ses-" + ses + ", Epoch " + iFreq);