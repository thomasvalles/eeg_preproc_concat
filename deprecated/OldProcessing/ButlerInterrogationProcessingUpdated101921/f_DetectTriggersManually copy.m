function [peaksX] = f_DetectTriggersManually(data,srate, startEndTime, option)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Option = 1, use demeaned data
%Option = 2, use detrended data

% substract mean
% aveData = mean(data);
% dataDemeaned = data - aveData;


%YPost2 = detrendnonlin((data-mean(data)),2)';
if(option == 1)
    DataCorrected = (data-mean(data));
elseif(option == 2)
    DataCorrected = detrendnonlin((data-mean(data)),2)';
end
     
figure;
% plot((1:numel(data))/srate, data);
plot((1:numel(data))/srate, DataCorrected);


threshold = input('threshold for peak detection:');
if(threshold > 0)
    [~, x] = findpeaks(DataCorrected, 'minpeakheight', threshold);
else
    [~, x] = findpeaks(-DataCorrected, 'minpeakheight', -threshold);
end
idxSelectPeaks = x/srate >startEndTime(1) & x/srate < startEndTime(2);
peaksX = x(idxSelectPeaks);

end

