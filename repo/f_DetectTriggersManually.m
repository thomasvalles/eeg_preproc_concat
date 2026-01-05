function [peaksX] = f_DetectTriggersManually(data, srate, startEndTime, option)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Option = 1, use demeaned data
%Option = 2, use detrended data

% substract mean
% aveData = mean(data);
% dataDemeaned = data - aveData;


%YPost2 = detrendnonlin((data-mean(data)),2)';
if (option == 1)
    DataCorrected = (data - mean(data));
    %DataCorrected=data;
elseif (option == 2)
    DataCorrected = (detrendnonlin((data - mean(data)), 2)');
elseif (option == 3)
    DataCorrected = (data - mean(data));
elseif (option ==4)
    DataCorrected = (detrendnonlin((data - mean(data)), 2)');



end

figure;
% plot((1:numel(data))/srate, data);
plot((1:numel(data))/srate, DataCorrected);


thresh = input('threshold for peak detection:');
if option == 3 || option ==4 
%close all
    %DataCorrected = (detrendnonlin(data - mean(data)));
    %plot(DataCorrected)
    if thresh < 0
        suprathresh = find(DataCorrected < thresh);
    else
        suprathresh = find(DataCorrected > thresh);
    end
    endpoints = suprathresh(find(diff(suprathresh/srate) > 5));
    startpoints = [suprathresh(1), suprathresh(find(diff(suprathresh/srate) > 5)+1)];
    numtrains = numel(startpoints);
    endpoints(numtrains) = suprathresh(end);
    startpoints = startpoints;
    endpoints = endpoints + 34;

   
    x = uniquetol(suprathresh, 1/20*srate, 'DataScale', 1);
    pulses = [];
    for itrain = 1:numel(startpoints)
        numpulses(itrain) = numel(find(x <= endpoints(itrain) & x >= startpoints(itrain)));
        pulses = [pulses, find(x <= endpoints(itrain) & x >= startpoints(itrain))];
    end
    table(startpoints', endpoints', numpulses')
    remidx = find(numpulses < 37 | numpulses>41);
    startpoints(remidx) = [];
    endpoints(remidx) = [];

    numtrains = numel(startpoints);
    pulses_final = [];
    numpulses_final = [];
    for itrain = 1:numel(startpoints)
        numpulses_final(itrain) = numel(find(x <= endpoints(itrain) & x >= startpoints(itrain)));
        pulses_final = [pulses_final, x(x <= endpoints(itrain) & x >= startpoints(itrain))];
    end


    table(startpoints', endpoints', numpulses_final')
    x = pulses_final;

else

    if (thresh > 0)
        [~, x] = findpeaks(DataCorrected, 'minpeakheight', thresh);
    else
        [~, x] = findpeaks(-DataCorrected, 'minpeakheight', -thresh);
    end
end

idxSelectPeaks = x / srate > startEndTime(1) & x / srate < startEndTime(2);
peaksX = x(idxSelectPeaks);

end
