function [ EEGcleanedlab, str_ChanLabCleaned] = f_CleanEdfChanLabels(EEG)
% the edf channel labels come in the format: EEG F4-CPz
% this functions cleans them to leave only the central part
% input : EEG struct

EEGcleanedlab = EEG;
str_ChanLabelsOrig = {EEG.chanlocs.labels}';
str_ChanLabCleaned = cell(1,numel(str_ChanLabelsOrig));
for i = 1:numel(str_ChanLabelsOrig)
    tmp1 = strsplit(str_ChanLabelsOrig{i}, '-');
    tmp1 = tmp1{1};
    tmp2 = strsplit(tmp1, ' ');
    str_ChanLabCleaned{i} = tmp2{2};
    EEGcleanedlab.chanlocs(i).labels = tmp2{2};
end

end

