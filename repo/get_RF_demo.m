clear 
cd ~/Documents/RF_StepsMatlab/
EEG = pop_loadset('/Users/neuromodit/Documents/Interrogation/21-173/Step6_FasterOutput/21-173_2021-11-24_14-51-11_Merged_EpochData.set');
chanlocs = {EEG.chanlocs.labels};
subjid = 'SCC22-020-2022-06-28 GABOR';
work_dir = '~/Documents/Interrogation/SCC_foo/';
ChanSet = {{'Fpz'},{'AF4'}};
ChanSet_old = {{'FC3','AF3'},{'FT8','T8','TP8','CP6','P8'}};
ChanSet = {{'Fz'},{'FC2'},{'C3'}}
weights_old = load('/Users/neuromodit/Documents/RF_StepsMatlab/RegressionValues_old.mat').b;
weights_matched = load('/Users/neuromodit/Documents/RF_StepsMatlab/RegressionMatched_01-17-2023.mat').b';

% [weighted_vals_sorted,UniqueFreqs] = get_RF([work_dir, subjid],ChanSet,'repeats','median','smoothed',chanlocs,'matched',weights_matched);
% close all
allsub = readtable('~/Documents/RF_Bruteforce/newInterrogationSubjects&SCC.xlsx');
allsub = allsub(contains(allsub.TMSID,'001-'),:);
for isubj=1:height(allsub)
    close all
    cd ~/Documents/RF_StepsMatlab/

    subj = char(strrep(allsub(isubj,:).TMSID,'001-','SCC22-'))
[weighted_vals_sorted,UniqueFreqs,topRF(isubj,:)] = get_RF([work_dir, subj],ChanSet,0.3,'median','smoothed',chanlocs,'matched',[1,1,-1,0]');
title(subj,'FontSize',25)
xline(allsub.Group(isubj),'LineWidth',2,'Color','g')

saveas(gcf,['~/Documents/RF Paper 03-2023/Fz_FC2_C3/' subj '.jpg'])

end 

out = table(allsub.TMSID,allsub.Group,topRF)


out.Properties.VariableNames = {'ID','Chosen Freq 13-17','Matched Top 5, C3 1hz '}

matches = sum(any(abs(allsub.Group - topRF)<=.5,2))
out.match = any(abs(allsub.Group - topRF)<=.5,2);
