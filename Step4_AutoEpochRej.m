clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;%clear;
close all;
clc;
eeglab;

Step0_SetDirectories;

%set working directories
str_DirInput = [work_dir subjid filesep 'Step4_FasterOutput' suffix filesep]; %converting the orignal input file data to a string array
str_DirOutput = [work_dir subjid filesep 'Step5_VisualInspectionOutput' suffix filesep]; %converting the output file data to a string array and renaming it
str_TriggerDir = [work_dir subjid filesep 'Step3_BrushedTriggers' suffix filesep];
str_TriggerOut = [work_dir subjid filesep 'Step5_FinalTriggerFiles' suffix filesep];
mkdir(str_TriggerOut);
mkdir(str_DirOutput);
%move to input directory

tmp = dir([str_DirInput '*set']);
listFiles = {tmp.name};
listFiles = transpose(listFiles);
listFiles = listFiles(~contains(string(listFiles),"._"));
numSET = numel(listFiles);


%start at 3 to skip directories '.' and '..'
for i = 1:size(listFiles,1)
    disp(["File: " num2str(i) " of " num2str(size(listFiles,1)) ': ' listFiles{i}]);
    ALLEEG = [];
    % load set files
    filename = strcat(listFiles{i});
    %directory = [ str_DirInput listFiles(i) ];
    EEG = pop_loadset('filepath', str_DirInput,'filename', [filename]);
    EEG = eeg_checkset(EEG);
    
    %%Load trigger file
    tmp_str = strsplit(filename,'_EpochData.set');
    load([str_TriggerDir tmp_str{1} '_TriggersInterrog.mat']);
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    if EEG.nbchan >= 64
        eog_chan =  find(contains({EEG.chanlocs.labels}, 'EOG'));
        channels_to_keep = setdiff(1:size(EEG.data, 1), eog_chan);
        EEG = pop_select(EEG, 'channel', channels_to_keep);
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');

    trialrej = main_epoch_rej(EEG);
    v_LogRej = zeros(1, EEG.trials);
    v_LogRej(trialrej) = 1;

    OUTEEG=EEG;
    OUTEEG.setname = [listFiles(i)];
    EEG = pop_saveset( OUTEEG,'filepath',[str_DirOutput listFiles{i}]);
    
    RemovedFreqs{1} = trialrej;
    RemovedFreqs{2} = burstIPI(trialrej);
    RemovedFreqs{3} = burstStartEnd(trialrej,:);
    RemovedFreqs{4} = stimFreqMedian(trialrej);
    original_f_order = stimFreqMedian;
    save([str_TriggerOut tmp_str{1} '_TriggersInterrog.mat'], 'burstIPI','burstStartEnd','stimFreqMedian','RemovedFreqs','original_f_order','trialrej');
    Removed_matrix=[];
    
end

