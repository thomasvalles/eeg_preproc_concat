% p_postFASTER_VI

clear;
close all;
clc;
clear;
close all;
clc;
%subjid = 'SCC22-019-2022-06-27';
work_dir = '~/Documents/Interrogation/Butler/';
%set working directories
%str_DirInput = [work_dir subjid '/Step6_FasterOutput/']; %converting the orignal input file data to a string array
%str_DirOutput = [work_dir subjid '/Step7_VisualInspectionOutput/']; %converting the output file data to a string array and renaming it
%str_TriggerDir = [work_dir subjid '/Step3_BrushedTriggers/'];
%str_TriggerOut = [work_dir subjid '/Step7_FinalTriggerFiles/'];

str_DirInput = work_dir;
str_DirOutput = [work_dir 'Output/batch2/'];
str_TriggerDir = work_dir;
str_TriggerOut = [work_dir 'Output/batch2/'];
%move to input directory
cd(str_DirInput); %changing the current folder to a new folder

tmp = dir('*set');
listFiles = {tmp.name};
listFiles = transpose(listFiles);
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

    EEG = eeg_checkset(EEG, 'eventconsistency');
    set(0,'units','pixels');
    scrsz = get(0,'ScreenSize');
    
    command = ...
        [  '[EEGTMP LASTCOM] = eeg_eegrej_AW(EEG,eegplot2event(TMPREJ, -1));' ...
        'v_Rej = TMPREJ;'...
        'clear EEGTMP tmpcom;' ];
    
   %remove and events not of  type 'X'
    k = 1;
    while(k <= size(EEG.event,2))
        if(EEG.event(k).type ~= 'X')
            EEG.event(k) = [];
        else
            k = k + 1;
        end
    end
    
    icacomp = 1;
    eegplotoptions = { 'events', EEG.event };
    if ~isempty(EEG.chanlocs) & icacomp
        eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
    end;
    eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities -- eegplot()', ...
        'limits', [EEG.xmin EEG.xmax]*1000 , 'position', scrsz, 'command', command, eegplotoptions{:});
    f_UpdateScale(60);
    f_UpdateTimeWin(7);
    
    uiwait
    if(isempty(v_Rej))
        disp("Empty");
        OUTEEG = EEG;
        OUTEEG.setname = [listFiles(i) '_Post_VI'];
        EEG = pop_saveset( OUTEEG,'filepath',[str_DirOutput listFiles{i}]);
        RemovedFreqs = [];
        save([str_TriggerOut tmp_str{1} '_TriggersInterrog.mat'], 'burstIPI','burstStartEnd','stimFreqMedian','RemovedFreqs');
    else
        [v_RejSorted, idx] = sort(v_Rej(:,2));
        trialrej = v_Rej(idx,2)/EEG.srate;
        trialrej = uint8([trialrej/(size(EEG.data,2)/EEG.srate)]);
        v_LogRej = zeros(1, EEG.trials);
        v_LogRej(trialrej) = 1;
        OUTEEG = pop_rejepoch( EEG, trialrej, 0);
        OUTEEG.setname = [listFiles(i)];
        EEG = pop_saveset( OUTEEG,'filepath',[str_DirOutput listFiles{i}]);
        
        RemovedFreqs{1} = trialrej;
        RemovedFreqs{2} = burstIPI(trialrej);
        RemovedFreqs{3} = burstStartEnd(trialrej,:);
        RemovedFreqs{4} = stimFreqMedian(trialrej);
        burstIPI(trialrej) = [];
        burstStartEnd(trialrej,:) = [];
        stimFreqMedian(trialrej) = [];
        save([str_TriggerOut tmp_str{1} '_TriggersInterrog.mat'], 'burstIPI','burstStartEnd','stimFreqMedian','RemovedFreqs');

    end

    
    
end
 