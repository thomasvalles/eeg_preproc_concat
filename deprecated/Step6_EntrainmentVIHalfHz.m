% p_postFASTER_VI

clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;
close all;
clc;
eeglab;
% addpath('/Users/thomas/Documents/GitHub/RF_StepsMatlab/')
% subjid = 'RFUCLA-001_IP1';
% work_dir = '/Users/thomas/Documents/interrogation/';%'/Volumes/Files/Interrogation/';
% work_dir = '/Volumes/Files/Interrogation/5x5/SCC23-030/';

Step0_SetDirectories;
%set working directories
str_DirInput = [work_dir subjid filesep 'new' '/Step4_FasterOutput' suffix filesep] %converting the orignal input file data to a string array
str_DirOutput = [work_dir subjid  filesep 'new' '/Step5_VisualInspectionOutput' suffix filesep]; %converting the output file data to a string array and renaming it
str_TriggerDir = [work_dir subjid  filesep 'new' '/Step3_BrushedTriggers' suffix filesep];
str_TriggerOut = [work_dir subjid  filesep 'new' '/Step5_FinalTriggerFiles' suffix filesep];
mkdir(str_TriggerOut);
mkdir(str_DirOutput);
%move to input directory
cd(str_DirInput); %changing the current folder to a new folder

tmp = dir('*set');
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
    %load([str_TriggerDir tmp_str{1} '_TriggersInterrog.mat']);
    f = dir([str_TriggerDir '*_TriggersInterrog.mat']);
    load(fullfile(f.folder, f.name));

    %%tmp
    %stimFreqMedian_new = [8,8,8,11.5,11.5,11.5,4,4,4,9.5,9.5,9.5,13,13,13,7,7,7,10,10,10,15.5,15.5,15.5,14,14,14,4.5,4.5,4.5,17.5,17.5,17.5,16,16,16,11,11,11,12.4,12.4,12.4,17,17,17,16.5,16.5,16.5,9,9,9,NaN,8.5,8.5,14.5,14.5,14.5,13.5,13.5,13.5,7.5,7.5,7.5,12,12,12,5,5,5,15,15,15,10.5,10.5,18,18,6,6,6,6.5,6.5,6.5,5.5,5.5,5.5]
    %stimFreqMedian(1) = 15.0;

%     stimFreqMedian = round(stimFreqMedian*2)/2;

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
        trialrej = [];
        original_f_order = stimFreqMedian;
        save([str_TriggerOut tmp_str{1} '_TriggersInterrog.mat'], 'burstIPI','burstStartEnd','stimFreqMedian','RemovedFreqs','trialrej');
    
    else
        [v_RejSorted, idx] = sort(v_Rej(:,2));
        trialrej = v_Rej(idx,2)/EEG.srate;
        trialrej = uint8([trialrej/(size(EEG.data,2)/EEG.srate)]);
        v_LogRej = zeros(1, EEG.trials);
        v_LogRej(trialrej) = 1;
        
        %% note we aren't actually removing the epochs here, we're just marking them in the mat file 
        % this is different from earlier versions where the epochs were cut
        % out. I changed this to not take the epochs out because it made
        % things easier with subjs with repeated interrogations, but it's
        % important to make sure that you are aware of whether the epochs
        % have already been removed. Late 2022 on will have the
        % 'trialrej' variable which is the easiest way to tell that the
        % epochs have *not* yet been cut out yet.

        
        %OUTEEG = pop_rejepoch( EEG, trialrej, 0);
        OUTEEG=EEG;
        OUTEEG.setname = [listFiles(i)];
        EEG = pop_saveset( OUTEEG,'filepath',[str_DirOutput listFiles{i}]);
        
        RemovedFreqs{1} = trialrej;
        RemovedFreqs{2} = nan;%burstIPI(trialrej);
        RemovedFreqs{3} = burstStartEnd(trialrej,:);
        RemovedFreqs{4} = stimFreqMedian(trialrej);
        original_f_order = stimFreqMedian;
        %burstIPI(trialrej) = [];
        %burstStartEnd(trialrej,:) = [];
        %stimFreqMedian(trialrej) = [];
        save([str_TriggerOut tmp_str{1} '_TriggersInterrog.mat'], 'burstIPI','burstStartEnd','stimFreqMedian','RemovedFreqs','original_f_order','trialrej');
        Removed_matrix=[];

            
    end

    
    
end
 
