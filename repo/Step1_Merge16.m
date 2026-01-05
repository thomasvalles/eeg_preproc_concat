%MasterRenaming files

clear all;
close all;
clc;
eeglab;

%set woring directories
Step0_SetDirectories;
str_DirInput = [work_dir subjid '/Step1_InputCNT/'];
str_DirOutput = [work_dir subjid '/Step2_MergedSetFiles/'];

mkdir(str_DirOutput);

%move to input directory
cd (str_DirInput);

listCNT = dir('*.cnt');
listCNT = {listCNT.name};
numCNT = numel(listCNT);

%create lists of files for running
tmp = dir('*.cnt');
ListFiles = {tmp.name};
ListFiles = transpose(ListFiles);
tmpSubjID = zeros(numel(ListFiles, 5));
AllFiles = {};

if contains(subjid, 'BH')

    for i = 1:numel(ListFiles)
        tmpFile = strsplit(ListFiles{i,1}, '_');
        TMSID_init = strsplit(tmpFile{1}, '-');
        TMSID = TMSID_init{2};
        TMSID = str2num(TMSID(3:end));
      
        date = 0;
        time = 0;

        %add to matrix
        tmpSubjID(i,1) = TMSID;
        tmpSubjID(i,2) = 1;
        tmpSubjID(i,3) = 0;
        tmpSubjID(i,4) = date;
        tmpSubjID(i,5) = time;
    end
else 

    for i = 1:numel(ListFiles)
        tmpFile = strsplit(ListFiles{i,1}, '_');
        TMSID_init = strsplit(tmpFile{1}, '-');
        TMSID = str2num([ TMSID_init{1} TMSID_init{2}]);
        
        %create Session variable
        date_init = strsplit(tmpFile{3},'-');
        
        %datesplit = strsplit(tmpFile{4}, '.');
        date = str2double([date_init{1} date_init{2} date_init{3}]);
    
        time = strsplit(tmpFile{4},'.');
        time_init = time{1};
        timesplit = strsplit(time_init, '-');
        time = str2double([timesplit{1} timesplit{2} timesplit{3}]);
        disp([num2str(date) ' ' num2str(time)]);
        
        %add to matrix
        tmpSubjID(i,1) = TMSID;
        tmpSubjID(i,2) = 1;
        tmpSubjID(i,3) = 0;
        tmpSubjID(i,4) = date;
        tmpSubjID(i,5) = time;
    end
end
subjID = unique(tmpSubjID(:,[1, 4]), 'rows');
 

for i = 1:numel(subjID(:,1))
    counter = 0;
    ListConcatFiles_temp = {};
    for j = 1:numel(ListFiles(:,1))
        if(tmpSubjID(j, 1) == subjID(i,1) && tmpSubjID(j, 4) == subjID(i,2))
            ListConcatFiles_temp{numel(ListConcatFiles_temp)+1,1} = ListFiles{j,1};
        end
    end
    
    
   disp('-----------------------------------------------------------------');
   ListConcatFiles = {};
   for j = 1:numel(ListConcatFiles_temp)
       if(~isempty(ListConcatFiles_temp{j,1}))
           ListConcatFiles{counter+1,1} = ListConcatFiles_temp{j,1};
           counter = counter + 1;
           disp(ListConcatFiles{counter,1})
       end
   end
   disp('-----------------------------------------------------------------');
    
    % load set files 
    ALLEEG = struct([]);
    s_HighSR = 0;
    nFiles = numel(ListConcatFiles);
    v_SzIndiv = zeros(1, nFiles);
    
    
     %Load Files
    for iFile = 1:nFiles
        disp(['Loading ' ListConcatFiles{iFile}]);
        sample1 = 16000 * 60 * 5; % ignore first 5 minutes
        sample2 = sample1 + 16000 * 60 * 20; % first eeg holds 20 minutes of interrogation (hopefully)

        % look in the first 20 minutes for the first train, then use that
        % as the first sample
        EEG1 = pop_loadeep_v4(ListConcatFiles{iFile}, 'sample1', sample1, 'sample2', sample2); % anteego plugin for eeglab must be downloaded
        trigInterrog = [EEG1.event.latency]';
        [burstStartEnd, burstIPI, stimFreqMedian] = f_CheckTriggers(trigInterrog, 30, EEG1.srate);

        sample1 = burstStartEnd(1, 1) - EEG1.srate * 10; % 10 seconds before first train
        sample2 = sample1 + EEG1.srate * 60 * 20; % first eeg holds first 20 minutes
        sample3 = sample2 + 1 + EEG1.srate * 60 * 20; % second eeg will hold the next 20 minutes

        clear EEG1
        EEG1 = pop_loadeep_v4(ListConcatFiles{iFile}, 'sample1', sample1, 'sample2', sample2); % anteego plugin for eeglab must be downloaded
        disp('first half okay');
        EEG2 = pop_loadeep_v4(ListConcatFiles{iFile}, 'sample1', sample2 + 1, 'sample2', sample3);
        disp('second half okay');
        
        EEG1=pop_chanedit(EEG1, 'lookup',ElectrodeLocation);
        EEG2=pop_chanedit(EEG2, 'lookup',ElectrodeLocation);
        %[ALLEEG, EEG1, CURRENTSET] = eeg_store( ALLEEG, EEG1, 0 );
    end
    
    % minsrate = min([ALLEEG.srate]);
    % for iEEG = 1:size(ALLEEG,2)
    %     if ALLEEG(iEEG).srate ~= minsrate
    %         ALLEEG(iEEG) = pop_resample( ALLEEG(iEEG), minsrate);
    %         v_SzIndiv(iEEG) = size(ALLEEG(iEEG).data,2);
    %     end
    % end
    % 
    
    %merge setfiles
    %EEG1 = pop_mergeset( ALLEEG, 1:nFiles, 0);
    %check sampling rate
    % disp(['total set length (min): ' num2str((size(EEG1.data,2)/EEG1.srate)/60)]);
    % if size(EEG1.data,2) ~= sum(v_SzIndiv)
    %     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    %     error('merged data set ~= sum individual sets');
    % end
    %create string variables for creating filenames
    str_splitname = strsplit(ListConcatFiles{1,1}, '.');
    str_splitname = strsplit(str_splitname{1}, '_');
    %counterr = counterr + counter;
   
    %create output filename and export
    OutputFilename = [str_splitname{1} '_' str_splitname{3} '_' str_splitname{4} '_Merged_a.set'];
    disp(OutputFilename);
    % Save set file
    EEG1.setname = OutputFilename;
    [~] = pop_saveset( EEG1 , 'filename',OutputFilename,'filepath',str_DirOutput);

    %create output filename and export
    OutputFilename = [str_splitname{1} '_' str_splitname{3} '_' str_splitname{4} '_Merged_b.set'];
    disp(OutputFilename);
    % Save set file
    EEG2.setname = OutputFilename;
    [~] = pop_saveset( EEG2 , 'filename',OutputFilename,'filepath',str_DirOutput);

end
disp('Done with all files');
