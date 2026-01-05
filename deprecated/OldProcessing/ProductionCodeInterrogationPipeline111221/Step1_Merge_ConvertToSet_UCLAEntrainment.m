%MasterRenaming files

clear;
close all;
clc;

subjid = 'SCC-011_20220324';

%set woring directories
str_DirInput = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step1_InputCNT/'];
str_DirOutput = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step2_MergedSetFiles/'];
ElectrodeLocation = '/Users/andrewwilson/Documents/MATLAB/eeglab2021.1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';
%move to input directory
cd (str_DirInput);

%set target sampling rate
%s_TargetSR = 2048;

%read in data for name conversions
%[numdata, textdata, rawdata] = xlsread('/Users/LBBPIT/Documents/EEG/TMSEEG/InterrogationProtocol/TMSResearchEEGDates.xlsx');
listCNT = dir('*.cnt');
listCNT = {listCNT.name};
numCNT = numel(listCNT);

%create lists of files for running
tmp = dir('*cnt');
ListFiles = {tmp.name};
ListFiles = transpose(ListFiles);
tmpSubjID = zeros(numel(ListFiles, 5));
AllFiles = {};
for i = 1:numel(ListFiles)
    tmpFile = strsplit(ListFiles{i,1}, '_');
    TMSID_init = strsplit(tmpFile{1}, '-');
    TMSID = str2num([ TMSID_init{1} TMSID_init{2}]);
    %[a,b] = ismember(tmpFile{1}, textdata(:,2));
    %ResearchID = numdata(b-1,1);
    %time = strsplit(tmpFile{4}, '.');
    
    %create Session variable
    date_init = strsplit(tmpFile{3},'-');
    
    %datesplit = strsplit(tmpFile{4}, '.');
    date = str2double([date_init{1} date_init{2} date_init{3}]);
    
    %create Series variable
%     series_init = strsplit(tmpFile{4},'Series');
%     series_init_1 = strsplit(series_init{2}, '.cnt');
%     
%     %datesplit = strsplit(tmpFile{4}, '.');
%     series = str2double(series_init_1(1));
       
%     %create time variable
%     if(strcmp(tmpFile{2},'PRE'))
%         time = 1;
%     elseif(strcmp(tmpFile{2},'During'))
%         time = 2;
%     elseif(strcmp(tmpFile{2},'POST'))
%         time = 3; 
%     else
%         time = '';
%     end

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
subjID = unique(tmpSubjID(:,[1, 4]), 'rows');
 
%counterr = 1;
for i = 1:numel(subjID(:,1))

    counter = 0;
    %indicesmatrix = zeros(numel(ListFiles(:,1)), 1);
    ListConcatFiles_temp = {};
    for j = 1:numel(ListFiles(:,1))
        if(tmpSubjID(j, 1) == subjID(i,1) && tmpSubjID(j, 4) == subjID(i,2))
            %disp(ListFiles{j,1});
            %indicesmatrix(j,1) = 1;
            ListConcatFiles_temp{numel(ListConcatFiles_temp)+1,1} = ListFiles{j,1};
%             if(tmpSubjID(j,5) == 1)
%                 ListConcatFiles_temp{1,1} = ListFiles{j,1};
%              elseif(tmpSubjID(j,5) == 2)
%                 ListConcatFiles_temp{2,1} = ListFiles{j,1};
%             elseif(tmpSubjID(j,5) == 3)
%                 ListConcatFiles_temp{3,1} = ListFiles{j,1};
%             end
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
    
    
    %disp(['Running Subject: ' num2str(subjID(i,1)) '. There are ' num2str(counter) ' files to concatenate.']);
    
    % load set files 
    ALLEEG = struct([]);
    s_HighSR = 0;
    nFiles = numel(ListConcatFiles);
    v_SzIndiv = zeros(1, nFiles);
    
     %Load Files
    for iFile = 1:nFiles
        disp(['Loading ' ListConcatFiles{iFile}]);
        try EEG = pop_loadeep_v4(ListConcatFiles{iFile}); % anteego plugin for eeglab must be downloaded
        catch
            EEG = pop_biosig(ListConcatFiles{iFile}, 'channels', [1:64]); % for edf files
        end
     
   %  if sampling rates are different across eeg files
%         if EEG.srate ~= s_TargetSR
%             disp(['sampling rate = ' num2str(EEG.srate) ' for file: ']);
%             disp(ListConcatFiles{iFile});
%             disp('resampling...');
%             EEG = pop_resample( EEG, s_TargetSR);
%             s_HighSR = 1;
%         end
        if(EEG.nbchan > 64)
            EEG = pop_select(EEG,'channel', 1:64)
            %disp('!!!!!!!!EEG HAS MORE THAN 64 CHANNELS!!!!!!!!!');
        end


        if strcmp(EEG.comments(end-2:end), 'edf') && length(EEG.chanlocs(1).labels)>4
%             [ EEG, str_ChanLab ] = f_CleanEDFLabels( EEG );
            [ EEG, str_ChanLab ] = f_CleanEdfChanLabels( EEG );           
        end
        str_ChanLab = {EEG.chanlocs(:).labels};
        [a, b] = ismember('EOG', str_ChanLab);
        if (~isempty(b) && b ~= 0)
            EEG.chanlocs(b).labels = 'CPz';
            % replace 'EOG' with CPz label (both not present), because eeglab chan editing cannot handle EOGs
        end
        
        %%%MAKE SURE CPZ IS ALL ZEROS
        if(numel(find(EEG.data(32,:) ~= 0)) > 0)
            disp("WARNING: CPZ NOT ALL ZEROS. ZEROING CPZ");
            EEG.data(32,:) = 0;
        end
        
        
        
        v_SzIndiv(iFile) = size(EEG.data,2);
        % edit channel locations (standard locations)
        EEG=pop_chanedit(EEG, 'lookup',ElectrodeLocation);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    end
    
    minsrate = min([ALLEEG.srate]);
    for iEEG = 1:size(ALLEEG,2)
        if ALLEEG(iEEG).srate ~= minsrate
            %disp(['sampling rate = ' num2str(EEG.srate) ' for file: ']);
            %disp(ListConcatFiles{iFile});
            %disp('resampling...');
            ALLEEG(iEEG) = pop_resample( ALLEEG(iEEG), minsrate);
            v_SzIndiv(iEEG) = size(ALLEEG(iEEG).data,2);
        end
    end
    
    
    %merge setfiles
    EEG = pop_mergeset( ALLEEG, 1:nFiles, 0);
    %check sampling rate
    disp(['total set length (min): ' num2str((size(EEG.data,2)/EEG.srate)/60)]);
    if size(EEG.data,2) ~= sum(v_SzIndiv)
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        error('merged data set ~= sum individual sets');
    end
    %create string variables for creating filenames
    str_splitname = strsplit(ListConcatFiles{1,1}, '.');
    str_splitname = strsplit(str_splitname{1}, '_');
    %counterr = counterr + counter;
   
    %create output filename and export
    OutputFilename = [str_splitname{1} '_' str_splitname{3} '_' str_splitname{4} '_Merged.set'];
    disp(OutputFilename);
    % Save set file
    EEG.setname = OutputFilename;
    OUTEEG = pop_saveset( EEG , 'filename',OutputFilename,'filepath',str_DirOutput);
   
end
disp('Done with all files');