%MasterRenaming files

clear all
close all;
clc;
eeglab;

%set working directories
Step0_SetDirectories;
sourcedata = '/Volumes/Leuchter/BIDS_DATABASE/sourcedata';
rawdata = '/Volumes/Leuchter/BIDS_DATABASE/rawdata';

source_folders = dir(sourcedata);
source_folders = source_folders(~startsWith({source_folders.name}, '.'));

ifail = 1;
errs = struct([]);
% for each subject
for isub = 1:numel(source_folders)
    sub = source_folders(isub).name;
    
    % get the sessions
    session_folders = dir([source_folders(isub).folder filesep sub]);
    session_folders = session_folders(~startsWith({session_folders.name}, '.'));

    % for each session
    for ises = 1:numel(session_folders)
        try
            ses = session_folders(ises).name;
    
            disp("Starting " + sub + ", " + ses);
    
            % clear all eegs for each new session
            ALLEEG = struct([]);
    
            % formulate the rawdata directory
            rawdata_dir = [rawdata filesep sub filesep ses filesep 'eeg'];
            rawdata_files = dir(rawdata_dir);
            rawdata_files = rawdata_files(~startsWith({rawdata_files.name}, '.'));
    
            % if the rawdata directory is empty
            if isempty(rawdata_files)
                cnt_files = dir([session_folders(ises).folder filesep ses filesep 'eeg' filesep '*.cnt']);
    
                % for each cnt in the session
                for icnt = 1:numel(cnt_files)
                    EEG = pop_loadeep_v4(fullfile(cnt_files(icnt).folder, cnt_files(icnt).name)); % anteego plugin for eeglab must be downloaded
                    EEG = f_FixCPzAndEOG(EEG);
                    EEG = pop_chanedit(EEG, 'lookup',ElectrodeLocation);
                    [ALLEEG, ~, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
    
                % resample to min sampling rate if necessary
                minsrate = min([ALLEEG.srate]);
                min_channels = min([ALLEEG.nbchan]);
                for iEEG = 1:size(ALLEEG,2)
                    if ALLEEG(iEEG).srate ~= minsrate
                        ALLEEG(iEEG) = pop_resample(ALLEEG(iEEG), minsrate);
                    end
    
                    % maybe they added the eog between pre/stim?
                    if ALLEEG(iEEG).nbchan > min_channels
                        ALLEEG(iEEG) = pop_select(EEG, 'nochannel', {'VEOG'});
                        ALLEEG(iEEG) = eeg_checkset(ALLEEG(iEEG));
                    end
                end
    
                % merge files
                disp("Merging " + numel(cnt_files) + " files.");
                EEG = pop_mergeset(ALLEEG, 1:numel(cnt_files), 0);
    
                if ~isfolder(rawdata_dir)
                    mkdir(rawdata_dir);
                end
                % save
                out_name = strrep(cnt_files(1).name, '_run-01', '');
                pop_saveset(EEG, 'filename', out_name, 'filepath', rawdata_dir);
    
                disp("Finished " + sub + "_" + ses);
            end
        catch
            errs(ifail).sub = sub;
            errs(ifail).ses = ses;
            ifail = ifail + 1;
        end
    end
end

    