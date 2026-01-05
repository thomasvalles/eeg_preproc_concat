
% InterrogationTriggers
close all;
clc;
clear;



%% 1) generate figures with EEG data at electrode F3
%% 2) brush data manually
%% 3) plot & save cleaned/updated triggers
%% 4) inspect updated figures separately
subjid = 'SCC-011_20220324';

rawDataDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step4_ICLabelRun/'];
SaveDirFasterInput = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step5_FasterInput_PostManualICRemoval/'];
if ~exist(SaveDirFasterInput, 'dir')
    mkdir(SaveDirFasterInput);
end
files = dir('*.set');
files = {files.name}';
cd(rawDataDir)

patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '.');
    patID{iFile} = tmpsplit{1};
end
nPat = numel(files);
for iPat = 1:nPat
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));

%     flag = 1
%     comptoRej = [];
%     while(flag ~= 0)
%        flag = input('View[HeadPlots=1|Activation Scroll=2|Ready to reject=0]?');
%        if(flag == 1)
%            pop_topoplot(EEG,0);
% 
%        elseif(flag == 2)
%            pop_eegplot( EEG, 0, 1, 1);
%        elseif(f5lag == 0)
%            comptoRej = input('Which components would you like to reject. PUT BRACKETS AROUND VALUES?');
%        end
%     end    

    %pop_topoplot(EEG,0);
    EEG = eeg_checkset(EEG, 'ica');
    pop_eegplot( EEG, 0, 1, 1);
    set(gcf, 'Position', [0,0,1,1]);
    f_UpdateScale(20);
    
    comptoRej = input('Which components would you like to reject. PUT BRACKETS AROUND VALUES?');

    disp('Removing components')
    EEGout = pop_subcomp( EEG, comptoRej, 0);
    disp('Saving New EEG')
    EEG = pop_saveset( EEGout, 'filename', fullfile(SaveDirFasterInput, files{iPat}));
    close all
    clear EEG EEGout
end
