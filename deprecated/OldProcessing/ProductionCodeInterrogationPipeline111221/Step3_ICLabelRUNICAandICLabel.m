% InterrogationTriggers
close all;
clc;
clear;



%% 1) generate figures with EEG data at electrode F3
%% 2) brush data manually
%% 3) plot & save cleaned/updated triggers
%% 4) inspect updated figures separately
subjid = 'SCC-011_20220324';

rawDataDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step3_CutFiles/'];
SaveDirFasterInput = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step4_ICLabelRun/'];

%% brush interrogations
cd(rawDataDir);
files = dir('*.cnt');
files = {files.name}';

if isempty(files)
    files = dir('*.set');
    files = {files.name}';
end
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = [tmpsplit{1} '_' tmpsplit{2} '_' tmpsplit{3} '_' tmpsplit{4}];
end
nPat = numel(files);


%% cut EEG data around each interrogation
for iPat =1:nPat

    
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
    %idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
   % idxCPz = find(strcmp({EEG.chanlocs.labels}', 'CPz'));
    EEG_ica = pop_runica(EEG, 'icatype','runica', 'dataset',1,'options',{'extended' 1}); %Run ICA, manual component rejection is required after this
    
    
    %nICAStart = size(EEG_ica.icaweights,1);
    
    %EEG = pop_iclabel(EEG, 'Default');
    
    %EEG = pop_icflag(EEG, [0 0.9;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
   
    
   % EEG_ica = iclabel(EEG_ica);
    %[~, idxMaxLabel] = max(EEG_ica.etc.ic_classification.ICLabel.classifications,[], 2);
    %idx_reject = find(idxMaxLabel ~= 1);
    
    
    %OUTEEG = pop_subcomp( EEG_ica, idx_reject, 0);
    idx_reject =  find(EEG.reject.gcompreject == 1)
    nICAEnd = size(EEG_ica.icaweights,1);
    %disp(['Starting Number of ICA Weights = ' num2str(nICAStart)]);
    %disp(['Ending Number of ICA Weights = ' num2str(nICAEnd)]);

    
    EEG = pop_saveset( EEG_ica, 'filename', [SaveDirFasterInput files{iPat}]);
    clear EEG OUTEEG EEG_ica interpEpochData
    
    
end