% p_Entrainment_ASC_InputOutput_WholeBrain_4ACW
close all;
clc;
clear;

tic;
% load data
%% 1) detect all pulses
%% 2) cut data around TMS pulse (Baseline =  a) before interrogatiopn b) before train)
%% 3) cut artifact
saveFig = 1;

subjid = 'SCC-011_20220324';

% dirData = '/Volumes/storage/rTMS/Interrogations/Entrainment/';
dirData = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step7_VisualInspectionOutput/'];
matDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step7_FinalTriggerFiles/'];
saveDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step8_ASCOutput/'];

%dirData = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/VIOutput/'];
%atDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/FinalTriggers/'];
%saveDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/SCC/'];


% TOIlength = 400;
TOIlength = 1000;
HalfHz = 2;
% files = dir([dirData '*_EpochData']);
% files = {files.name}';
% patID = cellfun(@(x) x(1:5), files, 'UniformOutput',false); 
% nPat = numel(files);


cd(dirData)
files = dir('*_EpochData.set');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    %patID{iFile} = tmpsplit{1};
    patID{iFile} = [tmpsplit{1} '_' tmpsplit{2} '_' tmpsplit{3} '_' tmpsplit{4}];
end
nPat = numel(files);

% exportToPPTX('close'); % close the document
% exportToPPTX('new'); %initiate
allStimFreq = cell(1,nPat);
aSCPreAll = cell(1,nPat);
aSCPostAll = cell(1,nPat);
for iPat = 1:nPat
    disp(['Patient ' num2str(iPat)]);
    
    EEG = pop_loadset(fullfile(dirData, [files{iPat}]));
    if(EEG.nbchan ~= 64)
        %%Need to remove channels greater than 64
        EEG = pop_select(EEG,'channel', 1:64)
        
    end
    
    chanLabels = {EEG.chanlocs.labels}';
    idxSeed = find(ismember(chanLabels, 'F3'));
    srate = EEG.srate;
    if srate ~= 1000
        error('srate ~= 1000');
    end
    
    nFreq = size(EEG.data,3);
    nChan = size(EEG.data,1);
    v_TimeAxis = ((1:size(EEG.data,2))/srate)*1000-1000;
    
    idxM1M2 = ismember({EEG.chanlocs.labels}', {'M1', 'M2', 'F3'});
    chanToDo = 1:64;
    chanToDo(idxM1M2) = [];
    
    %% load stimulation frequencies
    load(fullfile(matDir, [patID{iPat} '_TriggersInterrog.mat']));
    [sorted, idxSorted] =  sort(stimFreqMedian, 'ascend');
    burstStartEndSorted = burstStartEnd(idxSorted,:);
    stimFreqMedianSorted = stimFreqMedian(idxSorted);
    
    if TOIlength == 400
        nFreqOut = 10;  % 10 for 400 ms, 19 for 1000 ms
    elseif TOIlength == 1000
        nFreqOut = 19;  % 10 for 400 ms, 19 for 1000 ms
        %nFreqOut = 1;
        %freqofInterest = 15;
    end
    
    
    aSCPre = nan(nFreqOut, length(stimFreqMedianSorted),nChan, nChan);
    aSCPost = nan(nFreqOut, length(stimFreqMedianSorted),nChan, nChan);
    for iChan1 = 1:nChan
        
        disp(['Channel ' num2str(iChan1)]);
        for iFreqStim = 1:length(stimFreqMedianSorted)
            
            TOIbsl = v_TimeAxis < 0 & v_TimeAxis >= -TOIlength; %-1000; 
            TOI = v_TimeAxis > 0 & v_TimeAxis <=  TOIlength; %1000;
            
            %% PRE
            parFFT.srate = srate;
            winlength = nnz(TOIbsl); % at least 4 seconds of data to have a reasonable resolution.
            fftlength = 2^nextpow2(winlength);
            s_FreqRes = srate/fftlength; % Should be around 0.5
            parFFT.fftlength = fftlength;
            parFFT.winlength = winlength;
            parFFT.overlap = 50;
            parFFT.freqres = s_FreqRes;
            parFFT.subjID = patID{iPat};
            parFFT.chanLabels = chanLabels;
            clear abs_power;
            
            [ abs_powerPre,  rel_powerPre] = f_PSD(EEG.data(:,(TOIbsl),idxSorted(iFreqStim)), parFFT);
            v_FreqAxis = rel_powerPre.freqRange;
            
            for iChan2 = iChan1+1:nChan
                for iFreqOut = 1:numel(v_FreqAxis)
                    clear aSC;
                    [aSC] = f_aSC_Hz(rel_powerPre.spectraREL, rel_powerPre.freqRange, ...
                        v_FreqAxis(iFreqOut), {chanLabels{iChan1}, chanLabels{iChan2}}, chanLabels, HalfHz);
                    aSCPre(iFreqOut, iFreqStim,iChan1, iChan2) = aSC;
                end
            end
            
           
%             [ abs_powerPre,  rel_powerPre] = f_PSD(EEG.data(:,(TOIbsl),idxSorted(iFreqStim)), parFFT);
%             v_FreqAxis = freqofInterest;
% 
%             for iChan2 = iChan1+1:nChan
%                 for iFreqOut = 1:numel(v_FreqAxis)
%                     clear aSC;
%                     [aSC] = f_aSC_Hz(rel_powerPre.spectraREL, rel_powerPre.freqRange, ...
%                         v_FreqAxis(iFreqOut), {chanLabels{iChan1}, chanLabels{iChan2}}, chanLabels, HalfHz);
%                     aSCPre(iFreqOut, iFreqStim,iChan1, iChan2) = aSC;
%                 end
%             end
%             
            %% POST
            parFFT.srate = srate;
            winlength = nnz(TOIbsl); % at least 4 seconds of data to have a reasonable resolution.
            fftlength = 2^nextpow2(winlength);
            s_FreqRes = srate/fftlength; % Should be around 0.5
            parFFT.fftlength = fftlength;
            parFFT.winlength = winlength;
            parFFT.overlap = 50;
            parFFT.freqres = s_FreqRes;
            parFFT.subjID = patID{iPat};
            parFFT.chanLabels = chanLabels;
            clear abs_power;
           
            [ abs_powerPost,  rel_powerPost] = f_PSD(EEG.data(:,(TOI),idxSorted(iFreqStim)), parFFT);
            v_FreqAxis = rel_powerPost.freqRange;
            for iChan2 = iChan1+1:nChan
                for iFreqOut = 1:numel(v_FreqAxis)
                    clear aSC;
                    [aSC] = f_aSC_Hz(rel_powerPost.spectraREL, rel_powerPost.freqRange, ...
                        v_FreqAxis(iFreqOut), {chanLabels{iChan1}, chanLabels{iChan2}}, chanLabels, HalfHz);
                    aSCPost(iFreqOut, iFreqStim, iChan1, iChan2) = aSC;
                end
            end

%             [ abs_powerPost,  rel_powerPost] = f_PSD(EEG.data(:,(TOI),idxSorted(iFreqStim)), parFFT);
%             v_FreqAxis = freqofInterest;
%             
%             for iChan2 = iChan1+1:nChan
%                 for iFreqOut = 1:numel(v_FreqAxis)
%                     clear aSC;
%                     [aSC] = f_aSC_Hz(rel_powerPost.spectraREL, rel_powerPost.freqRange, ...
%                         v_FreqAxis(iFreqOut), {chanLabels{iChan1}, chanLabels{iChan2}}, chanLabels, HalfHz);
%                     aSCPost(iFreqOut, iFreqStim, iChan1, iChan2) = aSC;
%                 end
%             end
        end
    end
    
    clear EEG;
    aSCPreAll = aSCPre;
    aSCPostAll = aSCPost;
    allStimFreq = stimFreqMedianSorted;
    clear inputOutput;
    
    save(fullfile(saveDir, ['Entrainment_AlphaSCWholeBrain_' patID{iPat} '_' num2str(TOIlength) 'ms']), 'aSCPreAll', 'aSCPostAll', 'allStimFreq',...
    'chanLabels', 'v_FreqAxis', 'v_TimeAxis', 'parFFT');
end


