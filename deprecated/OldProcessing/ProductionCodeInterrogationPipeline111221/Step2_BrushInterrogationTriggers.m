% InterrogationTriggers
close all;
clc;
clear;



%% 1) generate figures with EEG data at electrode F3
%% 2) brush data manually
%% 3) plot & save cleaned/updated triggers
%% 4) inspect updated figures separately
subjid = 'SCC-011_20220324';

rawDataDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step2_MergedSetFiles/'];
saveDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step3_BrushedTriggers/'];
SaveDirFasterInput = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/SCCPipeline/' subjid '/Step3_CutFiles/'];
ElectrodeLocation = '/Users/andrewwilson/Documents/MATLAB/eeglab2021.1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% duration of the artifact to be cut
pre_sec = 0.005;
post_sec = 0.020;
msArtifact = 0.05;
msBSL = 0.5;
cd(saveDir)
files = dir('*_F3.fig');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{1};
end
nPat = numel(files);

%% generate fig files if they don't exist yet
if nPat == 0
    cd(rawDataDir);
    files = dir('*.cnt');
    files = {files.name}';
    
    if isempty(files)
        files = dir('*.set');
        files = {files.name}';
    end
    patID = [];
    for iFile = 1:numel(files)
        tmpsplit = strsplit(files{iFile}, '.');
        patID{iFile} = tmpsplit{1};
    end
    nPat = numel(files);
    
    for iPat = 1:nPat
        
        try
            EEG = pop_loadeep_v4(fullfile(rawDataDir, files{iPat})); % anteego plugin for eeglab must be downloaded
        catch
            EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
        end
        idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
        srate = EEG.srate;
        label = {EEG.chanlocs.labels}';
        %         if srate ~= 2000
        %             error('sampling rate not 2000');
        %         end
        %if(EEG.nbchan ~= 64)
        %   error('EEG has more than 64 channels')
        %end
        F3 = EEG.data(idxF3,:);
        idxValidTrigger = strcmpi({EEG.event.type}', '0001');
        if nnz(idxValidTrigger) == 0
            idxValidTrigger = strcmpi({EEG.event.type}', '01');
        end
        
        if nnz(idxValidTrigger) == 0
            idxValidTrigger = strcmpi({EEG.event.type}', '1');
        end
        
        %if nnz(idxValidTrigger) == 0
        %   error('no triggers');
        %end
        
        pulseTriggers = [EEG.event.latency]';
        shiftedTriggers = pulseTriggers - (17 * (EEG.srate/1000));
        save(fullfile(saveDir, [patID{iPat}, '_TriggersShifted']), 'shiftedTriggers', 'srate');
        
        %% important! the x-axis of the figure needs to be in seconds!
        figure('units', 'normalized', 'position', [0 0.5 1 0.5]);
        plot((1:numel(F3))/srate, F3);
        title(patID{iPat});
        savefig(gcf, [saveDir, patID{iPat}, '_F3']);
        clear EEG F3;
        close all;
    end
end


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
    tmpsplit = strsplit(files{iFile}, '.');
    patID{iFile} = tmpsplit{1};
end
nPat = numel(files);
for iPat = 1:nPat
    
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
    idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
    F3 = EEG.data(idxF3,:);
    
    fprintf('patient %d/%d \n', iPat, nPat);
    uiopen([saveDir patID{iPat} '_F3.fig'],1);
    load(fullfile(saveDir, [patID{iPat} '_TriggersShifted']));
    
    %----------------------------------------------------------------------
    % plot & brush data to create new segments for each stim freq
    %----------------------------------------------------------------------
    
    flag = 0;
    brushX = [];
    while flag == 0
        
        brush on
        
        disp('brush the interrogation, skip 10 Hz test, TBS & Tx');
        flag = input('done brushing? [0/1]');
        figure(gcf);
        
        
        [ brushX] = f_GetBrush;
        
        if flag == 0 && isempty(brushX)
            
            disp('Please redo brushing, skip 10 Hz test, TBS & Tx');
            flag = input('done brushing? [0/1]');
            figure(gcf); brush on;
            [ brushX] = f_GetBrush;
        end
        
    end
    %----------------------------------------------------------------------
    % remove triggers outside the selection
    %----------------------------------------------------------------------
    
    idxTrig2Remove = shiftedTriggers/srate < brushX(1) | shiftedTriggers/srate > brushX(end);
    trigInterrog = shiftedTriggers;
    trigInterrog(idxTrig2Remove) = []; % removes triggers outside of selected area
    flag_redo = 0;
    while(flag_redo == 0)
        
        if (numel(trigInterrog) < 1000)
            disp('WARNING: < 1000 triggers');
            option = input('Which option would you like to use for trigger detection[1|2]?');
            [trigInterrog] = f_DetectTriggersManually(F3,srate, [brushX(1) brushX(end)], option);
            figure;
            plot((1:numel(F3))/srate, F3);
            hold on;
            stem(trigInterrog/srate, ones(1,numel(trigInterrog))*1000, 'r')
            title(patID{iPat});
        end
        %----------------------------------------------------------------------
        % Remove all random and non group-of-40 pulses
        %----------------------------------------------------------------------
        
        minPulsePerBurst = 38;
        [burstStartEnd,burstIPI, stimFreqMedian] = f_CheckTriggers(trigInterrog, minPulsePerBurst, srate);
        save(fullfile(saveDir, [patID{iPat} '_TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian')
        
        
        %% save mat fig for visual inspection
        close all;
        uiopen([saveDir patID{iPat} '_F3.fig'],1);
        hold on;
        figure(gcf)
        stem(burstStartEnd/srate, ones(size(burstStartEnd))*-20000, 'r', 'marker', 'none', ...
            'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
        disp([num2str(size(burstStartEnd,1)) ' Pulse Trains Detected'])
        
        flag_redo = input('Done with this file(1) or Redo threshold(0)? [0/1]');
        if(flag_redo == 0)
            trigInterrog = [];
        else
            pptx = exportToPPTX();
            for i = 1:size(burstStartEnd,1)
                temp_indices = burstStartEnd(i,:);
                hFig = figure('units', 'normalized', 'position', [0 0 1 1], 'color', 'w');
                hold on;
                plot(F3(temp_indices(1) - srate:temp_indices(2) + srate))
                %plot(F3)
                
                
                temp_indices = temp_indices - temp_indices(1)+srate;
                stem((temp_indices), ones(size(temp_indices))*-20000, 'r', 'marker', 'none', ...
                    'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
                title(['Burst: ' num2str(i) ', ' num2str(stimFreqMedian(1,i)) ' Hz'])
                pptx.addSlide();
                pptx.addPicture(gcf,'Scale', 'maxfixed');
                close gcf;
            end
            pptx.save([saveDir patID{iPat} '_InterrogationMarkers']); % save to file
        end
    end
    savefig(gcf, [saveDir, patID{iPat}, '_TriggersInterrog' ]);
    
    close all;
    clear F3;
end


%% cut EEG data around each interrogation
for iPat = 1:nPat
    load(fullfile(saveDir, [patID{iPat} '_TriggersInterrog'])) % loads: 'burstStartEnd', 'burstIPI', 'stimFreqMedian');
    
    %% sort stimulation frequencies
    %[sorted, idxSorted] =  sort(stimFreqMedian, 'ascend');
    %burstStartEndSorted = burstStartEnd(idxSorted,:);
    %stimFreqMedianSorted = stimFreqMedian(idxSorted);
    
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
    idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
    srate = EEG.srate;
    label = {EEG.chanlocs.labels}';
    %     if srate ~= 2000
    %         error('sampling rate not 2000');
    %     end
    
    
    nChan = size(EEG.data,1);
    interpEpochData = zeros(nChan,3*srate,size(burstStartEnd,1));
    for iChan = 1:nChan
        
        disp(['Channel ' num2str(iChan)]);
        % 1) demean data for interpolation
        EEGdata = EEG.data(iChan,:);
        idxF3 = strcmp({EEG.chanlocs.labels}', 'F3');
        
        for iFreq = 1:size(burstStartEnd,1)
            
            if(srate == 2000)
                % concatenate 1s before & 2s after pulse
                idxArt2Remove = [burstStartEnd(iFreq,1)-srate*pre_sec burstStartEnd(iFreq,2)+srate*post_sec]; % the pulse to be remove
                idxWholeEpoch = [burstStartEnd(iFreq,1)-srate+1 burstStartEnd(iFreq,2)+srate*2-1]; % whole new epoch period
            elseif(srate == 2048)
                % concatenate 1s before & 2s after pulse
                idxArt2Remove = [burstStartEnd(iFreq,1)-srate*pre_sec burstStartEnd(iFreq,2)+srate*post_sec]; % the pulse to be remove
                idxWholeEpoch = [burstStartEnd(iFreq,1)-srate burstStartEnd(iFreq,2)+srate*2-1]; % whole new epoch period
            end
            YPre = EEGdata(idxWholeEpoch(1):idxArt2Remove(1)) - mean(EEGdata(idxWholeEpoch(1):idxArt2Remove(1)));
            YPost = EEGdata(idxArt2Remove(2):idxWholeEpoch(2)) - mean(EEGdata(idxArt2Remove(2):idxWholeEpoch(2)));
            
            % 2) detrend for interpolation (a) non-linear (b) linear between 0 & 100 samples
            YPre2 = detrendnonlin(YPre,2)';
            YPost2 = detrendnonlin(YPost,2)';
            detrend100 = detrend(YPost2(1:100));
            YPost21 = [detrend100 YPost2(101:end)];
            %
            %                         % deviation as SD
            %
            %             if std(YPost21(1,1:srate*msArtifact)) > 1.5*std(YPre2(1,1:srate*msBSL))
            %                 signal = YPost21(1,1:srate*msArtifact);
            %                 [noisySig, noise] = f_addGN(signal);
            %                 YPost21 = [noisySig YPost21(1,srate*msArtifact:end)];
            %             end
            
            % interpolation
            nInterpol = pre_sec*srate+post_sec*srate;
            nPrePost = 2;
            [X, Y, XX, YY] = f_CutInterpolArtifact(YPre2, YPost21, nInterpol, nPrePost);
            
            interpolData = [YPre2, YY, YPost21];
            
            % deviation as SD
            
            if std(interpolData(1,srate:srate+srate*msArtifact-1)) > 1.5*std(interpolData(1,1:srate*msBSL))
                signal = interpolData(1,srate:srate+srate*msArtifact-1);
                [noisySig, noise] = f_addGN(signal);
                if(srate == 2000)
                    interpolData = [interpolData(1,1:srate) noisySig interpolData(1,srate+srate*msArtifact+1:end)];
                else
                    interpolData = [interpolData(1,1:srate) noisySig interpolData(1,srate+srate*msArtifact:end)];
                end
            end
            
            interpEpochData(iChan,:,iFreq) = interpolData;
        end
    end
    
    clear EEG;
    root = strsplit(files{iPat}, '.');
    setname = fullfile(SaveDirFasterInput, [root{1}, '_EpochData.set']);
    [ EEGout ] = f_CreateEeglabStruct( interpEpochData, srate,label, setname );
    EEGout=pop_chanedit(EEGout, 'lookup',ElectrodeLocation);
    EEG = pop_saveset( EEGout, 'filename', setname);
    
    
end


%--------------------------------------------------------------------------
function [noisySig, noise] = f_addGN(signal)
% original version from FileExchange "addGaussianNoise.m"

reqSNR = 2;

a = 50;
b = 70;
f = round((b-a).*rand(1,1) + a); % random frequency selection between 50 & 70 hz

% generate sine wave at 50/60/70 Hz
lengthSig = length(signal);
srate = 1000;
t = [1:lengthSig]/srate;
signalGamma = sin(2*pi*f*t);
% figure, plot(signalGamma);

%% add white noise
sigEner = norm(signalGamma(:))^2;                    % energy of the signal
noiseEner = sigEner/(10^(reqSNR/10));        % energy of noise to be added
noiseVar = noiseEner/(length(signalGamma(:))-1);     % variance of noise to be added
noiseStd = sqrt(noiseVar);                   % std. deviation of noise to be added
noise = noiseStd*randn(size(signalGamma));           % noise
noisySig = signalGamma+noise;                        % noisy signal
end
