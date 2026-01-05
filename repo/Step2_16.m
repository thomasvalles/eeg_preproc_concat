% InterrogationTriggers
close all;
clc;
clear;
% change this script so it's -- if no events, then add the events using new
% scripts140



%% 1) generate figures with EEG data at electrode F3
%% 2) brush data manually
%% 3) plot & save cleaned/updated triggers
%% 4) inspect updated figures separately
eeglab;

Step0_SetDirectories;
chan_trig = 'F3';
rawDataDir = [work_dir subjid '/Step2_MergedSetFiles/'];
saveDir = [work_dir subjid '/Step3_BrushedTriggers/'];
SaveDirFasterInput = [work_dir subjid '/Step3_CutFiles/'];
mkdir(SaveDirFasterInput);

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Determine device used based on subject id. 
if contains(subjid, 'BH')
    device = "Magstim";
else
    device = "Magventure";
end


n_trains = 140;
pulses_per_train = 40;
cd(saveDir)
files = dir('*_F3.fig');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{1};
end
nPat = numel(files);

cd(rawDataDir);
files = dir('*.set');
files = {files.name}';

%% 

f1 = '/Volumes/My Passport/interrogation/BH0686_IP2/Step2_MergedSetFiles/sub-BH0686_task-INTER_eeg_Merged_a.set';
f2 = '/Volumes/My Passport/interrogation/BH0686_IP2/Step2_MergedSetFiles/sub-BH0686_task-INTER_eeg_Merged_b.set';
EEG1 = pop_loadset(f1); %(files{1});
%EEG1 = f_FixCPzAndEOG(EEG1);

EEG2 = pop_loadset(f2); %(files{2});
%EEG2 = f_FixCPzAndEOG(EEG2);

EEG = pop_mergeset(EEG1, EEG2);
clear EEG1 EEG2

pulseTriggers = [EEG.event.latency]';
shiftedTriggers=pulseTriggers; 

% downsample first
% all_data = zeros(iChan, 4800001);
% for iChan = 1:size(EEG.data, 1)
%     disp(iChan);
%     all_data(iChan, :) = downsample(EEG.data(iChan, :), 8);
% end
label = {EEG.chanlocs.labels}';
setname = fullfile(SaveDirFasterInput, [subjid, 'down_first.set']);
% [ EEGout ] = f_CreateEeglabStruct( all_data, 2000, label, setname );
% EEGout.trials = 1;
EEGout = pop_resample(EEG, 2000);
EEGout = pop_saveset( EEGout, 'filename', setname);


%% make cuts based on 16khz triggers
srate = EEG.srate;
minPulsePerBurst = 30;
[burstStartEnd, burstIPI, stimFreqMedian] = f_CheckTriggers(shiftedTriggers, minPulsePerBurst, EEG.srate);

post_secs = [0.020; 0.015; 0.010];
pre_sec = 0;
for iCut = 1:numel(post_secs)
    disp(iCut);
    post_sec = post_secs(iCut);
    srate_target = 2000;
    interpEpochData = zeros(EEG.nbchan, 2*srate_target - (post_sec*srate_target), size(burstStartEnd,1));
    for iFreq = 1:size(burstStartEnd, 1)
        burst_start = burstStartEnd(iFreq, 1);
        burst_end = burstStartEnd(iFreq, 2);
        % concatenate 2s before & 2s after pulse
        idxArt2Remove = [burst_start-srate*pre_sec burst_end+srate*post_sec]; % the pulse to be remove
        idxWholeEpoch = [burst_start-2*srate+1 burst_end+srate*2 - 1]; % whole new epoch period
        
        %interpEpochData(:, :, iFreq) = EEG.data(:, idxArt2Remove(2):idxWholeEpoch(2));
        tmp = EEG.data(:, idxArt2Remove(2):idxWholeEpoch(2));
        for iChan = 1:size(EEG.data, 1)
            interpEpochData(iChan, :, iFreq) = downsample(tmp(iChan,:), 8);
        end

    end
    label = {EEG.chanlocs.labels}';
    setname = fullfile(SaveDirFasterInput, [subjid, '_true_downsampled_cut-' num2str(1000*post_secs(iCut)) '.set']);
    [ EEGout ] = f_CreateEeglabStruct( interpEpochData, srate,label, setname );
    EEGout.trials = size(burstStartEnd, 1);
    %EEGout = pop_resample(EEGout, 2000);
    EEGout = pop_saveset( EEGout, 'filename', setname);
end

%% detrend on downsampled data
% now that we've already made the cut, can set burst_start = burst_end
burst_start = 4000;
burst_end = 4000;
pre_sec = 0.005;
post_sec = 0.015;
srate = EEG.srate;
interpEpochData = zeros(EEG.nbchan,3*srate,size(burstStartEnd,1));
parpool(4)
for iChan = 1:size(EEG.data, 1)
    disp(['Channel: ' num2str(iChan)])
    parfor iFreq = 1:size(burstStartEnd, 1)
        EEGdata = EEG.data(iChan, :, iFreq);
        % need to reshape
        [interpolData] = f_FixArtifact(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, device);

        % deviation as SD
        pre_sample = srate-(srate*pre_sec);
        post_sample = srate+srate*post_sec;

        % add Gaussian noise to the part we interpolate with a line 
        if iChan ~= 32 % exclude CPz
            signal = interpolData(1,pre_sample:post_sample-1);
            [noisySig, noise] = f_addGN(signal, srate);
            if(srate == 2000 || srate == 16000)
                interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
            else
                interpolData = [interpolData(1,1:pre_sample) noisySig interpolData(1,post_sample-1:end)];
            end
            interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
        end
        
        interpEpochData(iChan,:,iFreq) = interpolData;
    end
end
delete(gcp);
label = {EEG.chanlocs.labels}';
setname = fullfile(SaveDirFasterInput, [subjid, '_EpochData.set']);
[ EEGout ] = f_CreateEeglabStruct( interpEpochData, srate,label, setname );

EEGout=pop_chanedit(EEGout, 'lookup',ElectrodeLocation);
EEG = pop_saveset( EEGout, 'filename', setname);
    

%%
nPat = 1;
for iPat = 1:nPat
   
    idxF3 = find(strcmp({EEG.chanlocs.labels}', chan_trig));
    F3 = EEG.data(idxF3,:);
    
    trigInterrog = shiftedTriggers;
    srate = EEG.srate;
  
    flag_enough_trains = 0;
    start_cutoff = 1/2;
    attempts = 0;
    while (flag_enough_trains == 0) && (attempts < 5)
        flag_needs_manual_detection = 1;
        if numel(trigInterrog) > 1000 % if there are at least some triggers, see if we can find enough trains.
            minPulsePerBurst = 30;
            [burstStartEnd, burstIPI, stimFreqMedian] = f_CheckTriggers(trigInterrog, minPulsePerBurst, srate);

            if numel(stimFreqMedian) >= n_trains % if we have enough trains
                flag_needs_manual_detection = 0; % don't need manual detection
                burstStartEnd = burstStartEnd((end + 1)- n_trains:end, :); % leave out any extras from beginning
                burstIPI = burstIPI((end + 1) - n_trains:end);
                stimFreqMedian = stimFreqMedian((end + 1) - n_trains:end);
                flag_enough_trains = 1;

            end

            figure; % plot the trains
            %plot((1:numel(F3))/srate, F3);
            hold on;
            % stem(burstStartEnd/srate, ones(size(burstStartEnd))*-20000, 'r', 'marker', 'none', ...
            %     'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
            disp([num2str(size(burstStartEnd,1)) ' Pulse Trains Detected']);
            title([subjid '. ' num2str(numel(stimFreqMedian)) ' Pulse Trains Detected']);
        end

        if flag_needs_manual_detection % there weren't enough triggers or we couldn't find enough trains.
            option = 3;
            disp('Not enough trains detected. Using manual detection.')

            % don't bother with brushing, just check the whole file
            [trigInterrog] = f_DetectTriggersManuallyAuto(F3,srate, [1 size(EEG.data, 2)], n_trains, pulses_per_train, start_cutoff, option);

        end

        if(flag_enough_trains == 0)
            attempts = attempts + 1;
            start_cutoff = start_cutoff * 3/4;
      
        else
            disp('');
            % pptx = exportToPPTX();
            % for i = 1:size(burstStartEnd,1)
            %     temp_indices = burstStartEnd(i,:);
            %     hFig = figure('units', 'normalized', 'position', [0 0 1 1], 'color', 'w');
            %     hold on;
            %     plot(F3(temp_indices(1) - srate:temp_indices(2) + srate))
            % 
            %     temp_indices = temp_indices - temp_indices(1)+srate;
            %     stem((temp_indices), ones(size(temp_indices))*-20000, 'r', 'marker', 'none', ...
            %         'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
            %     title(['Burst: ' num2str(i) ', ' num2str(stimFreqMedian(1,i)) ' Hz'])
            %     pptx.addSlide();
            %     pptx.addPicture(gcf,'Scale', 'maxfixed');
            %     close gcf;
            % end
            % pptx.save([saveDir patID{iPat} '_InterrogationMarkers']); % save to file
        end
    end

    save(fullfile(saveDir, [subjid '_TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian')        
    savefig(gcf, [saveDir, subjid, '_TriggersInterrog' ]);
    close all;
    clear F3;
end


%% cut original file 
%% cut EEG data around each interrogation
for iPat = 1:nPat
    load(fullfile(saveDir, [patID{iPat} '_TriggersInterrog'])) % loads: 'burstStartEnd', 'burstIPI', 'stimFreqMedian');

    
    % EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
   
    idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
    srate = EEG.srate;
    label = {EEG.chanlocs.labels}';
    
    nChan = size(EEG.data,1);
    interpEpochData = zeros(nChan,3*srate,size(burstStartEnd,1));
    for iChan = 1:nChan
        
        disp(['Channel ' num2str(iChan)]);
        EEGdata = EEG.data(iChan,:);
        idxF3 = strcmp({EEG.chanlocs.labels}', 'F3');
        
        parfor iFreq = 1:size(burstStartEnd,1)
            
            burst_start = burstStartEnd(iFreq, 1);
            burst_end = burstStartEnd(iFreq, 2);

            % duration of the artifact to be cut
            pre_sec = 0.005;
            post_sec = 0.010;

            [interpolData, post_sec] = f_FixArtifact(EEGdata, pre_sec, post_sec, burst_start, burst_end, srate, device);
  
            % deviation as SD
            pre_sample = srate-(srate*pre_sec);
            post_sample = srate+srate*post_sec;

            % add Gaussian noise to the part we interpolate with a line 
            if iChan ~= 32 % exclude CPz
                signal = interpolData(1,pre_sample:post_sample-1);
                [noisySig, noise] = f_addGN(signal, srate);
                if(srate == 2000 || srate == 16000)
                    interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
                else
                    interpolData = [interpolData(1,1:pre_sample) noisySig interpolData(1,post_sample-1:end)];
                end
                interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
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
function [noisySig, noise] = f_addGN(signal, srate)
% original version from FileExchange "addGaussianNoise.m"

reqSNR = 2;

a = 50;
b = 70;
f = round((b-a).*rand(1,1) + a); % random frequency selection between 50 & 70 hz

% generate sine wave at 50/60/70 Hz
lengthSig = length(signal);
t = [1:lengthSig]/srate;
signalGamma = sin(2*pi*f*t);

%% add white noise
sigEner = norm(signalGamma(:))^2;                    % energy of the signal
noiseEner = sigEner/(10^(reqSNR/10));        % energy of noise to be added
noiseVar = noiseEner/(length(signalGamma(:))-1);     % variance of noise to be added
noiseStd = sqrt(noiseVar);                   % std. deviation of noise to be added
noise = noiseStd*randn(size(signalGamma));           % noise
noisySig = signalGamma+noise+signal;                        % noisy signal
end
