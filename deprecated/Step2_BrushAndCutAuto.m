% InterrogationTriggers
close all;
clc;
clear;
% change this script so it's -- if no events, then add the events using new
% scripts


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

% duration of the artifact to be cut
pre_sec = 0.005;
post_sec = 0.020;
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

%% generate fig files if they don't exist yet
if nPat >= 0
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

        idxF3 = find(strcmp({EEG.chanlocs.labels}', chan_trig));
        srate = EEG.srate;
        label = {EEG.chanlocs.labels}';

        F3 = EEG.data(idxF3,:);
 
        if(numel(EEG.event) > 0)
            idxValidTrigger = strcmpi({EEG.event.type}', '0001');
            if nnz(idxValidTrigger) == 0
                idxValidTrigger = strcmpi({EEG.event.type}', '01');
            end
            
            if nnz(idxValidTrigger) == 0
                idxValidTrigger = strcmpi({EEG.event.type}', '1');
            end
         
            pulseTriggers = [EEG.event.latency]';
            shiftedTriggers=pulseTriggers; % magventure      

            save(fullfile(saveDir, [patID{iPat}, '_TriggersShifted']), 'shiftedTriggers', 'srate');
        end
        %% important! the x-axis of the figure needs to be in seconds!
        figure('units', 'normalized', 'position', [0 0.5 1 0.5]);
        plot((1:numel(F3))/srate, F3);
        title(patID{iPat});
        savefig(gcf, [saveDir, patID{iPat}, '_' chan_trig]);
        clear EEG F3;
        close all;
    end
end


cd(rawDataDir);
for iPat = 1:nPat
    
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));

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
            plot((1:numel(F3))/srate, F3);
            hold on;
            stem(burstStartEnd/srate, ones(size(burstStartEnd))*-20000, 'r', 'marker', 'none', ...
                'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
            disp([num2str(size(burstStartEnd,1)) ' Pulse Trains Detected']);
            title([patID{iPat} '. ' num2str(numel(stimFreqMedian)) ' Pulse Trains Detected']);
        end

        if flag_needs_manual_detection % there weren't enough triggers or we couldn't find enough trains.
            option = 3;
            disp('Not enough trains detected. Using manual detection')

            % don't bother with brushing, just check the whole file
            [trigInterrog] = f_DetectTriggersManuallyAuto(F3,srate, [1 size(EEG.data, 2)], n_trains, pulses_per_train, start_cutoff, option);

        end

        if(flag_enough_trains == 0)
            attempts = attempts + 1;
            start_cutoff = start_cutoff * 3/4;
      
        else
            pptx = exportToPPTX();
            for i = 1:size(burstStartEnd,1)
                temp_indices = burstStartEnd(i,:);
                hFig = figure('units', 'normalized', 'position', [0 0 1 1], 'color', 'w');
                hold on;
                plot(F3(temp_indices(1) - srate:temp_indices(2) + srate))

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

    save(fullfile(saveDir, [patID{iPat} '_TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian')        
    savefig(gcf, [saveDir, patID{iPat}, '_TriggersInterrog' ]);
    close all;
    clear F3;
end


%% cut original file 
%% cut EEG data around each interrogation
for iPat = 1:nPat
    load(fullfile(saveDir, [patID{iPat} '_TriggersInterrog'])) % loads: 'burstStartEnd', 'burstIPI', 'stimFreqMedian');

    
    EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
   
    idxF3 = find(strcmp({EEG.chanlocs.labels}', 'F3'));
    srate = EEG.srate;
    label = {EEG.chanlocs.labels}';
    
    nChan = size(EEG.data,1);
    interpEpochData = zeros(nChan,3*srate,size(burstStartEnd,1));
    std_dev_count = 0;
    slope_count = 0;
    full_log_count = 0;
    failed_count = 0;
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
            gap_tolerance = -1;

            try % first try log detrend on quad detrend up to max std dev
                [YPost21, end_sample, stds] = f_DetrendLogisticStdMax(YPost2);
                gap = YPost21(end_sample+1) - YPost21(end_sample);

                if abs(gap) > gap_tolerance
                    YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;
                end
                std_dev_count = std_dev_count + 1;
            catch 
                try % try log detrend on quad detrend up to slope change (usually shorter than prev.)
                    [YPost21, end_sample] = f_DetrendLogisticSlopeChange(YPost2);
                    gap = YPost21(end_sample+1) - YPost21(end_sample);
    
                    if abs(gap) > gap_tolerance
                        YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;
                    end
                    
                    slope_count = slope_count + 1;
                catch
                    try % try log detrend on all of post
                        YPost21 = f_DetrendLogFull(YPost);
                        full_log_count = full_log_count + 1;
                    catch
                        [YPost21, end_sample, stds] = f_DetrendLinearStdMax(YPost2);
                        gap = YPost21(end_sample+1) - YPost21(end_sample);
    
                        if abs(gap) > gap_tolerance
                            YPost21(end_sample + 1:end) = YPost21(end_sample + 1:end) - gap;
                        end

                        failed_count = failed_count + 1;
                    end
                end
            end


            YPost21 = reshape(YPost21, 1, length(YPost21));

            nInterpol = pre_sec*srate+post_sec*srate;
            nPrePost = 2;
            [X, Y, XX, YY] = f_CutInterpolArtifact(YPre2, YPost21, nInterpol, nPrePost);
            
            interpolData = [YPre2, YY, YPost21];
            
            % deviation as SD
            pre_sample = srate-(srate*pre_sec);
            post_sample = srate+srate*post_sec;

            % add Gaussian noise to the part we interpolate with a line 
            if iChan ~= 32 % exclude CPz
                signal = interpolData(1,pre_sample:post_sample-1);
                [noisySig, noise] = f_addGN(signal, srate);
                if(srate == 2000)
                    interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
                else
                    interpolData = [interpolData(1,1:pre_sample) noisySig interpolData(1,post_sample-1:end)];
                end
                interpolData = [interpolData(1,1:pre_sample-1) noisySig interpolData(1,post_sample:end)];
            end
            
            interpEpochData(iChan,:,iFreq) = interpolData;
        end
    end
    
    counts = [std_dev_count slope_count full_log_count failed_count];
    save("counts.mat", "counts");
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
