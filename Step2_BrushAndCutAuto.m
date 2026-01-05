% InterrogationTriggers

clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;
close all;
clc;
%clear;
% change this script so it's -- if no events, then add the events using new
% scripts140


parpool(4);



%% 1) generate figures with EEG data at electrode F3
%% 2) brush data manually
%% 3) plot & save cleaned/updated triggers
%% 4) inspect updated figures separately
%eeglab;


eeglab
Step0_SetDirectories;
chan_trig = 'F3';
rawDataDir = [work_dir subjid filesep 'Step2_MergedSetFiles' suffix filesep];
saveDir = [work_dir subjid filesep 'Step3_BrushedTriggers' suffix filesep];
SaveDirFasterInput = [work_dir subjid filesep 'Step3_CutFiles' suffix filesep];
mkdir(SaveDirFasterInput);

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Determine device used based on subject id. 
if contains(subjid, 'LA')
    device = "Magventure";
else
    device = "Magstim";
end

% with the 5-20hz interrogation, there are 16 frequencies * 10 trains
%n_trains = 160;
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
            shiftedTriggers=pulseTriggers; 
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

    % these two vars are for when we don't have triggers and have to use
    % peak detection to id pulses. start cutoff is proportion of the max
    % peak to use as a threshold. If this doesn't work, we'll cut it back
    % until we hopefully find the right number of trains. We will try this
    % process at most 3 times. If we still can't find all the trains, give
    % up and use what we have found. On all of the files I've done, we
    % either have triggers or can find the pulses "manually" without
    % cutting back the threshold
    start_cutoff = 1/2;
    attempts = 0;
    while (flag_enough_trains == 0) && (attempts <= 3)
        flag_needs_manual_detection = 1;
        if numel(trigInterrog) > 1000 % if there are at least some triggers, see if we can find enough trains.
            minPulsePerBurst = 30;
            [burstStartEnd, burstIPI, stimFreqMedian] = f_CheckTriggers(trigInterrog, minPulsePerBurst, srate);

            train_options = [64; 75; 160];%[39; 75; 80; 108; 140; 160];
            n_found = numel(stimFreqMedian);
            n_train_tol = 0;
            f_dist_2_ntrains = @(x) (abs(x - n_found));
            [M, ind] = min(f_dist_2_ntrains(train_options));
            target_trains = train_options(ind);
            % we have 108, 140, and 160 train versions of the interrogation.
            % stop when we've gotten close to any of these numbers.
            if M <= n_train_tol | n_found > 160 % if we are close to one of the fixed numbers
                flag_needs_manual_detection = 0; % don't need manual detection
                burstStartEnd = burstStartEnd((end + 1)- min(target_trains, n_found):end, :); % leave out any extras from beginning
                burstIPI = burstIPI((end + 1) - min(target_trains, n_found):end);
                stimFreqMedian = stimFreqMedian((end + 1) - min(target_trains, n_found):end);
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
            disp('Not enough trains detected. Using manual detection.')
            n_trains = 160;
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

    save(fullfile(saveDir, [patID{iPat} '_TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian')        
    savefig(gcf, [saveDir, patID{iPat}, '_TriggersInterrog' ]);
    close all;
    clear F3;

    if ~ismember(size(stimFreqMedian, 2), train_options)
        disp("WARNING: Irregular number of trains found. Found " + num2str(size(stimFreqMedian, 2)) + " trains")
        proceed = input("Do you wish to proceed (0/1)?");
        if proceed == 0
            error("Irregular number of trains");
        end
    end
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
    for iChan = 1:nChan
        
        disp(['Channel ' num2str(iChan)]);
        EEGdata = EEG.data(iChan,:);
        idxF3 = strcmp({EEG.chanlocs.labels}', 'F3');
        
        parfor iFreq = 1:size(burstStartEnd,1)
            rng(iFreq, "twister");
            burst_start = burstStartEnd(iFreq, 1);
            burst_end = burstStartEnd(iFreq, 2);

            % duration of the artifact to be cut
            pre_sec = 0.005;
            post_sec = 0.020;

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
    
    clear EEG;
    root = strsplit(files{iPat}, '.');
    setname = fullfile(SaveDirFasterInput, [root{1}, '_EpochData.set']);
    [ EEGout ] = f_CreateEeglabStruct( interpEpochData, srate,label, setname );
    EEGout=pop_chanedit(EEGout, 'lookup',ElectrodeLocation);
    EEG = pop_saveset( EEGout, 'filename', setname);
    
    
end
delete(gcp);
%--------------------------------------------------------------------------
