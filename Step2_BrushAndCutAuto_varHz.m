% InterrogationTriggers

clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;
close all;
clc;
%clear;
% change this script so it's -- if no events, then add the events using new
% scripts140


% parpool(4);



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
nPat = 1;
iPat = 1;
%% generate fig files if they don't exist yet
if nPat >= 0
    % nagivate to raw data directory, and try to find the cnt/set file
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
        
    % open it
    try
        EEG = pop_loadeep_v4(fullfile(rawDataDir, files{iPat})); % anteego plugin for eeglab must be downloaded
    catch
        EEG = pop_loadset(fullfile(rawDataDir, files{iPat}));
    end

    % get the signal at F3
    idxF3 = find(strcmp({EEG.chanlocs.labels}', chan_trig));
    srate = EEG.srate;
    label = {EEG.chanlocs.labels}';

    F3 = EEG.data(idxF3,:);

    % if there are any event markers find them
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
    else
        shiftedTriggers = [];
    end

    % some additional quality control. We have noticed in the past some
    % single sample blips in these files that we would want to correct

    % highpass filter to smooth out signal didn't work (peak shifted indices by 1)
    % zero-phase filter made the peak less prominent. Trying a moving
    % median approach to detrend    
    DataCorrected = F3 - mean(F3);
    DataCorrected = detrendnonlin(DataCorrected, 2);

    thresh = 100 * std(DataCorrected);
    blip_inds = find((DataCorrected > thresh) | (DataCorrected < -thresh));

    % if there are blips, replace with linear interpolation from previous 5 points
    % to next 5 points
    if numel(blip_inds) > 0
        for ichan = 1:size(EEG.data, 1)
            for bad_ind = blip_inds
                if bad_ind > 0 && bad_ind < numel(F3)
                    EEG.data(ichan, (bad_ind-5:bad_ind+5)) = linspace(EEG.data(ichan, bad_ind - 5), EEG.data(ichan, bad_ind + 5), 11);
                end
            end
        end
    end

    F3 = EEG.data(idxF3,:);

    % save a png of the F3 signal.
    figure('units', 'normalized', 'position', [0 0.5 1 0.5]);
    plot((1:numel(DataCorrected))/srate, DataCorrected);
    hold on
    yline(thresh, 'color', 'r');
    title(patID{iPat});
    saveas(gcf, [saveDir, patID{iPat}, '_' chan_trig, '.png']);
    close all

    %% FIND PULSE TRAINS
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
    while (flag_enough_trains == 0) && (attempts <= 5)
        flag_needs_manual_detection = 1;
        if numel(trigInterrog) > 10 % if there are at least some triggers, see if we can find enough trains.
            minPulsePerBurst = 40 - 2;
            [burstStartEnd, burstIPI, stimFreqMedian] = f_CheckTriggers(trigInterrog, minPulsePerBurst, srate);

            % these are some common "expected" numbers of trains
            % We have them listed here because so we know when to stop
            % looking for trains. If we get close to one of these numbers,
            % but are slightly over (for tolerance trains) we can remove
            % them
            train_options = [71; 75; 80; 108; 140; 160];
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

           
        end

        if flag_needs_manual_detection % there weren't enough triggers or we couldn't find enough trains.
            option = 3;
            disp('Not enough trains detected. Using manual detection.')
            
            % this input in the function contributes to the number of peaks
            % used to estimate the threshold. here, were looking at the top
            % (30 * 40) / 2 peaks (/ 2 is used in the function itself).
            % It's okay if there are fewer than 60 trains.
            n_trains = 30;
            % don't bother with brushing, just check the whole file
            [trigInterrog] = f_DetectTriggersManuallyAuto_varHz(F3,srate, [1 size(EEG.data, 2)], n_trains, pulses_per_train, start_cutoff, option);

        end

        if(flag_enough_trains == 0)
            attempts = attempts + 1;
            start_cutoff = start_cutoff * 3/4;
      
        end
    end

    figure; % plot the trains
    plot((1:numel(F3))/srate, F3);
    hold on;
    stem(burstStartEnd/srate, ones(size(burstStartEnd))*-20000, 'r', 'marker', 'none', ...
        'linewidth', 1, 'basevalue', 20000, 'showbaseline', 'off');
    disp([num2str(size(burstStartEnd,1)) ' Pulse Trains Detected']);
    title([patID{iPat} '. ' num2str(numel(stimFreqMedian)) ' Pulse Trains Detected']);
    saveas(gcf, [saveDir, patID{iPat}, '_MarkedTrains', '.png']);

    save(fullfile(saveDir, [patID{iPat} '_TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian', 'blip_inds')        
    close all;

    % don't ask if the number of trains is okay. we expect a lot of
    % variaion
    % if ~ismember(size(stimFreqMedian, 2), train_options)
    %     disp("WARNING: Irregular number of trains found. Found " + num2str(size(stimFreqMedian, 2)) + " trains")
    %     proceed = input("Do you wish to proceed (0/1)?");
    %     if proceed == 0
    %         error("Irregular number of trains");
    %     end
    % end
    
    disp('');
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

    %% CUT AROUND TRAINS AND DETREND
    nChan = size(EEG.data,1);
    interpEpochData = zeros(nChan,3*srate,size(burstStartEnd,1));
    for iChan = 1:nChan
        
        disp(['Channel ' num2str(iChan)]);
        EEGdata = EEG.data(iChan,:);
        idxF3 = strcmp({EEG.chanlocs.labels}', 'F3');
        
        for iFreq = 1:size(burstStartEnd,1)
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
