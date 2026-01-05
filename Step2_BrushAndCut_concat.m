close all;
clear all;
eeglab
        
Step0_SetDirectories;
device = 'Magstim';

%% load the eeg
input_dir = '/Volumes/leuchter/BIDS_DATABASE/rawdata';
sub = 'SCC24035';
ses = '01a01';
output_dir = ['/Volumes/leuchter/BIDS_DATABASE/derivatives/int_rs_concat/sub-' sub filesep 'ses-' ses filesep 'Step3_CutFiles'];
if ~isfolder(output_dir)
    mkdir(output_dir);
end

set_files = dir([input_dir filesep 'sub-' sub filesep 'ses-' ses filesep 'eeg' filesep '*.set']);
EEG = pop_loadset(fullfile(set_files.folder, set_files.name));


%% request user to select the segments
segments = struct([]);
more_segments = 1;
isegments = 1;
while more_segments
    figure;
    plot(EEG.data(5, :));
    hold on
    % get the non-TMS events
    events = EEG.event((~strcmp({EEG.event.type}, '0001')) & (~strcmp({EEG.event.type}, '1')));
    
    T = struct2table(events);  % Convert struct to table
    T_sorted = sortrows(T, 'latency');  % Sort table by 'latency'
    events = table2struct(T_sorted)';  % Convert back to struct
    
    % add the vertical lines marking these events
    for ievent = 1:size(events, 2)
        xline(events(ievent).latency, 'color', 'r', 'Label', events(ievent).type);
        %text(events(ievent).latency, mean(EEG.data(5, :)), num2str(ievent), 'FontSize', 16)
    end

    title('Please select segment and name it');

    brush on
    figure(gcf);
    segment_name = input('Please select segment and name it', 's');
    [brushX] = f_GetBrush;
    segments(isegments).name = segment_name;
    segments(isegments).start = min(brushX);
    segments(isegments).end = max(brushX);
    close
    more_segments = input('More segments? (0/1)');
    isegments = isegments + 1;
end
% sort them by the start time
[~, idx] = sort([segments.start]);
segments = segments(idx);

% select the interrogation segment. will split the pre/post segments later
interrogation_segment = segments(strcmp({segments.name}, 'Interrogation'));
EEG_int = pop_select(EEG, 'point', [interrogation_segment.start interrogation_segment.end]);

%% process the interrogation segment
pulses_per_train = 40;
idxF3 = 5;
srate = EEG.srate;
label = {EEG.chanlocs.labels}';
F3 = EEG_int.data(idxF3,:);

% if there are any event markers find them
if(numel(EEG_int.event) > 0)
    idxValidTrigger = strcmpi({EEG_int.event.type}', '0001');
    if nnz(idxValidTrigger) == 0
        idxValidTrigger = strcmpi({EEG_int.event.type}', '01');
    end
    
    if nnz(idxValidTrigger) == 0
        idxValidTrigger = strcmpi({EEG_int.event.type}', '1');
    end
    

    pulseTriggers = [EEG_int.event.latency]';
    shiftedTriggers=pulseTriggers; 
    save(fullfile(output_dir, 'TriggersShifted'), 'shiftedTriggers', 'srate');
else
    shiftedTriggers = [];
end

trigInterrog = shiftedTriggers;

% sometimes the triggers are slightly shifted? If so just determine them by
% pulse height
close;
figure;
plot(EEG_int.data(5, :))
if ~isempty(trigInterrog)
    yLimits = ylim;
    line([trigInterrog(:) trigInterrog(:)]', ...
         repmat(yLimits, numel(trigInterrog), 1)', ...
         'Color', 'r');
end

triggers_ok = input('Do the triggers look okay (0/1)?');
if ~triggers_ok
    trigInterrog = [];
end

 %% FIND PULSE TRAINS  
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
        train_options = [51; 71; 75; 76; 80; 108; 140; 160];
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
        [trigInterrog] = f_DetectTriggersManuallyAuto_varHz(F3,srate, [1 size(EEG_int.data, 2)], n_trains, pulses_per_train, start_cutoff, option);

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
title([sub '. ' num2str(numel(stimFreqMedian)) ' Pulse Trains Detected']);
saveas(gcf, [output_dir, filesep, sub, '_MarkedTrains', '.png']);

save(fullfile(output_dir, ['TriggersInterrog.mat']), 'burstStartEnd', 'burstIPI', 'stimFreqMedian', 'segments');        
close all;


%% CUT AROUND TRAINS AND DETREND
nChan = size(EEG_int.data,1);
interpEpochData = zeros(nChan,4*srate,size(burstStartEnd,1));
for iChan = 1:size(EEG_int.data, 1)
    
    disp(['Channel ' num2str(iChan)]);
    EEGdata = EEG_int.data(iChan,:);
    idxF3 = strcmp({EEG_int.chanlocs.labels}', 'F3');
    
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

% stitch beginning of each epoch to end of previous
for iepoch = 2:size(interpEpochData, 3)
    for iChan = 1:size(interpEpochData, 1)
        difference = interpEpochData(iChan, 1, iepoch) - interpEpochData(iChan, end, iepoch-1);
        interpEpochData(iChan, :, iepoch) = interpEpochData(iChan, :, iepoch) - difference;
    end
end

% save the interrogation segment separately
setname = fullfile(output_dir, 'Interrogation_EpochData.set');
[ EEG_int_out ] = f_CreateEeglabStruct( interpEpochData, srate,label, setname );
EEG_int_out = pop_chanedit(EEG_int_out, 'lookup',ElectrodeLocation);
EEG_int_out.run = find(strcmp({segments.name}, 'Interrogation'));
EEG_int_out = pop_saveset(EEG_int_out, 'filename', setname);

% ensure no overlap between  segments
idx_interrogation = find(string({segments.name}) == 'Interrogation');

% if there is a segment before the interrogation
if idx_interrogation > 1
    % this is the "true" timing of the first pulse
    % there is only one segment before the interrogation.
    first_pulse_start = burstStartEnd(1, 1) + segments(idx_interrogation).start;
    % interrogation segment will start 4 seconds before the first pulse
    int_segment_start = first_pulse_start - srate * 4;
    % end of pre rest is right before the beginning of interrogation if there was
    % an overlap
    segments(idx_interrogation - 1).end = min(segments(idx_interrogation - 1).end, int_segment_start - 1);
end

% if there is a segment after the interrogation
if idx_interrogation < numel(segments)
    % this is the "true" timing of the last pulse
    last_pulse_end = burstStartEnd(end, 2) + segments(idx_interrogation).start;
    % interrogation segment will end 4 seconds after the last pulse
    int_segment_end = last_pulse_end + srate * 4;
    % start of post rest is 4 seconds after the last pulse if there was an
    % overlap
    segments(idx_interrogation + 1).start = max(segments(idx_interrogation).start, int_segment_end + 1);
end
% epoch the resting state segments


% epoch_length = 4;

% save segments separately
for iseg = 1:numel(segments)
    if ~strcmp(segments(iseg).name, 'Interrogation')
        EEG_tmp = pop_select(EEG, 'point', [segments(iseg).start segments(iseg).end]);
        setname = fullfile(output_dir, [segments(iseg).name '_EpochData.set']);
        % EEG_tmp = eeg_regepochs(EEG_tmp,epoch_length);
        EEG_tmp.run = iseg;
        EEG_tmp = pop_chanedit(EEG_tmp, 'lookup', ElectrodeLocation);
        EEG_tmp = pop_saveset(EEG_tmp, 'filename', setname);
        clear EEG_tmp
    end
end
