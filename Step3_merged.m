eeglab
        
clear all 
close all
Step0_SetDirectories;

%% load the eeg
sub = 'SCC24035';
ses = '01a01';
target_SR = 1000;
epoch_length = 4; % seconds
input_dir = ['/Volumes/leuchter/BIDS_DATABASE/derivatives/int_rs_concat/sub-' sub filesep 'ses-' ses filesep 'Step3_CutFiles'];
output_dir = ['/Volumes/leuchter/BIDS_DATABASE/derivatives/int_rs_concat/sub-' sub filesep 'ses-' ses filesep 'Step4_PipelineOutput'];
if ~isfolder(output_dir)
    mkdir(output_dir);
end

set_files = dir([input_dir filesep '*.set']);
ALLEEG = struct([]);
for ifile = 1: numel(set_files)
    EEG = pop_loadset(fullfile(set_files(ifile).folder, set_files(ifile).name));
    EEG = pop_resample(EEG, target_SR);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',50,'plotfreqz',0);
    [ALLEEG, ~, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end

% make a table for bad channels
n_rows = 64 * numel(set_files);
T = table('Size', [n_rows, 2], ...
          'VariableTypes', {'string', 'int32'}, ...
          'VariableNames', {'segment_name', 'channel'});

% sort by run number
runs = [ALLEEG.run];
[~, order] = sort(runs);
ALLEEG = ALLEEG(order);

%% identify and interpolate the bad channels by segment
row_idx = 1;
for ieeg = 1:numel(ALLEEG)

    % Detect bad channels
    [~, bad_chans, ~] = clean_channels_RS( ...
        eeg_epoch2continuous(ALLEEG(ieeg)), ...
        .80, [], [], 0.3, [], [] );

    bad_chans = find(bad_chans);

    % Log table entries
    for chan_idx = 1:numel(bad_chans)
        T.segment_name(row_idx) = string(erase(ALLEEG(ieeg).filename, '_EpochData.set'));
        T.channel(row_idx)      = bad_chans(chan_idx);
        row_idx = row_idx + 1;
    end

    % Remove bad channels FROM THIS DATASET
    EEG_remchan = pop_select(ALLEEG(ieeg), 'nochannel', bad_chans);

    % Interpolate back to original channel locations FROM THIS DATASET
    tmp = pop_interp(EEG_remchan, ALLEEG(ieeg).chanlocs, 'spherical');

    tmp = eeg_epoch2continuous(tmp);

    tmp.event(end+1).type = "end_run-" + tmp.run;
    tmp.event(end).latency = tmp.pnts + 0.5;   % place event after end
    tmp = eeg_checkset(tmp, 'eventconsistency');
    % Store back into ALLEEG properly
    [ALLEEG, tmp] = eeg_store(ALLEEG, tmp, ieeg);

end


T = rmmissing(T);
writetable(T, [output_dir filesep 'bad_chans.csv']);

%% merge and re-reference
merged = pop_mergeset(ALLEEG, 1:length(ALLEEG), 0);
merged = pop_reref(merged,[]);

%% split again
is_endrun = startsWith({merged.event.type}, 'end_run-');
endrun_events = find(is_endrun);
end_lat = round([merged.event(endrun_events).latency]);
end_lat = sort(end_lat);
segment_starts = [1, end_lat(1:end-1) + 1];
segment_ends = end_lat;
segments = [segment_starts(:), segment_ends(:)];

ALLEEG = [];  % reset container
EEGOUT = merged; % template

for k = 1:size(segments,1)
    EEG_tmp = pop_select(merged, 'point', segments(k,:));

    remove_idx = find(strcmpi({EEG_tmp.event.type}, 'boundary'));
    if ~isempty(remove_idx)
        EEG_tmp.event(remove_idx) = [];
        EEG_tmp = eeg_checkset(EEG_tmp, 'eventconsistency');
    end

    EEG_tmp = eeg_regepochs(EEG_tmp, epoch_length);
    EEG_tmp.setname = sprintf('Segment_%d', k);
    [ALLEEG, EEGOUT] = eeg_store(ALLEEG, EEG_tmp, k);
end

%% reject noisy epochs based on amplitude z-score and high-freq amp z-score
% TODO: consistency of rejected epochs in interrogation segment
for ieeg = 1:numel(ALLEEG)
    [rmepochs, rejE]= zscore_epoch_rej(ALLEEG(ieeg),3); 
    tmp = pop_select(ALLEEG(ieeg),'notrial',rmepochs);

    [rmepoch_f] = rej_noisy(tmp, 15, 49.5, 3);
    tmp = pop_select(tmp, 'notrial', rmepoch_f);

    [ALLEEG, tmp] = eeg_store(ALLEEG, tmp, ieeg);
end