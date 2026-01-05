%% Header
%Authors: Luke Acuff (Adapted from FASTER pipeline)
%Description: Modified FASTER Epoch Rejection
%Usage: Change inputs and Run
%Version 6/26/2024

%% Input
% clear, clc, close all
% data_path = 'C:\Users\WAcuff\OneDrive - Care New England\Documents\Projects\RFREQ\interrogation';
% eeglab_path = 'C:\Users\WAcuff\OneDrive - Care New England\Documents\MATLAB\eeglab2024.0';

function trialrej = main_epoch_rej(EEG)

    %channels of interest FP1, FPz, F2, F3 and their nearest neighbors
    %i.e. Fp1, Fpz, Fp2, F3, Fz, F4, FC2, AF7, AF3, AF4, AF8, F5, F1, F2, FC3
    chans = [1, 2, 3, 5, 6, 7, 11, 33, 34, 35, 36, 37, 38, 39, 41];
    
    %% Settings
    HF_low = 31;
    HF_high = 56;
    max_remove = 0.15;
    epoch_time_cutoff = 2;
    
    %% Base Thresholds
    % Values were selected based on observed values in 'bad' epochs as identified by a human rater.
    Rdev = 1; % Epoch Mean Channel Deviation from the Grand Channel Mean (Inspired by FASTER)
    Rvar = 80; % Mean Channel Variance (Inspired by FASTER)
    Ramp = 30; % Mean Channel Amplitude Range (Inspired by FASTER)
    RMamp = 100; % Max Channel Amplitude Range 
    Rhf = 0.7; % Mean Channel Average HF Power (Pwelch Maxhold, so each frequency bin is taken from max window within epoch before averaging)
    
    
    %% Load, Gen, and Save Scores
    %eeg_files = dir(['C:\Users\WAcuff\OneDrive - Care New England\Documents\Projects\RFREQ\interrogation' filesep '*' filesep 'Step4_FasterOutput' filesep '*.set']);
    epoch_table = {};
    
    
    %EEG = pop_loadset(filename); %load
    EEG = pop_epoch(EEG,{},[0 epoch_time_cutoff]); %remove last 1/3rd of data
    raw_scores = get_raw_scores(EEG, chans, HF_low, HF_high);
    
    
    %% Scoring Algo
    
    % Only relevant channels within the first two seconds are considered in each epoch.
    % All epochs which exceed at least one of the five thresholds are considered for rejection.
    % Epochs under consideration for rejection are scored and sorted by the proportion
    % each threshold was exceeded by. If more than 15% of all epochs are under
    % consideration for rejection, only the worst 15% of all epochs (determined by
    % score) are removed to avoid exceeding greater than 15% epoch removal.
    
    thresh_list = [Rdev,Rvar,Ramp,RMamp,Rhf];
    
    % create and process weight table
    weights_table = raw_scores;
    for i = 1:5
        weights_table(:,i) = weights_table(:,i)/thresh_list(i) - 1;
    end
    weights_table(weights_table(:,1:5) < 0) = 0;
    weights_table(:,7) = sum(weights_table(:,1:5),2);
    weights_table(weights_table(:,7) == 0,:) = [];
    weights_table = sortrows(weights_table,7,'descend');
    
    if length(weights_table) > length(raw_scores)*max_remove
        epochs_to_remove = weights_table(1:round(length(raw_scores)*max_remove),:);
    else
        epochs_to_remove = weights_table;
    end
    
    fileid = 1;
    % Save recommendation to file
    epoch_table(fileid,1) = {"test"};
    epoch_table(fileid,2:length(epochs_to_remove(:,6)) + 1) = num2cell(epochs_to_remove(:,6));
    trialrej = sort(uint8(epochs_to_remove(:, 6)));
end
