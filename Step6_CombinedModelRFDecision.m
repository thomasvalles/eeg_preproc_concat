clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;
clc
close all

eeglab;

Step0_SetDirectories;
%%Set Nearest neighbor pairs
load(NN_ANT); % nearest neighbors list

%Load example eeg with same structure. Use this for channel locations
EEG = pop_loadset(sample_eeg);
dirData = [work_dir, subjid, filesep 'Step6_ASCOutput' suffix filesep];
% %Where to save output mat files
saveDir = [work_dir, subjid, filesep 'Step7_FrequencyOutputDecision' suffix filesep];
mkdir(saveDir);

%% load ID
%go to data directory and get patient IDs for all files in that directory.
cd(dirData)


% get the scc and power features. f_GetFeatures averages over the band and
% smooths over neighboring channels
[chans_by_trials_scc, allStimFreq] = f_GetFeatures("SCC", "matched", NN_ANT);
[chans_by_trials_power, ~] = f_GetFeatures("Power", [1, 8], NN_ANT);

% hard coded fix for BH0645
allStimFreq(allStimFreq == 14.75) = 14.5;

% the training set was normalized when deriving the coefficients. need to
% make the same transformation according to the training set mean/std
% deviation for each channel
% training_set_stats = readtable("merged_86_stats.csv");
% power_stats_sorted = zeros(EEG.nbchan, 2);
% scc_stats_sorted = zeros(EEG.nbchan, 2);
% chan_labels = {EEG.chanlocs.labels};
% for iChan = 1:numel(chan_labels)
%     power_row = training_set_stats(strcmp(training_set_stats.Var1, [chan_labels{iChan} '[1, 8]']), 2:3);
%     scc_row = training_set_stats(strcmp(training_set_stats.Var1, [chan_labels{iChan} 'matched']), 2:3);
% 
%     power_row = [power_row{1, 1} power_row{1, 2}];
%     scc_row = [scc_row{1, 1} scc_row{1, 2}];
% 
%     power_stats_sorted(iChan, :) = power_row;
%     scc_stats_sorted(iChan, :) = scc_row;
% end
% 
% chans_by_trials_power{1} = (chans_by_trials_power{1} - power_stats_sorted(:, 1)) ./ power_stats_sorted(:, 2);
% chans_by_trials_scc{1} = (chans_by_trials_scc{1} - scc_stats_sorted(:, 1)) ./ scc_stats_sorted(:, 2);

% get the combined coefficients and names
load('RegressionValues_SCC.mat');
ROI_SCC = {};
ROI_power = {};
coef_rows = b(contains(b.vars, "[1, 8]") | contains(b.vars, "matched"), 1);
for r = 1:numel(coef_rows)
    var_name_parts = split(coef_rows{r, 1}, ["[1, 8]", "matched"]);
    if contains(coef_rows{r, 1}, "[1, 8]")
        ROI_power{end+1} = var_name_parts{1};
    else
        ROI_SCC{end+1} = var_name_parts{1};
    end
end


weighted_vals = ones(size(chans_by_trials_scc{1},2), 1)*b{strcmp(b.vars, 'Intercept'), 'coefs'};

% accumulate SCC weights
for iSite = 1:numel(ROI_SCC)
    ind_channel = find(ismember({EEG.chanlocs.labels}, ROI_SCC{iSite}));
    weighted_vals = weighted_vals + chans_by_trials_scc{1}(ind_channel, :)'*b{strcmp(b.vars, [ROI_SCC{iSite} 'matched']), 'coefs'};
end

% accumulate power weights
for iSite = 1:numel(ROI_power)
    ind_channel = find(ismember({EEG.chanlocs.labels}, ROI_power{iSite}));
    weighted_vals = weighted_vals + chans_by_trials_power{1}(ind_channel, :)'*b{strcmp(b.vars, [ROI_power{iSite} '[1, 8]']), 'coefs'};
end


sorted_freqs = sort(unique(allStimFreq));
weighted_vals_sorted = zeros(numel(sorted_freqs), sum(allStimFreq == min(allStimFreq)));
for iFreq = 1:numel(sorted_freqs)
    tmp_vals = weighted_vals(allStimFreq == sorted_freqs(iFreq))';
    tmp_vals = tmp_vals(1:min(size(weighted_vals_sorted, 2),size(tmp_vals, 2)));
    to_inject = NaN(1, size(weighted_vals_sorted, 2));
    to_inject(1:size(tmp_vals, 2)) = tmp_vals;
    weighted_vals_sorted(iFreq, :) = to_inject;
end
min_freq = min(sorted_freqs);
max_freq = max(sorted_freqs);    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create boxplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1500, 800]);

missing_none = sum(~isnan(weighted_vals_sorted), 2) >= size(weighted_vals_sorted, 2) * 0.75;
missing_some = sum(~isnan(weighted_vals_sorted), 2) >= size(weighted_vals_sorted, 2) * 0.5 & sum(~isnan(weighted_vals_sorted), 2) <= size(weighted_vals_sorted, 2) * 0.75;
missing_many = sum(~isnan(weighted_vals_sorted), 2) < size(weighted_vals_sorted, 2) * 0.5;
subplot(2, 1, 1)
boxplot(weighted_vals_sorted(missing_none, :)', 'positions', sorted_freqs(missing_none), 'labels', sorted_freqs(missing_none), 'colors', 'b', 'Widths', 0.2);
hold on
boxplot(weighted_vals_sorted(missing_some, :)', 'positions', sorted_freqs(missing_some), 'labels', sorted_freqs(missing_some), 'colors', 'm', 'Widths', 0.2);
boxplot(weighted_vals_sorted(missing_many, :)', 'positions', sorted_freqs(missing_many), 'labels', sorted_freqs(missing_many), 'colors', 'r', 'Widths', 0.2);

yline(b{strcmp(b.vars, 'Intercept'), 'coefs'}, 'LineWidth', 0.5)
%title('SCC Model', 'FontSize', 18);


intercept = b{strcmp(b.vars, 'Intercept'), 'coefs'};
medians = nanmedian(weighted_vals_sorted, 2);
[sorted_values, sorted_indices] = sort(medians, 'descend');
sorted_indices = sorted_indices(~isnan(sorted_values));
sorted_values = sorted_values(~isnan(sorted_values));

if length(sorted_freqs) > 1
    spacing = diff(sorted_freqs);
    spacing = 1 / spacing(1);
    rf_maxes = sorted_freqs(sorted_indices(1:3));%(sorted_indices(1:3) - 1) / spacing + (min(allStimFreq));
    rf_mins = flip(sorted_freqs(sorted_indices(end-2:end)));%(flip(sorted_indices(end - 2:end)) - 1)/ spacing + (min(allStimFreq)) ;
else
    spacing = 1;
    rf_maxes = [sorted_freqs(1); sorted_freqs(1); sorted_freqs(1)];
    rf_mins = [sorted_freqs(1); sorted_freqs(1); sorted_freqs(1)];
end

xticks(min_freq:0.5:max_freq)

labs = cellstr(string(min_freq:0.5:max_freq));
labs{1 + (rf_maxes(1) - min_freq) / 0.5} = [labs{1 + (rf_maxes(1) - min_freq) / 0.5} '*'];
labs{1 + (rf_mins(1) - min_freq) / 0.5} = [labs{1 + (rf_mins(1) - min_freq) / 0.5} '/'];

xticklabels(labs)


xlabel(['Top 3 medians: ', mat2str(rf_maxes), '    Bottom 3 medians: ', mat2str(rf_mins)], 'FontSize', 15);
ylabel('Predicted outcome', 'FontSize', 15);
axis([min_freq - 0.5,max_freq + 0.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])

save([saveDir subjid(1:10) '_' subjid(12:end) '_medians.mat'],'weighted_vals_sorted','sorted_freqs')

subplot(2,1,2)
plot(sorted_freqs,medians,'LineWidth',2)
    title('Regression Medians Order','FontSize',20)

axis([min_freq - 0.5,max_freq + 0.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])
xticks(min_freq:0.5:max_freq)
xticklabels(min_freq:0.5:max_freq)

hold on
scatter(rf_maxes(1),medians(sorted_indices(1)),100,'filled','g')
scatter(rf_mins(1),medians(sorted_indices(end)),100,'filled','r')
text(rf_maxes(1), sorted_values(1), num2str(round(sorted_values(1), 3)));
text(rf_mins(1), sorted_values(end), num2str(round(sorted_values(end), 3)));
yline(b{strcmp(b.vars, 'Intercept'), 'coefs'}, 'LineWidth', 0.5)

rfs = array2table(weighted_vals_sorted);
rfs=  cat(2,table(sorted_freqs'),rfs);

%% Clearly indicate rf max and rf min in a single csv %%
summary_table = table();
summary_table(:, 'RF Max') = {rf_maxes(1)};
summary_table(:, 'RF Min') = {rf_mins(1)};
writetable(summary_table, [saveDir 'results.csv'])

%% Get the full ranking of the spectrum %%
ranking = (sorted_indices - 1) / spacing + (min(allStimFreq));

full_table = table();
full_table.Rank = (1:numel(ranking))';
full_table.Frequency = ranking;
writetable(full_table, [saveDir 'full_ranking.csv'])

%% Get a table of the weighted values, make the std err ranking system log %%
[T_max, ranking_data_max] = make_wv_table(weighted_vals_sorted, intercept, sorted_freqs);
[T_min, ranking_data_min] = make_wv_table(-weighted_vals_sorted, -intercept, sorted_freqs);

% For maxes:
log_str = sprintf("New ranking method for maxes-\n\n");
log_str = log_str + sprintf("RF_max by median: %f\n", ranking_data_max.best_freq);
log_str = log_str + sprintf("Max median: %f\n", ranking_data_max.median_of_best);
log_str = log_str + sprintf("Frequencies above intercept: %s\n", mat2str(ranking_data_max.above_intercept));
log_str = log_str + sprintf("Frequencies within one standard error of max median: %s\n", mat2str(ranking_data_max.within_std_err_mx));
log_str = log_str + sprintf("Candidate maxes: %s\n", mat2str(ranking_data_max.candidate_freqs));

if ~isempty(ranking_data_max.final_order)
    log_str = log_str + sprintf("Candidate maxes sorted by smallest std err: %s\n\n\n", mat2str(ranking_data_max.final_order));
else
    log_str = log_str + sprintf("No candidate maxes\n\n\n");
end

% For mins:
log_str = log_str + sprintf("New ranking method for mins-\n\n");
log_str = log_str + sprintf("RF_min by median: %f\n", ranking_data_min.best_freq);
log_str = log_str + sprintf("Min median: %f\n", -ranking_data_min.median_of_best);
log_str = log_str + sprintf("Frequencies below intercept: %s\n", mat2str(ranking_data_min.above_intercept));
log_str = log_str + sprintf("Frequencies within one standard error of min median: %s\n", mat2str(ranking_data_min.within_std_err_mx));
log_str = log_str + sprintf("Candidate mins: %s\n", mat2str(ranking_data_min.candidate_freqs));

if ~isempty(ranking_data_min.final_order)
    log_str = log_str + sprintf("Candidate mins sorted by smallest std err: %s", mat2str(ranking_data_min.final_order));
else
    log_str = log_str + "No candidate mins";
end

log_file = fopen([saveDir 'std_err_ranking_system.txt'], 'w');
fprintf(log_file, '%s\n', log_str);
fclose(log_file);

writetable(T_max, [saveDir 'weighted_values.csv'])
if length(sorted_freqs) > 1
    sgtitle("RF Max: " + rf_maxes(1) + ". RF Min: " + rf_mins(1), 'fontsize', 18);
else
    sgtitle("Treated: " + sorted_freqs(1) + "Hz");
end

saveas(gca,[saveDir subjid(1:10) '_' subjid(12:end) '_Boxplots_SCC_model.jpg']);

