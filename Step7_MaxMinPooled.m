%Step0_SetDirectories;


clear all
% subjids = ['sub-LA0001'; 'sub-LA0002'; 'sub-LA0003'; 'sub-LA0004'; 'sub-LA0005'; 'sub-LA0006'; 'sub-LA0007';];
%subjids = [subjids; ['sub-BH0645'; 'sub-BH0655'; 'sub-BH0663'; 'sub-BH0486'; 'sub-BH0669'; 'sub-BH0686'; 'sub-BH0688'; 'sub-BH0691']];

subjid = 'sub-LA0011';
% parent folder where you have Steps1-7 output for the three interrogations
parent_dir = ['/Volumes/My Passport/RFREQ_BIDS_DATABASE/derivatives/interrogation/' subjid filesep];
tmp_name = split(parent_dir, filesep);
sub_id = tmp_name{end - 1};

out_dir = [parent_dir 'final_frequency_order'];
mkdir(out_dir);


% get all the weighted_values.csv tables
wv_dirs = dir([parent_dir '*/Step7_FrequencyOutputDecision/*medians.mat']);
if isempty(wv_dirs)
    wv_dirs = dir([parent_dir '*/Step9_FrequencyOutputDecision/*medians.mat']);
end

if size(wv_dirs, 1) ~= 3
    error('Number of Step7 folders not exactly 3. Please ensure that the parent dir has exactly 3 Step7 folders');
end
% to get the intercept of the model
load('RegressionValues_SCC.mat');

wvs = [];
% collect the weighted values across the three days
for i_wv = 1:size(wv_dirs, 1)
    clear weighted_vals_sorted sorted_freqs
    load(fullfile(wv_dirs(i_wv).folder, wv_dirs(i_wv).name));
    target_spacing = median(diff(sorted_freqs));
    
    min_freq = min(sorted_freqs);
    max_freq = max(sorted_freqs);

    expected_freqs = min_freq:target_spacing:max_freq;

    % Step 2: Find the missing frequencies
    [~, missing_idx] = setdiff(expected_freqs, sorted_freqs);
    missing_freqs = expected_freqs(missing_idx);
    
    % Step 3: Combine sorted_freqs with missing_freqs and sort
    all_freqs = sort([sorted_freqs, missing_freqs]);
    
    % Step 4: Initialize new weighted_vals_sorted with NaNs
    new_weighted_vals_sorted = NaN(length(all_freqs), size(weighted_vals_sorted, 2));
    
    % Step 5: Populate new_weighted_vals_sorted with the existing values
    [~, existing_idx] = ismember(sorted_freqs, all_freqs);
    new_weighted_vals_sorted(existing_idx, :) = weighted_vals_sorted;
    
    % Output the updated arrays
    sorted_freqs = all_freqs;
    weighted_vals_sorted = new_weighted_vals_sorted;

    wvs = [wvs weighted_vals_sorted];
end

try
    sorted_freqs(1);
catch
    sorted_freqs = UniqueFreqs;
end
% min_freq = min(sorted_freqs);
% max_freq = max(sorted_freqs);
% make pooled boxplots
gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1500, 800]);

missing_none = sum(~isnan(wvs), 2) >= size(wvs, 2) * 0.75;
missing_some = sum(~isnan(wvs), 2) >= size(wvs, 2) * 0.5 & sum(~isnan(wvs), 2) <= size(wvs, 2) * 0.75;
missing_many = sum(~isnan(wvs), 2) < size(wvs, 2) * 0.5;
subplot(2, 1, 1)
boxplot(wvs(missing_none, :)', 'positions', sorted_freqs(missing_none), 'labels', sorted_freqs(missing_none), 'colors', 'b', 'Widths', 0.2);
hold on
boxplot(wvs(missing_some, :)', 'positions', sorted_freqs(missing_some), 'labels', sorted_freqs(missing_some), 'colors', 'm', 'Widths', 0.2);
boxplot(wvs(missing_many, :)', 'positions', sorted_freqs(missing_many), 'labels', sorted_freqs(missing_many), 'colors', 'r', 'Widths', 0.2);

yline(b{strcmp(b.vars, 'Intercept'), 'coefs'}, 'LineWidth', 0.5)
%title('SCC model, all three interrogations pooled', 'FontSize', 18);

% plot the medians
intercept = b{strcmp(b.vars, 'Intercept'), 'coefs'};
medians = nanmedian(wvs, 2);
[sorted_values, sorted_indices] = sort(medians, 'descend');
spacing = diff(sorted_freqs);
spacing = 1 / spacing(1);
rf_maxes = (sorted_indices(1:3) - 1) / spacing + (min(sorted_freqs));
rf_mins = (flip(sorted_indices(end - 2:end)) - 1)/ spacing + (min(sorted_freqs)) ;

xticks(min_freq:0.5:max_freq)

labs = cellstr(string(min_freq:0.5:max_freq));
labs{1 + (rf_maxes(1) - min_freq) / 0.5} = [labs{1 + (rf_maxes(1) - min_freq) / 0.5} '*'];
labs{1 + (rf_mins(1) - min_freq) / 0.5} = [labs{1 + (rf_mins(1) - min_freq) / 0.5} '/'];

xticklabels(labs)

xlabel(['Top three medians: ', mat2str(rf_maxes), '    Bottom three medians: ', mat2str(rf_mins)], 'FontSize', 15);
ylabel('Predicted outcome', 'FontSize', 15);
axis([min_freq - 0.5,max_freq + 0.5,1.2*min(wvs,[],'all'),1.2*max(wvs,[],'all')])

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



% we will use the std. error ranking method for the pooled data too.
[T_max, ranking_data_max] = make_wv_table(wvs, intercept, sorted_freqs);
[T_min, ranking_data_min] = make_wv_table(-wvs, -intercept, sorted_freqs);


% abandoned std err ranking method
% if ~isempty(ranking_data_max.final_order)
%     rf_max = ranking_data_max.final_order(1);
% else
% % if theres somehow no max according to the std error ranking system
% % e.g. all freqs below intercept
%     rf_max = sortrows(T_max, 'median', 'descend').Freq(1);
% end
% 
% if ~isempty(ranking_data_min.final_order)
%     rf_min = ranking_data_min.final_order(1);
% else
%     rf_min = sortrows(T_min, 'median', 'descend').Freq(1);
% end

rf_max = rf_maxes(1);
rf_min = rf_mins(1);

rf_max_text = "RF max";
rf_min_text = "RF min";

% % if rf_max == 10, use the second best max. - we abandoned this. just use
% 10hz twice
% if rf_max == 10
% 
%     rf_max_text = "Second max";
%     % if there are other candidates according to the std err
%     if numel(ranking_data_max) > 1
%         rf_max = ranking_data_max(2);
% 
%     % otherwise, 10 is the standout max, must choose the next best median
%     else % 
%         rf_max = sortrows(T_max(T_max.Freq ~= 10, :), 'median', 'descend').Freq(1);
%     end
% end
% 
% % if rf_max == 10, use the second best max.
% if rf_min == 10
% 
%     rf_min_text = "Second min";
%     % if there are other candidates according to the std err
%     if numel(ranking_data_min) > 1
%         rf_min = ranking_data_min(2);
% 
%     % otherwise, 10 is the standout max, must choose the next best median
%     else % 
%         rf_min = sortrows(T_min(T_min.Freq ~= 10, :), 'median', 'descend').Freq(1);
%     end
% end

% make the randomized order

tx_freqs = [10; rf_max; rf_min];
tx_types = ["10Hz"; rf_max_text; rf_min_text];


r_permutation = [1; 2; 3];%randperm(3);
tx_order = tx_freqs(r_permutation);
t_blinded = table();
t_blinded.Frequency = tx_order;

t_unblinded = table();
t_unblinded.Frequency = tx_order;
t_unblinded.Tx_type = tx_types(r_permutation);


writetable(t_blinded, [out_dir filesep sub_id '_blinded_order.csv']);
writetable(t_unblinded, [out_dir filesep sub_id '_unblinded_order.csv']);
save([out_dir filesep sub_id '_weighted_vals.mat'],'wvs','sorted_freqs');


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

log_file = fopen([out_dir filesep sub_id '_std_err_ranking_system.txt'], 'w');
fprintf(log_file, '%s\n', log_str);
fclose(log_file);
writetable(T_max, [out_dir filesep 'weighted_values.csv'])

sgtitle("Pooled interrogations- RF Max: " + rf_maxes(1) + ". RF Min: " + rf_mins(1), 'fontsize', 18);
saveas(gca,[out_dir filesep sub_id '_pooled_boxplots.jpg']);
