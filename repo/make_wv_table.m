function [T, ranking_data] = make_wv_table(weighted_vals_sorted, intercept, sorted_freqs)
% Input:
%   weighted_vals_sorted: Matrix of weighted values calculated by the
%   model. Should by n_freqs x n_trains 
%
%   intercept:     The intercept of the SCC model
%
%   sorted_freqs: The frequencies sorted in increasing order(e.g. 5:20)
%
% Output:
%   T:             weighted_vals as a table, with median and std_err colunns 
%
%   ranking_data.best_freq: maximum freq based on median
%
%   ranking_data.median_of_best: max median
%
%   ranking_data.above_intercept: freqs above intercept
%
%   ranking_data.within_std_err_mx: freqs within one std err of max median
%
%   ranking_data.candidate_freqs: intersection of above_intercept and
%   within_std_err_mx
%
%   ranking_data.final_order: candidate freqs sorted by std error (smaller
%   is better)



    % set up table with median, std err columns
    Freq = sorted_freqs';%(min_freq:min_freq + size(weighted_vals_sorted, 1) - 1)';
    T = array2table(weighted_vals_sorted);
    T.Properties.VariableNames = arrayfun(@num2str, 1:size(weighted_vals_sorted, 2), 'UniformOutput', false);
    T = addvars(T, Freq, 'Before', 1);
    T.median = nanmedian(weighted_vals_sorted, 2);

    % std_err = std / sqrt(n)
    T.std_err = nanstd(weighted_vals_sorted, 0, 2) ./ sqrt(sum(~isnan(weighted_vals_sorted), 2));
    
    T = T(~isnan(T.median), :);
    % sort descending so we get largest on top
    maxes = sortrows(T, 'median', 'descend').Freq;

    % find freqs which have median above intercept, and whose medians are
    % within a std. error of the best median
    median_of_max = T{T.Freq == maxes(1), 'median'};
    std_err_of_max = T{T.Freq == maxes(1), 'std_err'};
    above_intercept = T{T.median > intercept, 'Freq'};
    within_std_err_mx = T{abs(T.median - median_of_max) < std_err_of_max, 'Freq'};
    candidate_freqs = intersect(above_intercept, within_std_err_mx);
    
    % sort the candidates by std err
    candidate_table = T(ismember(T.Freq, candidate_freqs), :);

    % sort ascending so that we get the smallest on top
    final_order = sortrows(candidate_table, 'std_err').Freq; 

    ranking_data.best_freq = maxes(1);
    ranking_data.median_of_best = median_of_max;
    ranking_data.above_intercept = above_intercept;
    ranking_data.within_std_err_mx = within_std_err_mx;
    ranking_data.candidate_freqs = candidate_freqs;
    ranking_data.final_order = final_order;
   
end 