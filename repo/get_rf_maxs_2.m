function [rf_maxs] = get_rf_maxs_2(wv, intercept)
    int_adj = intercept - 0;
    %wv = wv(:, 1:4);
    median_ = median(wv, 2, "omitmissing");
    min_ = min(wv, [], 2, "omitmissing");
    max_ = max(wv, [], 2, "omitmissing");
    num_above = sum(wv > int_adj, 2);
    %sos = meds.^2 + mins.^2;
    %Freq = linspace(5, 18, 27)';
    Freq = linspace(5, 18, size(wv, 1))';
    freq_stats = table(Freq, num_above, median_, min_, max_);

    ordered_freqs = sort_w_ties(freq_stats);
    rf_maxs = flipud(ordered_freqs.Freq);
    rf_maxs = rf_maxs(1:3);
end