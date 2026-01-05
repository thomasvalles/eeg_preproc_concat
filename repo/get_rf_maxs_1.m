function [rf_maxs] = get_rf_maxs_1(wv, intercept)
    int_adj = intercept - 1;
    wv = wv(:, 1:4);
    meds = median(wv, 2, "omitmissing");
    mins = min(wv, [], 2, "omitmissing");
    maxs = max(wv, [], 2, "omitmissing");
    nans = sum(isnan(wv),2);
    all_above = all(wv > int_adj, 2);
    three_above = sum(wv > int_adj, 2) >= 3;
    score = zeros(size(wv, 1), 1);
    num_above = sum(wv > int_adj, 2);
    %sos = meds.^2 + mins.^2;
    freqs = linspace(5, 18, 27)';
    freq_stats = table(freqs, meds, mins, maxs, all_above, three_above, num_above, nans, score);
   
    freq_stats.score(freq_stats.nans >= 2) = -100 + freq_stats.score(freq_stats.nans >= 2);
    freq_stats.score(freq_stats.nans == 1) = -3 + freq_stats.score(freq_stats.nans == 1);

    tmp = sortrows(freq_stats, 'meds', 'descend');
    largest_meds = tmp.freqs(1:5);
    freq_stats.score(freq_stats.freqs == largest_meds(1)) = 8 + freq_stats.score(freq_stats.freqs == largest_meds(1));
    freq_stats.score(freq_stats.freqs == largest_meds(2)) = 6 + freq_stats.score(freq_stats.freqs == largest_meds(2));
    freq_stats.score(freq_stats.freqs == largest_meds(3)) = 4 + freq_stats.score(freq_stats.freqs == largest_meds(3));
    freq_stats.score(freq_stats.freqs == largest_meds(4)) = 2 + freq_stats.score(freq_stats.freqs == largest_meds(4));
    freq_stats.score(freq_stats.freqs == largest_meds(5)) = 1 + freq_stats.score(freq_stats.freqs == largest_meds(5));

    tmp = sortrows(freq_stats, 'mins', 'descend');
    largest_mins = tmp.freqs(1:5);
    freq_stats.score(freq_stats.freqs == largest_mins(1)) = 3 + freq_stats.score(freq_stats.freqs == largest_mins(1));
    freq_stats.score(freq_stats.freqs == largest_mins(2)) = 2 + freq_stats.score(freq_stats.freqs == largest_mins(2));
    freq_stats.score(freq_stats.freqs == largest_mins(3)) = 1 + freq_stats.score(freq_stats.freqs == largest_mins(3));
    freq_stats.score(freq_stats.mins <= int_adj - 5) = -3 + freq_stats.score(freq_stats.mins <= int_adj - 5);

    tmp = sortrows(freq_stats, 'maxs', 'descend');
    largest_maxs = tmp.freqs(1:5);
    freq_stats.score(freq_stats.freqs == largest_maxs(1)) = 3 + freq_stats.score(freq_stats.freqs == largest_maxs(1));
    freq_stats.score(freq_stats.freqs == largest_maxs(2)) = 2 + freq_stats.score(freq_stats.freqs == largest_maxs(2));
    freq_stats.score(freq_stats.freqs == largest_maxs(3)) = 1 + freq_stats.score(freq_stats.freqs == largest_maxs(3));
    
    freq_stats.score(freq_stats.all_above) = 4 + freq_stats.score(freq_stats.all_above);
    freq_stats.score(freq_stats.three_above) = 2 + freq_stats.score(freq_stats.three_above);

    freq_stats.score(freq_stats.freqs >= 15) = 1 + freq_stats.score(freq_stats.freqs >= 15);

    ordered_freqs = sortrows(freq_stats, {'score' 'meds'}, 'descend');
    rf_maxs = ordered_freqs.freqs(1:3);
end
