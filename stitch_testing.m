nepochs = 108;
nsamples = 8000;
nrandom = 107;
nsample_diffs = 10;

between_epoch_diffs = zeros(nepochs-1, 1);
for i = 2:nepochs
    between_epoch_diffs(i - 1) = abs(interpEpochData(5, 1, i) - interpEpochData(5, end, i-1));
end

rand_epochs = randi(nepochs, nrandom);
rand_times = nsample_diffs + randi(nsamples, nrandom);
random_diffs = zeros(nrandom, nsample_diffs);

for i = 1:nrandom
    for j = 1:nsample_diffs
        random_diffs(i, j) = abs(interpEpochData(5, rand_times(i), rand_epochs(i)) - interpEpochData(5, rand_times(i) - j, rand_epochs(i)));
    end
end

boxplot([random_diffs between_epoch_diffs], 'Labels',{'1 sample', '2 sample', '3 sample', ...
    '4 sample', '5 sample', '6 sample', '7 sample', '8 sample', '9 sample', '10 sample', 'Between-epoch'});

[p, H] = ranksum(random_diffs(:, 1), between_epoch_diffs)
[p, H] = ranksum(random_diffs(:, 5), between_epoch_diffs)
[p, H] = ranksum(random_diffs(:, 10), between_epoch_diffs)

ylabel('Difference between samples');

%%
figure;
for i = 1:8
    subplot(2, 4, i);
    x = -20:19;
    ep = 1 + randi(nepochs-1);
    plot(x, [interpEpochData(5, end-19:end, ep-1)'; interpEpochData(5, 1:20, ep)'], 'Marker', 'o')
    xline(-0.5, 'Color', 'r')
    title("Epochs " + num2str(ep - 1) + "-" + ep);
    xlabel('Samples')
end