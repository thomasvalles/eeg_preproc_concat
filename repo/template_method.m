train = 20;
pulses_per_train = 40;

% how long the template should be (arbitrary)
template_length = 0.020 * srate;

% we always ignore 20ms after the pulse
buffer = 0.020 * srate;

sig = EEG.data(5,  burstStartEnd(train, 1) - 10 : burstStartEnd(train, 2));
x = burstStartEnd(train, 1) - 10 : burstStartEnd(train, 2);
[psor,lsor] = findpeaks(sig,x,'SortStr','descend');

% the top 40 peaks will represent the pulses
psor = psor(1:40);
lsor = lsor(1:40);

% remove the last pulse
[~, last_pulse_ind] = max(lsor);
lsor(last_pulse_ind) = [];
psor(last_pulse_ind) = [];

figure
plot(x, sig);
hold on
scatter(lsor, psor, "filled", "Color", "r");

template_tail = zeros(pulses_per_train - 1, template_length);
idx_matrix = lsor(:) + buffer + (0:template_length-1) - min(x);
template = sig(idx_matrix);
template = mean(template, 1);

figure
plot(template);

