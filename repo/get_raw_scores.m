function raw_scores = get_raw_scores(EEG, eeg_chans, HF_low, HF_high)
    list_properties = [];
    means = mean(EEG.data(eeg_chans,:),2);
    % 1 Epoch's mean deviation from channel means.
    for u = 1:size(EEG.data,3)
	    list_properties(u,1) = mean(abs(squeeze(mean(EEG.data(eeg_chans,:,u),2)) - means));
    end
    
    % 2 Epoch variance
    list_properties(:,2) = mean(squeeze(var(EEG.data(eeg_chans,:,:),0,2)));
    
    % 3 Mean Max amplitude difference
    for t = eeg_chans
	    for u = 1:size(EEG.data,3)
		    ampdiffs(t,u) = max(EEG.data(t,:,u)) - min(EEG.data(t,:,u));
	    end
    end
    list_properties(:,3) = mean(ampdiffs,1);

    % 4 Max Channel Max amplitude difference
    list_properties(:,4) = max(ampdiffs);

    % 5 HF Max Score
    max_freq = zeros(length(EEG.epoch),1);
    for i = 1:length(EEG.epoch)
        [spectra,freqs] = pwelch(double(EEG.data(eeg_chans,:,i)'),[],125,1000,EEG.srate,'maxhold');
        max_freq(i) = mean(mean(spectra(HF_low:HF_high,:)));
    end
    list_properties(:,5) = max_freq;

    % 6 Epoch
    list_properties(:,6) = 1:length(max_freq);%length(ampdiffs);

    raw_scores = list_properties;
end