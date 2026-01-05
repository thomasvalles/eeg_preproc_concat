function [EEG] = f_FixCPzAndEOG(EEG)
    if(EEG.nbchan > 64)
        if EEG.chanlocs(65).labels == "EOG"
            EEG = pop_select(EEG,'channel', 1:65);
        else
            EEG = pop_select(EEG,'channel', 1:64);
        end
    end

  
    chan_32_data = EEG.data(32, :); % get what's currently in chan 32
    chan_32_data_diffs = diff(chan_32_data);
    non_zero_diffs = chan_32_data_diffs(chan_32_data_diffs ~= 0);

    if (EEG.chanlocs(32).labels == "EOG") 
         % if channel 32 is actually EOG and holds non-constant data (at least 5% of differences are nonzero)
        if (numel(non_zero_diffs) / numel(chan_32_data_diffs) >= 0.05)
            EEG = pop_select(EEG, 'nochannel', 32); % remove channel 32,
            cpz_data = zeros(1, EEG.pnts);
            EEG = f_InsertChannel(EEG, 'CPz', cpz_data, 32); % fill in channel 32 with 0's
            EEG = f_InsertChannel(EEG, 'VEOG', chan_32_data, size(EEG.data,1) + 1); % insert eog data at the end.
        else % otherwise the label is EOG but it doesn't actually have EOG data. Just rename
            EEG.chanlocs(32).labels = 'CPz';
        end
        % if cpz/eog is completely missing (fewer than 64 channels)
    elseif ((EEG.chanlocs(32).labels ~= "CPz") && (size(EEG.data, 1) < 64))
        cpz_data = zeros(1, EEG.pnts);
        EEG = f_InsertChannel(EEG, 'CPz', cpz_data, 32); % fill in channel 32 with 0's
    end
   

    if strcmp(EEG.comments(end-2:end), 'edf') && length(EEG.chanlocs(1).labels)>4
        [ EEG, str_ChanLab ] = f_CleanEdfChanLabels( EEG );           
    end
    str_ChanLab = {EEG.chanlocs(:).labels};
        
    %%%MAKE SURE CPZ IS ALL ZEROS
    if(numel(find(EEG.data(32,:) ~= 0)) > 0)
        disp("WARNING: CPZ NOT ALL ZEROS. ZEROING CPZ");
        EEG.data(32,:) = 0;
    end

    EEG.nbchan = size(EEG.data, 1);
end