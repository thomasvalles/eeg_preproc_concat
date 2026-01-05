function [chans_by_trials, allStimFreq] = f_GetFeatures(scc_or_power, band, NN_ANT)
% f_GetFeatures - For each channel and each train, get the SCC/power
% values, averaged over an output band, smoothed (median) over neighboring
% channels
%
% Input:
%   scc_or_power:   "SCC" or "Power"
%   band:           The output band of interest. For now, either "matched"
%                   or [1, 8]
%   NN_ANT:         The nearest neighbors list                 
%
% Output:
%   chans_by_trials: array with dimensions n_chans x n_trials (64 x 140)
%                    containing the SCC/power values
%   allStimFreq:     ordered stimulation frequencies (140 x 1)
%
% T. Valles 2024

  
    if scc_or_power == "SCC"
        files = dir('SCC_F3*.mat');
    elseif scc_or_power == "Power"
        files =  dir('Power_forPy*.mat');
    end
    files = {files.name}';
    patID = [];
    for iFile = 1:numel(files)
        tmpsplit = strsplit(files{iFile}, '_');
        patID{iFile} = tmpsplit{3};
    end
    nPat = numel(files);
    %nPat = 1;
    %--------------------------------------------------------------------------
    
    AllPRE_vals = cell(1,nPat);
    AllPOST_vals = cell(1,nPat);
    AllStimFreq = cell(1,nPat);
    AllChanLabels = cell(64,nPat);
    
    EveryStimFreq = [];
    
    
    %%Load files in for the subset above.
    for iPat = 1:nPat
        load(files{iPat})
        if scc_or_power == "SCC"
            AllPRE_vals{iPat} = SCCPre;
            AllPOST_vals{iPat} = SCCPost;
        elseif scc_or_power == "Power"
            AllPRE_vals{iPat} = RelPre;
            AllPOST_vals{iPat} = RelPost;
        end
        AllStimFreq{iPat} = round(allStimFreq*4)/4;
    
        AllChanLabels(:,iPat) = chanLabels;
        EveryStimFreq = [EveryStimFreq allStimFreq];
    end
    
    
    %--------------------------------------------------------------------------
    %% Choosing Channels for analysis
    

    for iPat = 1:nPat
        if size(AllPRE_vals{iPat}, 1) ~= 95
            AllPRE_vals{iPat} = permute(AllPRE_vals{iPat}, [3, 2, 1]);
            AllPOST_vals{iPat} = permute(AllPOST_vals{iPat}, [3, 2, 1]);
        end
        
        pre_vals = AllPRE_vals{iPat};
        post_vals = AllPOST_vals{iPat};
        allStimFreq = AllStimFreq{iPat};
        chans_by_trials{iPat} = zeros(64, numel(allStimFreq));
        
        if strcmp(band, "matched")

            for iStim = 1:numel(allStimFreq)
                freq = allStimFreq(iStim);
                FreqofInterest = find(v_FreqAxis >= freq - 2 & v_FreqAxis <= freq + 2);
                pre_of_interest = pre_vals(FreqofInterest, :, iStim);
                post_of_interest = post_vals(FreqofInterest, :, iStim);
                chans_by_trials{iPat}(:, iStim) = nanmean(post_of_interest)' - nanmean(pre_of_interest)';
            end

        else
            FreqofInterest = find(v_FreqAxis >= min(band) & v_FreqAxis <= max(band));
            pre_of_interest = pre_vals(FreqofInterest, :, :);
            post_of_interest = post_vals(FreqofInterest, :, :);
            chans_by_trials{iPat} = nanmean(post_of_interest) - nanmean(pre_of_interest);
            chans_by_trials{iPat} = squeeze(chans_by_trials{iPat}(1, :, :));

        end
    end


    c_b_t_copy = chans_by_trials{iPat};
    for ichan = 1:64
        to_smooth_over = [ichan NN_ANT{ichan}];
        chans_by_trials{iPat}(ichan, :) = nanmedian(c_b_t_copy(to_smooth_over,:));
    end
    %standardize
end