clear
clc
close all



eeglab;
% addpath('/Users/thomas/Documents/GitHub/RF_StepsMatlab/')
% addpath('/Users/thomas/Documents/interrogation/')
% addpath('/Users/Thomas/Documents/rf/algorithmic')
% subjid = 'RFUCLA-001_IP1';
% work_dir = '/Users/thomas/Documents/interrogation/';%'/Volumes/Files/Interrogation/';
% work_dir = '/Volumes/Files/Interrogation/5x5/SCC23-030/';

set_rf_directories;
%%Set Nearest neighbor pairs
load(NN_ANT); % nearest neighbors list

%Load example eeg with same structure. Use this for channel locations
EEG = pop_loadset(sample_eeg);
dirData = [work_dir, subjid, '/Step8_ASCOutput/'];
% %Where to save output mat files
saveDir = [work_dir, subjid, '/Step9_FrequencyOutputDecision/'];
mkdir(saveDir);

%% load ID
%go to data directory and get patient IDs for all files in that directory.
cd(dirData)
files = dir('SCC_NoSort*.mat');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{3};
end
nPat = numel(files);
%nPat = 1;
%--------------------------------------------------------------------------

AllPREASC = cell(1,nPat);
AllPOSTASC = cell(1,nPat);
AllStimFreq = cell(1,nPat);
AllChanLabels = cell(64,nPat);

EveryStimFreq = [];
%%Load files in for the subset above.
for iPat = 1:nPat
    load(files{iPat})
    AllPREASC{iPat} = aSCPreAll;
    AllPOSTASC{iPat} = aSCPostAll;
    AllStimFreq{iPat} = round(allStimFreq*4)/4;

% AllStimFreq{iPat} = allStimFreq;
    AllChanLabels(:,iPat) = chanLabels;
    EveryStimFreq = [EveryStimFreq allStimFreq];
    clear aSCPostAll aSCPreAll allStimFreq chanLabels
end


%--------------------------------------------------------------------------
%% Choosing Channels for analysis

%define variables
AnalysisPRE = cell(1,nPat);
AnalysisPOST = cell(1,nPat);
seedChanIdx = [5];
%indices to remove
idxM1M2 = ismember({EEG.chanlocs.labels}', {'M1', 'M2', 'F3'});

%%Chan to do = 1:64 for all channels; choose channels to analyze
chanToDo = 1:64;
chanToDo(idxM1M2) = [];

%Calculate and store all ASC pre, post for all subjects
EveryDiffASC = [];
for iPat = 1:nPat
    temp_matPre = [];
    temp_matPost = [];
    chan1_counter = 1;
    for iChan1 = seedChanIdx
        for iChan2 = chanToDo
            %account for data structures being i = x; j = x+1;
            if(iChan2 < iChan1)
                temp_matPre(:,:,chan1_counter, iChan2) = AllPREASC{iPat}(:,:,iChan2,iChan1);
                temp_matPost(:,:,chan1_counter, iChan2) = AllPOSTASC{iPat}(:,:,iChan2,iChan1);
            else
                temp_matPre(:,:,chan1_counter, iChan2) = AllPREASC{iPat}(:,:,iChan1,iChan2);
                temp_matPost(:,:,chan1_counter, iChan2) = AllPOSTASC{iPat}(:,:,iChan1,iChan2);
            end
        end
        chan1_counter = chan1_counter + 1;
    end
    
    %reduce dimensions and average for multiple seed channels but keep 3
    %dimensiodfsdns
    AnalysisPRE_NoAveraging{iPat} = squeeze(nanmean(temp_matPre,3));
    AnalysisPOST_NoAveraging{iPat} = squeeze(nanmean(temp_matPost,3));

    %%For later use; Diff ASC
    clear temp_matPost temp_matPre
end


%All
%%%%LOOP THROUGH OLD AND NEW PIPELINES
for iMethod = 1
    if(iMethod == 1)
        ChanSet = {{'Fpz'},{'AF4'}};
    else
        ChanSet = {{'FC3','AF3'},{'FT8','T8','TP8','CP6','P8'}};
    end
    
    
    %Loop through and find indices of channels
    ChansToUse = {};
    ChanString = '';
    for ilistChans = 1:numel(ChanSet)
        TmpChans = [];
        for j = 1:numel(ChanSet{ilistChans})
            TmpChans(j) = find(ismember({EEG.chanlocs.labels}, ChanSet{ilistChans}{j}));
            ChanString = strcat(ChanString, ChanSet{ilistChans}{j});
            ChanString = strcat(ChanString, ',');
        end
        
        %%%FIND CHANNELS TO SMOOTH%%%
        if(iMethod == 1)
            TmpChans = [TmpChans NN_ANT{TmpChans}];
            ChansToUse{ilistChans} = TmpChans;
        else
            ChansToUse{ilistChans} = TmpChans;
        end
        
    end
    
    %% CORRECTING FOR MAGSTIM ERROR - no longer necessary bc magventure
%     AllStimFreq{1}(1) = AllStimFreq{1}(2);
    
    %Set structure for pre and Post data
    AllPostData{iPat,1} = AnalysisPOST_NoAveraging{1};
    AllPreData{iPat,1} = AnalysisPRE_NoAveraging{1};
    
    if((exist('RemovedFreqs')) && (numel(RemovedFreqs) > 0))
        %Set nan values based on what was removed during VI
        AllPreData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;
        AllPostData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;    
    end
    %Set window for averaging. 0.1 to get like frequencies
    window = 0.1;
    UpdatedPre1toPost1_sorted = [];
    
    %%Run on Session 1
    UniqueFreqs = unique(AllStimFreq{1});
    UniqueVals_sorted = nan(numel(UniqueFreqs), 6,size(AllPostData{iPat,1},3));
   
    %%%%%Run for all
    for iFreq = 1:numel(UniqueFreqs) % look at the frequencies in order .. 
        tmp_freqs = AllStimFreq{1};
        window_avg = find(tmp_freqs >= UniqueFreqs(iFreq) - window & tmp_freqs <= UniqueFreqs(iFreq) + window);
        
        %Set frquencies of interest based on method
        if(iMethod == 1)
            FreqofInterest = find(v_FreqAxis >= UniqueFreqs(iFreq) - 2 & v_FreqAxis <= UniqueFreqs(iFreq) + 2);
        else
            FreqofInterest = find(v_FreqAxis >= 13 & v_FreqAxis <= 17);
        end
        
        %Average over frequencies of interest. First axis so that we have a
        %difference for each stimulation frequency at each channel 
        Pre1toPost1 = mean(AllPostData{iPat,1}(FreqofInterest,:,:) - AllPreData{iPat,1}(FreqofInterest,:,:),1);
        %Set uniquevalues
        % Uniquevals_sorted gives you change in SCC for each trial  of each stim freq.     
        UniqueVals_sorted(iFreq,1:numel(window_avg),:) = squeeze(Pre1toPost1(1,window_avg,:));   
    end
    
    
    %Average over channels based on method
    FinalValsPre1toPost1_sorted = [];
    for iSite = 1:numel(ChanSet)
        if(iMethod == 1)
            %% for each trial, take the median value for the ROI. (axis=3).
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmedian(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        else
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmean(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        end
    end
    
    %Load regression coefficients values from model
    if(iMethod == 1)
        load(regression_matched)
        ROI1 = nanmedian(FinalValsPre1toPost1_sorted(:,:,1),2);
        ROI2 = nanmedian(FinalValsPre1toPost1_sorted(:,:,2),2);
        % median at each ROI gets multiplied by the regression coefficient
        weighted_diag_regression = b(end) + b(1) * ROI1 + b(2) * ROI2;
    else
        load(regression_old)
    end
    %Loop through and multiply regression values to matrix
    weighted_vals_sorted = ones(size(FinalValsPre1toPost1_sorted,1),size(FinalValsPre1toPost1_sorted,2))*b(end);
    for iSite = 1:numel(ChanSet)
        weighted_vals_sorted = weighted_vals_sorted + FinalValsPre1toPost1_sorted(:,:,iSite)*b(iSite);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Create boxplots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(iMethod == 1)
        gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1500, 800]);
    end
    allthree = sum(~isnan(weighted_vals_sorted), 2) >=4;
    onlytwo = sum(~isnan(weighted_vals_sorted), 2) == 3;
    onlyone = sum(~isnan(weighted_vals_sorted), 2) <= 2;
    subplot(2, 1, iMethod)
    boxplot(weighted_vals_sorted(allthree, :)', 'positions', UniqueFreqs(allthree), 'labels', UniqueFreqs(allthree), 'colors', 'b', 'Widths', 0.2);
    hold on
    boxplot(weighted_vals_sorted(onlytwo, :)', 'positions', UniqueFreqs(onlytwo), 'labels', UniqueFreqs(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
    scatter(UniqueFreqs(onlyone), weighted_vals_sorted(onlyone, :), 'k', 'filled');
    yline(b(end), 'LineWidth', 0.5)
    
    if(iMethod == 1)
        title('WEIGHTED: Post - Pre, Matched Method', 'FontSize', 18);
    else
        title('WEIGHTED: Post - Pre, 13-17Hz', 'FontSize', 18);
    end
    
    b = [18.8565; 9.5724; 21.5616];
    intercept = b(3);
    wv = weighted_vals_sorted(:, 1:4);
    m1_maxs = get_rf_maxs_1(wv, intercept);
    m2_maxs = get_rf_maxs_2(wv, intercept);
    m1_mins = get_rf_maxs_1(-wv, -intercept);
    m2_mins = get_rf_maxs_2(-wv, -intercept);
  

    
    xticks(5:0.5:18)

    labs = cellstr(string(5:0.5:18));
    labs{1 + (m1_maxs(1) - 5) / 0.5} = [labs{1 + (m1_maxs(1) - 5) / 0.5} '+'];
    labs{1 + (m2_maxs(1) - 5) / 0.5} = [labs{1 + (m2_maxs(1) - 5) / 0.5} '*'];
    labs{1 + (m1_mins(1) - 5) / 0.5} = [labs{1 + (m1_mins(1) - 5) / 0.5} '-'];
    labs{1 + (m2_mins(1) - 5) / 0.5} = [labs{1 + (m2_mins(1) - 5) / 0.5} '/'];
    
    %set(gca, 'XTickLabel', labs, 'TickLabelInterpreter', 'latex');
    xticklabels(labs)
    
    txt = ['+ : Method 1 max' newline '* : Method 2 max' newline '- : Method 1 min' newline '/ : Method 2 min'] ;
    yl = ylim;
    text(17,yl(2) - 5,txt, 'FontSize', 10)

    xlabel(['Method 1 Maxs: ', mat2str(m1_maxs), '    Method 2 Maxs: ', mat2str(m2_maxs), '    Method 1 Mins: ', mat2str(m1_mins), '    Method 2 Mins: ', mat2str(m2_mins)], 'FontSize', 15);
    ylabel('∆Weighted SCC', 'FontSize', 15);
    axis([4.5,18.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])
    if iMethod==2
    save([saveDir subjid '_1317_medians.mat'],'weighted_vals_sorted','UniqueFreqs')
    end 
     if iMethod==1
    save([saveDir subjid '_matched_medians.mat'],'weighted_vals_sorted','UniqueFreqs')
    end 
    subplot(2,1,2)
    plot(UniqueFreqs,weighted_diag_regression,'LineWidth',2)
        title('Regression Medians Order','FontSize',20)

    axis([4.5,18.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])
    xticks(5:0.5:18)
    xticklabels(5:0.5:18)
    [s,idx] = sort(weighted_diag_regression,"descend");
    hold on
    scatter(UniqueFreqs(idx(1:3)),weighted_diag_regression(idx(1:3)),100,'filled','r')
    yline(b(end), 'LineWidth', 0.5)

    
end
saveas(gca,[saveDir subjid '_Boxplots.jpg']);
rfs = array2table(weighted_vals_sorted);
rfs=  cat(2,table(UniqueFreqs'),rfs);


% close all;
