function [weighted_vals_sorted,UniqueFreqs,tops] = get_RF(subj_path,ChanSet,window,median_or_mean,smoothed,chanlocs,output,weights)
load('~/Documents/NN_ANT.mat'); % nearest neighbors list
chancat= [];
for i=1:numel(ChanSet)
chancat = [chancat;strjoin(string(ChanSet{i}))];

end

table([chancat;'Intercept'],weights)



dirData = [char(subj_path), '/Step8_ASCOutput/'];
saveDir = [char(subj_path), '/Step9_FrequencyOutputDecision/'];

cd(dirData)

if strcmp(smoothed, 'smoothed')
    smoothed = true;
else
    smoothed = false;
end
if strcmp(window,'repeats')
    repeats = true;
    window = 0.1;

else 
    repeats = false;
end
if strcmp(output,'matched')
    matched=true;
   
else
    matched=false;
end
table(matched,repeats,smoothed)
%% get the SCC .mat file
files = dir('*.mat');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{3};
end
nPat = numel(files);

disp(patID)



%Calculate and store all ASC pre, post for all subjects

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
    AllStimFreq{iPat} = allStimFreq;
    AllChanLabels(:,iPat) = chanLabels;
    EveryStimFreq = [EveryStimFreq allStimFreq];
    clear aSCPostAll aSCPreAll allStimFreq chanLabels
end





seedChanIdx = [5];
%indices to remove
idxM1M2 = ismember(chanlocs', {'M1', 'M2', 'F3'});

%%Chan to do = 1:64 for all channels; choose channels to analyze
chanToDo = 1:64;
chanToDo(idxM1M2) = [];

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

    %% find indices of channels
    ChansToUse = {};
    ChanString = '';
    for ilistChans = 1:numel(ChanSet)
        TmpChans = [];
        for j = 1:numel(ChanSet{ilistChans})
            TmpChans(j) = find(ismember(chanlocs, ChanSet{ilistChans}{j}));
            ChanString = strcat(ChanString, ChanSet{ilistChans}{j});
            ChanString = strcat(ChanString, ',');
        end
        
        %%%FIND CHANNELS TO SMOOTH%%%
        if(smoothed==true)
            TmpChans = [TmpChans NN_ANT{TmpChans}];
            ChansToUse{ilistChans} = TmpChans;
        elseif(smoothed==false)
            ChansToUse{ilistChans} = TmpChans;
        end
        
    end
    disp(ChansToUse)
    
    
    
    %Set structure for pre and Post data
    AllPostData{iPat,1} = AnalysisPOST_NoAveraging{1};
    AllPreData{iPat,1} = AnalysisPRE_NoAveraging{1};
    
    if((exist('RemovedFreqs')) && (numel(RemovedFreqs) > 0))
        %Set nan values based on what was removed during VI
        AllPreData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;
        AllPostData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;    
    end
    %Set window for averaging. 0.1 to get like frequencies
%     if repeats==true
%     window = 0.1;
%     end
    UpdatedPre1toPost1_sorted = [];
    
    %%Run on Session 1
    UniqueFreqs = unique(AllStimFreq{1});
    UniqueVals_sorted = nan(numel(UniqueFreqs), 6,size(AllPostData{iPat,1},3));
   
    
    
    
    
    
    for iFreq = 1:numel(UniqueFreqs)
        tmp_freqs = AllStimFreq{1};
        window_avg = find(tmp_freqs >= UniqueFreqs(iFreq) - window & tmp_freqs <= UniqueFreqs(iFreq) + window);
        
        %Set frquencies of interest based on method
        if(matched)
            FreqofInterest = find(v_FreqAxis >= UniqueFreqs(iFreq) - 2 & v_FreqAxis <= UniqueFreqs(iFreq) + 2);
        else
            FreqofInterest = find(v_FreqAxis >= output(1) & v_FreqAxis <= output(2));
        end
        
        %Average over frequencies of interest
        Pre1toPost1 = mean(AllPostData{iPat,1}(FreqofInterest,:,:) - AllPreData{iPat,1}(FreqofInterest,:,:),1);
        %Set uniquevalues
        UniqueVals_sorted(iFreq,1:numel(window_avg),:) = squeeze(Pre1toPost1(1,window_avg,:));        
    end
    
    
    
        FinalValsPre1toPost1_sorted = [];
    for iSite = 1:numel(ChanSet)
        if(strcmp(median_or_mean,'median'))
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmedian(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        elseif(strcmp(median_or_mean,'mean'))
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmean(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        end
    end
    
        %Loop through and multiply regression values to matrix
    weighted_vals_sorted = ones(size(FinalValsPre1toPost1_sorted,1),size(FinalValsPre1toPost1_sorted,2))*weights(end);
     weighted_diag_regression =[0];
    for iSite = 1:numel(ChanSet)
        weighted_vals_sorted = weighted_vals_sorted + FinalValsPre1toPost1_sorted(:,:,iSite)*weights(iSite);
        weighted_diag_regression = weighted_diag_regression + weights(iSite) * nanmedian(FinalValsPre1toPost1_sorted(:,:,iSite),2);
    end
    weighted_diag_regression = weighted_diag_regression + weights(end);
%             ROI1 = nanmedian(FinalValsPre1toPost1_sorted(:,:,1),2);
%         ROI2 = nanmedian(FinalValsPre1toPost1_sorted(:,:,2),2);
        
        
        
%         weighted_diag_regression = weights(end) + weights(1) * ROI1 + weights(2) * ROI2;
            [s,idx] = sort(nanmedian(weighted_diag_regression,2),"descend");

         tops = UniqueFreqs(idx(1:5));

    
    if repeats==true
    gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1500, 800]);
    subplot(2,1,1)
    allthree = sum(~isnan(weighted_vals_sorted), 2) >=3;
    onlytwo = sum(~isnan(weighted_vals_sorted), 2) == 2;
    onlyone = sum(~isnan(weighted_vals_sorted), 2) == 1;
    boxplot(weighted_vals_sorted(allthree, :)', 'positions', UniqueFreqs(allthree), 'labels', UniqueFreqs(allthree), 'colors', 'b', 'Widths', 0.2);
    hold on
    boxplot(weighted_vals_sorted(onlytwo, :)', 'positions', UniqueFreqs(onlytwo), 'labels', UniqueFreqs(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
    scatter(UniqueFreqs(onlyone), weighted_vals_sorted(onlyone, :), 'k', 'filled');
    yline(weights(3), 'LineWidth', 0.5)
    xticks(5:0.5:18)
    xticklabels(5:0.5:18)
    % make this appropriate
    title("WEIGHTED: Post - Pre: " + strjoin(string(output)) + " , " + strjoin(chancat), 'FontSize', 18);

    
    xlabel('Stim Frequency', 'FontSize', 15);
    ylabel('∆Weighted SCC', 'FontSize', 15);
    axis([4.5,18.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])
    
    subplot(2,1,2)
        plot(UniqueFreqs,weighted_diag_regression,'LineWidth',2)
        title('Regression Medians Order','FontSize',20)

    axis([4.5,18.5,1.2*min(weighted_vals_sorted,[],'all'),1.2*max(weighted_vals_sorted,[],'all')])
    xticks(5:0.5:18)
    xticklabels(5:0.5:18)
    [s,idx] = sort(weighted_diag_regression,"descend");
    hold on
    scatter(UniqueFreqs(idx(1:3)),weighted_diag_regression(idx(1:3)),100,'filled','r')
    yline(weights(end), 'LineWidth', 0.5)
    else
%         disp('test')
        if(strcmp(median_or_mean,'median'))
            plot(UniqueFreqs,nanmedian(weighted_diag_regression,2))
        else
            plot(UniqueFreqs,nanmean(weighted_diag_regression,2))
        end
        title("WEIGHTED: Post - Pre: " + strjoin(string(output)) + " , " + strjoin(chancat), 'FontSize', 18);
        hold on
        scatter(UniqueFreqs(idx(1:5)),weighted_diag_regression(idx(1:5)),100,'filled')

        xlabel('Input Frequency','FontSize',20)
        ylabel('Connectivity Change','FontSize',20)

        
    end
   
end
