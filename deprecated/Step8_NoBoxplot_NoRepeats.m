clear

clc
close all
load('~/Documents/NN_ANT.mat'); % nearest neighbors list
indiv_subj = false
if indiv_subj
    subjid = 'SCC22-004';
    work_dir = '~/Documents/Interrogation/';

    dirData = ['~/Documents/Interrogation/' subjid '/Step8_ASCOutput/'];
    saveDir = [work_dir, subjid, '/Step9_FrequencyOutputDecision/'];
else 
     dirData = '/Users/neuromodit/Downloads/Interrogation_Cpz0/SCCOutputGabor/';
%Where to save output mat files
 saveDir = '/Users/neuromodit/Downloads/Interrogation_Cpz0/SCCOutputGabor/Out/';
end 
%Load example eeg with same structure. Use this for channel locations
EEG = pop_loadset('/Users/neuromodit/Documents/Interrogation/OLD_SCC22-020-2022-06-28/Step7_VisualInspectionOutput/22-020_2022-06-28_16-19-14_Merged_EpochData.set')
 



%dirData = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/SCC/Do/'];
%Where to save output mat files
%saveDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/FrequencySelection/'];

filename_txfreq = readtable('~/Documents/GrantJan2023/filename_txfreq.csv','Delimiter',',');
%go to data directory and get patient IDs for all files in that directory.
cd(dirData)
files = dir('*.mat');
files = {files.name}';
patID = [];
for iFile = numel(files)
    close all
        clearvars -except indiv_subj patID files dirData saveDir EEG NN_ANT iFile RF_CTZ p r tops_matched_vals tops_matched txfreq filename_txfreq tx_val treated_at_matchedRF freqAntiRF

    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{3};
txfreq(iFile) = filename_txfreq((string(files(iFile)) == string(filename_txfreq.Filename)),:).TxFreq;
iFile
nPat = 1;%numel(files);

%--------------------------------------------------------------------------

AllPREASC = cell(1,nPat);
AllPOSTASC = cell(1,nPat);
AllStimFreq = cell(1,nPat);
AllChanLabels = cell(64,nPat);

EveryStimFreq = [];
%%Load files in for the subset above.
for iPat = 1
    load(files{iFile})
    AllPREASC{iPat} = aSCPreAll;
    
    AllPOSTASC{iPat} = aSCPostAll;
    AllStimFreq{iPat} = allStimFreq;

    if exist('RemovedFreqs')
    if(numel(RemovedFreqs) > 0)
        %Set nan values based on what was removed during VI
        AllPREASC{iPat}(:,RemovedFreqs{1}',:) = NaN;
        AllPOSTASC{iPat}(:,RemovedFreqs{1}',:) = NaN;    
        
    end
    end 
    [AllStimFreq{iPat},idx] = sort(AllStimFreq{iPat});
    AllPREASC{iPat} = AllPREASC{iPat}(:,idx,:,:);
    AllPOSTASC{iPat} = AllPOSTASC{iPat}(:,idx,:,:);

    AllChanLabels(:,iPat) = chanLabels;
    EveryStimFreq = [EveryStimFreq AllStimFreq{iPat}];
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
for iPat = 1
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
%         %%UPDATED
%          AnalysisPRE_NoAveraging{iPat} = squeeze(temp_matPre);
%          AnalysisPOST_NoAveraging{iPat} = squeeze(temp_matPost);
    
    %reduce to average ASC pre and post
%     temp_matPre = nanmean(squeeze(nanmean(temp_matPre,3)),3);
%     temp_matPost = nanmean(squeeze(nanmean(temp_matPost,3)),3);
%     
%     AnalysisPRE{iPat} = temp_matPre;
%     AnalysisPOST{iPat} = temp_matPost;
    
    %%For later use; Diff ASC
%     EveryDiffASC = [EveryDiffASC AnalysisPOST{iPat}-AnalysisPRE{iPat}];
    clear temp_matPost temp_matPre


%All

for iMethod = 1
    if(iMethod == 1)
        ChanSet = {{'Fpz'},{'AF4'}};
    elseif(iMethod==2)
        ChanSet = {{'FC3','AF3'},{'FT8','T8','TP8','CP6','P8'}};
    elseif(iMethod==3)
        ChanSet = {{'C3'}};
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
    
    %% CORRECTING FOR MAGSTIM ERROR
%     AllStimFreq{1}(1) = AllStimFreq{1}(2);
    
    %Set structure for pre and Post data
    AllPostData{iPat,1} = AnalysisPOST_NoAveraging{iPat};
    AllPreData{iPat,1} = AnalysisPRE_NoAveraging{iPat};
    
%     if(numel(RemovedFreqs) > 0)
%         %Set nan values based on what was removed during VI
%         AllPreData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;
%         AllPostData{iPat,1}(:,RemovedFreqs{1}',:) = NaN;    
%     end
    %Set window for averaging. 0.1 to get like frequencies
    if iMethod==1
        window=0.3; %.5 forhalf hz
    elseif iMethod==2
        window=0.5;
    elseif iMethod==3
        window = 0.5;
    end
%     window = 0.5;
    UpdatedPre1toPost1_sorted = [];
    
    %%Run on Session 1
    UniqueFreqs = unique(AllStimFreq{iPat});
    UniqueVals_sorted = nan(numel(AllStimFreq{iPat}), 11,size(AllPostData{iPat,1},3));
   
    %%%%%Run for all
    for iFreq = 1:numel(AllStimFreq{iPat})
        tmp_freqs = AllStimFreq{iPat};
        window_avg = find(tmp_freqs >= tmp_freqs(iFreq) - window & tmp_freqs <= tmp_freqs(iFreq) + window);
        
        %Set frquencies of interest based on method
        if(iMethod == 1 )
            FreqofInterest = find(v_FreqAxis >= tmp_freqs(iFreq) - 2 & v_FreqAxis <= tmp_freqs(iFreq) + 2);
        elseif(iMethod==2)
            FreqofInterest = find(v_FreqAxis >= 13 & v_FreqAxis <= 17);
        elseif(iMethod==3)
            FreqofInterest = find(v_FreqAxis >= tmp_freqs(iFreq) - .5 & v_FreqAxis <= tmp_freqs(iFreq) + .5);

        end
        
        %Average over frequencies of interest
        Pre1toPost1 = mean(AllPostData{iPat,1}(FreqofInterest,:,:) - AllPreData{iPat,1}(FreqofInterest,:,:),1);
        %Set uniquevalues
        UniqueVals_sorted(iFreq,1:numel(window_avg),:) = squeeze(Pre1toPost1(1,window_avg,:));        
    end
    
    
    %Average over channels based on method
    FinalValsPre1toPost1_sorted = [];
    for iSite = 1:numel(ChanSet)
        if(iMethod == 1)
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmedian(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        elseif(iMethod ==2 || iMethod==3)
            FinalValsPre1toPost1_sorted(:,:,iSite) = nanmean(UniqueVals_sorted(:,:,ChansToUse{iSite}),3);
        
        end
    end
    
    %Load b values from regressions
    if(iMethod == 1)
        load('/Users/neuromodit/Documents/RF_StepsMatlab/RegressionMatched_01-17-2023.mat')
        ROI1 = nanmedian(FinalValsPre1toPost1_sorted(:,:,1),2);
        ROI2 = nanmedian(FinalValsPre1toPost1_sorted(:,:,2),2);

    elseif(iMethod ==2)
        
        load('/Users/neuromodit/Documents/RF_StepsMatlab/RegressionValues_old.mat')
        ROI1 = nanmean(FinalValsPre1toPost1_sorted(:,:,1),2);
        ROI2 = nanmean(FinalValsPre1toPost1_sorted(:,:,2),2);
        
        
    elseif(iMethod ==3)
        b = [-30,0,0];
        ROI1 = nanmean(FinalValsPre1toPost1_sorted(:,:,1),2);
        ROI2 = zeros(size(ROI1));
%         ROI2 = nanmean(FinalValsPre1toPost1_sorted(:,:,2),2);
    end
    %Loop through and multiply regression values to matrix
%     weighted_vals_sorted = ones(size(FinalValsPre1toPost1_sorted,1),size(FinalValsPre1toPost1_sorted,2))*b(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating new MeanDiags with weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
freqTopPeak = [];
rank = [];
for iPat = 1
%     tmp_local = MeanDiags_Local;
%     tmp_RightTP = MeanDiags_RightTP;
    weighted_vals_sorted = ones(size(FinalValsPre1toPost1_sorted,1),size(FinalValsPre1toPost1_sorted,2))*b(end);
    for iSite = 1:numel(ChanSet)
        weighted_vals_sorted = weighted_vals_sorted + FinalValsPre1toPost1_sorted(:,:,iSite)*b(iSite);
    end
        if iMethod==1
    weighted_diag = nanmedian(weighted_vals_sorted,2);
    elseif iMethod==2
            weighted_diag = nanmean(weighted_vals_sorted,2);

    end 
    %Create weighted spectrum
%     weighted_diag = b(end) + b(1) * ROI1 + b(2) * ROI2;

    %find peak
    
    [sortedVals, idx_s] = sort(weighted_diag, 'desc');
    sortedStim = AllStimFreq{iPat}(idx_s);
    if iMethod==2
    [y,x] = findpeaks(weighted_diag, 'MinPeakHeight',-1);
    [~, idxSorted]= sort(y, 'descend');
    idxTopPeak = x(idxSorted);
    freqTopPeak(iPat) = AllStimFreq{iPat}(idxTopPeak(1));
    elseif iMethod==1 || iMethod==3
         freqTopPeak(iPat) = sortedStim(1);
         freqAntiRF(iFile) = sortedStim(end);

         idxTopPeak = idx_s(1);
    end
    
end


if iMethod==2
%%Individual plots AverageInertia vs RankAVerage
num_items = 15;
colors = [];
top_freqs = [];
closest = [];
for iPat = 1
    hFig = figure('units', 'normalized', 'position', [0 0 1 0.5], 'color', 'w');
    
    
    %find peak
    [y,x] = findpeaks(weighted_diag, 'MinPeakHeight',-1);
    [~, idxSorted]= sort(y, 'descend');
    idxTopPeak = x(idxSorted);
    
    temp_stim = sortedStim(1,1:num_items);
    
    %%KMeans Clustering 1D
    idx_k = {};
    C_k = {};
    average_inertia = {};
    av = [];
    average_rank = {};
    allowed_dist = 0.51;
    for iCluster = 1
        rng('default')
        rng(1); % For reproducibility
        [idx_k{iCluster},C_k{iCluster}] = kmeans(temp_stim',iCluster);
        temp_iCluster = NaN;
        for clust = 1:iCluster
            idx_i = find(idx_k{iCluster} == clust);
            
            freqs = temp_stim(idx_i);
            [sortedFreqs,s] = sort(freqs);
            distances = abs(sortedFreqs(1,1:end-1) - sortedFreqs(1,2:end));
            
            if(numel(freqs) > 3 && numel(find(distances > allowed_dist)))
                idx_k{iCluster}(idx_i) = NaN;
                new_freqs = [];
                flag = 0;
                count_1 = 1;
                count_2 = 1;
                for iFreq = 1:numel(sortedFreqs)-1
                    if(sortedFreqs(iFreq+1) - sortedFreqs(iFreq) < allowed_dist)
                        new_freqs{count_1}(1,count_2) = sortedFreqs(iFreq);
                        new_freqs{count_1}(2,count_2) = s(iFreq);
                        new_freqs{count_1}(1,count_2+1) = sortedFreqs(iFreq+1);
                        new_freqs{count_1}(2,count_2+1) = s(iFreq+1);
                        count_2 = count_2 + 1;
                    else
                        if(numel(new_freqs) ~= 0)
                            count_1 = count_1 + 1;
                            count_2 = 1;
                        end
                            
                    end
                end
                if(isnan(max(idx_k{iCluster})))
                    temp_iCluster = 0;
                else
                    temp_iCluster = max(idx_k{iCluster});
                end
                if(numel(new_freqs) == 1)
                    stop = 1;
                    for iClust = 1:stop
                        temp_iCluster = temp_iCluster + 1;
                        for jClust = 1:size(new_freqs{iClust},2)
                            idx_k{iCluster}(idx_i(new_freqs{iClust}(2,jClust))) = clust;
                            disp('');
                        end
                    end
                else
                    stop = numel(new_freqs)
                    for iClust = 1:stop
                        temp_iCluster = temp_iCluster + 1;
                        for jClust = 1:size(new_freqs{iClust},2)
                            idx_k{iCluster}(idx_i(new_freqs{iClust}(2,jClust))) = temp_iCluster;
                        end
                    end
                end
            end
        end
        
        %%Using Inertia for color KMeans
        inertia = zeros(max(idx_k{iCluster}),4);
        clusterToKeep{iCluster} = [];
        
        steps_with_nan = unique(idx_k{iCluster});
        steps = steps_with_nan(~isnan(steps_with_nan));
        for i=1:numel(steps)
            idx_i = find(idx_k{iCluster} == steps(i));
            C_k{iCluster}(steps(i)) = median(temp_stim(idx_i)); %mean(temp_stim(idx_i));%median(temp_stim(idx_i));%
            freqs = temp_stim(idx_i);
            
            [sortedFreqs,s] = sort(freqs);
            distances = abs(sortedFreqs(1,1:end-1) - sortedFreqs(1,2:end));
        
            %If only one value, then remove
            if(numel(freqs) < 3)
                disp(['REMOVE CLUSTER ' num2str(i)]);
                inertia(steps(i),1) = NaN;
                inertia(steps(i),2) = NaN;
                inertia(steps(i),3) = NaN;
                C_k{iCluster}(steps(i)) = NaN;
            elseif(min(distances) > allowed_dist)
                disp(['REMOVE CLUSTER ' num2str(i)]);
                inertia(steps(i),1) = NaN;
                inertia(steps(i),2) = NaN;
                inertia(steps(i),3) = NaN;
                C_k{iCluster}(steps(i)) = NaN;
            else
                clusterToKeep{iCluster}(1,size(clusterToKeep{iCluster},2)+1) = i;
                clusterToKeep{iCluster}(2,size(clusterToKeep{iCluster},2)) = steps(i);
                for j = 1:numel(idx_i)
                    %disp(temp_stim(idx_i(j)));
                    inertia(steps(i),1) = inertia(steps(i),1) + abs(temp_stim(idx_i(j)) - C_k{iCluster}(steps(i)));
                    inertia(steps(i),3) = inertia(steps(i),3) + idx_i(j);
                end
                inertia(steps(i),2) = numel(idx_i);
            end
        end
        
        inertia(clusterToKeep{iCluster}(1,:),4) = inertia(clusterToKeep{iCluster}(1,:),1)./inertia(clusterToKeep{iCluster}(1,:),2);
        
        av(iCluster,1) = mean(inertia(clusterToKeep{iCluster}(1,:),1)./inertia(clusterToKeep{iCluster}(1,:),2));
        average_inertia{iCluster} = inertia(:,1)./inertia(:,2);
        average_rank{iCluster} = inertia(:,3)./inertia(:,2);
        
    end
    idx_nonzero = find(av(:,1) ~= 0);
    NClusters = 1;%min(find(av(:) == min(av(idx_nonzero))));
    average_inertia = average_inertia{NClusters};
    average_rank = average_rank{NClusters};
    idx_k = idx_k{NClusters};
    idx_keep_k = find(sum(idx_k == clusterToKeep{NClusters}(2,:),2));
    idx_k = idx_k(idx_keep_k);
    Clusters = unique(idx_k);   
    
    temp_stim = temp_stim(idx_keep_k');
    
    C_k = C_k{NClusters};
    
    %%Plot 1: Plotting spectrum
    s1 = subplot(1, 3, 1);
    %title(['Peak Conn.: ' num2str(freqTopPeak(iPat))], 'FontSize', 16)
    plot(AllStimFreq{iPat},weighted_diag, 'linewidth', 2);
    hold on;
    plot(freqTopPeak(iPat), weighted_diag(idxTopPeak(1)), '.', 'Markersize', 30);
    ylim([-30 70]);
    
    %%Plot 2: Plotting Dots
    subplot(1, 3, 2);
    hold on;
    [idx_c, colorsort] = sort(average_inertia, 'asc');
    cmap = parula(max(idx_k));
    cmap_items = jet(num_items);
    colormap jet
    legend_text = {};
    
    %%PLOT KMEANS CLUSTERS
    colors = [];
    for i=1:numel(Clusters)
        %scatter(iPat,C(i), 100);
        idx_i = find(idx_k == Clusters(i));
        colorspot = find(Clusters(i) == colorsort);
        
        colors(i,:) = cmap(colorspot,:);
        h(i) = scatter(ones(1,numel(idx_i)), temp_stim(idx_i),80, 'filled', 'MarkerFaceColor', cmap(colorspot,:));
        %scatter(iPat, temp_stim(idx_i(j,1)),40, 'filled', 'MarkerFaceColor', cmap(idx_i(j,1),:));
        legend_text{i} = [num2str(C_k(Clusters(i))) ' N=' num2str(numel(idx_i))];
        
        for j = 1:numel(idx_i)
            if(idx_i(j,1) == 1)
                title([patID{iFile} ' Top Ranked Frequency: ' num2str(temp_stim(idx_i(j,1))) 'Hz'], 'FontSize', 16);
            end
            scatter(1.1, temp_stim(idx_i(j,1)),50, 'filled', 'MarkerFaceColor', cmap_items(idx_i(j,1),:));
        end
    end
    %Plot all original
    tmp_stim_plot = sortedStim(1,1:num_items);
    scatter(1.05*ones(1,numel(tmp_stim_plot)), tmp_stim_plot,50, 'x');
    xticks([1 1.05 1.1])
    xticklabels({'Final Clust.','Top 15 Freqs','Rank'})
    xtickangle(45)
    
    %xlabel('Subject', 'FontSize', 16);
    ylabel('Frequency', 'FontSize', 16);
    xlim([0.9 1.2]);
    for i=1:numel(Clusters)
        scatter(1,C_k(Clusters(i)), 100, 'MarkerEdgeColor', 'k');
    end
    legend(h(1:numel(Clusters)),legend_text);
    ylim([0 20]);
    colorbar;
       
    %%Plot 3: Average Inertia vs Rank
    subplot(1, 3, 3);
    hold on;

    min_dtz = 1000000;
    min_idx = 0;
    for i = 1:numel(Clusters)        
        plot([numel(idx_k) - inertia(Clusters(i),2)],average_rank(Clusters(i)), '.','MarkerSize',50, 'color',colors(i,:) );%,100, 'filled', 'MarkerFaceColor', colors(i,:));

        grid on;
        p12 = average_rank(Clusters(i));
        p13 = num_items - inertia(Clusters(i),2);
        p22 = 0;
        p23 = 0;
        if(sqrt((p22-p12)^2 + (p23-p13)^2) < min_dtz)
            min_idx = Clusters(i);
            min_dtz = sqrt( (p22-p12)^2 + (p23-p13)^2);
        end
    end
%     tops_1317{iFile} = C_k(Clusters);
    xlabel('Total Items - Number in Cluster', 'FontSize', 16)
    ylabel('Av. Rank', 'FontSize', 16)
    legend(legend_text);
    title(['Closest to Zero = ' num2str(C_k(min_idx))], 'FontSize', 16);
    
    %Plot Vertical line on subplot 1 with "closest value"
    axes(s1)
    title(['Peak Conn.: ' num2str(freqTopPeak(iPat)) 'Hz'], 'FontSize', 16)
    plot([C_k(min_idx) C_k(min_idx)],[min(weighted_diag) max(weighted_diag)], 'linewidth', 1.5, 'color', 'r');
    xlabel('Stimulation Frequency', 'FontSize', 16);
    ylabel('Weighted Conn. Value', 'FontSize', 16);
    RF_CTZ(iFile,iMethod) = C_k(min_idx);

%     close all
     plot(AllStimFreq{iPat},weighted_diag, 'linewidth', 2);
    hold on;
    plot(freqTopPeak(iPat), weighted_diag(idxTopPeak(1)), '.', 'Markersize', 30);
    ylim([-30 70]);       
        title(['13-17 Weighted, RF: ' char(string(C_k(min_idx)))])

    saveas(gcf, [saveDir patID{iFile} char(iMethod) '_1317_FreqDecisionMedian.jpg']);

end

elseif iMethod==1 || iMethod ==3 
        %%Plot 1: Plotting spectrum
%     s1 = subplot(1, 3, 1);
figure    
%title(['Peak Conn.: ' num2str(freqTopPeak(iPat))], 'FontSize', 16)
    plot(AllStimFreq{iPat},weighted_diag, 'linewidth', 2);
    hold on;
    plot(freqTopPeak(iPat), weighted_diag(idxTopPeak(1)), '.', 'Markersize', 30);
    ylim([-30 70]);
    RF_CTZ(iFile,iMethod) = freqTopPeak(iPat);
        rf_f = (find(weighted_diag==max(weighted_diag)));
        anti_rf_f = (find(weighted_diag==min(weighted_diag)));

%     hold on
tmp_top_stims = AllStimFreq{1}(idx_s(1));
istim = 1;
top_idx = idx_s(1);
while numel(tmp_top_stims)<3
    istim = istim+1;
    if all(abs(AllStimFreq{1}(idx_s(istim))-tmp_top_stims)>.5)
        tmp_top_stims = [tmp_top_stims,AllStimFreq{1}(idx_s(istim))];
        top_idx = [top_idx,idx_s(istim)];
    end
end
    
tops_matched{iFile} = tmp_top_stims;
tops_matched_vals{iFile} = weighted_diag(top_idx)';

tx_val{iFile} = nanmedian(weighted_diag(find(abs(AllStimFreq{1}-txfreq(iFile)) == min(abs(AllStimFreq{1}-txfreq(iFile))))));
    plot(AllStimFreq{1},weighted_diag,'LineWidth',2)
    hold on
    scatter(AllStimFreq{1}(rf_f(1)),weighted_diag(rf_f(1)),'filled','r','SizeData',150)

    title(['RF: ' char(string(AllStimFreq{1}(rf_f(1)))) ' , AntiRF: ' char(string(AllStimFreq{1}(anti_rf_f(1))))])
   
    saveas(gcf, [saveDir patID{iFile} char(iMethod) '_basicmatched_FreqDecisionMedian.jpg']);
    matched = weighted_diag;
    if any(abs(txfreq(iFile)-tops_matched{iFile})<=.5)
        treated_at_matchedRF(iFile) = 1;
    else
        treated_at_matchedRF(iFile) = 0;
    end
% close all

end 


end
[rtmp,ptmp] = corrcoef(matched,weighted_diag);
p(iFile) = ptmp(1,2); r(iFile) = rtmp(1,2);
% close all
end
end
% 
if ~indiv_subj
    out_table = table(files,RF_CTZ(:,1),freqAntiRF',treated_at_matchedRF');
    out_table.Properties.VariableNames = {'Filename','Matched RF','Anti-RF, Matched','TreatedAtRF'};
    writetable(out_table,'~/Documents/GrantJan2023/RF_CTZ_withAnti_half.csv');
%     close all;

    conn_table = table(files,treated_at_matchedRF',cell2mat(tx_val)');
    conn_table(conn_table.Var2==1,:).Var3
    conn_table.Properties.VariableNames={'files','treated_at_matchedRF','connectivity'};
    writetable(conn_table,'~/Documents/GrantJan2023/connectivity_rf_vs_nonrf1_half.csv')
end
% colors = colororder;
% color(1,:) = colors(2,:);
% color(2,:) = colors(1,:);
% clf
% h1 = boxplot([conn_table(conn_table.treated_at_matchedRF==1,:).connectivity;conn_table(conn_table.treated_at_matchedRF==0,:).connectivity],[repmat(0,height(conn_table(conn_table.treated_at_matchedRF==1,:)),1);repmat(1,height(conn_table(conn_table.treated_at_matchedRF==0,:)),1)],'Labels',{'RF','Non RF'})
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.8);
% end

