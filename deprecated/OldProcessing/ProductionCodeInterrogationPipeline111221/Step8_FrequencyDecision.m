clear
clc
% close all
figure
subjid = '18-072';
addpath '/Users/neuromodit/Documents/RF_StepsMatlab/deprecated/OldProcessing/ProductionCodeInterrogationPipeline111221'

%Load example eeg with same structure. Use this for channel locations
EEG = pop_loadset('/Users/neuromodit/Documents/Interrogation/OLD_SCC22-020-2022-06-28/Step7_VisualInspectionOutput/22-020_2022-06-28_16-19-14_Merged_EpochData.set')

dirData = ['~/Documents/Interrogation/' subjid '/Step8_ASCOutput/'];
% %Where to save output mat files
saveDir = ['~/Documents/Interrogation/' subjid '/Step9_FrequencyOutputDecision/'];

%dirData = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/SCC/Do/'];
%Where to save output mat files
%saveDir = ['/Users/andrewwilson/Documents/EEG/TMSEEG/InterrogationProtocol/MethodsAnalysis/Methods56/ButlerFiles/FrequencySelection/'];


%go to data directory and get patient IDs for all files in that directory.
cd(dirData)
files = dir('*.mat');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{3};
end
nPat = numel(files);

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
%         %%UPDATED
%          AnalysisPRE_NoAveraging{iPat} = squeeze(temp_matPre);
%          AnalysisPOST_NoAveraging{iPat} = squeeze(temp_matPost);
    
    %reduce to average ASC pre and post
    temp_matPre = nanmean(squeeze(nanmean(temp_matPre,3)),3);
    temp_matPost = nanmean(squeeze(nanmean(temp_matPost,3)),3);
    
    AnalysisPRE{iPat} = temp_matPre;
    AnalysisPOST{iPat} = temp_matPost;
    
    %%For later use; Diff ASC
    EveryDiffASC = [EveryDiffASC AnalysisPOST{iPat}-AnalysisPRE{iPat}];
    clear temp_matPost temp_matPre
end

%All
ChanSet = {{'FT8','T8','TP8','CP6','P8'},{'FC3','AF3'}};%{{'Fp1'},{'Fpz'},{'Fp2'},{'F7'},{'Fz'},{'F4'},{'F8'},{'FC5'},{'FC1'},{'FC2'},{'FC6'},{'T7'},{'C3'},{'Cz'},{'C4'},{'T8'},{'CP5'},{'CP1'},{'CP2'},{'CP6'},{'P7'},{'P3'},{'Pz'},{'P4'},{'P8'},{'POz'},{'O1'},{'O2'},{'CPz'},{'AF7'},{'AF3'},{'AF4'},{'AF8'},{'F5'},{'F1'},{'F2'},{'F6'},{'FC3'},{'FCz'},{'FC4'},{'C5'},{'C1'},{'C2'},{'C6'},{'CP3'},{'CP4'},{'P5'},{'P1'},{'P2'},{'P6'},{'PO5'},{'PO3'},{'PO4'},{'PO6'},{'FT7'},{'FT8'},{'TP7'},{'TP8'},{'PO7'},{'PO8'},{'Oz'}};
Output = {};
count_output = 1;
SCCWindow = [13 17];
for k = 1
    indices = [];
    listChans = nchoosek(ChanSet,k);
    for ilistChans = 1:size(listChans,1)
        ChansToUse = [listChans{ilistChans,:}];
        Channels = [];
        
        for j = 1:size(ChansToUse,2)
            Channels(j) = find(ismember({EEG.chanlocs(:).labels}, ChansToUse{j}));
        end
        %AllF3toCP2 = {};
        for iPat = 1:nPat
            %aSCDIFF = [];
            aSCDIFF = AnalysisPOST_NoAveraging{iPat} - AnalysisPRE_NoAveraging{iPat};
            tmp_SCCMat = mean(aSCDIFF(:,:,Channels),3);
            MeanDiags{ilistChans} = f_returnMeanDiag_OLDSCC(tmp_SCCMat,AllStimFreq{iPat}, v_FreqAxis, SCCWindow);
            %Updated
            %MeanDiags{ilistChans} = f_returnAverageSCC_Updated(tmp_SCCMat,AllStimFreq{iPat}, v_FreqAxis, SCCWindow);

        end


    end
end

MeanDiags_RightTP = MeanDiags{1};
MeanDiags_Local = MeanDiags{2};
        load('/Users/neuromodit/Documents/RF_StepsMatlab/deprecated/OldProcessing/ProductionCodeInterrogationPipeline111221/RegressionValues.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating new MeanDiags with weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
freqTopPeak = [];
rank = [];
for iPat = 1:nPat
    tmp_local = MeanDiags_Local;
    tmp_RightTP = MeanDiags_RightTP;
    
    %Create weighted spectrum
    weighted_diag = tmp_local*b(1) + tmp_RightTP*b(2) + ones(size(tmp_local,1),size(tmp_local,2))*b(3);
    %find peak
    
    [sortedVals, idx_s] = sort(weighted_diag, 'desc');
    sortedStim = AllStimFreq{iPat}(idx_s);
    
    [y,x] = findpeaks(weighted_diag, 'MinPeakHeight',-1);
    [~, idxSorted]= sort(y, 'descend');
    idxTopPeak = x(idxSorted);
    freqTopPeak(iPat) = AllStimFreq{iPat}(idxTopPeak(1));
end

%%Individual plots AverageInertia vs RankAVerage
num_items = 15;
colors = [];
top_freqs = [];
closest = [];
for iPat = 1:nPat
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
                title([patID{iPat} ' Top Ranked Frequency: ' num2str(temp_stim(idx_i(j,1))) 'Hz'], 'FontSize', 16);
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
%     saveas(gcf, [saveDir patID{iPat} '_FreqDecisionMedian.png']);
end
% close all;
