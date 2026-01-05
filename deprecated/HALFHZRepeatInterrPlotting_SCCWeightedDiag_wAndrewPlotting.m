clear
clc
close all

subjid = 'SCC22-020-2022-06-28 GABOR';
work_dir = '~/Documents/Interrogation/';
%Load example eeg with same structure. Use this for channel locations
EEG = pop_loadset('/Users/neuromodit/Documents/Interrogation/SCC22-020-2022-06-28/Step7_VisualInspectionOutput/22-020_2022-06-28_16-19-14_Merged_EpochData.set')
dirData = [work_dir, subjid, '/Step8_ASCOutput/'];
% %Where to save output mat files
saveDir = [work_dir, subjid, '/Step9_FrequencyOutputDecision/'];

% %Where to save output mat files%
%saveDir = [work_dir, subjid, '/Step9_FrequencyOutputDecision/'];

%go to data directory and get patient IDs for all files in that directory.
cd(dirData)
%load('rem_train_nums.mat')

files = dir('Entrain*.mat');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = tmpsplit{3};
end
nPat = numel(files);
pptx = exportToPPTX();
%--------------------------------------------------------------------------

AllPREASC = cell(1, nPat);
AllPOSTASC = cell(1, nPat);
AllStimFreq = cell(1, nPat);
AllChanLabels = cell(64, nPat);

EveryStimFreq = [];
%%Load files in for the subset above.
for iPat = 1:nPat
    load(files{iPat})
    AllPREASC{iPat} = aSCPreAll;
    AllPOSTASC{iPat} = aSCPostAll;
    AllStimFreq{iPat} = allStimFreq;
    AllChanLabels(:, iPat) = chanLabels;
    EveryStimFreq = [EveryStimFreq, allStimFreq];
    clear aSCPostAll aSCPreAll allStimFreq chanLabels
end
AllStimFreq = {round(AllStimFreq{1}*2) / 2};
%just for this subj

%--------------------------------------------------------------------------

%% Choosing Channels for analysis

%define variables
AnalysisPRE = cell(1, nPat);
AnalysisPOST = cell(1, nPat);
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
            if (iChan2 < iChan1)
                temp_matPre(:, :, chan1_counter, iChan2) = AllPREASC{iPat}(:, :, iChan2, iChan1);
                temp_matPost(:, :, chan1_counter, iChan2) = AllPOSTASC{iPat}(:, :, iChan2, iChan1);
            else
                temp_matPre(:, :, chan1_counter, iChan2) = AllPREASC{iPat}(:, :, iChan1, iChan2);
                temp_matPost(:, :, chan1_counter, iChan2) = AllPOSTASC{iPat}(:, :, iChan1, iChan2);
            end
        end
        chan1_counter = chan1_counter + 1;
    end

    %reduce dimensions and average for multiple seed channels but keep 3
    %dimensiodfsdns
    AnalysisPRE_NoAveraging{iPat} = squeeze(nanmean(temp_matPre, 3));
    AnalysisPOST_NoAveraging{iPat} = squeeze(nanmean(temp_matPost, 3));

    %reduce to average ASC pre and post
    temp_matPre = nanmean(squeeze(nanmean(temp_matPre, 3)), 3);
    temp_matPost = nanmean(squeeze(nanmean(temp_matPost, 3)), 3);

    AnalysisPRE{iPat} = temp_matPre;
    AnalysisPOST{iPat} = temp_matPost;

    %%For later use; Diff ASC
    %iEveryDiffASC = [EveryDiffASC, AnalysisPOST{iPat} - AnalysisPRE{iPat}];
    clear temp_matPost temp_matPre
end

%All
ChanSet = {{'FT8', 'T8', 'TP8', 'CP6', 'P8'}, {'FC3', 'AF3'}}; %{{'Fp1'},{'Fpz'},{'Fp2'},{'F7'},{'Fz'},{'F4'},{'F8'},{'FC5'},{'FC1'},{'FC2'},{'FC6'},{'T7'},{'C3'},{'Cz'},{'C4'},{'T8'},{'CP5'},{'CP1'},{'CP2'},{'CP6'},{'P7'},{'P3'},{'Pz'},{'P4'},{'P8'},{'POz'},{'O1'},{'O2'},{'CPz'},{'AF7'},{'AF3'},{'AF4'},{'AF8'},{'F5'},{'F1'},{'F2'},{'F6'},{'FC3'},{'FCz'},{'FC4'},{'C5'},{'C1'},{'C2'},{'C6'},{'CP3'},{'CP4'},{'P5'},{'P1'},{'P2'},{'P6'},{'PO5'},{'PO3'},{'PO4'},{'PO6'},{'FT7'},{'FT8'},{'TP7'},{'TP8'},{'PO7'},{'PO8'},{'Oz'}};
for iCompSet = 3
    Output = {};
    count_output = 1;
    SCCWindow = [13, 17];
    SCCMat = [];
    for k = 1
        indices = [];
        listChans = nchoosek(ChanSet, k);
        for ilistChans = 1:size(listChans, 1)
            ChansToUse = [listChans{ilistChans, :}];
            Channels = [];

            for j = 1:size(ChansToUse, 2)
                Channels(j) = find(ismember({EEG.chanlocs(:).labels}, ChansToUse{j}));
            end
            %AllF3toCP2 = {};
            for iPat = 1:nPat
                %aSCDIFF = [];
                if (iCompSet == 1)
                    aSCDIFF = AnalysisPRE_NoAveraging{iPat};
                    CompText = 'Pre';
                elseif (iCompSet == 2)
                    aSCDIFF = AnalysisPOST_NoAveraging{iPat}; %subj, trains, chans
                    CompText = 'Post';
                else
                    aSCDIFF = AnalysisPOST_NoAveraging{iPat} - AnalysisPRE_NoAveraging{iPat};
                    CompText = 'Post - Pre';
                end
                tmp_SCCMat = mean(aSCDIFF(:, :, Channels), 3);
                freq_window = find(v_FreqAxis >= SCCWindow(1) & v_FreqAxis <= SCCWindow(2));
                SCCMat(:, :, ilistChans, iPat) = mean(tmp_SCCMat(freq_window, :), 1);

                Pre = AnalysisPRE_NoAveraging{iPat};
                Pre_c(ilistChans, :, :) = mean(Pre(:, :, Channels), 3);
                Pre_f(ilistChans, :, :) = mean(Pre_c(ilistChans, freq_window, :), 2);
                Post = AnalysisPOST_NoAveraging{iPat};
                Post_c(ilistChans, :, :) = mean(Post(:, :, Channels), 3);
                Post_f(ilistChans, :, :) = mean(Post_c(ilistChans, freq_window, :), 2);
                %MeanDiags{iPat,ilistChans} = f_returnMeanDiagSCC(tmp_SCCMat,AllStimFreq{iPat}, v_FreqAxis, SCCWindow);
            end

        end


        %MeanDiags_RightTP = MeanDiags{1};
        %MeanDiags_Local = MeanDiags{2};
        load(['~/Downloads/ProductionCodeInterrogationPipeline111221/', 'RegressionValues.mat'])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Creating new MeanDiags with weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %figure;
        freqTopPeak = [];
        rank = [];
        for iPat = 1:nPat
            %tmp_local = MeanDiags{iPat,2};
            %tmp_RightTP = MeanDiags{iPat,1};


            tmp_localSCC = SCCMat(:, :, 2, iPat);
            tmp_RightTPSCC = SCCMat(:, :, 1, iPat);
            Pre_local = Pre_f(2, :);
            Post_local = Post_f(2, :);
            Pre_tp = Pre_f(1, :);
            Post_tp = Post_f(1, :);
            Pre_weighted = b(1) * Pre_local + b(2) * Pre_tp + b(3);
            Post_weighted = b(1) * Post_local + b(2) * Post_tp + b(3);
            %Create weighted spectrum
            weighted_valsSCC = tmp_localSCC * b(1) + tmp_RightTPSCC * b(2) + ones(size(tmp_localSCC, 1), size(tmp_localSCC, 2)) * b(3);


            %Create weighted spectrum
            %weighted_diag = tmp_local*b(1) + tmp_RightTP*b(2) + ones(size(tmp_local,1),size(tmp_local,2))*b(3);
            %find peak

            [sortedVals, idx_s] = sort(weighted_valsSCC, 'desc');
            sortedStim = AllStimFreq{iPat}(idx_s);

            [y, x] = findpeaks(weighted_valsSCC, 'MinPeakHeight', -1);
            [~, idxSorted] = sort(y, 'descend');
            idxTopPeak = x(idxSorted);
            freqTopPeak(iPat) = AllStimFreq{iPat}(idxTopPeak(1));
            weighted_medians = [];
            tmp_stdev = [];
            f_list = [4:0.5:18];

            weighted_vals = NaN(length(f_list), 3);
            local_vals = NaN(length(f_list), 3);
            RTP_vals = NaN(length(f_list), 3);
            Diff_median_pre = NaN(length(f_list), 3);
            Diff_first_pre = NaN(length(f_list), 3);
            Pre_weighted_all = NaN(length(f_list), 3);
            Post_weighted_all = NaN(length(f_list), 3);

            %% all values
            for i = 1:length(f_list)
                idx_i = find(AllStimFreq{1} == f_list(i));
                tmp_weighted = weighted_valsSCC(1, idx_i);
                tmp_local = tmp_localSCC(1, idx_i) * b(1);
                tmp_RTP = tmp_RightTPSCC(1, idx_i) * b(2);
                tmp_pre = Pre_weighted(1, idx_i);
                tmp_post = Post_weighted(1, idx_i);

                weighted_vals(i, 1:length(tmp_weighted)) = tmp_weighted;
                local_vals(i, 1:length(tmp_weighted)) = tmp_local;
                RTP_vals(i, 1:length(tmp_weighted)) = tmp_RTP;

                Pre_weighted_all(i, 1:length(tmp_weighted)) = tmp_pre;
                Post_weighted_all(i, 1:length(tmp_weighted)) = tmp_post;

                post_stdev(i) = nanstd(tmp_post);
                tmp_stdev(i) = nanstd(tmp_weighted);

            end

            %% now exclude the vi ones before calculating things
            rem_train_nums = logical(rem_train_nums);
            weighted_vals(rem_train_nums) = NaN;
            local_vals(rem_train_nums) = NaN;
            RTP_vals(rem_train_nums) = NaN;
            Pre_weighted_all(rem_train_nums) = NaN;
            Post_weighted_all(rem_train_nums) = NaN;
            tmp_stdev(sum(rem_train_nums, 2) > 1) = NaN;
            weighted_medians = nanmedian(weighted_vals, 2);
            local_medians = nanmedian(local_vals, 2);
            RTP_medians = nanmedian(RTP_vals, 2);


            %[p,tbl,stats] = anova1(tmp_mat,[],'off');
            %Fs(i) = tbl{2,5};
            %Ps(i) = tbl{2,6};


            first_pre = Pre_weighted_all(:, 1);


            median_pre = nanmedian(Pre_weighted_all, 2);
            %first_post(i) = Post_weighted(idx_i(1));
            %median_post(i) = median(Post_weighted(idx_i));

            Diff_median_pre = Post_weighted_all - median_pre;
            Diff_first_pre = Post_weighted_all - first_pre;

            %%


%             figure;
%             subplot(1, 2, 2)
%             hold on;
%             data1 = Pre_weighted_all(:, 2) - Pre_weighted_all(:, 1);
%             data2 = Pre_weighted_all(:, 3) - Pre_weighted_all(:, 2);
%             scatter(data1, data2, 80, 'filled');
%             lsline;
%             idx_notnan = find(~isnan(data1) & ~isnan(data2));
%             [r, p] = corr(data1(idx_notnan), data2(idx_notnan));
%             title(['r=', num2str(r), ',p=', num2str(p)], 'FontSize', 16);
%             xlabel('Second - First Pre', 'FontSize', 16)
%             ylabel('Third  - Second Pre', 'FontSize', 16)
% 
% 
%             subplot(1, 2, 1)
%             data = Pre_weighted_all;
%             %load carsmall MPG
%             %MPG(:,2)=MPG(:,1).*2;
%             %MPG(:,3)=MPG(:,1).*3;
%             boxplot(data);
%             hold on;
%             x = repmat(1:3, length(data), 1);
%             scatter(x(:), data(:), 'filled', 'MarkerFaceAlpha', 0.6', 'jitter', 'on', 'jitterAmount', 0.15);
%             title([patID{iPat}, ' PRE Data'], 'FontSize', 16)

            close all;
% 
%             %% Using median pre: Weighted Conn
%             allthree = sum(~isnan(Diff_median_pre), 2) == 3;
%             onlytwo = sum(~isnan(Diff_median_pre), 2) == 2;
%             onlyone = sum(~isnan(Diff_median_pre), 2) == 1;
%             gcf5 = figure('Renderer', 'painters', 'Position', [5, 5, 1300, 500]);
%             boxplot(Diff_median_pre(allthree, :)'+b(3), 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
%             hold on
%             boxplot(Diff_median_pre(onlytwo, :)'+b(3), 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
%             scatter(f_list(onlyone), Diff_median_pre(onlyone, :), 'k', 'filled');
%             yline(0, 'LineWidth', 0.5)
%             xticks(4:0.5:18)
%             xticklabels(4:0.5:18)
%             title('WEIGHTED: Post - Pre; (Using Median Pre)', 'FontSize', 18);
%             xlabel('Stim Frequency', 'FontSize', 15);
%             ylabel('∆Weighted SCC', 'FontSize', 15);
% 
%             %% using first pre : weighted conn
%             allthree = sum(~isnan(Diff_first_pre), 2) == 3;
%             onlytwo = sum(~isnan(Diff_first_pre), 2) == 2;
%             onlyone = sum(~isnan(Diff_first_pre), 2) == 1;
%             gcf6 = figure('Renderer', 'painters', 'Position', [5, 5, 1300, 500]);
%             boxplot(Diff_first_pre(allthree, :)', 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
%             hold on
%             boxplot(Diff_first_pre(onlytwo, :)', 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
%             scatter(f_list(onlyone), Diff_first_pre(onlyone, :), 'k', 'filled');
%             yline(0, 'LineWidth', 0.5)
%             xticks(4:0.5:18)
%             xticklabels(4:0.5:18)
%             title('WEIGHTED: Post - Pre; (Using First Pre)', 'FontSize', 18);
%             xlabel('Stim Frequency', 'FontSize', 15);
%             ylabel('∆Weighted SCC', 'FontSize', 15);

            %%


            f_list(isnan(weighted_medians)) = [];
            local_vals(isnan(weighted_medians), :) = [];
            RTP_vals(isnan(weighted_medians), :) = [];
            weighted_vals(isnan(weighted_medians), :) = [];
            Diff_median_pre(isnan(weighted_medians), :) = [];

            RTP_medians(isnan(weighted_medians)) = [];
            local_medians(isnan(weighted_medians)) = [];
            weighted_medians(isnan(weighted_medians)) = [];


            if (iCompSet == 3)
                allFigs{iPat} = figure('units', 'normalized', 'position', [0, 0, 1, 1], 'color', 'w');
            else
                set(0, 'CurrentFigure', allFigs{iPat})
            end


            %find peak
            [y, x] = findpeaks(weighted_medians, 'MinPeakHeight', -1);
            [~, idxSorted] = sort(y, 'descend');
            idxTopPeak = x(idxSorted);
            %make sure to find one with multiple points..
            peaklist_rank = 1;
            while sum(isnan(weighted_vals(idxTopPeak(peaklist_rank), :))) > 1
                peaklist_rank = peaklist_rank + 1;
            end
            top_peak(iPat) = f_list(idxTopPeak(peaklist_rank));

            %%Plot 1: Plotting spectrum
            %title(['Peak Conn.: ' num2str(freqTopPeak(iPat))], 'FontSize', 16)
            %             subplot(3, 1, iCompSet)


            plot(f_list, weighted_medians(1:end), 'linewidth', 2);
            hold on;
            if iCompSet == 3
                xline(top_peak);

            end
            %plot(top_peak(iPat), plot_medians(idxTopPeak(1)), '.', 'Markersize', 30);

            %axis([2,18,-50,50])
            %Fs = [];
            %             for i = 1:length(f_list)
            %                 idx_i = find(AllStimFreq{1} == f_list(i));
            %                 tmp_weighted = weighted_valsSCC(1,idx_i);
            %                 if sum(~isnan(tmp_weighted))>1
            %                 std_f(:,i) = nanstd(tmp_weighted);
            %                 else
            %                     std_f(:,i) = NaN;
            %                 end
            %                 %scatter(f_list(i)*ones(numel(tmp_mat)), tmp_mat, 'filled');
            %                 %[p,tbl,stats] = anova1(tmp_mat,[],'off');
            %                 %Fs(i) = tbl{2,5};
            %                 %Ps(i) = tbl{2,6};
            %             end
            %red is #1, y is #2, #3 is blue
            %             scatter(f_list, weighted_vals, [], [1, 0, 0; .9, .8, 0; 0, 0, 1], 'filled');

            title([patID{iPat}, ' ', CompText]);
            if (iCompSet == 3)

                %                 gcf2 = figure('Renderer', 'painters', 'Position', [5, 5, 1300, 500]);
                %                 hold on;
                %                 plot(f_list, weighted_medians(1:end), 'linewidth', 2);
                %                 scatter(f_list, weighted_vals, [], [1, 0, 0; .9, .85, 0; 0, 0, 1], 'filled');
                %                 axis([4, 18, -20 + min((weighted_vals), [], 'all'), 20 + max((weighted_vals), [], 'all')])
                %                 xline(top_peak);
                %                 title('Post - Pre', 'FontSize', 18);
                %                 xlabel('Stim Frequency', 'FontSize', 15);
                %                 ylabel('∆Weighted SCC', 'FontSize', 15);
                %                 legend("Median line", "Trial 1", "Trial 2", "Trial 3", "Peak Freq")
                %
                %                 %disp('')
                %                 saveas(gcf, [saveDir, patID{iPat}, 'POST-PRE_FreqDecisionMedian.png']);
                %                 saveas(allFigs{iPat}, [saveDir, patID{iPat}, 'threeplots_FreqDecisionMedian.png']);
                %
                %                 pptx.addSlide();
                %                 pptx.addPicture(gcf, 'Scale', 'maxfixed');

                %% Boxplots
weighted_vals = weighted_vals(:,1:3);
                %% Weighted Conn
                allthree = sum(~isnan(weighted_vals), 2) == 3;
                onlytwo = sum(~isnan(weighted_vals), 2) == 2;
                onlyone = sum(~isnan(weighted_vals), 2) == 1;
                gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1300, 1500]);
                subplot(3, 1, 3)
                boxplot(weighted_vals(allthree, :)', 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
                hold on
                boxplot(weighted_vals(onlytwo, :)', 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
                scatter(f_list(onlyone), weighted_vals(onlyone, :), 'k', 'filled');
                yline(0, 'LineWidth', 0.5)
                xticks(4:0.5:18)
                xticklabels(4:0.5:18)
                title('WEIGHTED: Post - Pre', 'FontSize', 18);
                xlabel('Stim Frequency', 'FontSize', 15);
                ylabel('∆Weighted SCC', 'FontSize', 15);
                %axis([3.5,18.5,1.2*min(weighted_vals,[],'all'),1.2*max(weighted_vals,[],'all')])
                %yline(20)
                
                
                % just post
                
                                %% Weighted Conn
%                 allthree = sum(~isnan(Post_weighted_all), 2) == 3;
%                 onlytwo = sum(~isnan(Post_weighted_all), 2) == 2;
%                 onlyone = sum(~isnan(Post_weighted_all), 2) == 1;
%                 gcf3 = figure('Renderer', 'painters', 'Position', [5, 5, 1500, 500]);
%                 %subplot(3, 1, 3)
%                 boxplot(Post_weighted_all(allthree, :)', 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
%                 hold on
%                 boxplot(Post_weighted_all(onlytwo, :)', 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
%                 scatter(f_list(onlyone), Post_weighted_all(onlyone, :), 'k', 'filled');
%                 yline(0, 'LineWidth', 0.5)
%                 xticks(4:0.5:18)
%                 xticklabels(4:0.5:18)
%                 title('WEIGHTED: Post', 'FontSize', 18);
%                 xlabel('Stim Frequency', 'FontSize', 15);
%                 ylabel('∆Weighted SCC', 'FontSize', 15);
%                 
%                 
                % saveas(gcf3,[saveDir patID{iPat} 'boxplots_FreqDecisionMedian.png']);

               % RTP_vals(:,4)=[];
%local_vals(:,4)=[];
                %% LOCAL
                %gcf4=figure('Renderer','painters','Position',[5, 5, 1300, 500]);
                subplot(3, 1, 2)
                % boxplot fig
                allthree = sum(~isnan(local_vals), 2) == 3;
                onlytwo = sum(~isnan(local_vals), 2) == 2;
                onlyone = sum(~isnan(local_vals), 2) == 1;

                boxplot(local_vals(allthree, :)'+b(3), 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
                hold on
                boxplot(local_vals(onlytwo, :)'+b(3), 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'Widths', 0.2);
                scatter(f_list(onlyone), local_vals(onlyone, :), 'k', 'filled');
                xticks(4:0.5:18)
                xticklabels(4:0.5:18)
                yline(0, 'LineWidth', 0.5)

                title('LOCAL: Post - Pre', 'FontSize', 18);
                xlabel('Stim Frequency', 'FontSize', 15);
                ylabel('∆Local SCC', 'FontSize', 15);
axis([3.5,18.5,1.2*min(weighted_vals,[],'all'),1.2*max(weighted_vals,[],'all')])
                %saveas(gcf4,[saveDir patID{iPat} '_boxplotsLOCAL.png']);

                %% RTP
                %gcf5=figure('Renderer','painters','Position',[5, 5, 1300, 500]);
                subplot(3, 1, 1)
                % boxplot fig
                allthree = sum(~isnan(RTP_vals), 2) == 3;
                onlytwo = sum(~isnan(RTP_vals), 2) == 2;
                onlyone = sum(~isnan(RTP_vals), 2) == 1;

                boxplot(RTP_vals(allthree, :)'+b(3), 'positions', f_list(allthree), 'labels', f_list(allthree), 'colors', 'b', 'Widths', 0.2);
                hold on
                boxplot(RTP_vals(onlytwo, :)'+b(3), 'positions', f_list(onlytwo), 'labels', f_list(onlytwo), 'colors', [1, 0, 0], 'widths', 0.2);
                scatter(f_list(onlyone), RTP_vals(onlyone, :), 'k', 'filled');
                yline(0, 'LineWidth', 0.5)
                xticks(4:0.5:18)
                xticklabels(4:0.5:18)
                axis([3.5,18.5,1.2*min(weighted_vals,[],'all'),1.2*max(weighted_vals,[],'all')])
                title('RTP: Post - Pre', 'FontSize', 18);
                xlabel('Stim Frequency', 'FontSize', 15);
                ylabel('∆RTP SCC', 'FontSize', 15);

                text(17, 100, '       ', 'BackgroundColor', 'b')
                text(17.5, 100, 'Three Trains')

                text(17, 80, '       ', 'BackgroundColor', 'r')
                text(17.5, 80, 'Two Trains')

                saveas(gcf3, [saveDir, patID{iPat}, '_boxplots_subplots.png']);

                %close(gcf);
            end
        end
    end
end
close all
% y= weighted_vals(allthree | onlytwo,:);
x = f_list(allthree | onlytwo);
%
% figure;
% highest_min = min(weighted_vals(allthree | onlytwo,:),[],2);
% smoothed_highest_min = halfhz_avg_inclusive(x,highest_min);
% plot(x,smoothed_highest_min,"Marker","*")
% title('Highest lower bound','FontSize',16)
% xticks(4:0.5:18)
%
% figure
% mdn_div_std = weighted_medians(allthree | onlytwo,:)./tmp_stdev(allthree | onlytwo)';
% plot(x,halfhz_avg_inclusive(x,mdn_div_std),"Marker","*")
% title('Median/stdev','FontSize',16)
% xticks(4:0.5:18)

% figure
 over0_has1_3 = ~isnan(weighted_vals(:, 1)) & ~isnan(weighted_vals(:, 3)) & ~any(weighted_vals < 0, 2);
 mdn_smooth = halfhz_avg_inclusive(f_list, weighted_medians);
% post_mdnsmooth = halfhz_avg_inclusive(f_list, Post_weighted_all);
% x_abovethresh = f_list(over0_has1_3); %freqs to look at
% y_thresh = weighted_vals(over0_has1_3, :);
% plot(f_list, mdn_smooth)
% hold on
% scatter(x_abovethresh, mdn_smooth(over0_has1_3), 'filled')
% 
% title('Min 0, Max Median', 'FontSize', 16)
% xticks(4:0.5:18)

%% CTZ

%% increasing resonance

% building?
% of_interest = [f_list', weighted_vals];
% of_interest = of_interest(over0_has1_3, :);
% delta_1_3 = of_interest(:, 4) - of_interest(:, 2);
% med_of_interest = mdn_smooth(over0_has1_3);
% med_norm = (med_of_interest / max(med_of_interest))';
% delta_norm = delta_1_3 / max(delta_1_3);
% dist = delta_norm .* abs(delta_norm) + med_norm .* abs(med_norm)
% best_idx = find(dist == max(dist));
% best_f = of_interest(best_idx, 1);
% 
% scatter(delta_1_3, med_of_interest, 'filled')
% hold on
% scatter(delta_1_3(best_idx), med_of_interest(best_idx), 'filled', 'SizeData', 150)
% labelpoints(delta_1_3, med_of_interest, char(string(of_interest(:, 1))), 'FontSize', 12, 'buffer', 0.15)
% xlabel('Trial 1->Trial 3 Increase', 'FontSize', 16)
% ylabel('Median value', 'FontSize', 16)
close all

%% Std deviation version
std_smooth = halfhz_avg_inclusive(f_list, tmp_stdev);
bool_interest = allthree & ~any(weighted_vals < 0, 2)
of_interest = [f_list', weighted_vals];
of_interest = of_interest(bool_interest, :);
med_of_interest = mdn_smooth(bool_interest);
stdev_interest = std_smooth(bool_interest);
med_norm = (med_of_interest / max(med_of_interest))';
std_norm = (1 - (stdev_interest / max(stdev_interest)))';
dist = std_norm .* abs(std_norm) + med_norm .* abs(med_norm)
best_idx = find(dist == max(dist));
best_f = of_interest(best_idx, 1);

scatter(std_norm, med_of_interest, 'filled')
hold on
scatter(std_norm(best_idx), med_of_interest(best_idx), 'filled', 'SizeData', 150)
labelpoints(std_norm, med_of_interest, char(string(of_interest(:, 1))), 'FontSize', 12, 'buffer', 0.15)
legend({'', 'RF Selection'})
xlabel('1-norm stdev', 'FontSize', 16)
ylabel('Median value', 'FontSize', 16)
% just post
% std_smooth = halfhz_avg_inclusive(f_list, post_stdev);
% bool_interest = allthree & ~any(Post_weighted_all < 40, 2)
% of_interest = [f_list', Post_weighted_all];
% of_interest = of_interest(bool_interest, :);
% med_of_interest = post_mdnsmooth(bool_interest);
% stdev_interest = std_smooth(bool_interest);
% med_norm = (med_of_interest / max(med_of_interest))';
% std_norm = (1 - (stdev_interest / max(stdev_interest)))';
% dist = std_norm .* abs(std_norm) + med_norm .* abs(med_norm)
% best_idx = find(dist == max(dist));
% best_f = of_interest(best_idx, 1);
% 
% scatter(std_norm, med_of_interest, 'filled')
% hold on
% scatter(std_norm(best_idx), med_of_interest(best_idx), 'filled', 'SizeData', 150)
% labelpoints(std_norm, med_of_interest, char(string(of_interest(:, 1))), 'FontSize', 12, 'buffer', 0.15)
% legend({'', 'RF Selection'})
% xlabel('1-norm stdev', 'FontSize', 16)
% ylabel('Median value', 'FontSize', 16)


% plot(f_list(allthree | onlytwo),(weighted_medians(allthree | onlytwo,:).*abs(weighted_medians(allthree | onlytwo,:)))./tmp_stdev(allthree | onlytwo)')
% title('Median/stdev','FontSize',16)
% xticks(4:0.5:18)
