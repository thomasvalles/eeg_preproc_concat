%% BE CAREFUL OF SORTING, changed 11.18.22
% p_Entrainment_ASC_InputOutput_WholeBrain_4ACW

close all;
clc;
clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;%clear;
tic;
% load data

saveFig = 1;

eeglab;

Step0_SetDirectories;
dirData = [work_dir subjid filesep 'Step4_FasterOutput' suffix filesep];
matDir = [work_dir subjid filesep 'Step5_FinalTriggerFiles' suffix filesep];
saveDir = [work_dir subjid filesep 'Step6_ASCOutput' suffix filesep];

if ~isdir(saveDir)
    mkdir(saveDir)
end

TOIlength = 1000;
HalfHz = 2;


cd(dirData)
files = dir('*.set');
files = {files.name}';
patID = [];
for iFile = 1:numel(files)
    tmpsplit = strsplit(files{iFile}, '_');
    patID{iFile} = files{iFile}(1:end - 14);
end
nPat = numel(files);


allStimFreq = cell(1, nPat);
aSCPreAll = cell(1, nPat);
aSCPostAll = cell(1, nPat);
for iPat = 1:nPat
    rel_powerPost = struct();
    rel_powerPre = struct();

    disp(['Patient ', num2str(iPat)]);

    EEG = pop_loadset(fullfile(dirData, [files{iPat}]));
    if (EEG.nbchan ~= 64)
        %%Need to remove channels greater than 64
        EEG = pop_select(EEG, 'channel', 1:64)

    end

    chanLabels = {EEG.chanlocs.labels}';
    idxSeed = find(ismember(chanLabels, 'F3'));
    srate = EEG.srate;
    if srate ~= 1000
        error('srate ~= 1000');
    end

    nFreq = size(EEG.data, 3);
    nChan = size(EEG.data, 1);
    v_TimeAxis = ((1:size(EEG.data, 2)) / srate) * 1000 - 1000;

    idxM1M2 = ismember({EEG.chanlocs.labels}', {'M1', 'M2', 'F3'});
    chanToDo = 1:64;
    chanToDo(idxM1M2) = [];

    %% load stimulation frequencies
    load(fullfile(matDir, [patID{iPat}, '_TriggersInterrog.mat']));

    %     [sorted, idxSorted] =  sort(stimFreqMedian, 'ascend');

    %     stimFreqMedian = round(stimFreqMedian*4)/4

    % first magstim pulses are a little too fast. make sure the first freq
    % is the same as the second
    stimFreqMedian(1) = stimFreqMedian(2);
    stimFreqMedianSorted = stimFreqMedian;

    burstStartEndSorted = burstStartEnd;
    if TOIlength == 400
        nFreqOut = 10; % 10 for 400 ms, 19 for 1000 ms
    elseif TOIlength == 1000
        nFreqOut = 19; % 10 for 400 ms, 19 for 1000 ms
    end


    aSCPre = nan(nFreqOut, length(stimFreqMedianSorted), nChan, nChan);
    aSCPost = nan(nFreqOut, length(stimFreqMedianSorted), nChan, nChan);

    %disp(['Channel ' num2str(iChan1)]);
    for iFreqStim = 1:length(stimFreqMedianSorted)
        disp(['Stim', num2str(iFreqStim)]);
        TOIbsl = v_TimeAxis < 0 & v_TimeAxis >= -TOIlength; %-1000;
        TOI = v_TimeAxis > 0 & v_TimeAxis <= TOIlength; %1000;

        fmin = 1;
        fmax = 20;
        s_StDevCycles = 3;
        s_FreqStep = 0.2;
        s_FreqSeg = (fmax - fmin) / s_FreqStep;
        alg = 'abs';
        clear m_TF_PRE m_TF_POST

        %%         returns PSD
        %         [ m_TF_PRE, m_TF_POST, v_FreqAxis, outpar ] = f_TF_SCCEntrainment( EEG.data(:,:,iFreqStim), TOIbsl,TOI,fmin, fmax, s_StDevCycles,  srate, alg, s_FreqSeg);
        % add for incorporating power to not calc sep
        [m_TF_PRE, m_TF_POST, rel_powerPre(iFreqStim).spectraABS, rel_powerPost(iFreqStim).spectraABS, v_FreqAxis, outpar] = f_TF_SCCEntrainment_Power(EEG.data(:, :, iFreqStim), TOIbsl, TOI, fmin, fmax, s_StDevCycles, srate, alg, s_FreqSeg);
        
        rel_powerPre(iFreqStim).spectraREL = m_TF_PRE';
        rel_powerPost(iFreqStim).spectraREL = m_TF_POST';

        rel_powerPre(iFreqStim).spectraABS = rel_powerPre(iFreqStim).spectraABS';
        rel_powerPost(iFreqStim).spectraABS = rel_powerPost(iFreqStim).spectraABS';

        %% do the correlations
        for iChan1 = 1:nChan
            for iFreqOut = 1:numel(v_FreqAxis)
                clear aSC;

                aSCPre(iFreqOut, iFreqStim, iChan1, :) = f_aSC_SeedtoAll(m_TF_PRE, v_FreqAxis, v_FreqAxis(iFreqOut), iChan1, chanLabels, HalfHz);
                aSCPost(iFreqOut, iFreqStim, iChan1, :) = f_aSC_SeedtoAll(m_TF_POST, v_FreqAxis, v_FreqAxis(iFreqOut), iChan1, chanLabels, HalfHz);
            end
        end

    end

    allStimFreq = stimFreqMedianSorted;

    RelPre = [];
    RelPost = [];
    AbsPre = [];
    AbsPost = [];
    for stimidx = 1:numel(allStimFreq)


        if ismember(stimidx, trialrej)
            rel_powerPre(stimidx).spectraREL = NaN(size(rel_powerPre(stimidx).spectraREL));
            rel_powerPre(stimidx).spectraABS = NaN(size(rel_powerPre(stimidx).spectraABS));
            rel_powerPost(stimidx).spectraREL = NaN(size(rel_powerPost(stimidx).spectraREL));
            rel_powerPost(stimidx).spectraABS = NaN(size(rel_powerPost(stimidx).spectraABS));

            aSCPre(:, stimidx, :, :) = NaN(size(aSCPre(:, stimidx, :, :)));
            aSCPost(:, stimidx, :, :) = NaN(size(aSCPre(:, stimidx, :, :)));
        end


        RelPre(stimidx, :, :) = rel_powerPre(stimidx).spectraREL;
        RelPost(stimidx, :, :) = rel_powerPost(stimidx).spectraREL;
        AbsPre(stimidx, :, :) = rel_powerPre(stimidx).spectraABS;
        AbsPost(stimidx, :, :) = rel_powerPost(stimidx).spectraABS;


    end


    save(fullfile(saveDir, ['Power_entrain_GABOR_', patID{iPat}, '_', num2str(TOIlength), 'ms']), 'rel_powerPre', ...
        'rel_powerPost', 'allStimFreq', 'chanLabels', 'v_FreqAxis', 'v_TimeAxis', 'RemovedFreqs');

    save(fullfile(saveDir, ['Power_forPy_', patID{iPat}, '_', num2str(TOIlength), 'ms']), 'RelPre', ...
        'RelPost', 'AbsPre', 'AbsPost', 'allStimFreq', 'chanLabels', 'v_FreqAxis', 'v_TimeAxis', 'RemovedFreqs');


    clear EEG;
    aSCPreAll = aSCPre;
    aSCPostAll = aSCPost;
    clear inputOutput;

    %% saving scc data
    
                seed =5 ; %f3
            SCCPre = NaN(size(aSCPreAll,1),size(aSCPreAll,4),size(aSCPreAll,2));
            SCCPost = NaN(size(aSCPreAll,1),size(aSCPreAll,4),size(aSCPreAll,2));

        for ichan=1:64
        if ichan > seed
           SCCPre(:, ichan, :) = aSCPreAll(:, :, seed, ichan);
            SCCPost(:, ichan, :) = aSCPostAll(:, :, seed, ichan);
        else
             SCCPre(:, ichan, :) = aSCPreAll(:, :, ichan, seed);
            SCCPost(:, ichan, :) = aSCPostAll(:, :, ichan, seed);
        end 
        end 

    save(fullfile(saveDir, ['SCC_NoSort_', patID{iPat}, '_', num2str(TOIlength), 'ms_UpdatedRel']), 'aSCPreAll', 'aSCPostAll', 'allStimFreq', ...
        'chanLabels', 'v_FreqAxis', 'stimFreqMedian', 'RemovedFreqs');
                save(fullfile(saveDir, ['SCC_F3_', patID{iPat}, '_', num2str(TOIlength), 'ms_UpdatedRel']), 'SCCPre', 'SCCPost', 'allStimFreq', ...
                'chanLabels', 'v_FreqAxis', 'stimFreqMedian', 'RemovedFreqs');
end
