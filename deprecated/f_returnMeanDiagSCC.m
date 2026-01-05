function [allPatMeanDiag,allstdDiag,allminDiag,allmaxDiag] = f_returnMeanDiagSCC(AllF3toCP2, AllStimFreq, v_FreqAxis, window,stim_idx)
% p_RF_SCC

nPat = size(AllF3toCP2,3);

allPatMeanDiag = cell(1,nPat);
for iPat = 1:nPat
    matrixPat = AllF3toCP2;
    clear stimFreq
    stimFreq = AllStimFreq;
    meanDiag = zeros(1,numel(stimFreq));
    for iInter = stim_idx
        currFreq = stimFreq(iInter);
        idxStim = stimFreq >= currFreq-0.5 & stimFreq <= currFreq+0.5;
        idxResp = find(v_FreqAxis >= window(1) & v_FreqAxis <= window(2));
        stdDiag(iInter) = nanstd(nanmean(matrixPat(idxResp,idxStim),1));
        meanDiag(iInter) = nanmean(nanmean(matrixPat(idxResp, idxStim)));
        minDiag(iInter) = min(nanmean(matrixPat(idxResp,idxStim),1));
        maxDiag(iInter) = max(nanmean(matrixPat(idxResp,idxStim),1));
    end
   
    allPatMeanDiag = meanDiag;
    allstdDiag = stdDiag;
    allminDiag = minDiag;
    allmaxDiag = maxDiag;
end

end