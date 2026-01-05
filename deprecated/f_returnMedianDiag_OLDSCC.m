function [allPatMeanDiag] = f_returnMedianDiag_OLDSCC(AllF3toCP2, AllStimFreq, v_FreqAxis, window)
% p_RF_SCC

nPat = size(AllF3toCP2,3);

allPatMeanDiag = cell(1,nPat);
for iPat = 1:nPat
    matrixPat = AllF3toCP2;
    clear stimFreq
    stimFreq = AllStimFreq;
    meanDiag = zeros(1,numel(stimFreq));
    for iInter = 1:numel(stimFreq)
        currFreq = stimFreq(iInter);
        idxStim = stimFreq >= currFreq-0.25 & stimFreq <= currFreq+0.25;
        idxResp = find(v_FreqAxis >= window(1) & v_FreqAxis <= window(2));
        
        meanDiag(iInter) = nanmedian(nanmean(matrixPat(idxResp, idxStim)));
    end
   
    allPatMeanDiag = meanDiag;

end

end