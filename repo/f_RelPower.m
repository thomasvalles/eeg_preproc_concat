function rel_power = f_RelPower(abs_power, par_relpower)

% parameter inputs:
%     par_relpower.fmin = rp_fmin;  % minimal freq for power spectrum
%     par_relpower.fmax = rp_fmax; % maximal freq for power spectrum
%     par_relpower.overlapBands = v_OverlapBands; % overlapping bands for
%     peak detection
%     par_relpower.subjID = str_SubjID{iSubj}; % subject ID
%     par_relpower.alg = 'pwelch';  % spectrum method
%     par_relpower.dirDB = str_PathDB; % path to database
%     par_relpower.dirFigOut = str_DirFigOut; % path for saving figures
%     par_relpower.peakCol = 'krg';  % colors of peaks triangles

% input from absolute power calculation
%     abs_power.cond = str_Condition;
%     abs_power.spect =  m_AllSpect;
%     abs_power.freq =  v_Freq;
%     abs_power.parFFT = parFFT;
%     abs_power.remvdChan = v_ChanIdx2Rem;
%     abs_power.chanLbls = str_ChanLab;
%     abs_power.chanLblsRmvd = str_Chan2Rem;



%----------------------------------------------------------------------
% relative power calculation
%----------------------------------------------------------------------

v_Freq = abs_power.freq;
fmin = par_relpower.fmin;
fmax = par_relpower.fmax;



m_CurrSpect = squeeze(mean(abs_power.spect,3));
if size(m_CurrSpect,1) > size(m_CurrSpect,2)
    m_CurrSpect = m_CurrSpect';
end

str_ChanLab = abs_power.chanLbls;


% normalize power
idx_fmin = find(v_Freq>=fmin, 1, 'first');
idx_fmax =  find(v_Freq>=fmax, 1, 'first');
idx_FreqRange = idx_fmin:idx_fmax;

v_FreqRange = v_Freq(idx_FreqRange);
v_SpectRange = m_CurrSpect(:,idx_FreqRange);
v_SpectRangePct = (v_SpectRange./repmat(sum(v_SpectRange,2), 1,size(v_SpectRange,2)))*100;

%--------------------------------------------------------------------------%
% write output
%--------------------------------------------------------------------------%

rel_power.subjID = par_relpower.subjID;
rel_power.spectRange = v_SpectRange;
rel_power.freq = v_Freq;
rel_power.spectraABS = m_CurrSpect;
rel_power.spectraREL = v_SpectRangePct;
rel_power.chanLbls =  str_ChanLab;
rel_power.freqRange = v_FreqRange;
rel_power.par = par_relpower;



end