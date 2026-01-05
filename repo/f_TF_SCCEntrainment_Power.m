function [ m_TFPRE,m_TFPOST, m_TFPRE_abs,m_TFPOST_abs,v_FreqAxis, outpar ] = f_TF_SCCEntrainment_Power( v_Sig, TOIBL,TOIPost, fmin, fmax, s_StDevCycles,  srate, alg, s_FreqSeg) % 'abs', 'zs', 'zsn', '01'
% algorithm input: 
% 'abs' = absolute
% 'zs' = zscore
% 'zsn' = zscore normalized by integrated value
% '01' = normalized between 0 & 1


%% example values from MV
% s_MinFreqHz = 4;
% s_MaxFreqHz = 550; 
% s_FreqSeg = 1000; % Frequency resolution between min & max
% s_StDevCycles = 3;
% s_Magnitudes = 1;
% s_SquaredMag = 0;
% s_MakeBandAve = 0;
% s_Phases = 0;
% s_TimeStep = [];
% s_DyadicScale = 0;


%% define values
if isempty(s_FreqSeg)
    s_FreqSeg = (fmax - fmin)/0.5; % frequency resolution
end
% s_FreqSeg = ceil( (fmax - fmin)/2);
% s_StDevCycles = 3;
s_Magnitudes = 1;
s_SquaredMag = 0;
s_MakeBandAve = 0;
s_Phases = 0;
s_TimeStep = [];
s_DyadicScale = 0;
s_PadPercent = 25;

outpar.freqres = s_FreqSeg;
outpar.timeres = s_StDevCycles;
outpar.magnitude = s_Magnitudes;
outpar.sqmag = s_SquaredMag;
outpar.average = s_MakeBandAve;
outpar.phase = s_Phases; 
outpar.timestep = s_TimeStep;
outpar.dyad = s_DyadicScale;


%%Run for all channels separately
v_PSD_PRE = [];
v_PSD_POST = [];
v_Freq = [];
for iChan = 1:size(v_Sig,1)
    %% calculation
    % timeaxis in sec
    [m_GaborWT, v_TimeAxis, v_FreqAxis] = ...
        f_GaborAWTransformMatlab(...
        v_Sig(iChan,:), ...
        srate, ...
        fmin, ...
        fmax, ...
        s_FreqSeg, ...
        s_StDevCycles, ...
        s_Magnitudes, ...
        s_SquaredMag, ...
        s_MakeBandAve, ...
        s_Phases, ...
        s_TimeStep, ...
        s_DyadicScale,...
        s_PadPercent); % , ...
    %    s_DyadicScale);
    
    %average over TOI
    
    %%ABSOLUTE%%%%%%%%%%%%%%%
    v_PSD_PRE_abs(:,iChan) = mean(m_GaborWT(:,TOIBL),2);
    v_PSD_POST_abs(:,iChan) = mean(m_GaborWT(:,TOIPost),2);
    %%ABSOLUTE%%%%%%%%%%%%%%%
    
    
    %%Relative%%%%%%%%%%%%%%%
    v_PSD_PRE(:,iChan) = mean(m_GaborWT(:,TOIBL),2)/sum(mean(m_GaborWT(:,TOIBL),2))*100;
    v_PSD_POST(:,iChan) = mean(m_GaborWT(:,TOIPost),2)/sum(mean(m_GaborWT(:,TOIPost),2))*100;
    %%Relative%%%%%%%%%%%%%%%

end

m_TFPRE = v_PSD_PRE;
m_TFPOST = v_PSD_POST;

m_TFPRE_abs = v_PSD_PRE_abs;
m_TFPOST_abs = v_PSD_POST_abs;

end

