function [v_SigPadded, s_PadSz] = f_Padding(m_DataMatSubS, s_PadPercent)
% input: reshaped data line vector 

    s_LenData = length(m_DataMatSubS);
    s_PadSz = floor(s_LenData/100*s_PadPercent);
    
    v_PadArrayPre = m_DataMatSubS(:,s_PadSz:-1:1);  
    v_PadArrayPost =m_DataMatSubS(:,end-s_PadSz-1:end); 
    v_SigPadded = [v_PadArrayPre, m_DataMatSubS, v_PadArrayPost]; 
   
end

