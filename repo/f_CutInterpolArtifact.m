function [X,Y,XX, YY] = f_CutInterpolArtifact(YPre, YPost, nInterpol, nPrePost)

% nPrePost - how many Y data points per side
% nInterpol - gap size

XX = nPrePost+1:nPrePost+1+nInterpol-1;
X = [1:nPrePost XX(end)+1:XX(end)+nPrePost+1];
Y = [YPre(end-(nPrePost-1):end) YPost(1:nPrePost+1)];
XX = nPrePost+1:nPrePost+1+nInterpol-1;
% YY = interp1(X,Y,XX, 'pchip'); % spchip : hape-preserving piecewise cubic interpolation
YY = interp1(X,Y,XX, 'linear');



% % XX = 2:2+nInterpol-1;
% % X = [1 XX(end)+1:XX(end)+1];
% % Y = [YPre(end) YPost(1)];
% % YY = interp1(X,Y,XX, 'cubic');




end

