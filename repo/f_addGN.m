function [noisySig, noise] = f_addGN(signal, srate)
% original version from FileExchange "addGaussianNoise.m"

reqSNR = 2;

a = 50;
b = 70;
f = round((b-a).*rand(1,1) + a); % random frequency selection between 50 & 70 hz

% generate sine wave at 50/60/70 Hz
lengthSig = length(signal);
t = [1:lengthSig]/srate;
signalGamma = sin(2*pi*f*t);

%% add white noise
sigEner = norm(signalGamma(:))^2;                    % energy of the signal
noiseEner = sigEner/(10^(reqSNR/10));        % energy of noise to be added
noiseVar = noiseEner/(length(signalGamma(:))-1);     % variance of noise to be added
noiseStd = sqrt(noiseVar);                   % std. deviation of noise to be added
noise = noiseStd*randn(size(signalGamma));           % noise
noisySig = signalGamma+noise+signal;                        % noisy signal
end
