function [burstStartEnd,burstIPI, stimFreqMedian] = f_CheckTriggers(triggers, minPulsePerBurst, srate)

% keeps triggers only marking start and end of the burst
% minimal number of pulses per burst = minPulsePerBurst, otherwise burst removed
triggers_bk = triggers;

IPI = diff(triggers);
countBurst = 1;
count = 1;
burstStartEnd = [];
countPulses = 1;
currIPI = [];
burst2delete = [];
currPulse = [];
consecutive_tolerance = srate / 500;
while count <= numel(IPI)-1
    
    currIPI(countPulses) = IPI(count);
    startIPI = IPI(count);
    startTrigger = triggers(count);
    count  = count + 1;
    countPulses = countPulses + 1;
    
    % if current IPI is similar to the previous IPI
    if IPI(count) >= startIPI-consecutive_tolerance && IPI(count) <= startIPI+consecutive_tolerance 
        nextBurst = 0;

        % while in the same burst
        while nextBurst ~= 1 & count <= numel(IPI)

            % if current IPI is similar to previous IPI
            if IPI(count) >= startIPI-consecutive_tolerance && IPI(count) <= startIPI+consecutive_tolerance  
                currIPI(countPulses) = IPI(count);
                count = count + 1; % move to the next pulse
                countPulses = countPulses + 1; % increment the counter

            % current IPI is not similar to the previous
            else
                if countPulses > 9 & countPulses < minPulsePerBurst & abs(IPI(count+1)-IPI(count+2))>3 % at least 10 pulses with the same IPI which are unlikely to be random
                    triggers(count+1) = [];
                    IPI = diff(triggers); % removes current pulse if IPI is different, but more than 2 pulses in a row & recalculates IPI
                else % start of the next burst
                    countPulses = 1;
                    nextBurst = 1;
                end
                
            end
        end
        
    else
        countPulses = 1;
        continue;
    end
    
    burstStartEnd(countBurst,1) = startTrigger;
    burstStartEnd(countBurst,2) = triggers(count);
    burstIPI{countBurst} = currIPI;
    if numel(currIPI) <= minPulsePerBurst-1
        burst2delete = [burst2delete countBurst];
        % disp(['sample: ' num2str(triggers(count))])
        % disp(['numel(currIPI): ' num2str(numel(currIPI))])
    end
    currIPI = [];
    countBurst = countBurst + 1;
    count = count + 1;
    
end
% remove bursts that have ~= 40 pulses
if ~isempty(burst2delete)
    burstStartEnd(burst2delete,:) = [];
    burstIPI(burst2delete) = [];
end


%% estimate stimulation frequency for each burst

stimFreqMost = zeros(1,numel(burstIPI));
stimFreqMean = zeros(1,numel(burstIPI));
stimFreqMedian = zeros(1,numel(burstIPI));
for iBurst = 1:numel(burstIPI)
    [a,b,c] = unique(burstIPI{iBurst});
    nHz = zeros(1,numel(a));
    for uniq = 1:numel(a)
        nHz(uniq) = nnz(c == uniq);
    end
    [m, idx] = max(nHz);
    stimFreqMost(iBurst) = round(srate./a(idx),1);
    stimFreqMean(iBurst) = round(srate/mean(burstIPI{iBurst}),1);
    stimFreqMedian(iBurst) = round(srate/median(burstIPI{iBurst}),1);
end

 % figure;
 % plot(stimFreqMost, 'linewidth',2);
 % hold on;
 % plot(stimFreqMedian, 'linewidth',2)
 % plot(stimFreqMean, 'linewidth',2)

end