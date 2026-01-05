%% load eeg Post-FASTER
%% copy to the VI folder

%% load out EEG from Cpz0 folder, do the epochs # match? (Faster epochs - numel(RemovedFreqs{1}) == Faster epochs)

%% load original triggers

%% load post-VI triggers

%% use RemovedFreqs{1} --> Save new trigger out file



%% Step1 - make new "final trigger" file
rootdir = '/Volumes/Files/Interrogation/Training_OldSubjs/'
cd(rootdir)
dirs = string({dir().name});
subjs = dirs(~contains(dirs, '.'));
redotrig=1;
if redotrig ==1 
for isubj = 1:numel(subjs)
    subj = subjs(isubj)
    cd(string(rootdir)+subj)

    dates = string({dir().name});
    if length(dates) > 2
        dates = dates(3:end); % get rid of . ..
        for idate = 1:numel(dates)
            cd(string(rootdir)+subj)

            date = dates(idate)
            if isdir(date)
                cd(date)
                step3trigger_file = dir(['Step3_BrushedTriggers/', '*Interrog.mat']);
                step3trigger_file = string({step3trigger_file.name});
                step3trigger_file = step3trigger_file(~contains(step3trigger_file, '._'));

                old_finaltrigger_file = dir(['Step7_FinalTriggerFiles_EpochsRemoved/', '*Interrog.mat']);
                old_finaltrigger_file = string({old_finaltrigger_file.name});
                old_finaltrigger_file = old_finaltrigger_file(~contains(old_finaltrigger_file, '._'));
                if ~isempty(old_finaltrigger_file) && ~isempty(step3trigger_file)
                    old_finaltrigger = load(['Step7_FinalTriggerFiles_EpochsRemoved/', char(old_finaltrigger_file)]);
                    step3trigger = load(['Step3_BrushedTriggers/', char(step3trigger_file)]);

                    % check correct sizes..
                    if numel(old_finaltrigger.RemovedFreqs) > 0
                        if ~(length(old_finaltrigger.RemovedFreqs{1}) + length(old_finaltrigger.stimFreqMedian) == length(step3trigger.stimFreqMedian))
                            if ~(old_finaltrigger.stimFreqMedian == step3trigger.stimFreqMedian) % could also be equal if newer file
                                disp("WRONG SIZE TRIGGERS: "+string(date))
                            end
                        end
                    else
                        if ~(numel(old_finaltrigger.stimFreqMedian) == numel(step3trigger.stimFreqMedian))
                            error(char("WRONG SIZE TRIGGERS: "+string(date)))
                        end
                    end

                    % if old step3 same as old step7 (ie export shouldn't be any
                    % different)
                    if (numel(old_finaltrigger.stimFreqMedian) == numel(step3trigger.stimFreqMedian))
                        disp(['New same as old: ', date'])
                    end


                    burstIPI = step3trigger.burstIPI;
                    burstStartEnd = step3trigger.burstStartEnd;
                    stimFreqMedian = step3trigger.stimFreqMedian;
                    RemovedFreqs = old_finaltrigger.RemovedFreqs;
                    if numel(old_finaltrigger.RemovedFreqs) > 0

                        trialrej = RemovedFreqs{1};
                    else
                        trialrej = [];
                    end


                    if ~isdir('Step7_FinalTriggerFiles')
                        mkdir('Step7_FinalTriggerFiles')
                    end
                    save(['Step7_FinalTriggerFiles/', char(old_finaltrigger_file)], 'burstIPI', 'burstStartEnd', 'stimFreqMedian', 'RemovedFreqs', 'trialrej')


                    clearvars step3trigger old_finaltrigger old_finaltrigger_file step3trigger_file trialrej burstIPI burstStartEnd stimFreqMedian RemovedFreqs


                end
            end
        end
    end
end
end 
%% step 2 - in the terminal move the trigger files to trigger dump for testing them..
%% step 3 - check whether old "final" eeg  + # of removed frequencies = full size 
%(eg we're removing the same # of frequencies as in the final)

final_eegs = string({dir('~/Downloads/Interrogation_Cpz0/FinalEEG/RF/*.set').name});
badeeg_match =[];
for ieeg = 1:numel(final_eegs)
EEG = pop_loadset(['~/Downloads/Interrogation_Cpz0/FinalEEG/RF/' char(final_eegs(ieeg))]);
subj = final_eegs(ieeg).split('_');

try 
    datetime(subj(2));
    date = subj(2);
catch 
    datetime(subj(3));
    date = subj(3);
end 
% date = subj(2);
subj = subj(1);
% subj = subj.replace('001-','SCC22-')


triggers=  string({dir('/Volumes/Files/Interrogation/Training_OldSubjs/FinalTrigger_Dump/*.mat').name});

new_trig = load(['/Volumes/Files/Interrogation/Training_OldSubjs/FinalTrigger_Dump/' char(triggers(contains(triggers,subj) & contains(triggers,date)))]);

% is the old "output eeg" (with epochs removed) plus the number of the
% epochs removed according to new trig == to the total number of epochs 
% if ~isfield(new_trig,'trialrej')
% if ~isempty(new_trig.RemovedFreqs)
    out_bool = (numel(new_trig.RemovedFreqs{1}) + size(EEG.data,3) == numel(new_trig.stimFreqMedian));
% disp(subj + date + ":   " + out_bool)
% else
if out_bool==0
    out_bool = (size(EEG.data,3) == numel(new_trig.stimFreqMedian));
% disp(subj + date + ":   " + )
end 
% else 
%      out_bool = (size(EEG.data,3) == numel(new_trig.stimFreqMedian));
% end 
if out_bool==0
    error(['mismatch' char(subj)])
end


    


new_eegs = string({dir('/Volumes/Files/Interrogation/Training_OldSubjs/FasterOut_Dump/*.set').name});
new_eeg = pop_loadset(['/Volumes/Files/Interrogation/Training_OldSubjs/FasterOut_Dump/' char(new_eegs(contains(new_eegs,subj) & contains(new_eegs,date)))]);

iepoch = 1;
while ~isempty(new_trig.RemovedFreqs) && ismember(iepoch,new_trig.RemovedFreqs{1}) 
    
    iepoch=iepoch+1;

end
same_epoch_data = all(all(new_eeg.data(:,:,iepoch) == EEG.data(:,:,1)));
if (~same_epoch_data )
    badeeg_match = [badeeg_match,string(subj)]
%     error(['Epochs are off for ' char(subj) char(date)])
end
end 
%% step 4 -- check that final EEG epoch 1 (or 2... so on) is equal to the faster output eeg (no cpz issues)
% the files would be different if i have a bad faster ouput file 





