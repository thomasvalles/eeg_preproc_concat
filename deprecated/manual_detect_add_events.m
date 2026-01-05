
close all
EEG= pop_loadset('/Users/neuromodit/Documents/Interrogation/SCC22-022-2022-06-29/Step2_MergedSetFiles/00-022_2022-06-29_14-13-55_Merged.set')
clearvars -except EEG
electrode = 'FCz';
e_idx = find(string({(EEG.chanlocs.labels)})==electrode)
data = EEG.data(e_idx,:);
DataCorrected = abs(detrendnonlin((data-mean(data)),2)') ;
plot(DataCorrected)
thresh= input('enter thresh:')
if thresh<0
suprathresh = find(DataCorrected<-thresh);
else 
    suprathresh = find(DataCorrected>thresh);
end

endpoints = suprathresh(find(diff(suprathresh/2000)>5));
startpoints = [suprathresh(1),suprathresh(find(diff(suprathresh/2000)>5)+1)] ;
numtrains = numel(startpoints);
endpoints(numtrains) = suprathresh(end);
startpoints =startpoints;
endpoints=endpoints+34;

% 
% checkthis = suprathresh(1);
% badidx=[]
% for idx= 2:numel(suprathresh)
% if suprathresh(idx)-checkthis<1/19*2000
%     badidx=[badidx idx];
% end
% 
% checkthis=suprathresh(idx);
% end
% 
% 
x = uniquetol(suprathresh,1/20*2000,'DataScale',1);
pulses=[];
for itrain = 1:numel(startpoints)
    numpulses(itrain) = numel(find(x<=endpoints(itrain) & x>=startpoints(itrain)));
    pulses = [pulses,find(x<=endpoints(itrain) & x>=startpoints(itrain))]
end
table(startpoints',endpoints',numpulses')
remidx = find(numpulses<39);
startpoints(remidx)=[];
endpoints(remidx)=[];

numtrains = numel(startpoints);
pulses_final=[];
numpulses_final=[];
for itrain = 1:numel(startpoints)
    numpulses_final(itrain) = numel(find(x<=endpoints(itrain) & x>=startpoints(itrain)));
    pulses_final = [pulses_final,x(x<=endpoints(itrain) & x>=startpoints(itrain))];
end

EEG.event=[];
for ipulse = 1:numel(pulses_final)
         EEG.event(ipulse).latency=pulses_final(ipulse);
     EEG.event(ipulse).type='0001';
end
table(startpoints',endpoints',numpulses_final')

% for itrain=1:numtrains
% startpulse=startpoints(itrain);
% endpulse=endpoints(itrain);
% init_event= numel(EEG.event);
% for ipulse =1:40
%    eventnum = init_event+ipulse;
%     EEG.event(eventnum).latency=startpulse + (ipulse-1)*(endpulse-startpulse)/39;
%     EEG.event(eventnum).type='0001';
% end
% 
% end
pop_saveset(EEG,'/Users/neuromodit/Documents/Interrogation/SCC22-022-2022-06-29/Step2_MergedSetFiles/00-022_2022-06-29_14-13-55_Merged.set')