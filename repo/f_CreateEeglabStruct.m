function [ EEG ] = f_CreateEeglabStruct( data, fs,label, setname )
% mimicks the format of the eeglab EEG struct


EEG.data = data;
EEG.srate = fs;
EEG.setname = setname;
EEG.filename= '';
EEG.filepath = '';
EEG.subject = '';
EEG.group =  '';
EEG.condition = '';
EEG.session = [];
for i = 1:size(data,1)
    EEG.chanlocs(i).labels = label{i};
    EEG.chanlocs(i).ref = '';
    EEG.chanlocs(i).theta = [];
    EEG.chanlocs(i).radius = [];
    EEG.chanlocs(i).X = [];
    EEG.chanlocs(i).Y = [];
    EEG.chanlocs(i).Z = [];
    EEG.chanlocs(i).sph_theta = [];
    EEG.chanlocs(i).sph_phi = [];
    EEG.chanlocs(i).sph_radius = [];
    EEG.chanlocs(i).type = '';
    EEG.chanlocs(i).urchan = [];
end

EEG.comments = setname;
EEG.nbchan = size(data,1);
EEG.trials = 1;
EEG.pnts = size(data,2);
EEG.xmin = 0;
EEG.xmax = size(data,2)/fs;
EEG.times = linspace(EEG.xmin, EEG.xmax, numel(EEG.pnts));
EEG.icaact  = [];
EEG.icawinv  = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icachansind = [];
EEG.urchanlocs = [];
EEG.chaninfo.plotrad = [];
EEG.chaninfo.shrink = [];
EEG.chaninfo.nosedir = '+X';
EEG.chaninfo.nodatchans = [];
EEG.chaninfo.icachansind = [];
EEG.ref = '';
EEG.event = [];
EEG.urevent = [];
EEG.eventdescription = {};
EEG.epoch = [];
EEG.epochdescription = {};
EEG.reject = [];
EEG.stats = [];
EEG.spectdata = [];
EEG.specticaact = [];
EEG.splinefile = '';
EEG.icasplinefile = '';
EEG.dipfit = [];
EEG.history = '';
EEG.saved = 'no';
EEG.etc = [];
eeg_checkset(EEG);

end

