% clearvars -except subjids subjid sessions i_subject i_session work_dir suffix;
% eeglab;
% Step0_SetDirectories;
% disp(subjid);
L = load('Step3_FasterInput.mat');

ElectrodeLocation = '/Users/thomas/Documents/GitHub/rf_matlab/repo/standard-10-5-cap385.elp';
L.option_wrapper.options.job_filename = 'Step3_FasterInput.eegjob';
L.option_wrapper.options.file_options.channel_locations = ElectrodeLocation;
L.option_wrapper.options.current_file =[];
L.option_wrapper.options.file_options.current_file = [];
L.option_wrapper.options.file_options.resume=0;
L.option_wrapper.option.file_options.make_subdirectories=0;
% 
data_dir = [work_dir subjid filesep  'Step3_CutFiles' suffix filesep];
disp(data_dir)
cd(data_dir);
files = dir('*.set');
files = {files.name};

EEG = pop_loadset(fullfile(data_dir, files{1}));

% if we have EOG, use it
if EEG.nbchan > 64
    L.option_wrapper.options.ica_options.EOG_channels = [1,2,3,find(strcmp({EEG.chanlocs.labels}', 'VEOG'))];
else
    L.option_wrapper.options.ica_options.EOG_channels = [1,2,3];
end

L.option_wrapper.options.file_options.folder_name = [work_dir subjid filesep 'Step3_CutFiles' suffix filesep];

mkdir([work_dir subjid filesep 'Step4_FasterOutput' suffix filesep]);
L.option_wrapper.options.file_options.output_folder_name = [work_dir subjid filesep 'Step4_FasterOutput' suffix filesep];
% L.option_wrapper.options.channel_options.channel_rejection_on = 0;
% L.option_wrapper.options.channel_options.bad_channels = [32];

rng(1234, "twister");
FASTER(L.option_wrapper)

