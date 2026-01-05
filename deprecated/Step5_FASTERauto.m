eeglab;
%addpath('/Users/thomas/Documents/GitHub/RF_StepsMatlab/')
set_rf_directories;

L = load('Step5_FasterInput.mat');
L.option_wrapper.options.job_filename = 'Step5_FasterInput.eegjob';
L.option_wrapper.options.file_options.channel_locations = ElectrodeLocation;
L.option_wrapper.options.current_file =[]
L.option_wrapper.options.file_options.current_file = []
L.option_wrapper.options.file_options.resume=0;
L.option_wrapper.options.file_options.make_subdirectories=0;


%subjid = 'RFUCLA-001_IP1';
%work_dir = '/Users/thomas/Documents/interrogation/';%'/Volumes/Files/Interrogation/';
% work_dir = '/Volumes/Files/Interrogation/5x5/SCC23-030/';
L.option_wrapper.options.file_options.folder_name = [work_dir subjid '/Step5_FasterInput_PostManualICRemoval/']

mkdir([work_dir subjid '/Step6_FasterOutput/']);
L.option_wrapper.options.file_options.output_folder_name = [work_dir subjid '/Step6_FasterOutput/']
FASTER(L.option_wrapper)
