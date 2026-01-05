addpath('/Users/Thomas/Documents/Github/rf_matlab/');

suffix = '';
subjids = {'sub-0540'; 'sub-0592'};%'sub-0347'; 'sub-0413'; 'sub-0472'; 
sessions = {'ses-01'};

for i_subject = 1:size(subjids, 1)
    for i_session = 1:size(sessions, 2)
         %if (i_subject > 4) || ((i_subject == 4) && (i_session == 3))
            %subjid = [subjids(i_subject, :) filesep sessions(i_session, :)];
            subjid = [subjids{i_subject}];% filesep sessions{i_session}];
            %work_dir = '/Volumes/My Passport/alzheimers/derivatives/interrogation_processing/';
            work_dir = '/Volumes/My Passport/OldVarHzData/derivatives/interrogation/';
            %workdir = ['/Volumes/My Passport/tb0s_interrogation/' subjid];
            % info = dir([work_dir subjid filesep 'Step4_FasterOutput/*.set']);
            % info2 = dir([work_dir subjid filesep 'Step3_CutFiles_hfl/*.set']);
            % info3 = dir([work_dir subjid filesep 'Step4_FasterOutput_hfl/*.set']);
            % if (~isempty(info2dele) && isempty(info3))
            %     Step3_FASTERauto;
            % end
            % Step1_Merge;
            % Step2_BrushAndCutAuto_varHz;
            % Step3_FASTERauto;
            % Step4_AutoEpochRej;
            % Step5_Entrainment_ASC_InputOutput_WholeBrain_Gabor;
                            Step6_CombinedModelRFDecision;

            % try
            %     Step6_CombinedModelRFDecision;
            % catch
            %    disp("Failed step 6: " + subjid);
            % end
        %end

    end
end


%Step6_EntrainmentVIHalfHz