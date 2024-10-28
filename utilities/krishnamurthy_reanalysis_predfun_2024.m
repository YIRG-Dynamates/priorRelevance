function generated_responses = krishnamurthy_reanalysis_predfun_2024(subj_nr,cond_nr,M,P,D,simple_flag,full_prior_flag,HR,sd_exp)
%Create predictions for a single subject (subj_nr) and condition number (cond_nr) 

%Add paths
addpath('BCPmodel2024');

%Load behavioral data
data_path = 'final_behav_eye_data_Krishnamurthy';
load(fullfile(data_path,'subj_nrs.mat'),'subj_nrs');
subj_ID = subj_nrs{subj_nr};
load(fullfile(data_path,subj_ID,[subj_ID '_behav_eye_data.mat']),'trials_cell','trl_cond_nrs');

%Select data of the requested experimental condition (1 (sd_exp = 10 deg) or 2 (sd_exp = 20 deg), not the training trials == 3)   
idx_rel_trl_nrs = find(trl_cond_nrs == cond_nr);
num_trials = numel(idx_rel_trl_nrs);

%Create the input data
input_data = [];
input_data.trials_cell = cell(num_trials,1);
for j=1:num_trials
    input_data.trials_cell{j}.x_true = trials_cell{idx_rel_trl_nrs(j)}.x;   %Including the last A-only stimulus!
end
input_data.trl_cond_nrs = trl_cond_nrs(idx_rel_trl_nrs);
input_data.trl_cond_nrs(:) = 1; %Only one condition

%Set the number of "responses" to generate as a character vector
input_data.responses = '1';

%Set options structure
options_struct = [];

pruning_methods = {'mixture_var','max_weight','forget_oldest'};
decision_functions = {'model_averaging','model_selection'};

options_struct.model_settings.memory_capacity = M;                      
options_struct.model_settings.pruning_method = pruning_methods{P};
options_struct.model_settings.decision_fun = decision_functions{D};
options_struct.model_settings.approximate_algorithm = simple_flag;
options_struct.model_settings.full_prior_flag = full_prior_flag;

options_struct.param_settings.sd_exp = sd_exp;                              %standard deviation of experimental noise
options_struct.param_settings.cp_hazard_rate = HR;                          %hazard rate
options_struct.param_settings.sd_motor = 0;                                 %standard deviation of motor noise
options_struct.param_settings.lapse_rate = 0;                               %lapse rate 

options_struct.fit_settings.gen_predictions = false;                        %Should we generate predictions with the fitted parameters? Logical or ...
                                                                            %Set fit_settings.gen_predictions to a vector of positive integers indicating the indices of trials you wish to compute predictions for
%Display settings
options_struct.disp_settings.overall = false;                               %Should we display overall predictions? (Logical)
options_struct.disp_settings.trials = false;    

%Start a timer
cStart = clock;  

%Generate responses for the data of this subject
BCPfitResults = BCPfitModel(input_data, options_struct);  

%Report computation time in command window
fprintf('Finished generating data for subject %i (%s),\nElapsed time (days hours:minutes:seconds) %s \n', ... 
                                              subj_nr, subj_ID, datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

%Save output variable
generated_responses = BCPfitResults.generated_responses;
                                          
end %[EoF]
