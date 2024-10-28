function krishnamurthy_reanalysis_fitfun_2024(subj_nr,cond_nr,fit_nr,simple_flag)
%Fit the data of a single subject (subj_nr) and condition number (cond_nr) 
%where "fit_nr" determines 'memory_n', 'pruning_method', 'decision_func',
%and "simple_flag" determines whether we use the simplified model.  

fit_memory_n       = [1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4];
fit_pruning_method = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
fit_decision_func  = [1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2];

memory_n = fit_memory_n(fit_nr);
pruning_method = fit_pruning_method(fit_nr);
decision_func = fit_decision_func(fit_nr);

%Add paths
addpath('bads-master');
addpath('BCPmodel2024');

%Load behavioral data
data_path = 'final_behav_eye_data_Krishnamurthy';
load(fullfile(data_path,'subj_nrs.mat'),'subj_nrs');
subj_ID = subj_nrs{subj_nr};
load(fullfile(data_path,subj_ID,[subj_ID '_behav_eye_data.mat']),'trials_cell','trl_cond_nrs');
num_subj = 29;

%Select data of the requested experimental condition (1 (sd_exp = 10 deg) or 2 (sd_exp = 20 deg), not the training trials == 3)   
idx_rel_trl_nrs = find(trl_cond_nrs == cond_nr);
num_trials = numel(idx_rel_trl_nrs);

%Create the input data
input_data = [];
input_data.trials_cell = cell(num_trials,1);
input_data.responses = cell(num_trials,1);
for j=1:num_trials
    input_data.trials_cell{j}.x_true = trials_cell{idx_rel_trl_nrs(j)}.x;
    input_data.trials_cell{j}.mu_true = trials_cell{idx_rel_trl_nrs(j)}.mu;
    input_data.trials_cell{j}.cp = trials_cell{idx_rel_trl_nrs(j)}.cp;
    
    input_data.trials_cell{j}.cp(1) = true;   %Start each trial with a changepoint
    cp = input_data.trials_cell{j}.cp;                     
    idx_cp = find(cp);
    input_data.trials_cell{j}.SAC = (1:numel(cp))-(idx_cp(cumsum(cp))-1);
    
    input_data.trials_cell{j}.i_vis = true(1,numel(trials_cell{idx_rel_trl_nrs(j)}.x));
    input_data.trials_cell{j}.i_vis(end) = false;
    
    input_data.responses{j}.x_pred = trials_cell{idx_rel_trl_nrs(j)}.prediction;
    input_data.responses{j}.x_est = trials_cell{idx_rel_trl_nrs(j)}.estimate;
end
input_data.trl_cond_nrs = trl_cond_nrs(idx_rel_trl_nrs);
input_data.trl_cond_nrs(:) = 1; %Only one condition

%Set options structure
options_struct = [];

pruning_methods = {'mixture_var','max_weight','forget_oldest'};
decision_functions = {'model_averaging','model_selection'};

options_struct.model_settings.memory_capacity = memory_n;                      
options_struct.model_settings.pruning_method = pruning_methods{pruning_method};
options_struct.model_settings.decision_fun = decision_functions{decision_func};

options_struct.model_settings.approximate_algorithm = simple_flag;

if (memory_n==1) && (pruning_method == 3) 
    %Special case -> oldest posterior is thrown out, so prior mean is equal to last stim only, i.e. CP is always assumed (HR not used)
    if simple_flag
        options_struct.fit_settings.fit_param_names = {'sd_motor','lapse_rate'};    
        options_struct.fit_settings.fit_param_nrs_per_cond = {1:2};
    else
        options_struct.fit_settings.fit_param_names = {'sd_motor','lapse_rate','sd_exp'};    %sd_exp is used to bias the intended response location towards the centre via the expectation of the truncated normal
        options_struct.fit_settings.fit_param_nrs_per_cond = {1:3};                          %(larger sd_exp = more bias to centre)
    end
else
    options_struct.fit_settings.fit_param_names = {'sd_motor','lapse_rate','sd_exp','cp_hazard_rate'};
    options_struct.fit_settings.fit_param_nrs_per_cond = {1:4};
end

options_struct.fit_settings.optim_tol_mesh = 1e-4;                          %BADS converges when the parameter values don't change more than optim_tol_mesh*(PUB-PLB) 
options_struct.fit_settings.optim_num_attempts = [4 4];                     %Number of BADS convergence attempts [MIN MAX]. The highest log-probability solution is chosen as best out of all converged solutions.
options_struct.fit_settings.optim_num_grid = 1000;                          %Number of samples that are randomly drawn from within plausible bounds to find suitable starting points for BADS

options_struct.fit_settings.gen_predictions = false;                        %Should we generate predictions with the fitted parameters? Logical or ...
                                                                            %Set fit_settings.gen_predictions to a vector of positive integers indicating the indices of trials you wish to compute predictions for
%Display settings
options_struct.disp_settings.overall = false;                               %Should we display overall predictions? (Logical)
options_struct.disp_settings.trials = false;    

%Fit the BCP model
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['%%% SUBJECT ' num2str(subj_nr) ' out of ' num2str(num_subj) ' %%%']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%Start a timer
cStart = clock;  

%Fit parameters to the data of this subject
BCPfitResults = BCPfitModel(input_data, options_struct);

%Save results
save_path = fullfile('fitted_data_Krishna2017_2024',subj_ID);
if ~exist(save_path,'dir')
    mkdir(save_path)
end
file_name = ['BCPfitResults_Krishna2017_2024_F' num2str(fit_nr) '_' subj_ID '_C' num2str(cond_nr) '_M' num2str(memory_n) '_P' num2str(pruning_method) '_D' num2str(decision_func) '_S' num2str(double(simple_flag)) '.mat'];
save(fullfile(save_path,file_name),'BCPfitResults','-v7.3');   

%Report computation time in command window
fprintf('Finished BADS fit and saved data for subject %i (%s),\nElapsed time (days hours:minutes:seconds) %s \n', ... 
                                              subj_nr, subj_ID, datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end %[EoF]
