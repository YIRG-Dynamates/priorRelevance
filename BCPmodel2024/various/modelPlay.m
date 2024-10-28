% This script demonstrates how the BCP model can be used to simulate
% responses and how to fit the model's parameters to (simulated) data.

clearvars;
close all;

%Add the model and its sub-functions to the path
model_path = fileparts(cd);                                                 %Assuming that we are in "various" folder now
addpath(genpath(model_path));    

%% Generate trials and conditions

stop_rate = 0.1;
hazard_rate = 0.15;
sd_exp = [10 20];

num_conds = 2;
num_trials_per_cond = 100;
trials_cell = cell(num_trials_per_cond,num_conds);  
for c=1:num_conds
    for j=1:num_trials_per_cond
        mu = rand(1)*180-90;         %First generative mean location (between -90 and +90)
        x = [];
        done = false;
        while ~done
            
            cp = rand(1) < hazard_rate;
            if cp ; mu = rand(1)*180-90; end
            x = [x, sd_exp(c)*randn(1)+mu];
            
            %randomly abort trial if long enough
            if (numel(x) > 8) && (rand(1) < stop_rate)
                done = true;
            end
        end
        trials_cell{j,c}.x_true = x;
    end
end
trl_cond_nrs = [1, 2] .* ones(num_trials_per_cond,1);

trials_cell = trials_cell(:);
trl_cond_nrs = trl_cond_nrs(:);
num_trials = numel(trials_cell);

%% Create input_data structure 

input_data.trials_cell = trials_cell(:);                                    %Ensure column vectors
input_data.trl_cond_nrs = trl_cond_nrs(:);                                    

%% 1. Produce some figures of predicted responses with the true parameters       

options_struct = [];

%Set some model specifics
options_struct.model_settings.approximate_algorithm = true;                 %Use the simplified approximate algorithm? This avoids the computation of normcdf differences and is therefore quite a lot faster..    

options_struct.model_settings.memory_capacity = 1;                          %How many truncated normal distributions (0 <= integer <= inf) are maximally used to summarize posteriors as new priors.   
options_struct.model_settings.pruning_method = 'mixture_var';               %Choose from {'mixture_var','run_length','max_weight','forget_oldest'}
options_struct.model_settings.decision_fun = 'model_averaging';             %Choose from {'model_averaging','model_selection','probability_matching'}
options_struct.model_settings.full_prior_flag = false;                      %Set to true if you want to include the option of a changepoint in the modelled prediction responses

%Set the parameter values
options_struct.param_settings.sd_exp = [10,20];                             %standard deviation of experimental noise for both conditions
options_struct.param_settings.cp_hazard_rate = 0.15;                        %hazard rate

options_struct.param_settings.sd_motor = 5;                                 %standard deviation of motor noise 
options_struct.param_settings.lapse_rate = 0.01;                            %lapse rate 

%Ask for predictions and figures (default = false)
options_struct.fit_settings.gen_predictions = true;                         
options_struct.disp_settings.overall = true;                                %Should we display overall predictions?             --> Not implemented yet! See plotAllTrials.m
options_struct.disp_settings.trials = true;                                 %Should we display predictions for each trial?      --> Not implemented yet! See plotOneTrial.m
   
%Call the model fitting function
BCPfitResults_1 = BCPfitModel(input_data,options_struct);  

%% 2. Simulate responses for one participant with the pre-set parameters

input_data.responses = '2';                                                 %Two simulated responses per trial 
BCPfitResults_2 = BCPfitModel(input_data,options_struct);

%Collect and merge the simulated responses for the same trials
input_data.responses = BCPfitResults_2.generated_responses(:,1);            %We'll use the simulated responses for the fits below
for j=1:num_trials
    input_data.responses{j,1}.x_pred = [input_data.responses{j,1}.x_pred, BCPfitResults_2.generated_responses{j,2}.x_pred];
end

%% 3. Call model to compute just a single log likelihood (LL)

options_struct.fit_settings.gen_predictions = false;                        %Don't create predictions (therefore also no figures), default = true 

%With the correct param settings
BCPfitResults_3a = BCPfitModel(input_data,options_struct);             
disp('LL with correct params: '); disp(BCPfitResults_3a.LL_total);

%Set different parameter settings and compute again
param_settings_backup = options_struct.param_settings;
options_struct.param_settings.sd_exp = [5,10];                              %standard deviation of experimental noise for both conditions
options_struct.param_settings.cp_hazard_rate = 0.3;                         %hazard rate

BCPfitResults_3b = BCPfitModel(input_data,options_struct);                  %The LL with wrong params should be lower than the LL with correct params
disp('LL with wrong params: '); disp(BCPfitResults_3b.LL_total);

%% 4. Fit parameters to the simulated dataset (i.e. see if we can recover the original parameters)

options_struct.fit_settings.fit_param_names = {'sd_exp','sd_exp','cp_hazard_rate','sd_motor','lapse_rate'};                       %The names of the parameters to fit 
options_struct.fit_settings.fit_param_nrs_per_cond = {[1 3 4 5],[2 3 4 5]}; %sd_exp is fit separately to each condition, while the other three params are shared across both conditions

options_struct.fit_settings.gen_predictions = true;                         %Generate predictions with the fitted parameters 
options_struct.disp_settings.overall = true;                                %Display overall results  
options_struct.disp_settings.trials = false;                                %Do not display results for each trial
   
%Fit the model
BCPfitResults_4 = BCPfitModel(input_data,options_struct);                   

disp('Comparison with the backed-up correct parameters:'); 
disp(param_settings_backup);
