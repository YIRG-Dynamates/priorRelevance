function generated_responses = genRespAllTrials(N,fitResults)
%Generate N responses for all trials.
%Output is a cell array of size [num_trials x N]

% Unpack
trials_cell = fitResults.data.trials_cell;
trl_cond_nrs = fitResults.data.trl_cond_nrs;
num_trials = numel(trials_cell);

param_settings = fitResults.settings.param_settings;
model_settings = fitResults.settings.model_settings;
fit_settings = fitResults.settings.fit_settings;

% Subdivide the param_settings struct per condition
param_settings_cond = divideParamsPerCond(param_settings,fit_settings);

% Preprocess trials_cell for all trials based on the parameter settings
trials_cell = preprocessTrials(trials_cell,trl_cond_nrs,param_settings_cond,model_settings,fit_settings);

% Loop through all the trials and collect the generated responses
generated_responses = cell(num_trials,N);
for j=1:num_trials 
    c = trl_cond_nrs(j);
    generated_responses(j,:) = genRespOneTrial(N,trials_cell{j},param_settings_cond{c},model_settings,fit_settings);
end

end %[EoF]
