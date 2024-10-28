function checkInputData(fitResults)
%Do some model-specific checks on the input_data and model settings S

input_data = fitResults.data;
S = fitResults.settings;

num_trials = length(input_data.trials_cell);

%% Check necessary items in "trials_cell"
    
assert(all(cellfun(@(C) isfield(C,'x_true'),input_data.trials_cell)),'Not all cells in input_data.trials_cell contain a field x_true with the true stimulus positions');

%% Check necessary items in "responses"
if ischar(input_data.responses)
    
    %%% Simulate N responses %%%

elseif isempty(input_data.responses)
    
    %%% Generate model predictions %%%
    if S.fit_settings.gen_predictions
        
    end
    
else
    
    %%% Compute LLs %%%
    
    %Check that sd_motor is not equal to zero
    if any(S.param_settings.sd_motor <= 0)
        error('sd_motor must be larger than zero to compute likelihoods');
    end
    
    %Check that relevant responses are present to compute the log-likelihood
    assert(all(cellfun(@(C) isfield(C,'x_pred'),input_data.responses)),'"x_pred" must be present for all trials, not just some. Use NaN for missing responses');
    x_pred_all = [];
    for j=1:num_trials
        x_pred_all = [x_pred_all; input_data.responses{j}.x_pred(:)];
    end
    assert(~isempty(x_pred_all) || ~all(isnan(x_pred_all)),'All "x_pred" responses are empty or NaN! Check the input data and try again');
    
    resp_range = S.model_settings.response_range;
    assert((min(x_pred_all) >= resp_range(1)) && (max(x_pred_all) <= resp_range(2)),'One or more provided input_data.responses exceed the response range as set in setDefaults.m')
    
    %%% Fit parameters to the responses? %%%
    if S.fit_settings.num_params > 0
        
    end
    
    %%% Compute predictions too? %%%
    if S.fit_settings.gen_predictions
        
    end
end

%% Check "trl_cond_nrs" against the expected number of conditions

uniq_cond_indices = unique(input_data.trl_cond_nrs(:));
assert(isequal(uniq_cond_indices',1:S.fit_settings.num_conds),'Check the "trl_cond_nrs": some expected condition indices are either not present or some entries exceed fit_settings.num_conds');

%% Check some other settings

%%% S.param_settings


%%% S.model_settings

%Check memory capacity
M = S.model_settings.memory_capacity;
assert((M>=0) && (floor(M)==M),'Memory capacity M must be an integer >= 0 or inf');
assert(~((M==0) && (~S.model_settings.full_prior_flag)),'Memory capacity M cannot be 0 if full_prior_flag is false');

%Check pruning method
assert(ismember(S.model_settings.pruning_method,{'mixture_var','run_length','max_weight','forget_oldest'}),'Unknown pruning method');

%Check decision function
assert(ismember(S.model_settings.decision_fun,{'model_averaging','model_selection','probability_matching'}),'Unknown decision function');

%Ensure that we can form a uniform prior distribution over generative mean mu
mean_range = diff(S.model_settings.generative_mean_range);
assert((mean_range > 0) && (mean_range < inf),'Check the range for the generative mean; it has to be larger than zero and smaller than inf');


%%% S.fit_settings


end %[EoF]
