function LL_trials = compLLallTrials(params,fitResults,ref_LL)
%Compute a column vector of log-likelihood values, one for each response.

% Check ref_LL input parameter (serves as a break from further computing LLs - used for speed-up during MCMC)   
if nargin < 3
    ref_LL = -inf;
end

% Unpack
trials_cell = fitResults.data.trials_cell;
trl_cond_nrs = fitResults.data.trl_cond_nrs;
responses = fitResults.data.responses;
num_trials = numel(trials_cell);

param_settings = fitResults.settings.param_settings;
model_settings = fitResults.settings.model_settings;
fit_settings = fitResults.settings.fit_settings;

% Back-transform the parameters and overwrite param_settings with the given values
params = transformParams(params,fit_settings.transforms,'trans2real');
param_settings = overwriteParams(params,param_settings,fit_settings);

% Subdivide the param_settings struct per condition
param_settings_cond = divideParamsPerCond(param_settings,fit_settings);

% Preprocess trials_cell for all trials based on the parameter settings
trials_cell = preprocessTrials(trials_cell,trl_cond_nrs,param_settings_cond,model_settings,fit_settings,responses);

% Loop through all the trials and compute the log-likelihood(s) for each
LL_total = 0;
LL_trials_cell = cell(num_trials,1);        %Using a cell-array ensures support for multiple responses per trial and for
for j=1:num_trials                          %alternative methods to compute the likelihood (see postprocessLLs function)
    
    c = trl_cond_nrs(j);
    LL_trials_cell{j} = compLLoneTrial(responses{j},trials_cell{j},param_settings_cond{c},model_settings,fit_settings);
    LL_total = postprocessLLs(LL_trials_cell,LL_total,j,trl_cond_nrs,param_settings_cond,model_settings,fit_settings);
    
    % Prematurely stop computing log-likelihoods if the current LL_total goes below some reference log-likelihood
    if LL_total < ref_LL       
        break;
    end
end

% Finally, compute the output vector LL_trials
[~,LL_trials] = postprocessLLs(LL_trials_cell,LL_total,num_trials,trl_cond_nrs,param_settings_cond,model_settings,fit_settings);

end %[EoF]
