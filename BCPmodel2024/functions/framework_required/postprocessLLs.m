function [LL_total,LL_trials] = postprocessLLs(LL_trials_cell,LL_total,trial_nr,trl_cond_nrs,param_settings_cond,model_settings,fit_settings)
% Use content in LL_trials_cell to post-process or compute log likelihoods

% Sum the log likelihoods
LL_total = LL_total + sum(LL_trials_cell{trial_nr}.LL(:));                  %sum over multiple responses

% At the end: Create one column vector of log-likelihoods (across trials and responses)
if nargout >= 2
    LL_trials = cell2mat(reshape(cellfun(@(x) x.LL(:),LL_trials_cell,'UniformOutput',false),[numel(LL_trials_cell),1]));

    % Check for NaNs in the log-likelihoods
    assert(~any(isnan(LL_trials)),'There are NaNs in the log-likelihoods, something is likely wrong. E.g. check that NaN responses are removed successfully in compLLoneTrial.m');
end
    
end %[EoF]
