function LL_trial_struct = compLLoneTrial(trial_responses,trial_data,param_settings,model_settings,fit_settings)
%Compute log likelihood for one or more responses on this particular trial.   

%Unpack some parameters
sd_motor = param_settings.sd_motor;
lapse_rate = param_settings.lapse_rate;

mu_range = model_settings.generative_mean_range;
resp_range = model_settings.response_range;

%Check if x_pred responses are present
resp_x_pred = [];
if isfield(trial_responses,'x_pred')
    i_NaN = isnan(trial_responses.x_pred);
    if ~all(i_NaN)                                                          %Note that all([]) == 1
        resp_x_pred = trial_responses.x_pred(~i_NaN);                       %Remove the NaN responses
        resp_x_pred = reshape(resp_x_pred,[1 numel(resp_x_pred)]);          %Ensure row vector
    end
end

%Check if there's any response for this trial
if isempty(resp_x_pred)
    LL_trial_struct.LL = [];
    return;
else
    num_resp = numel(resp_x_pred);
    LL_trial_struct.LL = zeros(num_resp,1);
end

%Generate model-predicted estimates for the last sound
mod_pred_struct = genEstimatesOneTrial(trial_data, param_settings, model_settings, fit_settings);

%Unpack structure
mod_x_pred = mod_pred_struct.x_pred(:);
mod_weights = mod_pred_struct.prob(:);
mod_weights = mod_weights ./ sum(mod_weights);

%Check for NaNs
assert(~any(isnan(mod_x_pred)),'There are NaNs in mod_x_pred, something is likely wrong');

%Compute the likelihood of each response given the model predictions and random motor noise  
lik_resps_model = sum(mod_weights.*exp(pdfTruncatedNormal(resp_x_pred,mod_x_pred,sd_motor,resp_range(1),resp_range(2),true)),1);                 %weighted average likelihood

%Compute the likelihood of each response given that a random response is given, plus random motor noise   
lik_resps_lapse = exp(pdfTruncatedConvUniformAndNormal(resp_x_pred,mu_range(1),mu_range(2),sd_motor,resp_range(1),resp_range(2),true));

%Correct the likelihoods for the lapse rate
lik_resps = lapse_rate.*lik_resps_lapse + (1-lapse_rate)*lik_resps_model; 

%Finally, compute log of the likelihoods and save those
LL_trial_struct.LL(1:num_resp,1) = log(lik_resps);     

end %[EoF]
