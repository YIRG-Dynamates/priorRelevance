function pred_trial_struct = genPredictionsOneTrial(trial_data,param_settings,model_settings,fit_settings,trial_responses)
%Generate predictions for one trial 

%Unpack
sd_motor = param_settings.sd_motor;
lapse_rate = param_settings.lapse_rate;

mu_range = model_settings.generative_mean_range;
resp_range = model_settings.response_range;

%Generate model-predicted estimates and predictions for the last sound
[mod_pred_struct,latent_vars] = genEstimatesOneTrial(trial_data, param_settings, model_settings, fit_settings);

%Unpack structure
mod_x_pred = mod_pred_struct.x_pred(:);
mod_weights = mod_pred_struct.prob(:);
mod_weights = mod_weights ./ sum(mod_weights);

%Check the response predictions
assert(~any(isnan(mod_x_pred)),'There are NaNs in mod_x_pred, something is likely wrong');

%Save in output structure
pred_trial_struct.x_pred = mod_x_pred;
pred_trial_struct.x_pred_prob = mod_weights;

pred_trial_struct.latent_vars = latent_vars;

%%%

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
    return;
end

%Compute the posterior probability that each response was a "lapse"
%I.e. it's a lapse if it cannot be explained by the model predictions and motor noise only
lik_resps_NoLapse = sum(mod_weights.*exp(pdfTruncatedNormal(resp_x_pred,mod_x_pred,sd_motor,resp_range(1),resp_range(2),true)),1);                 %weighted average likelihood
lik_resps_Lapse = exp(pdfTruncatedConvUniformAndNormal(resp_x_pred,mu_range(1),mu_range(2),sd_motor,resp_range(1),resp_range(2),true));
pred_trial_struct.prob_lapse_resp = lapse_rate.*lik_resps_Lapse ./ (lapse_rate.*lik_resps_Lapse + (1-lapse_rate).*lik_resps_NoLapse);

end %[EoF]
