function generated_responses_trial = genRespOneTrial(N,trial_data,param_settings,model_settings,fit_settings)
%Generate N responses for one particular trial. 

%Initialize output
generated_responses_trial = cell(1,N);

%Unpack some parameters
sd_motor = param_settings.sd_motor;
lapse_rate = param_settings.lapse_rate;

mu_range = model_settings.generative_mean_range;
resp_range = model_settings.response_range;

%Generate model-predicted estimates for the last sound
mod_pred_struct = genEstimatesOneTrial(trial_data, param_settings, model_settings, fit_settings);

%Unpack structure
mod_x_pred = mod_pred_struct.x_pred(:);
mod_weights = mod_pred_struct.prob(:);
mod_weights = mod_weights ./ sum(mod_weights);

%Loop over the responses to be generated
for j=1:N

    %Sample from the discrete probablity distribution:
    resp_x_pred = mod_x_pred(discretize(rand(1),[0; cumsum(mod_weights(:))]));
    
    %Replace lapses by a random draw from a uniform distribution on the generative mean range
    if rand(1) <= lapse_rate
        resp_x_pred = rand(1)*diff(mu_range)+mu_range(1);
    end
    
    %Add motor noise but take care of response range. I.e. sample from truncated normal distibution:  
    %https://en.wikipedia.org/wiki/Truncated_normal_distribution#Generating_values_from_the_truncated_normal_distribution
    if sd_motor > eps
        normalized_noise = norminv(normcdf((resp_range(1)-resp_x_pred)./sd_motor) + rand(1)*normcdf_diff((resp_range(1)-resp_x_pred)./sd_motor,(resp_range(2)-resp_x_pred)./sd_motor));
        resp_x_pred = resp_x_pred + sd_motor*normalized_noise;
    end
    
    %Ensure that the response falls within the response range (this should have been taken care of by the noise sampling procedure above)
    assert(~any((resp_x_pred < resp_range(1)) | (resp_x_pred > resp_range(2))),'One or more of the generated responses are outside of the response range; something is likely wrong');
    assert(~any(isnan(resp_x_pred)),'The generated responses contain NaNs; something is likely wrong');
    
    %Save responses in the cell-array
    generated_responses_trial{1,j}.x_pred = resp_x_pred;
end

end %[EoF]
