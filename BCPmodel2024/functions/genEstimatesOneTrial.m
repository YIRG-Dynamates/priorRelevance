function [pred_struct, latent_vars] = genEstimatesOneTrial(trial_struct, param_settings, model_settings, fit_settings)
% Generate generative mean estimates for the last stimulus in the trial. 
% Optionally, also compute latent variables for all/only the last stimuli.

% Unpack trial_struct
x_true = trial_struct.x_true;

% Gather number of stimuli in sequence
num_stim = length(x_true); 

% Gather the number of simulations for every stimulus
num_sim = fit_settings.num_sim_per_trial;

% Initialize latent_vars cell_array
if nargout >= 2
    latent_var_names = model_settings.latent_var_names;
    latent_vars = cell(1,num_stim);
end

% Initialize a structure for the prior on the generative mean
prior_mu.mean = ones(num_sim,0);
prior_mu.var = ones(num_sim,0);
prior_mu.runlength = ones(num_sim,0);
prior_mu.weight = ones(num_sim,0);
    
% Loop through all stimuli in the sequence
for t=1:num_stim
    
    % Temporarily set the hazard rate to be 1 (first stimulus only)
    if t==1
        true_cp_hazard_rate = param_settings.cp_hazard_rate;
        param_settings.cp_hazard_rate = 1;
    end
    
    % Default with minimal computations for speed
    if nargout == 1
        
        % Update the prior on every stimulus
        if t < num_stim
            prior_mu = processOneStimulus(prior_mu,x_true(t),param_settings,model_settings,fit_settings,[]);    
            
        % Compute estimates on the last stimulus only
        elseif t == num_stim
            [prior_mu,~,pred_struct] = processOneStimulus(prior_mu,x_true(t),param_settings,model_settings,fit_settings,[]);     %Latent variables are not computed because "latent_var_names = []" input argument
        end
        
    % Also compute latent variables
    elseif nargout >= 2
        
        if t < num_stim
            [prior_mu,latent_vars{t}] = processOneStimulus(prior_mu,x_true(t),param_settings,model_settings,fit_settings,latent_var_names);
        elseif t == num_stim
            [prior_mu,latent_vars{t},pred_struct] = processOneStimulus(prior_mu,x_true(t),param_settings,model_settings,fit_settings,latent_var_names);
        end
    end
    
    % Return the true hazard rate parameter
    if t==1
        param_settings.cp_hazard_rate = true_cp_hazard_rate;
    end
    
end

end % [EoF]
