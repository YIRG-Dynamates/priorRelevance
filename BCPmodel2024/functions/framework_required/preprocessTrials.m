function trials_cell = preprocessTrials(trials_cell,trl_cond_nrs,param_settings_cond,model_settings,fit_settings,responses)
% Preprocess trials_cell for all trials based on the parameter settings
% This could be useful in case of dependencies between trials

if nargin < 5
    responses = [];
end

%Perform some checks on the parameters
for c=1:fit_settings.num_conds
    
    sd_exp = param_settings_cond{c}.sd_exp;
    cp_hazard_rate = param_settings_cond{c}.cp_hazard_rate;
    sd_motor = param_settings_cond{c}.sd_motor;
    lapse_rate = param_settings_cond{c}.lapse_rate;

    assert(sd_exp > 0,'sd_exp must be greater than zero');
    assert((cp_hazard_rate >= 0) && (cp_hazard_rate <= 1),'cp_hazard_rate must be between 0 and 1');
    assert(sd_motor >= 0,'sd_motor must be greater than zero');
    assert((lapse_rate >= 0) && (lapse_rate <= 1),'lapse_rate must be between 0 and 1');
end


end %[EoF]
