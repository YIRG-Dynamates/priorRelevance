function predictions = genPredictionsAllTrials(fitResults,fitted_params,idx_trial)
% Generate predictions for all (or some) trials and optionally plot results

%Unpack settings
param_settings = fitResults.settings.param_settings;
model_settings = fitResults.settings.model_settings;
fit_settings = fitResults.settings.fit_settings;
disp_settings = fitResults.settings.disp_settings;

%Unpack data
trials_cell = fitResults.data.trials_cell;
trl_cond_nrs = fitResults.data.trl_cond_nrs;
num_trials = numel(trials_cell);

%Deal with responses (in case of plotting results against predictions)
if isempty(fitResults.data.responses)
    responses = cell(num_trials,1);
else
    responses = fitResults.data.responses;
end

%Set idx_trial (using fit_settings.gen_predictions) to indicate indexes of trials that you wish to compute predictions for
if nargin < 3
    if islogical(fit_settings.gen_predictions)      %Don't use fit_settings.gen_predictions select all trials for predictions (even if "false", because this function was somehow called anyway..)
        idx_trial = 1:num_trials;    
    elseif isempty(fit_settings.gen_predictions)    %Special case, predictions will not be made for any trial
        idx_trial = [];
    else                                            %Use fit_settings.gen_predictions to select the trials that you wish to compute predictions for
        idx_trial = unique(fit_settings.gen_predictions);
        assert(isnumeric(idx_trial) && isvector(idx_trial),'If not a logical, then "fit_settings.gen_predictions" should be a vector of postive integers');
        assert(all((floor(idx_trial)==idx_trial) & (idx_trial>0) & (max(idx_trial)<=num_trials)),'If not a logical, then "fit_settings.gen_predictions" should be a vector of postive integers');
    end
end

%Set idx_trial_disp (using disp_settings.trials) to indicate which trials should be displayed   
if islogical(disp_settings.trials)  
    if disp_settings.trials
        idx_trial_disp = idx_trial;    
    else
        idx_trial_disp = [];  
    end
elseif isempty(disp_settings.trials)                %Special case, no trials will be displayed
    idx_trial_disp = [];
else %Use disp_settings.trials to indicate indexes of trials that should be displayed
    idx_trial_disp = unique(disp_settings.trials);
    assert(isnumeric(idx_trial_disp) && isvector(idx_trial_disp),'If not a logical, then "disp_settings.trials" should be a vector of postive integers');
    assert(all((floor(idx_trial_disp)==idx_trial_disp) & (idx_trial_disp>0) & (max(idx_trial_disp)<=num_trials)),'If not a logical, then "disp_settings.trials" should be a vector of postive integers');
    if ~isempty(setdiff(idx_trial_disp,idx_trial))
        warning('Cannot display trials that we do not compute predictions for. Some of the requested trials can therefore not be displayed');
        idx_trial_disp = intersect(idx_trial_disp,idx_trial);
    end
end

%Overwrite param_settings with the fitted 'params' 
if nargin >= 2
    if ~isempty(fitted_params)
        param_settings = overwriteParams(fitted_params,param_settings,fit_settings);
    end
end         

%Subdivide the param_settings struct per condition
param_settings_cond = divideParamsPerCond(param_settings,fit_settings);

%Preprocess trials_cell for all trials based on the parameter settings
trials_cell = preprocessTrials(trials_cell,trl_cond_nrs,param_settings_cond,model_settings,fit_settings,responses);

%Loop through all the trials and collect the predicted responses
predictions = cell(num_trials,1);
for i=idx_trial 
    c = trl_cond_nrs(i);
    predictions{i} = genPredictionsOneTrial(trials_cell{i},param_settings_cond{c},model_settings,fit_settings,responses{i});
    
    %Plot predictions for this trial
    if ismember(i,idx_trial_disp)
        plotOneTrial(trials_cell{i},responses{i},predictions{i},c,param_settings_cond{c},model_settings,fit_settings);
    end
end

%Plot predictions overview for all trials
if disp_settings.overall && ~isempty(idx_trial)
    plotallTrials(trials_cell(idx_trial),trl_cond_nrs(idx_trial),predictions(idx_trial),responses(idx_trial),param_settings_cond,model_settings,fit_settings);
end

end %[EoF]