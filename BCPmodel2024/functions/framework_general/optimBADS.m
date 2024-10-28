function fit_struct = optimBADS(fitResults)
%Prepare a call to the BADS wrapper function to optimize the parameters

%% Start a timer
disp('Starting BADS fit ...');
cStart = clock;  

%% Unpack the transformed bounds
LB = fitResults.settings.fit_settings.bounds_transf.LB;
UB = fitResults.settings.fit_settings.bounds_transf.UB;
PLB = fitResults.settings.fit_settings.bounds_transf.PLB;
PUB = fitResults.settings.fit_settings.bounds_transf.PUB;

%% Create anonymous function for input to BADS
LL_fun = @(params) compLLallTrials(params,fitResults);                                      %Returns a log-likelihood vector (one LL value for each response)

MLE_or_MAP = fitResults.settings.fit_settings.optim_MLE_or_MAP;
if strcmp(MLE_or_MAP,'MLE')
    prob_fun = @(params) sum(LL_fun(params));                                               %Returns the total log-likelihood 
elseif strcmp(MLE_or_MAP,'MAP')
    prob_fun = @(params) sum(LL_fun(params))+compLogTrapPrior(params,LB,PLB,PUB,UB);        %Returns the total posterior probability (LL + logPrior)
else
    error('Unknown optimization objective method. Choose MLE or MAP.');
end
fit_fun = @(params) -prob_fun(params);                                                      %Returns the negative log probability

%% Set some options for BADS
BADSoptions = bads('defaults');
%BADSoptions.Display = 'off';

BADSoptions.TolMesh = fitResults.settings.fit_settings.optim_tol_mesh;                      %The mesh is scaled between the plausible bounds (PUB-PLB = 1). 
                                                                                            %When the values don't change more than TolMesh, then BADS calls it converged (BADS default = 1e-6).
num_attempts = fitResults.settings.fit_settings.optim_num_attempts;
num_grid_search = fitResults.settings.fit_settings.optim_num_grid;

%% Call BADS wrapper function
fit_struct = bads_wrapper(fit_fun,[],LB,UB,PLB,PUB,[],BADSoptions,num_attempts,num_grid_search);            

%% Collect the probabilities (prior, likelihood, posterior)
fittedValue = -fit_struct.fittedValue;                                                      %Note the change from negative log(prob) to positive log(prob)
fit_struct = rmfield(fit_struct,'fittedValue');
prob.logPrior = compLogTrapPrior(fit_struct.fittedParams,LB,PLB,PUB,UB);
if strcmp(MLE_or_MAP,'MLE')
    prob.logLikelihood = fittedValue;
    prob.logPosterior = fittedValue + prob.logPrior;
elseif strcmp(MLE_or_MAP,'MAP')
    prob.logLikelihood = fittedValue - prob.logPrior;
    prob.logPosterior = fittedValue;
end
fit_struct.prob = prob;

%Transform the fitted parameters back to the 'normal' domain
fit_struct.fittedParams = transformParams(fit_struct.fittedParams,fitResults.settings.fit_settings.transforms,'trans2real');
    
%Report computation time in command window
fprintf('Finished parameter optimization in (days hours:minutes:seconds) %s \n',datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end %[EoF]
