function fitResults = fitModelStart(input_data,options_struct)
% Starting function of David Meijer's general modelling framework
%
% INPUTS
% * input_data 
%   A struct with the following fields:
%
%   - trials_cell (required)
%   A cell array of size [nTrials x 1]. Each cell contains (a structure 
%   with) the relevant information about the trial's characteristics. 
%
%   - responses (optional)
%   A cell array of size [nTrials x 1]. Each cell contains (a structure
%   with) the relevant information about the trial's responses.  
%   
%   If "responses" is omitted or empty, then the model cannot be fit. 
%   Instead, model predictions are computed using fixed parameter values.
%
%   Alternatively, "responses" can be set as a character array / string of
%   a positive integer N. In this case, N random responses per trial are
%   generated using the model's inference rules. 
% 
%   - trl_cond_nrs (optional)
%   A numerical vector of size [nTrials x 1] containing the condition 
%   numbers of each trial. Leave empty or use 'ones(nTrials,1)' if your 
%   dataset has only one experimental condition. The use of multiple 
%   conditions allows users to fit the same parameter once for each 
%   experimental condition, while other parameters may be shared across 
%   conditions. For further information on this, please see the commentary 
%   below for "options_struct.fit_settings.fit_param_names" and
%   "options_struct.fit_settings.fit_param_nrs_per_cond".
%
% * options_struct (optional)
%   Struct with various fields that are used as options in the model. 
%   Basically, this structure allows the user to deviate from the default
%   settings that are found in "setDefaults.m" (replacement of defaults by
%   options_struct inputs happens in "setOptions.m"). All settings for 
%   'param_settings', 'model_settings', 'fit_settings', and 'disp_settings' 
%   can be changed. The following are especially important:
%
%   - "fit_settings.fit_param_names"
%   A cell array of size [1 x nParams2Fit] with the names of the parameters 
%   to fit (one per cell). E.g. {'param_1','param_2'}. Default = {};
%
%   All parameters must be set to some default value in "setDefaults.m".
%   They will not be fit if unmentioned in "fit_settings.fit_param_names"). 
%   If you wish to assign different fixed values to them, then add them as 
%   separate fields (e.g. options_struct.param_settings.param_3 = [10 20];
%   for two conditions with different values for param_3. 
%
%   Parameter bounds must also be set in "setDefaults.m" (hard lower = LB, 
%   plausible lower = PLB, plausible upper = PUB, hard upper = UB). These
%   will be used when fitting the parameters. Any of the defaults can be 
%   changed by adding them as vectors of size [1 x 4] in 'options_struct'. 
%   E.g. set options_struct.fit_settings.bounds.param_1 = [0, 1, 3, 10];. 
%   Beware that the order should be preserved: [LB < PLB < PUB < UB].
%
%   If numel(fit_settings.fit_param_names) equals zero but "responses"
%   are present, then no parameters will be fit but one log-likelihood will
%   be computed using all of the fixed parameter values.
%
%   - "fit_settings.fit_param_nrs_per_cond"
%   A cell array of size [1 x nConditions]. Each cell should contain a
%   vector of size [1 x nParams2FitInThisCondition]. The vector in cell j
%   contains the indices of the 'fit_param_names' that belong to the j'th
%   experimental condition. For example, for two conditions, one may set
%   "fit_param_names = {'param_1','param_1','param_2'}" and 
%   "fit_param_nrs_per_cond = {[1 3],[2 3]}". This would mean that two 
%   param_1's are fit separately, using trial's of either condition only, 
%   while a single param_2 is fit using trials from both conditions. 
%   Set "fit_param_nrs_per_cond = {[1:length(fit_param_names)]} if there is
%   only one experimental condition (all params belong to condition 1).
%  
%
% OUTPUT
%
% *fitResults
% A structure with various fields containing fitting results, generated
% responses or model predictions. Please explore the output by yourself.
%
%
% N.B. Parameters are optimized using the "Bayesian Adaptive Direct Search"
% algorithm by Luigi Acerbi and Wei Ji Ma (Advances in NIPS, 2017; please 
% see https://github.com/acerbilab/bads for further information).  
%
%
% This modelling framework was developed by David Meijer during his PhD at
% the University of Birmingham, United Kingdom, and as a postdoc at the 
% Acoustics Research Institute, Austrian Academy of Sciences, Wien, Austria
% Communication: MeijerDavid1@gmail.com
%
%
% Version: 02-11-2023 


%% Assess input arguments

%Ensure that input_data.trials_cell exists and is of the right type
if ~isfield(input_data, "trials_cell")                                   
    error('input_data.trials_cell is required as input argument');
else
    assert(~isempty(input_data.trials_cell),'input_data.trials_cell must be non-empty');
    assert(iscell(input_data.trials_cell),'input_data.trials_cell must a cell array');
    if ~iscolumn(input_data.trials_cell)
        warning('input_data.trials_cell will be treated as a column vector: i.e. input_data.trials_cell = input_data.trials_cell(:);');
        input_data.trials_cell = input_data.trials_cell(:);
    end
    num_trials = length(input_data.trials_cell);
end

%Check the optional field input_data.responses
if ~isfield(input_data, "responses") || isempty(input_data.responses)
    input_data.responses = [];                    %Empty: Don't fit anything, just generate predictions
elseif ~isempty(input_data.responses)
    if ischar(input_data.responses)               %Character: Generate 'gen_N_resp' input_data.responses for a hypothetical observer
        gen_N_resp = str2double(input_data.responses);
        assert((mod(gen_N_resp,1) == 0) && (gen_N_resp > 0),'The number of input_data.responses to generate is not a valid positive integer');
    else
        %Cell array with responses for each trial: Use input_data.responses to compute log-likelihoods and fit parameters 
        assert(iscell(input_data.responses),'input_data.responses mus be a cell array');
        assert(numel(input_data.responses) == num_trials,'Cell array for input_data.responses should have one cell per trial');
        if ~isequal(size(input_data.responses), size(input_data.trials_cell))
            warning('input_data.responses will be treated as a vector: i.e. input_data.responses = input_data.responses(:);');
            input_data.responses = input_data.responses(:); 
        end
    end
end

%Check the optional field input_data.trl_cond_nrs
if ~isfield(input_data, "trl_cond_nrs") || isempty(input_data.trl_cond_nrs)
    input_data.trl_cond_nrs = ones(num_trials,1);
else
    assert(length(input_data.trl_cond_nrs) == num_trials,'The number of "input_data.trl_cond_nrs" does not match the number of trials in input_data.trials_cell');
    if ~iscolumn(input_data.trl_cond_nrs)
        warning('input_data.trl_cond_nrs will be treated as a column vector: i.e. input_data.trl_cond_nrs = input_data.trl_cond_nrs(:);');
        input_data.trl_cond_nrs = input_data.trl_cond_nrs(:);
    end
    assert(all((mod(input_data.trl_cond_nrs,1) == 0) & (input_data.trl_cond_nrs > 0)),'input_data.trl_cond_nrs must exclusively consist of positive integers');
end

%If there is no options_struct, then use all defaults
if (nargin < 2) || isempty(options_struct)        
    options_struct = struct([]);
else
    assert(isstruct(options_struct),'The 2nd input argument "options_struct" must be a structure whose fields describe non-default settings');
end

%Evaluate 'options_struct' and set default settings
S = setOptions(options_struct);

%Initialize output structure
fitResults = [];
fitResults.data = input_data;
fitResults.settings = S;

%Model-specific checks on the input data
checkInputData(fitResults);

%Shuffle the random number generator and set it to the faster algorithm
rng('shuffle','simdTwister'); 
fitResults.rng_seed = rng;                            %Save the seed

%% Perform the action!

if ischar(input_data.responses) || isempty(input_data.responses)
	
    assert(fitResults.settings.fit_settings.num_params == 0,'Cannot fit parameters without input_data.responses in the input argument');
    
    %Generate gen_N_resp input_data.responses for a hypothetical observer    
    if ischar(input_data.responses)
        fitResults.generated_responses = genRespAllTrials(gen_N_resp,fitResults); 
    
    else
        %Generate (and optionally plot) model predictions using the given/default parameters 
        if S.fit_settings.gen_predictions
            fitResults.predictions = genPredictionsAllTrials(fitResults);
        else
            warning('There are no responses to fit and "fit_settings.gen_predictions" is false, so the program does nothing');
        end
    end
    
else
    
    %Relevant responses are present, but no parameters are requested to be fitted
    if S.fit_settings.num_params == 0
        
        %Compute a vector of log-likelihood values: one per response (using the given/fixed parameter values)   
        params = [];
        fitResults.LL_trials = compLLallTrials(params,fitResults);
        fitResults.LL_total = sum(fitResults.LL_trials);
        
        %Generate (and optionally plot) model predictions using the given/default parameters 
        if S.fit_settings.gen_predictions
            fitResults.predictions = genPredictionsAllTrials(fitResults);
        end
        
    %Default fitting behaviour    
    else 
    
        %Chech that BADS has been added to the path
        if ~exist('bads','file')
            error('Please ensure that BADS (https://github.com/acerbilab/bads) is added to the Matlab path');
        else
            try 
                bads_version = bads('version');
                disp(['BADS v.' bads_version ' was found on path: ' which('bads')]);
                disp('BADS (https://github.com/acerbilab/bads) will be used to optimize parameter values.'); disp(' ');
            catch
                error('BADS (https://github.com/acerbilab/bads) function exists but behaves in unexpected manner');
            end
        end

        %Optimize parameters to find MLE or MAP
        fitResults.fit = optimBADS(fitResults);
        
        %Report fitted params per condition in command window
        for c=1:fitResults.settings.fit_settings.num_conds
            fitted_params_cond = [];
            for j=1:numel(fitResults.settings.fit_settings.fit_param_nrs_per_cond{c})
                param_index = fitResults.settings.fit_settings.fit_param_nrs_per_cond{c}(j);
                fitted_params_cond.(fitResults.settings.fit_settings.fit_param_names{param_index}) = fitResults.fit.fittedParams(param_index);
            end
            disp(''); disp(['Fitted parameter values for condition ' num2str(c) ':']); disp(fitted_params_cond);
        end
        
        %Generate (and optionally plot) model predictions using the fitted parameters
        if S.fit_settings.gen_predictions
            fitResults.predictions = genPredictionsAllTrials(fitResults,fitResults.fit.fittedParams);
        end
    end
end

end %[EoF]
