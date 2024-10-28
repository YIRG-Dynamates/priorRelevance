function S = setDefaults()
% Define the defaults parameter settings and parameter bounds. 

%N.B. All stimuli locations are in degrees visual angle
%Left = 90° 
%Right = -90°
%Straight ahead = 0°

%Parameter values for params that can be fitted
S.param_settings.sd_exp = 10;                                               %standard deviation of experimental noise
S.param_settings.cp_hazard_rate = 0.15;                                     %hazard rate
S.param_settings.sd_motor = 5;                                              %standard deviation of motor noise
S.param_settings.lapse_rate = 0.01;                                         %lapse rate 

%Model settings
S.model_settings.approximate_algorithm = true;                              %Use the simplified approximate algorithm? This avoids the computation of normcdf differences and is therefore quite a lot faster..    

S.model_settings.memory_capacity = 1;                                       %How many truncated normal distributions (0 <= integer <= inf) are maximally used to summarize posteriors as new priors.   
S.model_settings.pruning_method = 'mixture_var';                            %Choose from {'mixture_var','run_length','max_weight','forget_oldest'}
S.model_settings.decision_fun = 'model_averaging';                          %Choose from {'model_averaging','model_selection','probability_matching'}
S.model_settings.full_prior_flag = false;                                   %Set to true if you want to include the option of a changepoint in the modelled prediction responses

S.model_settings.generative_mean_range = [-90, 90];                         %Restrict the range of the generative mean
S.model_settings.response_range = [-90, 90];                                %Restrict the response range
                                                                            
                                                                            %By default we compute all latent variables
S.model_settings.latent_var_names = {'mean_priorNoCP_mu'       , 'sd_priorNoCP_mu'        , 'weight_priorNoCP_mu'   , 'runlength_priorNoCP_mu', ...
                                     'reliability_priorNoCP_mu', 'surprisal_Shannon_x'    , 'relevance_priorNoCP_mu', ...  
                                     'mean_post_mu'            , 'sd_post_mu'             , 'weight_post_mu'        , 'runlength_post_mu'};                                     

%Fitting settings
S.fit_settings.num_sim_per_trial = 1;                                       %Number of Monte Carlo simulations (increase when there is stochasticity, e.g. when introducing sensory noise)

S.fit_settings.param_log = {'sd_exp','sd_motor'};                           %Fit these parameters in log space (lower bound is 0)
S.fit_settings.param_logit = {'cp_hazard_rate','lapse_rate'};               %Fit these parameters in logit space (bounded by 0 and 1)

S.fit_settings.fit_param_names = {};                                        %The names of the parameters to fit (e.g. {'cp_hazard_rate','sd_exp'}). By default we don't fit anything. Indices here determine the param numbers.
S.fit_settings.fit_param_nrs_per_cond = cell(1);                            %The number of cells defines the number of conditions. An integer c inside a vector of cell j means that parameter c belongs to condition j.
                                                                            %E.g. {[1 2]} means that the first two parameters belong to the first (and only) condition. 

S.fit_settings.optim_MLE_or_MAP = 'MLE';                                    %Optimize parameters to obtain MLE or MAP (using trapezoid priors, see below)    
S.fit_settings.optim_tol_mesh = 1e-6;                                       %BADS converges when the parameter values don't change more than optim_tol_mesh*(PUB-PLB) 
S.fit_settings.optim_num_attempts = [1 4];                                  %Number of BADS convergence attempts [MIN MAX]. The highest log-probability solution is chosen as best out of all converged solutions.
S.fit_settings.optim_num_grid = 1e3;                                        %Number of samples that are randomly drawn from within plausible bounds to find suitable starting points for BADS

S.fit_settings.gen_predictions = false;                                     %Should we generate predictions with the fitted parameters? Logical or ...
                                                                            %Set fit_settings.gen_predictions to a vector of positive integers indicating the indices of trials you wish to compute predictions for
%Display settings
S.disp_settings.overall = false;                                            %Should we display overall predictions? (Logical)
S.disp_settings.trials = false;                                             %Should we display predictions for each trial? Logical or ...
                                                                            %Set disp_settings.trials to a vector of positive integers indicating the indices of trials you wish to display predictions for
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter bounds %%% 
%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter bounds are only used for parameters that are fitted. They are ignored for fixed parameter values.  
% The bounds can sometimes be more restrictive than the default setting because of the log/logit transformations that are used for fitted parameters.   

% Please note that the parameter bounds also define the prior probability distributions. 
% Priors are defined as trapezoids with highest probability between the plausible bounds, and linearly decreasing on either side towards the hard bounds). 

% Be realistic when changing the parameter bounds. Extraordinarily small/large bounds do not help the fitting algorithms!
% In fact, BADS convergence depends on the plausible bounds settings. The algorithm converges when the parameter values don't change more than optim_tol_mesh*(PUB-PLB).

% Default bounds:                           [Hard Lower, Plausible Lower, Plausible Upper, Hard Upper]    
S.fit_settings.bounds.sd_exp =             [     1e-2           1              50             200      ];     
S.fit_settings.bounds.sd_motor =           [     1e-2           1              50             200      ];     
S.fit_settings.bounds.cp_hazard_rate =     [     1e-6           0.01             0.99           1-1e-6 ];     
S.fit_settings.bounds.lapse_rate =         [     1e-6           1e-3             0.25           1-1e-6 ];       %These are parameter values in the non-transformed space

end %[EOF]
