function [prior_mu_new,latent_vars,pred_struct] = processOneStimulus(prior_mu,x_true,param_settings,model_settings,fit_settings,latent_var_names)
% 1. Adjust prior weights according to hazard rate    
% 2. Optionally, predict the upcoming stimulus based on the prior
% 3. Compute the posterior on mu by integrating likelihood and prior 
% 4. Update the prior on mu for the next stimulus (posterior --> prior)
% 5. Optionally, collect some latent variables (output: latent_vars). 

%% 0. Unpack some useful settings

sd_exp = param_settings.sd_exp;
hazard_rate = param_settings.cp_hazard_rate;

simplified_flag = model_settings.approximate_algorithm;
mu_range = model_settings.generative_mean_range;

num_sim = fit_settings.num_sim_per_trial; 

%%  1. Adjust prior weights according to hazard rate    

prior_weights_wCP = [hazard_rate*ones(num_sim,1), (1-hazard_rate)*prior_mu.weight];         %The first column is for a changepoint

%% 2. Optionally, Predict the upcoming stimulus based on the prior

%if requested...(also compute them if latent variables are requested)
if nargout >= 2 
    
    %Ensure that the output variable exists
    pred_struct = [];
    
    %Compute the mean of the response range
    mu_range_mean = mu_range(1) + 0.5*diff(mu_range);

    %Special case for t=1 and M=0
    if size(prior_mu.mean,2) == 0

        pred_struct.x_pred = mu_range_mean*ones(num_sim,1);                     %Best guess without any sensory evidence is in middle of response range
        pred_struct.prob = ones(num_sim,1);

    else
        
        %Approximate algorithm, simply use TN's location parameters
        if simplified_flag
            prior_means = prior_mu.mean;
        else %Compute the TN's actual mean (this also ensures that prediction is within generative mean bounds)
            prior_means = TNmeanvar(mu_range(1),mu_range(2),prior_mu.mean,sqrt(prior_mu.var));
        end
        
        %Include the option of an upcoming changepoint?
        if model_settings.full_prior_flag
            prior_means = [mu_range_mean*ones(num_sim,1), prior_means];
            prior_weights = prior_weights_wCP;
        else %Ignore the option of a changepoint                                %This was found to fit data better in Nassar 2010
            prior_weights = prior_mu.weight;
        end

        %Select response based on decision function
        if strcmp(model_settings.decision_fun,'model_averaging')

            pred_struct.x_pred = sum(prior_weights.*prior_means,2);                             %weighted average
            pred_struct.prob = ones(num_sim,1);

        elseif strcmp(model_settings.decision_fun,'model_selection')

            [~,idx_max] = max(prior_weights,[],2);
            pred_struct.x_pred = prior_means(sub2ind(size(prior_means),(1:num_sim)',idx_max));  %select maximum
            pred_struct.prob = ones(num_sim,1);

        elseif strcmp(model_settings.decision_fun,'probability_matching')

            pred_struct.x_pred = prior_means;                                                   %probabilistically select all options
            pred_struct.prob = prior_weights;  
        end
        
        %Ensure that prediction responses are bound to within the range of the generative mean   
        if simplified_flag
            pred_struct.x_pred = min(max(mu_range(1),pred_struct.x_pred),mu_range(2));
        end
    end
end

%% 3. Compute the posterior on mu by integrating likelihood and prior

%Variance of the likelihood about mu is equal to the experimental noise
var_lik_mu = sd_exp.^2;

%Compute the posterior conditional on a changepoint  
mean_post_mu_cp = x_true*ones(num_sim,1);
var_post_mu_cp = var_lik_mu*ones(num_sim,1);
runlength_post_mu_cp = ones(num_sim,1);

%Special case for t=1, M=0, and hazard_rate=1
if (size(prior_mu.mean,2) == 0) || (hazard_rate == 1)
 
    mean_post_mu = mean_post_mu_cp;
    var_post_mu = var_post_mu_cp;
    runlength_post_mu = runlength_post_mu_cp;
    weight_post_mu = ones(num_sim,1);
    
else
    
    %Compute the posterior conditional on NO changepoint
    weight_priorVsLike_mu = var_lik_mu ./ (var_lik_mu + prior_mu.var);      %i.e. prior reliability "tau"
     
    mean_post_mu_noCP = weight_priorVsLike_mu.*prior_mu.mean + (1-weight_priorVsLike_mu).*x_true;
    var_post_mu_noCP = prior_mu.var.*weight_priorVsLike_mu;
    runlength_post_mu_noCP = prior_mu.runlength+1;
    
    %Concatenate the conditional posteriors (n.b. the first column is for a changepoint)
    mean_post_mu = [mean_post_mu_cp, mean_post_mu_noCP];
    var_post_mu = [var_post_mu_cp, var_post_mu_noCP];
    runlength_post_mu = [runlength_post_mu_cp, runlength_post_mu_noCP];
    
    %Also update the weights (i.e. posterior probability of the mixture components... 
    if simplified_flag
        %lik_x_cp = 1/diff(mu_range);
        %lik_x_noCP = normpdf(x_true,prior_mu.mean,sd_prior_x);
        LL_x_cp = -log(diff(mu_range))*ones(num_sim,1);
        LL_x_noCP = normlogpdf(x_true,prior_mu.mean,sqrt(prior_mu.var+sd_exp.^2));
    else
        %lik_x_cp = pdfConvUniformAndNormal(x_true,mu_range(1),mu_range(2),sd_exp);
        %lik_x_noCP = pdfConvTruncatedNormalAndNormal(x_true,mu_range(1),mu_range(2),prior_mu.mean,sqrt(prior_mu.var),sd_exp);
        LL_x_cp = pdfConvUniformAndNormal(x_true,mu_range(1),mu_range(2),sd_exp,true);
        LL_x_noCP = pdfConvTruncatedNormalAndNormal(x_true,mu_range(1),mu_range(2),prior_mu.mean,sqrt(prior_mu.var),sd_exp,true);
    end
    %weight_post_mu = prior_weights_wCP .* [lik_x_cp, lik_x_noCP];
    %weight_post_mu = weight_post_mu ./ sum(weight_post_mu,2);
    weight_post_mu = log(prior_weights_wCP) + [LL_x_cp, LL_x_noCP];
    weight_post_mu = exp(weight_post_mu - logsumexp(weight_post_mu,2));
end

%% 4. Update the prior on mu for the next stimulus (posterior --> prior)
% Impose memory constraints (prune the posterior nodes)

M = model_settings.memory_capacity;   

%Special case for likelihood-only model. I.e. always assume a changepoint.    
if M == 0
    
    prior_mu_new.mean = ones(num_sim,0);
    prior_mu_new.var = ones(num_sim,0);
    prior_mu_new.runlength = ones(num_sim,0);
    prior_mu_new.weight = ones(num_sim,0);

%Prune the nodes if they exceed the memory capacity    
elseif (M >= 1) && (size(mean_post_mu,2) > M)
    
    oldest_means = mean_post_mu(:,M:(M+1));
    oldest_vars = var_post_mu(:,M:(M+1));
    oldest_runlengths = runlength_post_mu(:,M:(M+1));
    oldest_weights = weight_post_mu(:,M:(M+1));
    
    %Normalize the weights of the two oldest nodes
    if any(sum(oldest_weights,2) == 0,'all')                                %Due to numerical limits, component weights can become so small that they are essentially zero.
        row_idx = find(sum(oldest_weights,2) == 0);                         %If all mixture components are essentially zero, then they are given equal weights here.      
        oldest_weights(row_idx,:) = ones(numel(row_idx),size(oldest_weights,2));              
    end
    oldest_weights = oldest_weights ./ sum(oldest_weights,2);
    
    %Similar to Nassar 2012 model (https://doi.org/10.1038/nn.3130)  
    %Merge oldest two nodes by weighted averaging means and combining the variances by computing the mixture variance   
    if strcmp(model_settings.pruning_method,'mixture_var')
        
        %The new prior mean and variance are equal to the mean and variance of the mixture distribution of the two nodes  
        [oldest_mean,oldest_var] = twoMomentsOfMixture(oldest_means,oldest_vars,oldest_weights,2);
        
        %Compute weighted average runlength  
        oldest_runlength = sum(oldest_weights.*oldest_runlengths,2);
        
    %Similar to Nassar 2010 model (https://doi.org/10.1523/JNEUROSCI.0822-10.2010)    
    %Merge oldest two nodes by weighted averaging means and runlengths, the variance is then based on that new runlength
    elseif strcmp(model_settings.pruning_method,'run_length')
        
        %Compute weighted average mean and runlength (described in Nassar 2010, Eq 18 and 21)
        oldest_mean = sum(oldest_weights.*oldest_means,2);
        oldest_runlength = sum(oldest_weights.*oldest_runlengths,2);
        
        %The new prior variance is based on the weighted average runlength 
        %Compare Eq. 15 in Nassar 2010. N.b. this is prior on generative mean mu, not the predictive distribution on x
        oldest_var = var_lik_mu ./ oldest_runlength;  
        
    %Similar to what was proposed by Adams & MacKay 2007 (https://doi.org/10.48550/arXiv.0710.3742)   
    %Discard the node with the lowest weight/relevance
    elseif strcmp(model_settings.pruning_method,'max_weight')
        
        [~,idx2keep] = max(oldest_weights,[],2);                            %Note that max returns "idx = 1" for "[~,idx] = max([1 1],[],2);". I.e. in case weights are equal, it selects the least old node.
        idx2keep = sub2ind(size(oldest_weights),(1:num_sim)',idx2keep);
        oldest_mean = oldest_means(idx2keep);
        oldest_var = oldest_vars(idx2keep);
        oldest_runlength = oldest_runlengths(idx2keep);
        
    %Similar to what was proposed by Skerritt-Davis & Elhilali, 2018 (https://doi.org/10.1371/journal.pcbi.1006162)  
    %Discard the oldest node
    elseif strcmp(model_settings.pruning_method,'forget_oldest')
        
        oldest_mean = mean_post_mu(:,M);
        oldest_var = var_post_mu(:,M);
        oldest_runlength = runlength_post_mu(:,M);
        
    end
    
    %Merge and save as the next prior
    prior_mu_new.mean = [mean_post_mu(:,1:(M-1)), oldest_mean];
    prior_mu_new.var = [var_post_mu(:,1:(M-1)), oldest_var];
    prior_mu_new.runlength = [runlength_post_mu(:,1:(M-1)), oldest_runlength];
    
    %Sum the weights of the last (merged) nodes
    prior_mu_new.weight = [weight_post_mu(:,1:(M-1)), sum(weight_post_mu(:,M:(M+1)),2)];

else
    
    %No memory constraint applied
    prior_mu_new.mean = mean_post_mu;
    prior_mu_new.var = var_post_mu;
    prior_mu_new.runlength = runlength_post_mu;
    prior_mu_new.weight = weight_post_mu;
end

%% 5. Optionally, collect some latent variables (output: latent_vars).

if nargout >= 2
    
    %Ensure that the output variable exists
    latent_vars = [];
    
    %PriorNoCP means
    if any(strcmp(latent_var_names,'mean_priorNoCP_mu'))
        if size(prior_mu.mean,2) == 0
            latent_vars.mean_priorNoCP_mu = mu_range_mean*ones(num_sim,1);                %Note that this is actually the prior conditional on a CP (noCP is not available)
        else
            latent_vars.mean_priorNoCP_mu = prior_mu.mean;
        end
    end
    
    %PriorNoCP sds
    if any(strcmp(latent_var_names,'sd_priorNoCP_mu'))
        if size(prior_mu.var,2) == 0
            latent_vars.sd_priorNoCP_mu = sqrt(diff(mu_range).^2/12)*ones(num_sim,1);     %SD of uniform distribution on mu_range
        else
            latent_vars.sd_priorNoCP_mu = sqrt(prior_mu.var);
        end
    end
    
    %PriorNoCP weights
    if any(strcmp(latent_var_names,'weight_priorNoCP_mu'))
        if size(prior_mu.weight,2) == 0
            latent_vars.weight_priorNoCP_mu = ones(num_sim,1);     
        else
            latent_vars.weight_priorNoCP_mu = prior_mu.weight;
        end
    end
    
    %PriorNoCP runlengths
    if any(strcmp(latent_var_names,'runlength_priorNoCP_mu'))
        if size(prior_mu.runlength,2) == 0
            latent_vars.runlength_priorNoCP_mu = zeros(num_sim,1);     
        else
            latent_vars.runlength_priorNoCP_mu = prior_mu.runlength;
        end
    end
    
    %Prior reliability, i.e. normalized weight of prior vs. likelihood
    if any(strcmp(latent_var_names,'reliability_priorNoCP_mu'))
        if (size(prior_mu.mean,2) == 0) || (hazard_rate == 1)
            latent_vars.reliability_priorNoCP_mu = zeros(num_sim,1);        %Dummy, no prior available
        else
            latent_vars.reliability_priorNoCP_mu = weight_priorVsLike_mu;
        end
    end
    
    %Surprisal (Shannon information in natural units of information: "nat")
    if any(strcmp(latent_var_names,'surprisal_Shannon_x'))
        if (size(prior_mu.mean,2) == 0) || (hazard_rate == 1)
            latent_vars.surprisal_Shannon_x = log(diff(mu_range))*ones(num_sim,1);
        else
            latent_vars.surprisal_Shannon_x = -logsumexp(log(prior_weights_wCP)+[LL_x_cp, LL_x_noCP],2);
        end
    end
    
    %Prior relevance, i.e. posterior probability of no changepoint
    if any(strcmp(latent_var_names,'relevance_priorNoCP_mu'))
        latent_vars.relevance_priorNoCP_mu = 1-weight_post_mu(:,1);
    end
    
    %Posterior mean, sd, runlength and weights of all nodes before pruning(!) 
    if any(strcmp(latent_var_names,'mean_post_mu'))
        latent_vars.mean_post_mu = mean_post_mu;
    end
    if any(strcmp(latent_var_names,'sd_post_mu'))
        latent_vars.sd_post_mu = sqrt(var_post_mu);
    end
    if any(strcmp(latent_var_names,'runlength_post_mu'))
        latent_vars.runlength_post_mu = runlength_post_mu;
    end
    if any(strcmp(latent_var_names,'weight_post_mu'))
        latent_vars.weight_post_mu = weight_post_mu;                        
    end
    
end %end of if-statement nargout >= 2 (optional output)

end %[EoF]
