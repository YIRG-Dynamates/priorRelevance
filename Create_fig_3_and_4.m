clc
clearvars;

%Add some functions to the path
addpath('utilities');

%Load the raw data
raw_data_folder = 'final_behav_eye_data_Krishnamurthy';
load(fullfile(raw_data_folder,'subj_nrs.mat'),'subj_nrs','num_subj');

%Also load the fitted parameter values
fitted_data_folder = 'fitted_data_Krishna2017_2024';
load(fullfile(fitted_data_folder,'fitted_params_2024.mat'),'fitted_sd_exp_all','fitted_cp_hazard_rate_all');

%Set true generative parameter values
hazard_rate = 0.15;
exp_noise = [10 20];

%Collect "subj_data" (incl. modelled data)
subj_data = cell(num_subj,1);
cond_nrs = cell(num_subj,1);
num_trials = nan(num_subj,1);

for j_subj=1:num_subj
    
    %Load the data for this subject
    load(fullfile(raw_data_folder,subj_nrs{j_subj},[subj_nrs{j_subj} '_behav_eye_data.mat']),'trials_cell','trl_cond_nrs');
    
    %Remove training trials
    trials_cell(trl_cond_nrs == 3) = [];
    trl_cond_nrs(trl_cond_nrs == 3) = [];
    
    num_trials(j_subj) = numel(trials_cell);
    cond_nrs{j_subj} = trl_cond_nrs;
    
    %Generate model responses for this subject, per condition
    num_mdls = 3;                                                           %Fully Bayesian observer (1), reduced Bayesian observer model without (2) and with (3) fitted parameters
    
    generated_responses = cell(num_trials(j_subj),num_mdls);
    for j_cond = 1:2
        
        i_cond = (cond_nrs{j_subj} == j_cond);
                
        %Loop over models and generate one response per trial (without lapses or response noise)  
        for j_mdl=1:num_mdls
            if j_mdl == 1
                M = inf;    %memory capacity
                P = 1;      %WAVG pruning (irrelevant, not used)
                D = 1;      %model averaging
                simple_flag = false;
                full_prior_flag = true;
            else
                M = 1;      %memory capacity
                P = 1;      %WAVG pruning 
                D = 1;      %model averaging
                simple_flag = true;
                full_prior_flag = false;
            end
            if j_mdl == 3
                HR = fitted_cp_hazard_rate_all(j_subj,j_cond,M,P,D,simple_flag+1);
                sd_exp = fitted_sd_exp_all(j_subj,j_cond,M,P,D,simple_flag+1);
            else
                HR = hazard_rate;
                sd_exp = exp_noise(j_cond); 
            end
            generated_responses(i_cond,j_mdl) = krishnamurthy_reanalysis_predfun_2024(j_subj,j_cond,M,P,D,simple_flag,full_prior_flag,HR,sd_exp);
        end
    end
    
    %Save relevant trial data
    subj_data{j_subj} = cell(num_trials(j_subj),1);
    for j_trial=1:num_trials(j_subj)
        
        %trial conditions
        num_stim = numel(trials_cell{j_trial}.x)-1;                         %Exclude last sound (i.e. save AV stim only)
        subj_data{j_subj}{j_trial}.x = trials_cell{j_trial}.x(1:num_stim);
        subj_data{j_subj}{j_trial}.mu = trials_cell{j_trial}.mu(1:num_stim);
        subj_data{j_subj}{j_trial}.cp = trials_cell{j_trial}.cp(1:num_stim);
        
        %Correct the CP and SAC levels such that each trial starts with CP=1 and SAC=1
        subj_data{j_subj}{j_trial}.cp(1) = true;
        cp_tmp = subj_data{j_subj}{j_trial}.cp; idx_cp = find(cp_tmp);
        subj_data{j_subj}{j_trial}.SAC = (1:numel(cp_tmp))-(idx_cp(cumsum(cp_tmp))-1);
        
        %Collect responses
        subj_data{j_subj}{j_trial}.resp_subj = trials_cell{j_trial}.prediction;
        subj_data{j_subj}{j_trial}.resp_mdls = nan(1,3);
        for j_mdl=1:num_mdls
            subj_data{j_subj}{j_trial}.resp_mdls(j_mdl) = generated_responses{j_trial,j_mdl}.x_pred;
        end
    end
    clear trials_cell trl_cond_nrs
end

% Extract some useful trial conditions, per experimental condition
resp_subj = cell(num_subj,2);
resp_mdls = cell(num_subj,2,num_mdls+2);    %We will add the naive- and omniscient observer models  

last_SAC = cell(num_subj,2);
last_stim = cell(num_subj,2);               
last_mu_true = cell(num_subj,2);
last_mu_omni = cell(num_subj,2);            
last_disp = cell(num_subj,2);

nBack_runlength = cell(num_subj,2,4);       %Going up to 4 stimuli back
nBack_stim = cell(num_subj,2,4);    
nBack_mu_true = cell(num_subj,2,4);
nBack_mu_omni = cell(num_subj,2,4);
nBack_disp = cell(num_subj,2,4);

for j_subj=1:num_subj
    for j_cond=1:2
        i_cond = (cond_nrs{j_subj} == j_cond);
        
        resp_subj{j_subj,j_cond} = cellfun(@(x) x.resp_subj,subj_data{j_subj}(i_cond)); 
        for j_mdl=1:num_mdls
            resp_mdls{j_subj,j_cond,j_mdl} = cellfun(@(x) x.resp_mdls(j_mdl),subj_data{j_subj}(i_cond)); 
        end
        
        last_SAC{j_subj,j_cond} = cellfun(@(x) x.SAC(end),subj_data{j_subj}(i_cond));
        last_stim{j_subj,j_cond} = cellfun(@(x) x.x(end),subj_data{j_subj}(i_cond));
        last_mu_true{j_subj,j_cond} = cellfun(@(x) x.mu(end),subj_data{j_subj}(i_cond));
        last_mu_omni{j_subj,j_cond} = cellfun(@(x) mean(x.x((end-x.SAC(end)+1):end)),subj_data{j_subj}(i_cond));
        last_disp{j_subj,j_cond} = cellfun(@(x) x.x(end-1)-x.x(end),subj_data{j_subj}(i_cond));                                                     %preceding stim minus last stim
        
        resp_mdls{j_subj,j_cond,num_mdls+1} = last_stim{j_subj,j_cond};         %Naive observer model
        resp_mdls{j_subj,j_cond,num_mdls+2} = last_mu_omni{j_subj,j_cond};      %Omniscient observer model
        
        for j_nBack = 1:4
            nBack_runlength{j_subj,j_cond,j_nBack} = cellfun(@(x) x.SAC(end-j_nBack),subj_data{j_subj}(i_cond));
            nBack_stim{j_subj,j_cond,j_nBack} = cellfun(@(x) x.x(end-j_nBack),subj_data{j_subj}(i_cond));
            nBack_mu_true{j_subj,j_cond,j_nBack} = cellfun(@(x) x.mu(end-j_nBack),subj_data{j_subj}(i_cond));
            nBack_mu_omni{j_subj,j_cond,j_nBack} = cellfun(@(x) mean(x.x((end-j_nBack-x.SAC(end-j_nBack)+1):(end-j_nBack))),subj_data{j_subj}(i_cond));
            nBack_disp{j_subj,j_cond,j_nBack} = cellfun(@(x) x.x(end-j_nBack-1)-x.x(end-j_nBack),subj_data{j_subj}(i_cond));                        %preceding stim minus stim@nBack
        end
    end
end

num_mdls = num_mdls+2;

%Quick check on the SAC distribution. 
%figure; histogram(min(6,cat(1,last_SAC{:})));
%There are relatively few trials with SAC > 5. So, we'll use a maximum SAC of 5. 


% Start to prepare plots
subj_color = [68,119,170]/255;

mdl_names = {'Full Bayes','Default Params','Fitted Params','Naïve','Omniscient'};
mdl_colors = {[204,187,68]/255,[170,51,119]/255, [204 77 119]/255, [34,136,51]/255, [204,187,68]/255};
mdl_lineStyles = {'-.','-','--','-','-'};

LW = 1.5; %LineWidth

cond_titles = {'Experimental Noise = 10°','Experimental Noise = 20°'};

n_boot = 100;   %Number of bootstraps to obtain robust median traces at the group-level


%% Prediction responses vs omniscient predictions

%For each participant, get trace for 50% weighted percentile (weights by Gaussian kernel). 
%Then compute group-level quartiles (Q1, Q2, Q3) for each grid point.

figure('Name','Prediction responses vs omniscient predictions','WindowState', 'maximized'); 

x_grid = -90:90;
sd = 5;

for j_cond=1:2
    
    %Compute group-average traces
    median_traces_subj = squeeze(computeGroupQtraces(last_mu_omni(:,j_cond),resp_subj(:,j_cond),x_grid,sd,.5,[.25 .5 .75]))';
    
    median_traces_mdls = cell(1,num_mdls);
    for j_mdl=[1 2:4]
        median_traces_mdls{j_mdl} = computeGroupQtraces(last_mu_omni(:,j_cond),resp_mdls(:,j_cond,j_mdl),x_grid,sd,.5);
    end
    
    %Do the plots
    subplot(1,2,j_cond); hold on; box on;
    hd1 = plot([-90,90],[-90,90],'k:');  %Diagonal
    h_subj = boundedLine_DM(x_grid,median_traces_subj,subj_color);
    h_mdls = nan(1,num_mdls); 
    for j_mdl=[1 2:4]
        h_mdls(j_mdl) = plot(x_grid,median_traces_mdls{j_mdl},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW);
    end
    xlim([-90 90]); ylim([-90 90]); title(cond_titles{j_cond});
    xlabel('Omniscient Prediction'); ylabel('Prediction Response');
    legend([h_subj,h_mdls([1 2:4])],[{'Participants'},mdl_names([1 2:4])],'location','northwest');
end

%Compute correlation coefficients and regression lines
rho = nan(num_subj,2);
slopes = nan(num_subj,2);
for j_cond=1:2
    for j_subj=1:num_subj
        rho_tmp = corr([last_mu_omni{j_subj,j_cond}, resp_subj{j_subj,j_cond}]);
        rho(j_subj,j_cond) = rho_tmp(2);
        i_include = abs(last_mu_omni{j_subj,j_cond}) < 100;
        slopes_tmp = regress(resp_subj{j_subj,j_cond}(i_include),[last_mu_omni{j_subj,j_cond}(i_include), ones(sum(i_include),1)]);
        slopes(j_subj,j_cond) = slopes_tmp(1);
    end
end
quantile(rho,[.25 .5 .75])
quantile(slopes,[.25 .5 .75])

%% Absolute error per SAC level

figure('Name','Absolute error per SAC level','WindowState', 'maximized'); 

absolute_error_subj_stats = nan(num_subj,10);

for j_cond=1:2

    %Get the individual traces
    absolute_error_subj = nan(5,num_subj);
    absolute_error_mdls = cell(1,num_mdls);
    for j_mdl=2:5
        absolute_error_mdls{j_mdl} = nan(5,num_subj);
    end
    for j_subj=1:num_subj

        SACs = min(5,last_SAC{j_subj,j_cond}); %1-5
        
        abs_errors_subj = abs(resp_subj{j_subj,j_cond}-last_mu_true{j_subj,j_cond});
        absolute_error_subj(:,j_subj) = accumarray(SACs,abs_errors_subj,[],@(x) quantile(x,.50));
        for j_mdl=2:5
            abs_errors_mdl = abs(resp_mdls{j_subj,j_cond,j_mdl}-last_mu_true{j_subj,j_cond});
            absolute_error_mdls{j_mdl}(:,j_subj) = accumarray(SACs,abs_errors_mdl,[],@(x) quantile(x,.50));
        end
    end
    
    %Save absolute errors for statistical analyses
    absolute_error_subj_stats(:,(j_cond-1)*5+(1:5)) = absolute_error_subj';
    
    %Get the group-average traces via bootstrapping
    medians_subj = nan(5,3,n_boot);
    medians_mdls = cell(1,num_mdls);
    for j_mdl=2:5
        medians_mdls{j_mdl} = nan(5,1,n_boot);
    end
    for b=1:n_boot
        idxB = randsample(num_subj,num_subj,true);                                                        
        medians_subj(:,:,b) = quantile(absolute_error_subj(:,idxB),[.25 .5 .75],2);  
        for j_mdl=2:5
            medians_mdls{j_mdl}(:,:,b) = quantile(absolute_error_mdls{j_mdl}(:,idxB),.5,2);         
        end
    end
    medians_subj = mean(medians_subj,3)';
    for j_mdl=2:5
        medians_mdls{j_mdl} = mean(medians_mdls{j_mdl},3)';
    end
    
    %Do the plots
    subplot(1,2,j_cond); hold on; box on; h_mdls = nan(1,num_mdls);
    boundedLine_DM(1:5,medians_subj,subj_color);
    h_subj = plot(1:5,medians_subj(2,:),'-o','Color',subj_color,'Linewidth',LW,'Marker','o','MarkerSize',8,'MarkerFaceColor',subj_color);
    for j_mdl=2:5
        h_mdls(j_mdl) = plot(1:5,medians_mdls{j_mdl},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW,'Marker','o','MarkerSize',8,'MarkerFaceColor',mdl_colors{j_mdl});  
    end
    xlim([1 5]); xticks(1:5); xticklabels({'1','2','3','4','>=5'}); ylim([0 16]); title(cond_titles{j_cond}); 
    xlabel('Stimulus After Changepoint (SAC)'); ylabel('Absolute Error');
    legend([h_subj,h_mdls(2:5)],[{'Participants'},mdl_names(2:5)],'location','southwest');
end

%Save the data for statistics
T = array2table(absolute_error_subj_stats,'VariableNames',{'SAC1_low','SAC2_low','SAC3_low','SAC4_low','SAC5_low','SAC1_high','SAC2_high','SAC3_high','SAC4_high','SAC5_high'});
writetable(T,fullfile('stats','abs_errors.csv'));


%% Normalized error per SAC level

figure('Name','Normalized error per SAC level','WindowState', 'maximized'); 

for j_cond=1:2
    
    %Get the individuals' traces
    normalized_bias_subj = nan(5,num_subj);
    normalized_bias_mdls = cell(1,num_mdls);
    for j_mdl=2:num_mdls
        normalized_bias_mdls{j_mdl} = nan(5,num_subj);
    end
    for j_subj=1:num_subj

        runlengths = min(5,last_SAC{j_subj,j_cond}); %1-5
        disparities = last_stim{j_subj,j_cond} - last_mu_true{j_subj,j_cond};
        i_disp0 = disparities == 0;

        directed_bias_subj = resp_subj{j_subj,j_cond}-last_mu_true{j_subj,j_cond};
        normalized_bias_subj_all = directed_bias_subj(~i_disp0) ./ disparities(~i_disp0);
        normalized_bias_subj(:,j_subj) = accumarray(runlengths(~i_disp0),normalized_bias_subj_all,[],@(x) quantile(x,.50));
        
        for j_mdl=2:num_mdls
            directed_bias_mdl = resp_mdls{j_subj,j_cond,j_mdl}-last_mu_true{j_subj,j_cond};
            normalized_bias_mdl_all = directed_bias_mdl(~i_disp0) ./ disparities(~i_disp0);
            normalized_bias_mdls{j_mdl}(:,j_subj) = accumarray(runlengths(~i_disp0),normalized_bias_mdl_all,[],@(x) quantile(x,.50));
        end
    end
    
    %Perform a t-test for each SAC level, that the normalized errors are different from 1
    p = nan(1,5);
    for j_SAC = 1:5
        [p(j_SAC),~,STATS] = signtest(normalized_bias_subj(j_SAC,:)',1);
    end
    disp(p*8);                                  %Bonferroni correction
    
    %Get the group-average traces via bootstrapping
    medians_subj = nan(5,3,n_boot);
    medians_mdls = cell(1,num_mdls);
    for j_mdl=2:num_mdls
        medians_mdls{j_mdl} = nan(5,1,n_boot);
    end
    for b=1:n_boot
        idxB = randsample(num_subj,num_subj,true);                                                        
        medians_subj(:,:,b) = quantile(normalized_bias_subj(:,idxB),[.25 .5 .75],2);  
        for j_mdl=2:num_mdls
            medians_mdls{j_mdl}(:,:,b) = quantile(normalized_bias_mdls{j_mdl}(:,idxB),.5,2);            
        end
    end
    medians_subj = mean(medians_subj,3)';
    for j_mdl=2:num_mdls
        medians_mdls{j_mdl} = mean(medians_mdls{j_mdl},3)';          
    end
    
    %Do the plots
    subplot(1,2,j_cond); hold on; box on; h_mdls = nan(1,num_mdls);
    boundedLine_DM(1:5,medians_subj,subj_color);
    h_subj = plot(1:5,medians_subj(2,:),'-o','Color',subj_color,'Linewidth',LW,'Marker','o','MarkerSize',8,'MarkerFaceColor',subj_color);
    for j_mdl=2:num_mdls
        h_mdls(j_mdl) = plot(1:5,medians_mdls{j_mdl},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW,'Marker','o','MarkerSize',8,'MarkerFaceColor',mdl_colors{j_mdl});          
    end
    plot([1 5],[0 0],'k:'); ylim([0, 1.1]);
    xlim([1 5]); xticks(1:5); xticklabels({'1','2','3','4','>=5'}); title(cond_titles{j_cond}); 
    xlabel('Stimulus After Changepoint (SAC)'); ylabel('Normalized Error');
    legend([h_subj,h_mdls(2:num_mdls)],[{'Participants'}, mdl_names(2:num_mdls)],'location','southwest');
end


%% Normalized bias for SAC 1

%Absolute disparities are roughly distributed as a triangular distribution: 
%disp=0:180; pdf=2/180-disp*2/180^2; cdf=disp*2/180-disp.^2/180^2; disp=180-sqrt(180^2*(1-cdf));

sd_cdf = 0.10;

disp2cdf = @(disp) disp*2/180-disp.^2/180^2;
cdf2disp = @(cdf) 180-sqrt(180^2*(1-cdf));

max_disp = 180;
max_cdf = disp2cdf(max_disp);

x_cdf = linspace(0,max_cdf,100);
x_disp = cdf2disp(x_cdf);

x_ticks_disp = [0 10 20 30 40 50 60 80 110 180];
x_ticks_cdf = disp2cdf(x_ticks_disp);
x_tick_labels = cellfun(@(x) num2str(x),num2cell(x_ticks_disp),'UniformOutput',false);

figure('Name','Normalized bias for SAC 1','WindowState', 'maximized'); 

mdl_lineStyles{num_mdls} = '--'; %tmp because of overlap

avgIQR_traces_subj = cell(1,2);
ind_values_subj_low_high = cell(2,2);

for j_cond=1:2 
    
    %Compute the normalized biases and transform the disparities to cdf values
    disp_cdf = cell(num_subj,1);
    norm_bias_subj = cell(num_subj,1);
    norm_bias_mdls = cell(num_subj,num_mdls);
    for j_subj=1:num_subj
        
        runlengths = min(5,last_SAC{j_subj,j_cond}); %1-5
        directed_disp = nBack_mu_omni{j_subj,j_cond,1} - last_stim{j_subj,j_cond};
        i_rel = (runlengths==1) & (directed_disp ~= 0);
        
        disp_cdf{j_subj,1} = disp2cdf(abs(directed_disp(i_rel)));
        
        directed_bias_subj = resp_subj{j_subj,j_cond}(i_rel)-last_stim{j_subj,j_cond}(i_rel);
        norm_bias_subj{j_subj,1} = directed_bias_subj ./ directed_disp(i_rel);
        
        for j_mdl=2:num_mdls
            directed_bias_mdl = resp_mdls{j_subj,j_cond,j_mdl}(i_rel)-last_stim{j_subj,j_cond}(i_rel);
            norm_bias_mdls{j_subj,j_mdl} = directed_bias_mdl ./ directed_disp(i_rel);
        end
    end
    
    %Compute the group-average traces
    [group_traces,ind_traces] = computeGroupQtraces(disp_cdf,norm_bias_subj,x_cdf,sd_cdf,.5,[.25 .5 .75]);
    avgIQR_traces_subj{j_cond} = squeeze(group_traces)';
    avgIQR_traces_mdls = cell(1,num_mdls);
    for j_mdl=2:num_mdls
        avgIQR_traces_mdls{j_mdl} = computeGroupQtraces(disp_cdf,norm_bias_mdls(:,j_mdl),x_cdf,sd_cdf,.5);
    end
    
    %Collect the values to compare across conditions
    [~,idx_low] = min(abs(disp2cdf(20)-x_cdf));
    [~,idx_high] = min(abs(disp2cdf(40)-x_cdf));
    ind_values_subj_low_high{1,j_cond} = squeeze(ind_traces(1,idx_low(1),:));
    ind_values_subj_low_high{2,j_cond} = squeeze(ind_traces(1,idx_high(1),:));
    
    %Do the plots
    subplot(1,2,j_cond); hold on; box on; h_mdls=nan(1,num_mdls);
    h_subj = boundedLine_DM(x_cdf,avgIQR_traces_subj{j_cond},subj_color);
    for j_mdl=2:num_mdls
        h_mdls(j_mdl) = plot(x_cdf,avgIQR_traces_mdls{j_mdl},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW); 
    end
    
    %Overwrite last handle for legend
    h_mdls(j_mdl) = plot(-x_cdf,avgIQR_traces_mdls{j_mdl},'-','Color',mdl_colors{j_mdl},'Linewidth',LW); 
    
    plot([0 max_cdf],[1 1],'k:'); YLIM = [-.05, 0.8]; ylim(YLIM);
    plot([1 1]*disp2cdf(10),YLIM,'k:'); plot([1 1]*disp2cdf(20),YLIM,'k:'); plot([1 1]*disp2cdf(30),YLIM,'k:'); plot([1 1]*disp2cdf(40),YLIM,'k:'); 
    plot([1 1]*disp2cdf(50),YLIM,'k:'); plot([1 1]*disp2cdf(60),YLIM,'k:'); plot([1 1]*disp2cdf(80),YLIM,'k:'); plot([1 1]*disp2cdf(110),YLIM,'k:');
    plot([0 max_cdf],[.2 .2],'k:'); plot([0 max_cdf],[.4 .4],'k:'); plot([0 max_cdf],[.6 .6],'k:'); plot([0 max_cdf],[.8 .8],'k:');
    xlim([0 max_cdf]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); title(cond_titles{j_cond}); 
    xlabel('Omniscient Prediction Error'); ylabel('Normalized Bias');
    
    if j_cond==2
        legend([h_subj,h_mdls(2:num_mdls)],[{'Participants'},mdl_names(2:num_mdls)],'location','northeast');
    end
end
mdl_lineStyles{num_mdls} = '-'; %undo tmp 

%Perform a t-test on condition differences for the low and high prediction errors separately 
disp('Normalized bias at low prediction errors (median, Q1, Q3])');
disp('low noise:'); disp(quantile(ind_values_subj_low_high{1,1},[.5 .25 .75]));
disp('high noise:'); disp(quantile(ind_values_subj_low_high{1,2},[.5 .25 .75]));
[p,~,STATS1] = signrank(ind_values_subj_low_high{1,1},ind_values_subj_low_high{1,2}); disp(p);

disp('Normalized bias at high prediction errors (median, Q1, Q3])');
disp('low noise:'); disp(quantile(ind_values_subj_low_high{2,1},[.5 .25 .75]));
disp('high noise:'); disp(quantile(ind_values_subj_low_high{2,2},[.5 .25 .75]));
[p,~,STATS2] = signrank(ind_values_subj_low_high{2,1},ind_values_subj_low_high{2,2}); disp(p);


%% Signed biases towards penultimate stim as a function of noCP disparity (SAC level 2, following large CPs only)

cdf2disp_RL1 = @(cdf) 180-sqrt(180^2*(1-cdf));
cut_off_disp_RL1 = cdf2disp_RL1(0.5);   %median split of CP disparities

%Absolute disparities are roughly distributed as a half normal distribution with width sqrt(2)*sd_exp 
sd_cdf = 0.20;

max_disp = [25,50];
x_ticks_disp = [0 10 20 30 40];
x_tick_labels = cellfun(@(x) num2str(x),num2cell(x_ticks_disp),'UniformOutput',false);

% figure('Name','Signed biases towards penultimate stim as a function of noCP disparity (SAC level 2, following large CPs only)','WindowState', 'maximized'); 

%Save the individual traces for statistics afterwards
group_traces_SAC2 = cell(1,2);
ind_traces_SAC2 = cell(2,2);
avgIQR_traces_mdls_SAC2 = cell(num_mdls,2); 

for j_cond=1:2 
    
    sd_exp = sqrt(2)*exp_noise(j_cond);
    disp2cdf = @(disp) 2*(normcdf(disp,0,sd_exp)-0.5);
    cdf2disp = @(cdf) norminv(cdf/2+0.5,0,sd_exp);
    max_cdf = disp2cdf(max_disp(j_cond));
    x_cdf = linspace(0,max_cdf,100);
    x_ticks_cdf = disp2cdf(x_ticks_disp);
    
    ind_traces_SAC2{2,j_cond} = cdf2disp(x_cdf);
    
    %Get the signed biases and transform the corresponding disparities to cdf values
    disp_cdf = cell(num_subj,1);
    directed_bias_subj = cell(num_subj,1);
    directed_bias_mdls = cell(num_subj,num_mdls);
    for j_subj=1:num_subj
        
        runlengths = min(5,last_SAC{j_subj,j_cond}); %1-5
        directed_disp = last_disp{j_subj,j_cond};
        abs_disp_CP = abs(nBack_mu_omni{j_subj,j_cond,2} - nBack_stim{j_subj,j_cond,1});
        i_rel = (runlengths==2) & (directed_disp ~= 0) & (abs_disp_CP > cut_off_disp_RL1);                         %runlengths==2 ! 
        disp_cdf{j_subj,1} = disp2cdf(abs(directed_disp(i_rel)));
        
        directed_bias_subj{j_subj,1} = sign(directed_disp(i_rel)).*(resp_subj{j_subj,j_cond}(i_rel)-last_stim{j_subj,j_cond}(i_rel));
        for j_mdl=2:num_mdls
            directed_bias_mdls{j_subj,j_mdl} = sign(directed_disp(i_rel)).*(resp_mdls{j_subj,j_cond,j_mdl}(i_rel)-last_stim{j_subj,j_cond}(i_rel));
        end
    end
    
    %Compute the individuals' IQR traces
    [group_traces_SAC2{j_cond},ind_traces] = computeGroupQtraces(disp_cdf,directed_bias_subj,x_cdf,sd_cdf,.5,[.25 .5 .75]);
    ind_traces_SAC2{1,j_cond} = squeeze(ind_traces)';
    avgIQR_traces_subj = squeeze(group_traces_SAC2{j_cond})';
    for j_mdl=2:num_mdls
        avgIQR_traces_mdls_SAC2{j_mdl,j_cond} = computeGroupQtraces(disp_cdf,directed_bias_mdls(:,j_mdl),x_cdf,sd_cdf,.5);
    end
    
%     %Do the plots
%     subplot(1,2,j_cond); hold on; box on; h_mdls = nan(1,num_mdls);
%     h_subj = boundedLine_DM(x_cdf,avgIQR_traces_subj,subj_color); 
%     for j_mdl=2:num_mdls
%         h_mdls(j_mdl) = plot(x_cdf,avgIQR_traces_mdls_SAC2{j_mdl,j_cond},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW);
%     end
%     YLIM = [-1,10]; ylim(YLIM);
%     plot([1 1]*disp2cdf(10),YLIM,'k:'); plot([1 1]*disp2cdf(20),YLIM,'k:'); plot([1 1]*disp2cdf(30),YLIM,'k:'); plot([1 1]*disp2cdf(40),YLIM,'k:');
%     plot([1 1]*disp2cdf(50),YLIM,'k:'); plot([1 1]*disp2cdf(60),YLIM,'k:'); plot([1 1]*disp2cdf(80),YLIM,'k:'); plot([1 1]*disp2cdf(110),YLIM,'k:');
%     xlim([0 max_cdf]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); title(cond_titles{j_cond}); 
%     xlabel('noCP disparity: abs(penultimate- minus last-stim)'); ylabel('Group median (Q1-Q3) of local median of signed biases towards penultimate stim from last stim');
%     legend([h_subj,h_mdls(2:num_mdls)],[{'Participants'},mdl_names(2:num_mdls)],'location','northwest');
end


% Signed bias towards omniscient prior as a function of last noCP disparity (SAC level 3, small preceding disparities and large CP-disparity only)

cdf2disp_RL1 = @(cdf) 180-sqrt(180^2*(1-cdf));
cut_off_disp_RL1 = cdf2disp_RL1(0.5);   %median split of CP disparities

%Absolute disparities are roughly distributed as a half normal distribution with width sqrt(2)*sd_exp 
sd_cdf = 0.20;

max_disp = [25,50];
x_ticks_disp = [0 10 20 30 40];
x_tick_labels = cellfun(@(x) num2str(x),num2cell(x_ticks_disp),'UniformOutput',false);

% figure('Name','Signed bias towards omniscient prior as a function of last noCP disparity (SAC level 3, small preceding disparities and large CP-disparity only)','WindowState', 'maximized'); 

%Save the individual traces for statistics afterwards
group_traces_SAC3 = cell(1,2);
ind_traces_SAC3 = cell(2,2);
avgIQR_traces_mdls_SAC3 = cell(num_mdls,2); 

for j_cond=1:2 
    
    %Absolute disparities are roughly distributed as a half normal distribution, where the width depends on the runlength:
    %RL=2: sqrt(exp_var/1 + exp_var = 2.00*exp_var) = sqrt(2/1)*sd_exp;
    %RL=3: sqrt(exp_var/2 + exp_var = 1.50*exp_var) = sqrt(3/2)*sd_exp;
    %RL=4: sqrt(exp_var/3 + exp_var = 1.33*exp_var) = sqrt(4/3)*sd_exp;
    %RL=5: sqrt(exp_var/4 + exp_var = 1.25*exp_var) = sqrt(5/4)*sd_exp;
    %RL=6: sqrt(exp_var/5 + exp_var = 1.20*exp_var) = sqrt(6/5)*sd_exp;
    sd_exp = sqrt(3/2)*exp_noise(j_cond);
    disp2cdf = @(disp) 2*(normcdf(disp,0,sd_exp)-0.5);
    cdf2disp = @(cdf) norminv(cdf/2+0.5,0,sd_exp);
    max_cdf = disp2cdf(max_disp(j_cond));
    x_cdf = linspace(0,max_cdf,100);
    x_ticks_cdf = disp2cdf(x_ticks_disp);
    
    cdf2disp_RL2 = @(cdf) norminv(cdf/2+0.5,0,sqrt(2)*exp_noise(j_cond));
    
    ind_traces_SAC3{2,j_cond} = cdf2disp(x_cdf);
    
    %Compute the signed biases
    disp_cdf = cell(num_subj,1);
    signed_bias_subj = cell(num_subj,1);
    signed_bias_mdls = cell(num_subj,num_mdls);
    for j_subj=1:num_subj
        
        runlengths = min(5,last_SAC{j_subj,j_cond}); %1-5
        
        disp_CP_RL1 = abs(nBack_mu_omni{j_subj,j_cond,3} - nBack_stim{j_subj,j_cond,2});        
        disp_noCP_RL2 = abs(nBack_stim{j_subj,j_cond,2} - nBack_stim{j_subj,j_cond,1});    
        directed_disp = nBack_mu_omni{j_subj,j_cond,1} - last_stim{j_subj,j_cond};
        
        i_rel = (runlengths==3) & (directed_disp ~= 0) & (disp_CP_RL1 > cut_off_disp_RL1) & (disp_noCP_RL2 < cdf2disp_RL2(2/3));                         %runlengths==3 !        
        disp_cdf{j_subj,1} = disp2cdf(abs(directed_disp(i_rel)));
        
        directed_bias_subj = resp_subj{j_subj,j_cond}(i_rel)-last_stim{j_subj,j_cond}(i_rel);
        signed_bias_subj{j_subj,1} = sign(directed_disp(i_rel)).*directed_bias_subj;
        
        for j_mdl=2:num_mdls
            directed_bias_mdl = resp_mdls{j_subj,j_cond,j_mdl}(i_rel)-last_stim{j_subj,j_cond}(i_rel);
            signed_bias_mdls{j_subj,j_mdl} = sign(directed_disp(i_rel)).*directed_bias_mdl;
        end
    end
    
    %Get the group-average traces
    [group_traces_SAC3{j_cond},ind_traces] = computeGroupQtraces(disp_cdf,signed_bias_subj,x_cdf,sd_cdf,.5,[.25 .5 .75]);
    ind_traces_SAC3{1,j_cond} = squeeze(ind_traces)';
    avgIQR_traces_subj = squeeze(group_traces_SAC3{j_cond})';
    for j_mdl=2:num_mdls
        avgIQR_traces_mdls_SAC3{j_mdl,j_cond} = computeGroupQtraces(disp_cdf,signed_bias_mdls(:,j_mdl),x_cdf,sd_cdf,.5);
    end
    
%     %Do the plots
%     subplot(1,2,j_cond); hold on; box on; h_mdls=nan(1,num_mdls);
%     h_subj = boundedLine_DM(x_cdf,avgIQR_traces_subj,subj_color); 
%     for j_mdl=2:num_mdls
%         h_mdls(j_mdl) = plot(x_cdf,avgIQR_traces_mdls_SAC3{j_mdl,j_cond},mdl_lineStyles{j_mdl},'Color',mdl_colors{j_mdl},'Linewidth',LW); 
%     end
%     YLIM = [-1, 10]; ylim(YLIM);
%     for j_tick=1:numel(x_ticks_disp)
%         plot([1 1]*x_ticks_cdf(j_tick),YLIM,'k:');
%     end 
%     xlim([0 max_cdf]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); title(cond_titles{j_cond}); 
%     xlabel('Last noCP disparity: abs(omniscient-prior minus last-stim)'); ylabel('Group median of median bias towards omniscient prior (1) from last stim (0)');
%     legend([h_subj,h_mdls(2:num_mdls)],[{'Participants'},mdl_names(2:num_mdls)],'location','northwest');
end


% Plot SAC 2 and SAC 3 on top of eachother

% Start to prepare plots
%subj_colors2 = {[102,153,204]/255; [0,68,136]/255};                        %light,dark (blue)
subj_colors2 = {[51,187,238]/255; [0,68,136]/255};                          %light,dark (blue)

mdl_colors2 = {[],[];...                                                    %empty (Full Bayes)
               [194,165,207]/255,[118 42 131]/255; ...                      %purple (Default)
               [238,153,170]/255,[153 68 85]/255; ...                       %red (Fitted, dashed)
               [34,136,51]/255,[34,136,51]/255; ...                         %green (naive)
               [238,204,102]/255,[153 119 0]/255}; ...                      %yellow (omniscient)

%mdl_names = {'Full Bayes','Default','Fitted','Naive','Omniscient'};        %Just a reminder
%mdl_lineStyles = {'--','-','--','-','-'};

figure('Name','Signed Bias for SAC 2 and SAC 3','WindowState', 'maximized'); 

for j_cond=1:2 
    
    sd_exp = sqrt(1+2/3)*exp_noise(j_cond);                                 %Halfway between N=1 and N=2, i.e. N=1.5, where, 2/3 is N^-1
    disp2cdf = @(disp) 2*(normcdf(disp,0,sd_exp)-0.5);
    cdf2disp = @(cdf) norminv(cdf/2+0.5,0,sd_exp);
    max_cdf = disp2cdf(max_disp(j_cond));
    x_cdf = linspace(0,max_cdf,100);
    x_ticks_cdf = disp2cdf(x_ticks_disp);
    
    %Do the plotting
    subplot(1,2,j_cond); hold on; box on; h_subj=nan(1,2); h_mdls=nan(1,num_mdls*2);
    
    %SAC 2
    x_cdf_SAC2 = disp2cdf(ind_traces_SAC2{2,j_cond});
    h_subj(1) = boundedLine_DM(x_cdf_SAC2,squeeze(group_traces_SAC2{j_cond})',subj_colors2{1}); 
    for j_mdl=2:num_mdls
        h_mdls(j_mdl) = plot(x_cdf_SAC2,avgIQR_traces_mdls_SAC2{j_mdl,j_cond},mdl_lineStyles{j_mdl},'Color',mdl_colors2{j_mdl,1},'Linewidth',LW); 
    end
    
    %SAC 3
    x_cdf_SAC3 = disp2cdf(ind_traces_SAC3{2,j_cond});
    h_subj(2) = boundedLine_DM(x_cdf_SAC3,squeeze(group_traces_SAC3{j_cond})',subj_colors2{2}); 
    for j_mdl=2:num_mdls
        h_mdls(num_mdls+j_mdl) = plot(x_cdf_SAC3,avgIQR_traces_mdls_SAC3{j_mdl,j_cond},mdl_lineStyles{j_mdl},'Color',mdl_colors2{j_mdl,2},'Linewidth',LW); 
    end
    
    YLIM = [-1, 10]; ylim(YLIM);
    for j_tick=1:numel(x_ticks_disp)
        plot([1 1]*x_ticks_cdf(j_tick),YLIM,'k:');
    end 
    xlim([0 max_cdf]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); title(cond_titles{j_cond}); 
    xlabel('Omniscient Prediction Error'); ylabel('Signed Bias');
    
    if j_cond==1
        legend([h_subj(2) h_subj(1),h_mdls([num_mdls+2, 2, num_mdls+3, 3, 4, num_mdls+5, 5])], ...
            {'Participants SAC 3','Participants SAC 2','Default Params SAC 3','Default Params SAC 2', ...
             'Fitted Params SAC 3','Fitted Params SAC 2','Naive','Omniscient SAC 3','Omniscient SAC 2'},'location','northwest');
    end
    
end

%% Check for differences between SAC2 and SAC 3 at peak and end of traces

peaks_disp = [10 20];

peaks_bias = cell(2,2);
for j=1:4; peaks_bias{j} = nan(29,1); end

for j_subj=1:num_subj
    [~,idx] = min(abs(ind_traces_SAC2{2,1}-peaks_disp(1)));
    peaks_bias{1,1}(j_subj) = ind_traces_SAC2{1,1}(j_subj,idx(1));          %SAC 2 low noise
    [~,idx] = min(abs(ind_traces_SAC2{2,2}-peaks_disp(2)));
    peaks_bias{1,2}(j_subj) = ind_traces_SAC2{1,2}(j_subj,idx(1));          %SAC 2 high noise
    [~,idx] = min(abs(ind_traces_SAC3{2,1}-peaks_disp(1)));
    peaks_bias{2,1}(j_subj) = ind_traces_SAC3{1,1}(j_subj,idx(1));          %SAC 3 low noise
    [~,idx] = min(abs(ind_traces_SAC3{2,2}-peaks_disp(2)));
    peaks_bias{2,2}(j_subj) = ind_traces_SAC3{1,2}(j_subj,idx(1));          %SAC 3 high noise
end

%Perform statistics on peak and on difference towards end
disp('Peak bias at low noise (median, Q1, Q3])');
disp('SAC2:'); disp(quantile(peaks_bias{1,1},[.5 .25 .75]));
disp('SAC3:'); disp(quantile(peaks_bias{2,1},[.5 .25 .75]));
[p,~,STATS1] = signrank(peaks_bias{1,1},peaks_bias{2,1}); disp(p);

disp('Peak bias at high noise (median, Q1, Q3])');
disp('SAC2:'); disp(quantile(peaks_bias{1,2},[.5 .25 .75]));
disp('SAC3:'); disp(quantile(peaks_bias{2,2},[.5 .25 .75]));
[p,~,STATS2] = signrank(peaks_bias{1,2},peaks_bias{2,2}); disp(p);


%% Simulate bias for reduced bayesian observer model with true params

H = 0.15;

sd_exp = 20;
sd_prior = sd_exp./sqrt(2);

pred_error = 0:180;

sd_pred = sqrt((sd_exp^2 + sd_prior^2));
prior_reliability = sd_exp^2 / sd_pred^2;

Q = log(((1/H-1).*180)./(sqrt(2*pi).*sd_pred))-.5.*pred_error.^2./sd_pred^2;
prior_relevance = 1 ./ (1 + exp(-Q));

norm_bias = prior_reliability.*prior_relevance;
bias = norm_bias.*pred_error;

figure;
subplot(2,1,1); box on; hold on; 
plot(pred_error,norm_bias,'m');
subplot(2,1,2); box on; hold on; 
plot(pred_error,bias,'m');

%peak bias is reached at:
%- 29 for SAC 2 and low noise ('b')
%- 26 for SAC 3 and low noise ('r')
%- 53 for SAC 2 and high noise ('c')
%- 47 for SAC 3 and high noise ('m')

%% Compute R-squared values for the models

mdl_nrs = 2:4;
R2_regular = cell(3,2);     %3 models (Default, Fitted, Naive) and 2 conditions
for j=1:6; R2_regular{j} = nan(num_subj,1); end

for j_cond=1:2
    for j_subj=1:num_subj
        SS_tot = sum((resp_subj{j_subj,j_cond}-mean(resp_subj{j_subj,j_cond})).^2);
        for j_mdl=1:3
            SS_res = sum((resp_subj{j_subj,j_cond}-resp_mdls{j_subj,j_cond,mdl_nrs(j_mdl)}).^2);
            R2_regular{j_mdl,j_cond}(j_subj) = 1-SS_res/SS_tot;
        end
    end
end

disp('regular R2 for Default model:');
disp('Low noise:'); disp(quantile(R2_regular{1,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_regular{1,2},[.5 .25 .75]));

disp('regular R2 for Fitted model:');
disp('Low noise:'); disp(quantile(R2_regular{2,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_regular{2,2},[.5 .25 .75]));

disp('regular R2 for Naive model:');
disp('Low noise:'); disp(quantile(R2_regular{3,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_regular{3,2},[.5 .25 .75]));

disp('Fitted vs Default (regular R2)');
disp('Low noise:'); [p,~,STATS1] = signrank(R2_regular{2,1},R2_regular{1,1}); disp(p);
disp('High noise:'); [p,~,STATS2] = signrank(R2_regular{2,2},R2_regular{1,2}); disp(p);

disp('Fitted vs Naive (regular R2)');
disp('Low noise:'); [p,~,STATS3] = signrank(R2_regular{2,1},R2_regular{3,1}); disp(p);
disp('High noise:'); [p,~,STATS4] = signrank(R2_regular{2,2},R2_regular{3,2}); disp(p);

mdl_nrs = 2:4;
R2_adjusted = cell(3,2);     %3 models (Default, Fitted, Naive) and 2 conditions
for j=1:3; R2_adjusted{j} = nan(num_subj,1); end

for j_cond=1:2
    for j_subj=1:num_subj
        errors_subj_tmp = resp_subj{j_subj,j_cond} - last_mu_true{j_subj,j_cond};
        SS_tot = sum((errors_subj_tmp-mean(errors_subj_tmp)).^2);
        for j_mdl=1:3
            errors_mdl_tmp = resp_mdls{j_subj,j_cond,mdl_nrs(j_mdl)} - last_mu_true{j_subj,j_cond};
            SS_res = sum((errors_subj_tmp-errors_mdl_tmp).^2);
            R2_adjusted{j_mdl,j_cond}(j_subj) = 1-SS_res/SS_tot;
        end
    end
end

disp('adjusted R2 for Default model:');
disp('Low noise:'); disp(quantile(R2_adjusted{1,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_adjusted{1,2},[.5 .25 .75]));

disp('adjusted R2 for Fitted model:');
disp('Low noise:'); disp(quantile(R2_adjusted{2,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_adjusted{2,2},[.5 .25 .75]));

disp('adjusted R2 for Naive model:');
disp('Low noise:'); disp(quantile(R2_adjusted{3,1},[.5 .25 .75]));
disp('High noise:'); disp(quantile(R2_adjusted{3,2},[.5 .25 .75]));

disp('Fitted vs Default (adjusted R2)');
disp('Low noise:'); p = signrank(R2_adjusted{2,1},R2_adjusted{1,1}); disp(p);
disp('High noise:'); p = signrank(R2_adjusted{2,2},R2_adjusted{1,2}); disp(p);

disp('Fitted vs Naive (adjusted R2)');
disp('Low noise:'); p = signrank(R2_adjusted{2,1},R2_adjusted{3,1}); disp(p);
disp('High noise:'); p = signrank(R2_adjusted{2,2},R2_adjusted{3,2}); disp(p);

