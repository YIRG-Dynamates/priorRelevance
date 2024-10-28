clearvars;

%Add paths
addpath('BCPmodel2024');

data_path = 'final_behav_eye_data_Krishnamurthy';
load(fullfile(data_path,'subj_nrs.mat'),'subj_nrs');

num_sds = 25;
num_HRs = 25;
SDs = exp(linspace(1,5.5,num_sds));
HRs = linspace(.02,.98,num_HRs);

num_conds = 2;

%% Compute the Surprisal surface values

surprisal = nan(num_sds,num_HRs,num_conds);

for j_SD = 1:num_sds
    
    disp(j_SD);
    
    for j_HR = 1:num_HRs
        
        disp(j_HR);
        
        surprisal_tmp1 = nan(29,2);             %Separately for the two conditions
        
        %Start a timer
        cStart = clock;  
        
        %Loop over subjects and conditions
        for subj_nr=1:29

            %Load behavioral data
            subj_ID = subj_nrs{subj_nr};
            load(fullfile(data_path,subj_ID,[subj_ID '_behav_eye_data.mat']),'trials_cell','trl_cond_nrs');

            for cond_nr=1:2

                %Select data of the requested experimental condition (1 (sd_exp = 10 deg) or 2 (sd_exp = 20 deg), not the training trials == 3)   
                idx_rel_trl_nrs = find(trl_cond_nrs == cond_nr);
                num_trials = numel(idx_rel_trl_nrs);

                %Create the input data
                input_data = [];
                input_data.trials_cell = cell(num_trials,1);
                for j=1:num_trials
                    input_data.trials_cell{j}.x_true = trials_cell{idx_rel_trl_nrs(j)}.x(1:(end-1));    %Do not include the last auditory stimulus
                end

                %Set options structure
                options_struct = [];

                options_struct.param_settings.sd_exp = SDs(j_SD);           %standard deviation of experimental noise
                options_struct.param_settings.cp_hazard_rate = HRs(j_HR);   %hazard rate

                options_struct.model_settings.memory_capacity = 1;
                options_struct.model_settings.pruning_method = 'mixture_var'; 
                options_struct.model_settings.decision_fun = 'model_averaging';             %Choose from {'model_averaging','model_selection','probability_matching'

                options_struct.fit_settings.gen_predictions = true;                         %Should we generate predictions with the fitted parameters? 
                options_struct.model_settings.latent_var_names = {'surprisal_Shannon_x'};   %Only compute surprisal as a latent variable
                
                %Display settings
                options_struct.disp_settings.overall = false;                               %Don't plot anything
                options_struct.disp_settings.trials = false;    

                %Create predictions for all trials
                BCPfitResults = BCPfitModel(input_data, options_struct);
                
                %Collect the surprisal values of each trial and sum them
                surprisal_tmp2 = nan(num_trials,1);
                for j_trial=1:num_trials
                    surprisal_tmp2(j_trial,1) = sum(cellfun(@(x) x.surprisal_Shannon_x,BCPfitResults.predictions{j_trial,1}.latent_vars));  %Sum over stimuli
                end
                surprisal_tmp1(subj_nr,cond_nr) = sum(surprisal_tmp2);                                                                      %Sum over trials
            end
        end
        
        %Also sum across subjects
        surprisal(j_SD,j_HR,:) = sum(surprisal_tmp1,1);
        
        %Display computation time
        disp(['Finished: SD ' num2str(j_SD) '/' num2str(num_sds), ' HR ' num2str(j_HR) '/' num2str(num_HRs) '. Elapsed time: ' datestr(etime(clock,cStart)/86400,'dd HH:MM:SS')]);
    end
end

%% Save the results

save('surprisal_surface_plot_data.mat','surprisal','SDs','HRs');
