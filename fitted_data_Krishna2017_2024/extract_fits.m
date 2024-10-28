clearvars;
clc;

fitted_data_folder = 'fitted_data_Krishna2017_2024';

raw_data_folder = 'final_behav_eye_data_Krishnamurthy';
load(fullfile(raw_data_folder,'subj_nrs.mat'),'subj_nrs','num_subj');

%% Collect fitted parameters and LLs of all fits

fit_memory_n       = [1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4];
fit_pruning_method = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
fit_decision_func  = [1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2];

DimNamesCell = {'Subj','C','M','P','D','S'};
LevelNamesCell = cell(1,6);
LevelNamesCell{1} = subj_nrs';
LevelNamesCell{2} = {'1','2'};
LevelNamesCell{3} = {'1','2','3','4'};
LevelNamesCell{4} = {'1','2','3'};
LevelNamesCell{5} = {'1','2'};
LevelNamesCell{6} = {'0','1'};

num_conds = 2;
num_mems = 4;
num_pruns = 3;
num_decfns = 2;
num_simpl = 2;

LL = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);
num_trials = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);
num_params = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);

sd_motor = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);
lapse_rate = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);
sd_exp = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);
cp_hazard_rate = nan(num_subj,num_conds,num_mems,num_pruns,num_decfns,num_simpl);

for j_subj=1:num_subj 
    
    disp(['Loading subject ' num2str(j_subj) ' out of ' num2str(num_subj)]);
    
    subj_ID = subj_nrs{j_subj};
    subj_path = fullfile(fitted_data_folder,subj_ID);

    for j_cond=1:num_conds
        for j_simpl=1:num_simpl
            for j_mem=1:num_mems

                i_mem = (fit_memory_n == j_mem);

                for j_prun=1:num_pruns

                    i_prun = (fit_pruning_method == j_prun);

                    for j_decfn=1:num_decfns

                        if (j_mem==1 && j_decfn==2); continue; end  

                        i_decfn = (fit_decision_func == j_decfn);
                        fit_nr = find(i_mem & i_prun & i_decfn);
                        assert(numel(fit_nr) == 1,'Found more than one fit that meets the criteria..');

                        file_name = ['BCPfitResults_Krishna2017_2024_F' num2str(fit_nr) '_' subj_ID '_C' num2str(j_cond) '_M' num2str(j_mem) '_P' num2str(j_prun) '_D' num2str(j_decfn) '_S' num2str(j_simpl-1) '.mat'];
                        load(fullfile(subj_path,file_name),'BCPfitResults');

                        LL(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = BCPfitResults.fit.prob.logLikelihood;
                        num_trials(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = numel(BCPfitResults.data.trials_cell);
                        if (j_mem==1 && j_prun==3)
                            if j_simpl == 2
                                num_params(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = 2;
                            else
                                num_params(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = 3;
                            end
                        else
                            num_params(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = 4;
                        end

                        sd_motor(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = BCPfitResults.fit.fittedParams(1);
                        lapse_rate(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = BCPfitResults.fit.fittedParams(2);
                        
                        if (j_mem==1 && j_prun==3 && j_simpl == 2); continue; end  

                        sd_exp(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = BCPfitResults.fit.fittedParams(3);
                        
                        if (j_mem==1 && j_prun==3 && j_simpl == 1); continue; end  
                        
                        cp_hazard_rate(j_subj,j_cond,j_mem,j_prun,j_decfn,j_simpl) = BCPfitResults.fit.fittedParams(4);
                    end
                end
            end
        end
    end
end

%Correct the LL according to the Aikaike/Bayesian information criterion, then sum over both conditions. 
AIC = -2*LL + 2.*num_params;
LL_AIC_all = squeeze(sum(-.5*AIC,2));
BIC = -2*LL + log(num_trials).*num_params;
LL_BIC_all = squeeze(sum(-.5*BIC,2));

%% Save the reuslts in a separate small file

fitted_sd_motor_all = sd_motor;
fitted_lapse_rate_all = lapse_rate;                   
fitted_sd_exp_all = sd_exp;                    
fitted_cp_hazard_rate_all = cp_hazard_rate;

save(fullfile(fitted_data_folder,'fitted_params_2024.mat'),'LL_AIC_all','LL_BIC_all','fitted_sd_motor_all','fitted_lapse_rate_all','fitted_sd_exp_all','fitted_cp_hazard_rate_all');
