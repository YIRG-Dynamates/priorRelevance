clearvars;
clc;

addpath(genpath('BCPmodel2024'));
addpath('utilities');

fitted_data_folder = 'fitted_data_Krishna2017_2024';
load(fullfile(fitted_data_folder,'fitted_params_2024.mat'),'LL_BIC_all');

num_subj = 29;

%% Fixed effects analysis of all models

lme = LL_BIC_all;       %1=num_subj(29), 2=num_mems(4), 3=num_pruns(3), 4=num_decfns(2), 5=num_simpl(2)

lme_P1_D1_S0 = squeeze(lme(:,:,1,1,1));      %4 models
lme_P2_D1_S0 = squeeze(lme(:,:,2,1,1));      %4 models
lme_P3_D1_S0 = squeeze(lme(:,:,3,1,1));      %4 models
lme_P1_D2_S0 = squeeze(lme(:,2:4,1,2,1));    %3 models - but we also share M1
lme_P2_D2_S0 = squeeze(lme(:,2:4,2,2,1));    %3 models - but we also share M1
lme_P3_D2_S0 = squeeze(lme(:,2:4,3,2,1));    %3 models - but we also share M1

lme_P1_D1_S1 = squeeze(lme(:,:,1,1,2));      %4 models
lme_P2_D1_S1 = squeeze(lme(:,:,2,1,2));      %4 models
lme_P3_D1_S1 = squeeze(lme(:,:,3,1,2));      %4 models
lme_P1_D2_S1 = squeeze(lme(:,2:4,1,2,2));    %3 models - but we also share M1
lme_P2_D2_S1 = squeeze(lme(:,2:4,2,2,2));    %3 models - but we also share M1
lme_P3_D2_S1 = squeeze(lme(:,2:4,3,2,2));    %3 models - but we also share M1

lme_all = [lme_P1_D1_S0,lme_P2_D1_S0,lme_P3_D1_S0,lme_P1_D2_S0,lme_P2_D2_S0,lme_P3_D2_S0, ...
           lme_P1_D1_S1,lme_P2_D1_S1,lme_P3_D1_S1,lme_P1_D2_S1,lme_P2_D2_S1,lme_P3_D2_S1];

num_models = size(lme_all,2);
ref_model_nr = 22;

%Bootstrap the mean differences with the reference model
rng('shuffle');
num_bootstrap = 10000;
diff_with_ref = lme_all-repmat(lme_all(:,ref_model_nr),[1 num_models]);     %Subtract the lme of the reference model from all other models
bootstrap_sum_diffs = nan(num_bootstrap,num_models);
for b=1:num_bootstrap
    idxB = randsample(num_subj,num_subj,true);                              %Pick at random nSubjects subjIDs (with replacement) 
    bootstrap_sum_diffs(b,:) = sum(diff_with_ref(idxB,:),1);                %Sum the differences (fixed effects analysis)
end
median_sum_diffs = quantile(bootstrap_sum_diffs,.5,1);                      %Compute the median across bootstraps

%Plot the results
colour_sets = [238,5,21; 181,101,29; 144,201,135; ...
               178,24,43; 123,80,25; 17,85,34]/255;
lineStyles2 = {'-','--'};

figure('Name','Model comparison','WindowState', 'maximized'); 

subplot(1,2,1); 
cla; hold on; box on; h = nan(1,12);
h(1) = plot(1:4,median_sum_diffs(1,21+(1:4)),['o' lineStyles2{1}],'Color',colour_sets(1,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(1,:));
h(2) = plot(1:4,median_sum_diffs(1,21+(5:8)),['o' lineStyles2{1}],'Color',colour_sets(2,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(2,:));
h(3) = plot(1:4,median_sum_diffs(1,21+(9:12)),['o' lineStyles2{1}],'Color',colour_sets(3,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(3,:));
h(4) = plot(1:4,[median_sum_diffs(1,21+(1)), median_sum_diffs(1,21+(13:15))],['o' lineStyles2{1}],'Color',colour_sets(4,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(4,:));
h(5) = plot(1:4,[median_sum_diffs(1,21+(5)), median_sum_diffs(1,21+(16:18))],['o' lineStyles2{1}],'Color',colour_sets(5,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(5,:));
h(6) = plot(1:4,[median_sum_diffs(1,21+(9)), median_sum_diffs(1,21+(19:21))],['o' lineStyles2{1}],'Color',colour_sets(6,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(6,:));
plot(1,median_sum_diffs(1,21+(1)),'o','Color',colour_sets(4,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(4,:),'MarkerEdgeColor',colour_sets(1,:));
plot(1,median_sum_diffs(1,21+(5)),'o','Color',colour_sets(5,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(5,:),'MarkerEdgeColor',colour_sets(2,:));
plot(1,median_sum_diffs(1,21+(9)),'o','Color',colour_sets(6,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(6,:),'MarkerEdgeColor',colour_sets(3,:));
%xlim([1 4]); ylim([-2000 200]); xticks(1:4); xlabel('Memory capacity'); ylabel('Summed lme difference with reduced Bayesian observer model');
%title('Fixed effects model comparison');
%legend(h,{'Model Avg DecFunc + wAvg Pruning','Model Avg DecFunc + maxP Pruning','Model Avg DecFunc + Last Pruning', ...
%          'Model Sel DecFunc + wAvg Pruning','Model Sel DecFunc + maxP Pruning','Model Sel DecFunc + Last Pruning'},'location','southeast');

%subplot(1,2,2); cla; hold on; box on; h = nan(1,6);
h(7) = plot(1:4,median_sum_diffs(1,1:4),['o' lineStyles2{2}],'Color',colour_sets(1,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(1,:));
h(8) = plot(1:4,median_sum_diffs(1,5:8),['o' lineStyles2{2}],'Color',colour_sets(2,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(2,:));
h(9) = plot(1:4,median_sum_diffs(1,9:12),['o' lineStyles2{2}],'Color',colour_sets(3,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(3,:));
h(10) = plot(1:4,[median_sum_diffs(1,1), median_sum_diffs(1,13:15)],['o' lineStyles2{2}],'Color',colour_sets(4,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(4,:));
h(11) = plot(1:4,[median_sum_diffs(1,5), median_sum_diffs(1,16:18)],['o' lineStyles2{2}],'Color',colour_sets(5,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(5,:));
h(12) = plot(1:4,[median_sum_diffs(1,9), median_sum_diffs(1,19:21)],['o' lineStyles2{2}],'Color',colour_sets(6,:),'LineWidth',0.5,'MarkerFaceColor',colour_sets(6,:));
plot(1,median_sum_diffs(1,1),'o','Color',colour_sets(4,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(4,:),'MarkerEdgeColor',colour_sets(1,:));
plot(1,median_sum_diffs(1,5),'o','Color',colour_sets(5,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(5,:),'MarkerEdgeColor',colour_sets(2,:));
plot(1,median_sum_diffs(1,9),'o','Color',colour_sets(6,:),'LineWidth',1.5,'MarkerFaceColor',colour_sets(6,:),'MarkerEdgeColor',colour_sets(3,:));
xlim([1 4]); ylim([-2000 200]); xticks(1:4); xlabel('Memory Capacity (M)'); ylabel('Group-level Summed LME Differences');
title('Fixed Effects Model Comparison');
legend(h(1:6),{'WAVG pruning + Model Averaging','MAXP pruning + Model Averaging','LAST pruning + Model Averaging', ...
               'WAVG pruning + Model Selection','MAXP pruning + Model Selection','LAST pruning + Model Selection'},'location','southeast');

%% Individual memory capacity (Model Avg DecFunc and Mixture Variance Pruning)

Mem_colours = [218,34,34; 238,102,119; 187 187 187]/255;

lme_P1_D1_S1_norm = lme_P1_D1_S1 - lme_P1_D1_S1(:,1);
mem1_wins = lme_P1_D1_S1(:,1) > lme_P1_D1_S1(:,2);

subplot(1,2,2); cla; hold on; box on;
for j_subj = 1:29
    if abs(lme_P1_D1_S1_norm(j_subj,2)) < 1
        h1 = plot(1:4,lme_P1_D1_S1_norm(j_subj,:),'o-','Color',Mem_colours(3,:),'LineWidth',1,'MarkerFaceColor',[.5 .5 .5]);
    elseif mem1_wins(j_subj)
        h2 = plot(1:4,lme_P1_D1_S1_norm(j_subj,:),'o-','Color',Mem_colours(2,:),'LineWidth',1,'MarkerFaceColor',Mem_colours(2,:));
    else
        h3 = plot(1:4,lme_P1_D1_S1_norm(j_subj,:),'o-','Color',Mem_colours(1,:),'LineWidth',1,'MarkerFaceColor',Mem_colours(1,:));
    end
end
xlim([1 4]); ylim([-8 4]); xticks(1:4); yticks(-8:4);
xlabel('Memory capacity'); ylabel('Individual LME differences');
title('Individual Memory Capacity');
legend([h3 h1 h2],{'Larger memory model preference','No clear preference','Limited memory model preference'},'location','southwest');


%% Check which of the Bayesian models is best

%addpath('E:\Matlab Toolboxes\spm12');
addpath('E:\Matlab Toolboxes\rbms_Acerbi\rbms_Acerbi');

%First do a model comparison of the M=1 models only: Three pruning functions (P = WAVG vs. PMAX vs. LAST) and two simplification options (S0 vs. S1)    
lme = LL_BIC_all(:,1,:,1,:);
lme = reshape(lme,[29 6]);
P_levels = reshape(repmat((1:3)',[1 2]),[1 6]);
S_levels = reshape(repmat(0:1,[3 1]),[1 6]);

factors = cell(1,2);
factors{1} = [P_levels==1; P_levels==2; P_levels==3];
factors{2} = [S_levels==0; S_levels==1];

[~,bms_fac] = bms_Acerbi(lme,factors);
disp('Pruning: '); disp(bms_fac{1}.pxp);        %1.0000    0.0000    0.0000     %Mixture_var (WAVG) wins - no evidence for simple switching or simple forget
disp('Simplified: '); disp(bms_fac{2}.pxp);     %0.2701    0.7299               %Simplified wins

%Second, do a model comparison with longer memory models only (M = 2:4)
lme = LL_BIC_all(:,2:4,:,:,:);
lme = reshape(lme,[29 36]);
M_levels = reshape(repmat((2:4)',[1 3 2 2]),[1 36]);
P_levels = reshape(repmat((1:3),[3 1 2 2]),[1 36]);
D_levels = reshape(repmat(permute(1:2,[1 3 2]),[3 3 1 2]),[1 36]);
S_levels = reshape(repmat(permute(0:1,[1 3 4 2]),[3 3 2 1]),[1 36]);

factors = cell(1,4);
factors{1} = [M_levels==2; M_levels==3; M_levels==4];
factors{2} = [P_levels==1; P_levels==2; P_levels==3];
factors{3} = [D_levels==1; D_levels==2];
factors{4} = [S_levels==0; S_levels==1];

[~,bms_fac] = bms_Acerbi(lme,factors);
disp('Memory: '); disp(bms_fac{1}.pxp);         %0.5912    0.2103    0.1987     %Memory 2 wins
disp('Pruning: '); disp(bms_fac{2}.pxp);        %0.9298    0.0406    0.0296     %Mixture wins
disp('Dec Func: '); disp(bms_fac{3}.pxp);       %0.9717    0.0283               %Model Averaging wins
disp('Simplified: '); disp(bms_fac{4}.pxp);     %0.0430    0.9570               %Simplified wins

%Third, select only the WAVG pruning function and Model Averaging decision function and Simplified models.    
%Then test the model memory (M = 1:4)
lme = LL_BIC_all(:,:,1,1,2);
bms = bms_Acerbi(lme);
disp(bms.pxp)                                   %[0.2465    0.2418    0.2471    0.2646]     %No clear winner  for memory

%Finally, compare the memory 1 and 2 models only
lme = LL_BIC_all(:,[1 2],1,1,2);
bms = bms_Acerbi(lme);
disp(bms.pxp)                                   %[0.3141    0.6859]         %No clear winner, but small preference for larger memory model (M2)


%% Also load the fitted parameter values - to compare the M1 and M2 reduced Bayesian observer models

fitted_data_folder = 'fitted_data_Krishna2017_2024';
load(fullfile(fitted_data_folder,'fitted_params_2024.mat'),'fitted_sd_exp_all','fitted_cp_hazard_rate_all');

fitted_logSDs_M1 = log(fitted_sd_exp_all(:,:,1,1,1,2));    %M=1, P=1, D=1, S=2 (i.e. simple = true, not false)
fitted_HRs_M1 = fitted_cp_hazard_rate_all(:,:,1,1,1,2);  

fitted_logSDs_M2 = log(fitted_sd_exp_all(:,:,2,1,1,2));    %M=2, P=1, D=1, S=2 (i.e. simple = true, not false)
fitted_HRs_M2 = fitted_cp_hazard_rate_all(:,:,2,1,1,2);  

rho = corr(fitted_logSDs_M1(:,1),fitted_logSDs_M2(:,1));
disp('Log SDs for low noise condition'); disp(rho);

rho = corr(fitted_logSDs_M1(:,2),fitted_logSDs_M2(:,2));
disp('Log SDs for high noise condition'); disp(rho);

rho = corr(fitted_HRs_M1(:,1),fitted_HRs_M2(:,1));
disp('Hazard rates for low noise condition'); disp(rho);

rho = corr(fitted_HRs_M1(:,2),fitted_HRs_M2(:,2));
disp('Hazard rates for high noise condition'); disp(rho);

%Plot the fitted parameters by participants' preference of M1 vs M2
cond_titles = {'Experimental Noise = 10°','Experimental Noise = 20°'};
HR2plot = .15:.15:.9;
SD2plot = [3 4.5 7 10 14 20 30 45 70 120 200];

figure; 
for j_cond=1:2
    subplot(1,2,j_cond); cla; hold on; box on;
    for j_subj = 1:29
        if abs(lme_P1_D1_S1_norm(j_subj,2)) < 1
            h1 = plot(fitted_HRs_M1(j_subj,j_cond),fitted_logSDs_M1(j_subj,j_cond),'o','Color',Mem_colours(3,:),'MarkerSize',10,'MarkerFaceColor',[.5 .5 .5]);
        elseif mem1_wins(j_subj)
            h2 = plot(fitted_HRs_M1(j_subj,j_cond),fitted_logSDs_M1(j_subj,j_cond),'o','Color',Mem_colours(2,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(2,:));
        else
            h3 = plot(fitted_HRs_M1(j_subj,j_cond),fitted_logSDs_M1(j_subj,j_cond),'o','Color',Mem_colours(1,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(1,:));
        end
    end
    xlim([0 1]); ylim([log(3) log(210)]);
    xticks(HR2plot); yticks(log(SD2plot)); yticklabels(cellfun(@(x) num2str(x),num2cell(SD2plot),'UniformOutput',false));
    xlabel('hazard rate H'); ylabel('sd exp (log-scale)'); title(cond_titles{j_cond});
end
