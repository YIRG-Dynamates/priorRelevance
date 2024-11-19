%% Load surface data and fitted params

clc
clearvars;

%Make sure to run computeSurprisal.m first!
load('surprisal_surface_plot_data.mat','surprisal','SDs','HRs');

%Also load the fitted parameter values
fitted_data_folder = 'fitted_data_Krishna2017_2024';
load(fullfile(fitted_data_folder,'fitted_params_2024.mat'),'fitted_sd_exp_all','fitted_cp_hazard_rate_all','fitted_lapse_rate_all','fitted_sd_motor_all','LL_BIC_all');

%% Plot the fitted parameters on top of the 2D surprisal surface

fitted_logSDs = log(fitted_sd_exp_all(:,:,1,1,1,2));    %M=1, P=1, D=1, S=2 (i.e. simple = true, not false)
fitted_HRs = fitted_cp_hazard_rate_all(:,:,1,1,1,2);  
fitted_lapse_rates = fitted_lapse_rate_all(:,:,1,1,1,2);
fitted_sd_motor = fitted_sd_motor_all(:,:,1,1,1,2);  

num_SDs = numel(SDs);
num_HRs = numel(HRs);
SDs_mat = repmat(SDs',[1 num_HRs]);
HRs_mat = repmat(HRs,[num_SDs 1]);

HR2plot = .15:.15:.9;
SD2plot = [3 4.5 7 10 14 20 30 45 70 120 200];

figure('Name','Fitted parameters on top of surprisal surface'); 
cond_titles = {'Experimental Noise = 10°','Experimental Noise = 20°'};
sd_exp_true = [10 20];
HR_true = 0.15;

wavg_color = [204 77 119]/255; 

%Use colour coding similar to figure 6B (i.e. M=1 vs. M=2 preference)
Mem_colours = [218,34,34; 238,102,119; 187 187 187]/255;

lme = LL_BIC_all;                                                           %1=num_subj(29), 2=num_mems(4), 3=num_pruns(3), 4=num_decfns(2), 5=num_simpl(2)
lme_M1 = lme(:,1,1,1,2);                                                    %M=1, P=1, D=1, S=2 (i.e. simple = true, not false)
lme_M2 = lme(:,2,1,1,2);                                                    %M=1, P=1, D=1, S=2 (i.e. simple = true, not false)
lme_diff = abs(lme_M1-lme_M2);
mem1_wins = lme_M1 > lme_M2;

surprisal2plot = cell(1,2);
for j_cond = 1:2
    
    %Note that surprisal values for SD=10 range between (1) very low for the true values, and (2) very high for anything slightly away from it. 
    %Surprisal values for SD=20 on the other hand, are more diffuse and moderate all over, higher for the true values, and lower for the wrong values.    
    %This is why we scale the values per condition, to see the colour differences well. 
    surprisal2plot{j_cond} = surprisal(:,:,j_cond);
    surprisal2plot{j_cond} = surprisal2plot{j_cond}-min(surprisal2plot{j_cond},[],'all');
    surprisal2plot{j_cond} = surprisal2plot{j_cond}/surprisal2plot{j_cond}(num_SDs,num_HRs);         
    surprisal2plot{j_cond} = min(surprisal2plot{j_cond},1.5);
    %surprisal2plot{j_cond} = 63*(surprisal2plot{j_cond} / max(surprisal2plot{j_cond},[],'all'))+1;      %Map to 1-64 range for direct colormap
    
    subplot(1,2,j_cond); cla; hold on; box on;
    s = pcolor(HRs_mat,log(SDs_mat),surprisal2plot{j_cond}); 
    s.FaceColor = 'interp'; set(s,'edgecolor','none'); %colorbar; %s.CDataMapping = 'direct';
    YLIM = [min(log(SDs)),max(log(SDs))]; plot([0 1],log(sd_exp_true(j_cond)*[1 1]),'k:'); plot(HR_true*[1 1],YLIM,'k:');
    
    %Find minimum surprisal indices in the sd_exp dimension for each value of the hazard rate and plot as a dashed line    
    [~,idx] = min(surprisal(:,:,j_cond));
    min_logSDs = smooth(log(SDs(idx)));
    plot(HRs,min_logSDs,'k--');
    
    %Plot markers for fitted subject parameters
    for j_subj = 1:29
        %h = plot(fitted_HRs(j_subj,j_cond),fitted_logSDs(j_subj,j_cond),'o','Color',wavg_color,'MarkerSize',10,'MarkerFaceColor',wavg_color);
        if lme_diff(j_subj) < 1
            h1 = plot(fitted_HRs(j_subj,j_cond),fitted_logSDs(j_subj,j_cond),'o','Color',Mem_colours(3,:),'MarkerSize',10,'MarkerFaceColor',[.5 .5 .5]);
        elseif mem1_wins(j_subj)
            h2 = plot(fitted_HRs(j_subj,j_cond),fitted_logSDs(j_subj,j_cond),'o','Color',Mem_colours(2,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(2,:));
        else
            h3 = plot(fitted_HRs(j_subj,j_cond),fitted_logSDs(j_subj,j_cond),'o','Color',Mem_colours(1,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(1,:));
        end
    end
    
    xlim([min(HRs),max(HRs)]); ylim(YLIM); 
    xticks(HR2plot); yticks(log(SD2plot)); yticklabels(cellfun(@(x) num2str(x),num2cell(SD2plot),'UniformOutput',false));
    xlabel('Hazard Rate'); ylabel('Experimental Noise'); title(cond_titles{j_cond});
    %legend(h,'wAvg','location','northwest');
    
    %Check the correlation between the parameters
    [rho,p] = corr([fitted_HRs(:,j_cond), fitted_logSDs(:,j_cond)]);
    disp(['Fitted parameter correlation for condition number = ' num2str(j_cond)]); 
    disp(['rho = ' num2str(rho(2))]); disp(['p = ' num2str(p(2))]);
end

%% Plot the fitted parameters as a correlation between sd_exp 10 and sd_exp 20

figure('Name','Fitted parameter correlations across noise conditions'); 

%Hazard rate
subplot(1,2,1); cla; hold on; box on;
YLIM = [min(HRs),max(HRs)]; plot(YLIM,YLIM,'k:');                   %diagonal
plot([0.15 0.15],YLIM,'k:'); plot(YLIM,[0.15 0.15],'k:');           %vertical and horizontal
for j_subj = 1:29
    if lme_diff(j_subj) < 1
        h1 = plot(fitted_HRs(j_subj,1),fitted_HRs(j_subj,2),'o','Color',Mem_colours(3,:),'MarkerSize',10,'MarkerFaceColor',[.5 .5 .5]);
    elseif mem1_wins(j_subj)
        h2 = plot(fitted_HRs(j_subj,1),fitted_HRs(j_subj,2),'o','Color',Mem_colours(2,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(2,:));
    else
        h3 = plot(fitted_HRs(j_subj,1),fitted_HRs(j_subj,2),'o','Color',Mem_colours(1,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(1,:));
    end
end
xlim(YLIM); ylim(YLIM); xticks(HR2plot); yticks(HR2plot); 
xlabel('Low Noise Condition'); ylabel('High Noise Condition'); title('Hazard Rates');
%legend([h3 h1 h2],{'Larger memory model preference','No clear preference','Limited memory model preference'},'location','southeast');

%Compute correlation coefficient for hazard rates
rho = corr(fitted_HRs);
disp('Fitted hazard rate correlation:'); disp(rho(2));

%Experimental noise
subplot(1,2,2); cla; hold on; box on;
YLIM = [min(log(SDs)),max(log(SDs))]; plot(YLIM,YLIM,'k:');         %diagonal      
%plot(log(SDs),log(2*SDs),'k--');                                   %double values
plot(log([10 10]),YLIM,'k:'); plot(YLIM,log([20 20]),'k:');         %vertical and horizontal
for j_subj = 1:29
    if lme_diff(j_subj) < 1
        h1 = plot(fitted_logSDs(j_subj,1),fitted_logSDs(j_subj,2),'o','Color',Mem_colours(3,:),'MarkerSize',10,'MarkerFaceColor',[.5 .5 .5]);
    elseif mem1_wins(j_subj)
        h2 = plot(fitted_logSDs(j_subj,1),fitted_logSDs(j_subj,2),'o','Color',Mem_colours(2,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(2,:));
    else
        h3 = plot(fitted_logSDs(j_subj,1),fitted_logSDs(j_subj,2),'o','Color',Mem_colours(1,:),'MarkerSize',10,'MarkerFaceColor',Mem_colours(1,:));
    end 
end
xlim(YLIM); ylim(YLIM); xticks(log(SD2plot)); yticks(log(SD2plot)); labels_tmp = cellfun(@(x) num2str(x),num2cell(SD2plot),'UniformOutput',false);
xticklabels(labels_tmp); yticklabels(labels_tmp);
xlabel('Low Noise Condition'); ylabel('High Noise Condition'); title('Experimental Noise (°)');
legend([h3 h1 h2],{'Larger memory model preference','No clear preference','Limited memory model preference'},'location','southeast');

%Compute statistics for experimental noise
disp('Fitted experimental noise parameters:');
disp('Low noise:'); disp(quantile(exp(fitted_logSDs(:,1)),[.5 .25 .75]));
disp('High noise:'); disp(quantile(exp(fitted_logSDs(:,2)),[.5 .25 .75]));
[p,~,STATS] = signrank(fitted_logSDs(:,1),fitted_logSDs(:,2)); disp(p);


%% Descriptive stats on other params

disp('Fitted lapse rate parameters:');
disp('Low noise:'); disp(quantile(fitted_lapse_rates(:,1),[.5 .25 .75]));
disp('High noise:'); disp(quantile(fitted_lapse_rates(:,2),[.5 .25 .75]));

disp('Fitted motor noise parameters:');
disp('Low noise:'); disp(quantile(fitted_sd_motor(:,1),[.5 .25 .75]));
disp('High noise:'); disp(quantile(fitted_sd_motor(:,2),[.5 .25 .75]));
