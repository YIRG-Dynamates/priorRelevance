clc
clearvars;

%Add the utilities folder
addpath('utilities'); 

%% Draw sequence of stimuli (Figure 2A)

HR = 0.15;
var_exp = 20^2;

n=10;

% sample_cp = @() rand(1)<HR;
% sample_mu = @() rand(1)*120-60;
% sample_x = @(mu) sqrt(var_exp)*randn(1)+mu;
% 
% x = nan(n,1);
% cp = zeros(n,1); cp([1 4 6]) = true;          %Hardcode, SAC 5 trial
% for j_stim=1:n
%     if cp(j_stim) %sample_cp()
%         mu = sample_mu();
%     end
%     x(j_stim) = sample_x(mu);
% end

x = [35.2046   72.9987   50.6596  -49.7294  -54.0150   6.5510   45.4038   17.0930   33.1000   9.8446];

%Define reduced Bayesian observer model
comp_tau = @(var_prior) var_exp / (var_exp + var_prior);                    %Eq 3

comp_Q = @(x,mu_prior,var_prior) log(((1/HR-1)*180)/(sqrt(2*pi)*sqrt(var_prior+var_exp)))-.5*(mu_prior-x)^2/(var_prior+var_exp);  
comp_PI = @(Q) 1/(1+exp(-Q));                                               %Eq 4

update_var = @(var_prior,mu_prior,x,tau,PI) (1-PI)*var_exp + PI*tau*var_prior + PI*(1-PI)*tau^2*(mu_prior-x)^2;
update_mu = @(mu_prior,x,tau,PI) x+PI*tau*(mu_prior-x);                     %Eq 5    

int_mu = @(mu_prior,x,tau) x+tau*(mu_prior-x); 
int_var = @(var_prior,tau) tau*var_prior; 

%Initialize model variables
tau = nan(n,1);
Q = nan(n,1);
PI = nan(n,1);
mu_prior = nan(n+1,1);
var_prior = nan(n+1,1);
mu_post_int = nan(n,1);
var_post_int = nan(n,1);

%Run the model
mu_prior(2) = x(1);
var_prior(2) = var_exp;
for j_stim=2:n
    tau(j_stim) = comp_tau(var_prior(j_stim));
    Q(j_stim) = comp_Q(x(j_stim),mu_prior(j_stim),var_prior(j_stim));
    PI(j_stim) = comp_PI(Q(j_stim));
    var_prior(j_stim+1) = update_var(var_prior(j_stim),mu_prior(j_stim),x(j_stim),tau(j_stim),PI(j_stim));
    mu_prior(j_stim+1) = update_mu(mu_prior(j_stim),x(j_stim),tau(j_stim),PI(j_stim));
    
    mu_post_int(j_stim) = int_mu(mu_prior(j_stim),x(j_stim),tau(j_stim)); 
    var_post_int(j_stim) = int_var(var_prior(j_stim),tau(j_stim)); 
end

%plot sequence
figure;
ax = subplot(3,4,[1 5 9]); cla; box on; hold on;
for j_stim=1:n
    plot(x(j_stim),j_stim,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
end
grey = [187 187 187]/255;
patch([-90 90 90 -90],[3.6 3.6 4.4 4.4],grey,'facealpha', 0.2, 'edgecolor', grey);
patch([-90 90 90 -90],[4.6 4.6 5.4 5.4],grey,'facealpha', 0.2, 'edgecolor', grey);
patch([-90 90 90 -90],[5.6 5.6 6.4 6.4],grey,'facealpha', 0.2, 'edgecolor', grey);
set(ax,'YDir','reverse');
xlim([-90, 90]); ylim([.5, n+.5]);
yticks(1:n); xticks([-90 -45 0 45 90]);


prior_color = [102 204 238]/255;        %cyan
lik_color = [34 134 51]/255;            %green
post_color = [221 170 51]/255;          %yellow    
reduce_color = [238 102 119]/255;       %red

LW = 1.5;
alpha = 0.4;

n_grid = 1000;
x_grid = linspace(-90,90,n_grid);
x_grid2 = [x_grid,flip(x_grid,2)];

%Zoom in on changepoint
j_stim = 4;
ax = subplot(3,4,2:4); cla; box on; hold on;

y_prior = normpdf(x_grid,mu_prior(j_stim),sqrt(var_prior(j_stim)));
y_lik = normpdf(x_grid,x(j_stim),sqrt(var_exp));
y_reduce = normpdf(x_grid,mu_prior(j_stim+1),sqrt(var_prior(j_stim+1)));
y_post_C1 = normpdf(x_grid,mu_post_int(j_stim),sqrt(var_post_int(j_stim)));
y_post = PI(j_stim)*y_post_C1 + (1-PI(j_stim))*y_lik;

patch(x_grid2,[y_prior,zeros(1,n_grid)],prior_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_prior,'-','Color',prior_color,'LineWidth',LW);

patch(x_grid2,[y_lik,zeros(1,n_grid)],lik_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_lik,'-','Color',lik_color,'LineWidth',LW);

patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);

plot(x_grid,y_reduce,'--','Color',reduce_color,'LineWidth',LW+1);
plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);

xlim([-90,90]); ylim([0 0.031]); yticks([]); xticks([-90 -45 0 45 90]);
title('Changepoint');

%Zoom in on NO changepoint
j_stim = 5;
ax = subplot(3,4,6:8); cla; box on; hold on;

y_prior = normpdf(x_grid,mu_prior(j_stim),sqrt(var_prior(j_stim)));
y_lik = normpdf(x_grid,x(j_stim),sqrt(var_exp));
y_reduce = normpdf(x_grid,mu_prior(j_stim+1),sqrt(var_prior(j_stim+1)));
y_post_C1 = normpdf(x_grid,mu_post_int(j_stim),sqrt(var_post_int(j_stim)));
y_post = PI(j_stim)*y_post_C1 + (1-PI(j_stim))*y_lik;

patch(x_grid2,[y_prior,zeros(1,n_grid)],prior_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_prior,'-','Color',prior_color,'LineWidth',LW);

patch(x_grid2,[y_lik,zeros(1,n_grid)],lik_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_lik,'-','Color',lik_color,'LineWidth',LW);

patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);

h5 = plot(x_grid,y_reduce,'--','Color',reduce_color,'LineWidth',LW+1);
h2 = plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);

xlim([-90,90]); ylim([0 0.031]); yticks([]); xticks([-90 -45 0 45 90]);
title('No Changepoint');

h1 = patch(181+x_grid2,[y_prior,zeros(1,n_grid)],prior_color,'facealpha', alpha, 'edgecolor', prior_color);
h3 = patch(181+x_grid2,[y_lik,zeros(1,n_grid)],lik_color,'facealpha', alpha, 'edgecolor', lik_color);
h4 = patch(181+x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', post_color);
legend([h1 h2 h3 h4 h5],{'Prior','Stimulus','Likelihood','Posterior','Reduced'})


%Zoom in on causal uncertainty
j_stim = 6;
ax = subplot(3,4,10:12); cla; box on; hold on;

y_prior = normpdf(x_grid,mu_prior(j_stim),sqrt(var_prior(j_stim)));
y_lik = normpdf(x_grid,x(j_stim),sqrt(var_exp));
y_reduce = normpdf(x_grid,mu_prior(j_stim+1),sqrt(var_prior(j_stim+1)));
y_post_C1 = normpdf(x_grid,mu_post_int(j_stim),sqrt(var_post_int(j_stim)));
y_post = PI(j_stim)*y_post_C1 + (1-PI(j_stim))*y_lik;

patch(x_grid2,[y_prior,zeros(1,n_grid)],prior_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_prior,'-','Color',prior_color,'LineWidth',LW);

patch(x_grid2,[y_lik,zeros(1,n_grid)],lik_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_lik,'-','Color',lik_color,'LineWidth',LW);

patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);
plot(x_grid,PI(j_stim)*y_post_C1,'--','Color',post_color,'LineWidth',LW);
plot(x_grid,(1-PI(j_stim))*y_lik,'--','Color',post_color,'LineWidth',LW);

plot(x_grid,y_reduce,'--','Color',reduce_color,'LineWidth',LW+1);
plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);

xlim([-90,90]); ylim([0 0.031]); yticks([]); xticks([-90 -45 0 45 90]);
title('Causal Uncertainty');


%% Figure 2B. Logit transformed prior relevance as function of prediction error

disp2cdf = @(disp) disp*2/180-disp.^2/180^2;
cdf2disp = @(cdf) 180-sqrt(180^2*(1-cdf));

n_grid = 1000;
max_disp = 180;
max_cdf = disp2cdf(max_disp);
x_cdf = linspace(0,max_cdf,n_grid);
x_disp = cdf2disp(x_cdf);

x_ticks_disp = [0 10 20 30 40 50 60 80 110 180];
x_ticks_cdf = disp2cdf(x_ticks_disp);
x_tick_labels = cellfun(@(x) num2str(x),num2cell(x_ticks_disp),'UniformOutput',false);

%Settings for 8 different lines
sd_exp = [10 10 10 10, 20 20 20 20];
var_exp = sd_exp.^2;

sd_prior = sd_exp./[1 1 sqrt(2) sqrt(2) 1 1 sqrt(2) sqrt(2)];
var_prior = sd_prior.^2;

tau = var_exp ./ (var_exp + var_prior); %prior reliability                  (SAC 2 2 3 3 2 2 3 3)

HR = [.15 .6 .15 .6 .15 .6 .15 .6];

colours = [[253 179 102]; [253 179 102]; ...      %10, SAC2 (light orange), HR true (solid) vs. too high (dashed)         
           [209 24 7]; [209 24 7]; ...        %10, SAC3 (dark orange), HR true (solid) vs. too high (dashed)
           [187 204 51]; [187 204 51]; ...      %20, SAC2 (light pear), HR true (solid) vs. too high (dashed)
           [110 110 0]; [120 120 0]]/255;       %20, SAC3 (dark olive), HR true (solid) vs. too high (dashed)

linestyles = {'-','--','-','--','-','--','-','--'};

line_labels = {'10° noise, \tau = 1/2, H = 0.15', '10° noise, \tau = 1/2, H = 0.60', ...
               '10° noise, \tau = 2/3, H = 0.15', '10° noise, \tau = 2/3, H = 0.60', ...
               '20° noise, \tau = 1/2, H = 0.15', '20° noise, \tau = 1/2, H = 0.60', ...
               '20° noise, \tau = 2/3, H = 0.15', '20° noise, \tau = 2/3, H = 0.60'};

comp_Q = @(disp,HR,var_prior,var_exp) log(((1./HR-1).*180)./(sqrt(2*pi).*sqrt(var_prior+var_exp)))-.5.*(disp).^2./(var_prior+var_exp); 
comp_PI = @(Q) 1./(1+exp(-Q));                                                                                              %Eq 4       
comp_sd = @(disp,var_prior,var_exp,tau,PI) sqrt((1-PI).*var_exp + PI.*tau.*var_prior + PI.*(1-PI).*tau.^2.*(disp).^2);      %Eq 5 


Q = nan(8,n_grid);
PI = nan(8,n_grid);
SD = nan(8,n_grid);
for j=1:8
    Q(j,:) = comp_Q(x_disp,HR(j),var_prior(j),var_exp(j));
    PI(j,:) = comp_PI(Q(j,:));
    SD(j,:) = comp_sd(x_disp,var_prior(j),var_exp(j),tau(j),PI(j,:));
end

logit = @(p) log(p./(1-p));
Q_ref = logit([.1 .9]);


figure;
LW = 1.5;

%Transformed prior relevance
subplot(3,1,1); cla; hold on; box on; h = nan(1,8);
for j=1:8
    h(j) = plot(x_cdf,Q(j,:),'LineStyle',linestyles{j},'Color',colours(j,:),'LineWidth',LW);
end
plot([0 1],[0 0],'k:','LineWidth',0.5);
plot([0 1],Q_ref(1)*[1 1],'k:','LineWidth',0.5); plot([0 1],Q_ref(2)*[1 1],'k:','LineWidth',0.5);
ylim([-5 5]); 
xlim([0 1]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); %xlabel('Prediction Error'); ylabel('Q (logit Prior Relevance)');
%legend(h,line_labels,'location','northeast');

%Normalized Bias
subplot(3,1,2); cla; hold on; box on; h = nan(1,8);
for j=1:8
    h(j) = plot(x_cdf,tau(j).*PI(j,:),'LineStyle',linestyles{j},'Color',colours(j,:),'LineWidth',LW);
end
plot([0 1],tau(1)*[1 1],'k:','LineWidth',0.5); plot([0 1],tau(3)*[1 1],'k:','LineWidth',0.5);
plot([0 1],[0 0],'k:','LineWidth',0.5);
ylim([-.05 0.8]); 
xlim([0 1]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); %xlabel('Prediction Error'); ylabel('Normalized Bias');
legend(h,line_labels,'location','northeast');

%Spatial Uncertainty
subplot(3,1,3); cla; hold on; box on; h = nan(1,8);
for j=1:8
    h(j) = plot(x_cdf,SD(j,:),'LineStyle',linestyles{j},'Color',colours(j,:),'LineWidth',LW);
end
plot([0 1],[10 10],'k:','LineWidth',0.5);
plot([0 1],[20 20],'k:','LineWidth',0.5);
ylim([0 30]); 
xlim([0 1]); xticks(x_ticks_cdf); xticklabels(x_tick_labels); %xlabel('Prediction Error'); ylabel('Spatial Uncertainty (SD°)');

