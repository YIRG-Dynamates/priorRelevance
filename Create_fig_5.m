clc
clearvars;

%Add the utilities folder
addpath('utilities'); 

HR = 0.15;
sd_exp = 20;
var_exp = sd_exp^2;

n=3;

x = [-54.0150   11.0710   mean([-54.0150   11.0710])] - 20;

%Define reduced Bayesian observer model
comp_tau = @(var_prior) var_exp ./ (var_exp + var_prior);                   %Eq 3

comp_Q = @(x,mu_prior,var_prior) log(((1./HR-1).*180)./(sqrt(2*pi).*sqrt(var_prior+var_exp)))-.5.*(mu_prior-x).^2./(var_prior+var_exp);  
comp_PI = @(Q) 1./(1+exp(-Q));                                              %Eq 4

update_var = @(var_prior,mu_prior,x,tau,PI) (1-PI).*var_exp + PI.*tau*var_prior + PI.*(1-PI).*tau.^2.*(mu_prior-x).^2;
update_mu = @(mu_prior,x,tau,PI) x+PI.*tau.*(mu_prior-x);                   %Eq 5    

%Define simplified Bayesian observer model
update_weights = @(w_post) [(1-HR)*w_post, HR];
lik_cp = 1/180;
comp_lik_noCP = @(x,mu_prior,var_prior) normpdf(x,mu_prior,sqrt(var_prior));
comp_weights = @(w_prior,lik_noCP) w_prior.*[lik_noCP lik_cp] / sum(w_prior.*[lik_noCP lik_cp],2);

int_mu = @(mu_prior,x,tau) x+tau.*(mu_prior-x); 
int_var = @(var_prior,tau) tau.*var_prior; 

%Initialize model variables
PI = zeros(n,1);
mu_prior = nan(n+1,1);
var_prior = nan(n+1,1);
mu_prior_full = cell(n+1,1);
var_prior_full = cell(n+1,1);
w_prior_full = cell(n+1,1);

%Run the models
mu_prior(2) = x(1);
var_prior(2) = var_exp;
mu_prior_full{2} = x(1);
var_prior_full{2} = var_exp;
w_prior_full{2} = 1;
for j_stim=2:n
    
    %Reduced
    tau = comp_tau(var_prior(j_stim));
    PI(j_stim) = comp_PI(comp_Q(x(j_stim),mu_prior(j_stim),var_prior(j_stim)));
    var_prior(j_stim+1) = update_var(var_prior(j_stim),mu_prior(j_stim),x(j_stim),tau,PI(j_stim));
    mu_prior(j_stim+1) = update_mu(mu_prior(j_stim),x(j_stim),tau,PI(j_stim));
    
    %Simplified full Bayes
    w_prior = update_weights(w_prior_full{j_stim});                         %Posterior weights to prior weights
    lik_noCP = comp_lik_noCP(x(j_stim),mu_prior_full{j_stim},var_prior_full{j_stim}+var_exp);
    w_prior_full{j_stim+1} = comp_weights(w_prior,lik_noCP);                %Prior weights to posterior weights
    
    tau_full = comp_tau(var_prior_full{j_stim});
    mu_prior_full{j_stim+1} = [int_mu(mu_prior_full{j_stim},x(j_stim),tau_full), x(j_stim)];
    var_prior_full{j_stim+1} = [int_var(var_prior_full{j_stim},tau_full), var_exp];
end

%Quick check on equality for the 2nd posteriors 
[expectation, variance] = twoMomentsOfMixture(mu_prior_full{3}, var_prior_full{3}, w_prior_full{3}, 2);
disp([expectation mu_prior(3)]);
disp([variance var_prior(3)]);

% plot sequence
LW = 1.5;
alpha = 0.4;

n_grid = 1000;
x_grid = linspace(-90,90,n_grid);
x_grid2 = [x_grid,flip(x_grid,2)];

prior_color = [102 204 238]/255;        %cyan
lik_color = [34 134 51]/255;            %green
post_color = [221 170 51]/255;          %yellow    
reduce_color = [238 102 119]/255;       %red

XLIM = [-90 45];
YLIM = [0 0.035];

figure('Name','Modelling Framework','WindowState', 'maximized');

for j_stim=1:n
    
    ax = subplot(3,2,2*(j_stim-1)+1); cla; box on; hold on;
    
    y_post = zeros(1,n_grid);
    for j_node=1:numel(w_prior_full{j_stim+1})
        y_post = y_post + w_prior_full{j_stim+1}(j_node).*normpdf(x_grid,mu_prior_full{j_stim+1}(j_node),sqrt(var_prior_full{j_stim+1}(j_node)));
    end
    patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
    for j_node=1:numel(w_prior_full{j_stim+1})
        y_post_node = w_prior_full{j_stim+1}(j_node).*normpdf(x_grid,mu_prior_full{j_stim+1}(j_node),sqrt(var_prior_full{j_stim+1}(j_node)));
        h3 = plot(x_grid,y_post_node,'--','Color',post_color,'LineWidth',LW);
    end
    plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);
    y_reduce = normpdf(x_grid,mu_prior(j_stim+1),sqrt(var_prior(j_stim+1)));
    h4 = plot(x_grid,y_reduce,'--','Color',reduce_color,'LineWidth',LW+1);
    h1 = plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
    xlim(XLIM); ylim(YLIM); yticks([]); xticks([-90 -45 0 45 90]);
    
    %ylabel(['Stimulus Number ' num2str(j_stim)]);
    
    if j_stim==1
        title('Memory Extension');
    else
        %title(['Prior Relevance = ' num2str(round(100*PI(j_stim))/100)]);
    end
    
    if j_stim==2
        h2 = patch(181+x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', post_color);
        legend([h1 h2 h3 h4],{'Stimulus','Posterior','Nodes','Reduced'});
    end
end

%%%%%%%%%%%%%%%%%%
%Decision function
%%%%%%%%%%%%%%%%%%

%Take posterior from 2nd stimulus, apply decision strategy: model averaging vs. model selection   
j_stim = 2;

model_avg_color = [238,5,21]/255; 
model_sel_color = [178,24,43]/255;

ax = subplot(3,2,4); cla; box on; hold on;

y_post = zeros(1,n_grid);
for j_node=1:numel(w_prior_full{j_stim+1})
    y_post = y_post + w_prior_full{j_stim+1}(j_node).*normpdf(x_grid,mu_prior_full{j_stim+1}(j_node),sqrt(var_prior_full{j_stim+1}(j_node)));
end
patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
for j_node=1:numel(w_prior_full{j_stim+1})
    y_post_node = w_prior_full{j_stim+1}(j_node).*normpdf(x_grid,mu_prior_full{j_stim+1}(j_node),sqrt(var_prior_full{j_stim+1}(j_node)));
    plot(x_grid,y_post_node,'--','Color',post_color,'LineWidth',LW);
end
plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);

%Add the expectations 
expectation = twoMomentsOfMixture(mu_prior_full{j_stim+1}, var_prior_full{j_stim+1}, w_prior_full{j_stim+1}, 2);
[~,idx] = min(abs(expectation-x_grid));
h2 = plot(expectation*[1 1],[0 y_post(idx)],'-','Color',model_avg_color,'LineWidth',LW+1);

y_post_node = w_prior_full{j_stim+1}(1).*normpdf(x_grid,mu_prior_full{j_stim+1}(1),sqrt(var_prior_full{j_stim+1}(1)));
[~,idx] = min(abs(mu_prior_full{j_stim+1}(1)-x_grid));
h3 = plot(mu_prior_full{j_stim+1}(1)*[1 1],[0 y_post_node(idx)],'-','Color',model_sel_color,'LineWidth',LW+1);

h1 = plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);

xlim(XLIM); ylim(YLIM); yticks([]); xticks([-90 -45 0 45 90]);
title('Decision strategy');
legend([h2 h3],{'Model Averaging','Model Selection'},'location','northeast');

%%%%%%%%
%Pruning
%%%%%%%%

%Take posterior from last stimulus, apply M=2 pruning: WAVG, PMAX, LAST
j_stim = 3;

%wavg_color = red/purple
%pmax_color = yellow
%last_color = green
colour_sets = [238,5,21; 181,101,29; 144,201,135; ...                       %Figure 6 colors
               178,24,43; 123,80,25; 17,85,34]/255;
pruning_colors = [mean([colour_sets(1,:); colour_sets(4,:)],1); ...         %1 = WAVG
                  mean([colour_sets(2,:); colour_sets(5,:)],1); ...         %2 = PMAX
                  mean([colour_sets(3,:); colour_sets(6,:)],1)];            %3 = LAST

w_prior_pruned = [sum(w_prior_full{j_stim+1}(1:2)) w_prior_full{j_stim+1}(3)];

[expectation, variance] = twoMomentsOfMixture(mu_prior_full{j_stim+1}(1:2), var_prior_full{j_stim+1}(1:2), w_prior_full{j_stim+1}(1:2), 2);
mu_WAVG = [expectation, mu_prior_full{j_stim+1}(3)];
var_WAVG = [variance, var_prior_full{j_stim+1}(3)];

mu_PMAX = mu_prior_full{j_stim+1}([1 3]);
var_PMAX = var_prior_full{j_stim+1}([1 3]);

mu_LAST = mu_prior_full{j_stim+1}([2 3]);
var_LAST = var_prior_full{j_stim+1}([2 3]);

mu_all = {mu_WAVG,mu_PMAX,mu_LAST};
var_all = {var_WAVG,var_PMAX,var_LAST};

%Plot the results 
%figure('Name','Pruning functions','WindowState', 'maximized');
ax = subplot(3,2,6); cla; box on; hold on; h = nan(1,3);

for j_pn=1:n
    
    y_prior_1 = w_prior_pruned(1).*normpdf(x_grid,mu_all{j_pn}(1),sqrt(var_all{j_pn}(1)));
    y_prior_2 = w_prior_pruned(2).*normpdf(x_grid,mu_all{j_pn}(2),sqrt(var_all{j_pn}(2)));
    y_prior = y_prior_1 + y_prior_2;
    
    patch(x_grid2,[y_prior,zeros(1,n_grid)],pruning_colors(j_pn,:),'facealpha', alpha, 'edgecolor', 'none');
    plot(x_grid,y_prior_1,'--','Color',pruning_colors(j_pn,:),'LineWidth',LW);
    plot(x_grid,y_prior_2,'--','Color',pruning_colors(j_pn,:),'LineWidth',LW);
    plot(x_grid,y_prior,'-','Color',pruning_colors(j_pn,:),'LineWidth',LW);
    
    h(j_pn) = patch(181+x_grid2,[y_post,zeros(1,n_grid)],pruning_colors(j_pn,:),'facealpha', alpha, 'edgecolor', pruning_colors(j_pn,:));
end

h1 = plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
   
xlim(XLIM); ylim(YLIM); yticks([]); xticks([-90 -45 0 45 90]);
title('Pruning function');

legend(h,{'Keep WAVG','Keep PMAX','Keep LAST'},'location','northeast');


%%%%%%%%
%Late truncation simplification
%%%%%%%%

j_stim = 1;

ax = subplot(3,2,2); cla; box on; hold on; 

y_post = normpdf(x_grid,mu_prior_full{j_stim+1}(1),sqrt(var_prior_full{j_stim+1}(1)));

patch(x_grid2,[y_post,zeros(1,n_grid)],post_color,'facealpha', alpha, 'edgecolor', 'none');
plot(x_grid,y_post,'-','Color',post_color,'LineWidth',LW);

h1 = plot(x(j_stim),0,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10);
xlim(XLIM); ylim(YLIM); yticks([]); xticks([-90 -45 0 45 90]);

[~,idx] = min(abs(x(j_stim)-x_grid));
h2 = plot(x(j_stim)*[1 1],[0 y_post(idx)],'k-','LineWidth',LW+.5);

expectation = TNmeanvar(-90,90,mu_prior_full{j_stim+1}(1),sqrt(var_prior_full{j_stim+1}(1)));
[~,idx] = min(abs(expectation-x_grid));
h3 = plot(expectation*[1 1],[0 y_post(idx)],'k--','LineWidth',LW+.5);

title('Late truncation simplification');
legend([h2 h3],{'With Simplification','Without Simplification'},'location','northeast');
