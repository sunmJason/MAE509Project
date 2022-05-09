%% Metropolis algorithm
% This code illustrates the implementation of the Metropolis algorithm for
% plastic deformation model, e: strain, T: stress
%
%     e = T/E + 0.002*(T/Ty)^3.5, where E and Ty are model parameters
%
% Here we are sampling the Young's modulus and yield stress
% parameters theta = [E,Ty]

clear all; close all; clc	

 
%% Define forward model
% plastic deformation model is defined as an anonymous function

forward = @(t,p) 10*3+ sum([p t].^2-10*cos(2*pi*[p t]),2);  %2*sin(-parameters(1)*t*0.5).^3+(cos(0.25*t.*sqrt(parameters(2)-parameters(1)^2/4))).^3;%T./parameters(1) + 0.002*(T./parameters(2)).^3.5 ;



%% Generate synthetic data from forward model

% true parameters
C = 0;
K = 0;

% true stress-strain points 
n = 501;
T_true = linspace(-5.16,5.16, n);
for i=1:n
    e_true(i) = forward(T_true(i),[C, K]);
end
% adding fixed data noise to the true strain
num_exprmt = 5; % number of stress-strain curves
data_stDev = 0.1;

e_exprmts=zeros(num_exprmt,length(e_true));
for i = 1:num_exprmt
    e_exprmts(i,:) = e_true + normrnd(0,data_stDev, size(e_true));
end
e_data_mean = mean(e_exprmts);
e_data_stDev = std(e_exprmts);

% plot synthetic data
fig = figure();
plot(T_true, e_true, '--bo', 'LineWidth', 2); hold on;
errorbar(T_true, e_data_mean, e_data_stDev, '-ro', 'LineWidth', 2);
% xlim([25,350])
xlabel('Stress, T')
ylabel('Strain, e')
legend('True model','Sythetic data','location','northwest')
set(gca,'FontSize',20)
set(gcf,'Position',[100 100 650 500])
saveas(fig,"P2_Synthetic_Data.eps",'epsc');


%% Define Gaussian prior

% % parameter1: C 
C_min  = -2000;          % lowerbound
C_max  = 5000;          % upperbound
C_prior_mean  = -2000+7000*rand   % prior mean
C_prior_stDev = 50*rand   % prior standard deviation

% parameter1: K
K_min  = -2000;          % lowerbound
K_max  = 5000;          % upperbound
K_prior_mean  = -2000+7000*rand   % prior mean
K_prior_stDev = 50*rand % prior standard deviation

% rearrange prior infos
prior_mean = [C_prior_mean K_prior_mean];
prior_var  = diag([C_prior_stDev K_prior_stDev].^2);
lobounds   = [C_min K_min];
upbounds   = [C_max K_max];


%% Define inputs of Metropolis: initial point and proposal variance

% start at prior means by default
MH_initial = [C_prior_mean  K_prior_mean];

% MH step size, assumed 1% of prior range here
proposal_variance = [(C_max-C_min)*0.05   (K_max-K_min)*0.05 ].^2;  


%% Define log-likelihood function

misfit = @(sample) (forward(T_true,sample) - e_data_mean).^2 ./ (e_data_stDev.^2) ;
loglike = @(sample) -sum( 0.5*misfit(sample) );


%% Construct a Metropolis chain of length N

N = 1000;                     
Nparam = length(MH_initial);

% initilize parameter samples
samples = zeros(N,Nparam)*nan; 

% initial sample at MH start
samples(1,:) = MH_initial;

loglike_current = loglike(samples(1,:));
logprior_current = log(  mvnpdf(samples(1,:), prior_mean, prior_var)  );

for i = 2:N
    
    % proposed a parameter sample
    % resample if the sample is not in prior bounds
    in_range = 0; 
    while ~in_range
        proposed_sample = samples(i-1,:) + mvnrnd( [0,0],  diag(proposal_variance)  );
        if sum(  (proposed_sample>upbounds) + (proposed_sample<lobounds)   ) == 0
            in_range = 1;
        end            
    end
    
    loglike_proposal = loglike(proposed_sample);
    logprior_proposal = log(  mvnpdf(proposed_sample, prior_mean, prior_var)  );
    
    
    % acceptance ratio
    accept =  exp(  (loglike_proposal + logprior_proposal) - (loglike_current + logprior_current)  );
    
    if  rand(1) < min(1,accept)
        % accept the proposed sample and update the current infos
        samples(i,:) = proposed_sample;
        loglike_current  = loglike_proposal;
        logprior_current = logprior_proposal;
            
    else
        % reject proposal, stay at the same sample
        samples(i,:) = samples(i-1,:);
        
    end  
end

% burn-in period: 20% of the samples
burn_in = round(N*0.2);


%% Plot MH results


% marginal paths of the chain
fig = figure();
title('Marginal paths')
subplot(1,2,1)
plot(samples(burn_in+1:end,1)); hold on;
plot(0,C, '*', 'Markersize', 15,'LineWidth',3);
set(gca,'FontSize',20)
xlabel('Chain iteration')
ylabel('C')

% set(gca,'FontSize',20)

subplot(1,2,2)
plot(samples(burn_in+1:end,2)); hold on;
plot(0,K, '*', 'Markersize', 15,'LineWidth',3);
xlabel('Chain iteration')
ylabel('K')
set(gca,'FontSize',20)
set(gcf,'Position',[100 100 900 500])
saveas(fig,"P2_Marginal_Paths.eps",'epsc');
% set(gca,'FontSize',20)

% marginal posterior distribution
fig = figure();
title('Marginal posterior distributions of parameters')
subplot(1,2,1)
histogram(samples(burn_in+1:end,1),50,'edgecolor','r','Normalization','pdf'); hold on;
plot(C,0,'*', 'Markersize', 15,'LineWidth',3);
Clocs = linspace(lobounds(1),upbounds(1),200);
plot(Clocs,normpdf(Clocs, prior_mean(1), sqrt(prior_var(1,1))  ),'LineWidth',5);
xlim([C_min, C_max])
xlabel('C')
legend('Posterior Samples','True','Prior','location','northwest')
set(gca,'FontSize',20)
% set(gca,'FontSize',20)

subplot(1,2,2)
histogram(samples(burn_in+1:end,2),50,'edgecolor','r','Normalization','pdf'); hold on;
plot(K,0,'*', 'Markersize', 15,'LineWidth',3);
Klocs = linspace(lobounds(2),upbounds(2),200);
plot(Klocs,normpdf(Klocs, prior_mean(2), sqrt(prior_var(2,2))  ),'LineWidth',5);
xlim([K_min, K_max])
xlabel('K')
legend('Posterior Samples','True','Prior','location','northwest')
% set(gca,'FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[100 100 900 500])
saveas(fig,"P2_Marginal_Posterior.eps",'epsc');
% print(gcf,'-depsc2',["P2_Marginal_Posterior" '.eps']);


% posterior samples
fig = figure();
plot(C, K, 'b*', 'Markersize', 15,'LineWidth',3)
hold on
plot(samples(burn_in+1:end,1), samples(burn_in+1:end,2), 'ro', 'Markersize', 7);
xlabel('C'); ylabel('K')
legend('True Parameters','Sample Points')
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
saveas(fig,"P2_Random_Walk_Sampling.eps",'epsc');
% print(gcf,'-depsc2',["P2_Random_Walk_Sampling" '.eps']);

% autocorelation
fig = figure();
nlags =30;
[ACF_C,lags,bounds] = autocorr(samples(1:end,1), nlags, 0);
[ACF_K,lags,bounds]= autocorr(samples(1:end,2), nlags, 0);
plot(lags,ACF_C,'bo-',lags,ACF_K,'r*-','linewidth',3);
ylabel('Autocorrelation');
xlabel('Lag');
grid on;
legend('C','K');
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
saveas(fig,"P2_Autocorelation.eps",'epsc');
% print(gcf,'-depsc2',["P2_Autocorelation" '.eps']);

% Covariance and Correlation Matrices
cov_matrix  = cov(samples);
corr_matrix = corr(samples);



%% Forward uncertianty propagation

fwds = [];
for i = burn_in+1:N
    fwds(end+1,:) = forward(T_true,samples(i,:));
end

fwd_mean = mean(fwds);
fwd_stDev  = std(fwds);
    
    
% plot data and model prediction
fig = figure();
errorbar(T_true, e_data_mean, e_data_stDev, 'ob', 'LineWidth', 1);hold on;
plot(T_true,fwd_mean+fwd_stDev, 'k-', 'LineWidth', 4);
plot(T_true,fwd_mean, 'r--', 'LineWidth', 3);
plot(T_true,fwd_mean-fwd_stDev, 'g', 'LineWidth', 2);

% xlim([25,350])
xlabel('Stress, T')
ylabel('Strain, e')
legend('Data','Model (upper)','Model (mean)', 'Model (lower)',...
    'location','northwest')
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
saveas(fig,"P2_Model_Prediction.eps",'epsc');
% print(gcf,'-depsc2',["P2_Model_Prediction" '.eps']);
