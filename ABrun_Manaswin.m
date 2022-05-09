% Example "A<->B"

% Model function: ABfun.m
% Sum of squares function: ABss.m
% Prior function ABprior.m
% Jacobian function: ABjac.m
clear;clc;
% select the method used
method = ['dram']; % 'mh','am','dr', or 'dram', see below

% data
clear model data params options
% data.tdata = [2 4 6 8 10]';
% data.ydata = [0.661 0.668 0.663 0.682 0.650]';

data.tdata = linspace(10, 66, 15); % x
data.ydata = [96.1 80.12 67.66 57.96 50.90 44.84 39.75 36.16 ...
    33.31 31.15 29.28 27.88 27.18 26.40 25.86]; % T
data.std = [0.2 0.5 0.8 0.45 0.32 0.15 0.7 0.65 ...
    0.54 0.48 0.84 0.56 0.74 0.36 0.75];

% LSQ fit and its Jacobian 
k0   = [-18 0]; %initial guess [phi, h]
[kopt, rss] = fminsearch(@ABss,k0,[],data);
xjac = jacob(@ABfun, data.tdata, kopt);      % numerical Jacobian
% J    = ABjac(data.tdata,kopt);             % analytical Jacobian
n=length(data.tdata); p = length(k0);
mse = rss/(n-p); % mean squared error

% J'*J is numerically singular, so increase the smallest singular value
[u,s,v]=svd(xjac); s=diag(s);s(s<1e-3)=1e-3;
% covariance matrix of the lsq estimates, used as initial proposal
cmat = v'*diag(1./s.^2)*v*mse;

% parameters for mcm
switch method
 case 'mh'
   nsimu    = 3000;
   drscale  = 0;
   adaptint = 0;
 case 'dr'
  nsimu    = 3000;
  drscale  = 2; 
  adaptint = 0;
 case 'am'
  nsimu    = 3000;
  drscale  = 0; 
  adaptint = 100;
 case 'dram'
  nsimu    = 10000;
  drscale  = 0; 
  adaptint = 0;
end

% create input arguments for the dramrun function

model.ssfun    = @ABss2;
model.priorfun = @ABprior;

params.par0    = [-18 0]; % initial parameter values
params.n       = n;    % number of observations
params.sigma2  = 0.5;  % prior for error variance sigma^2
params.n0      = 1;    % prior accuracy for sigma^2
params.bounds  = [-40 -5; 0 5];
	
params.parmu0   =  [-18 0];%params.bounds(1,:)+(params.bounds(2,:)-params.bounds(1,:)).*rand(1,2); %;            % prior mean of theta
params.parsig0  = [4 0.1];            % prior std of theta

options.nsimu    = nsimu;               % size of the chain
options.adaptint = adaptint;            % adaptation interval
options.drscale  = drscale;
options.qcov     = diag(params.parsig0.^2);%cmat.*2.38^2./p;      % initial proposal covariance 
% options.qcov     = [1 0.9;0.9 1];

% run the chain
[results,chain, s2chain] = dramrun(model,data,params,options);

%tau=iact(chain); % Integrated Autocorrelation Time

% figure(1);clf
% t = linspace(10, 66, 15);
% plot(data.tdata,data.ydata,'o',t,ABfun(t,kopt),'-')
% legend('data','LSQ estimate')
% 
% figure(2);clf
% plot(chain(:,1),chain(:,2),'.')
% xlabel('k_1'); ylabel('k_2');title('MCMC chain');

% add 95% ellipses of the proposal to the plot
%hold on;axis manual
%ellipse(kopt+[100 100],cmat*6,'Linewidth',1,'Color','black')
%ellipse(mean(chain)-[50,0],results.R'*results.R*6,'Linewidth',1,'Color','red')
%hold off; axis normal

figure(3);clf
subplot(2,1,1)
plot(chain(:,1),'.');ylabel('k_1');
title(sprintf('%s chain. Accepted %.1f%%',upper(method),results.accepted*100))
%title(sprintf('\\tau = %.1f',tau(1)));
subplot(2,1,2)
plot(chain(:,2),'.');ylabel('k_2');
%title(sprintf('\\tau = %.1f',tau(2)));
N = nsimu;
burn_in = round(N*0.2);

lobounds = params.bounds(1,:);
upbounds =params.bounds(2,:);
prior_mean = params.parmu0;
prior_var = params.parsig0;
figure()
title('Marginal posterior distributions of parameters')
subplot(1,2,1)
histogram(chain(burn_in+1:end,1),50,'edgecolor','k','Normalization','pdf'); hold on;
plot(params.parmu0(1),0,'*', 'Markersize', 15,'LineWidth',3);
Elocs = linspace(lobounds(1),upbounds(1),200);
plot(Elocs,normpdf(Elocs, prior_mean(1), prior_var(1)  ),'LineWidth',5);
% xlim([phi_min, phi_max])
xlabel('phi')
legend('Posterior Samples','True','Prior','location','northwest')
set(gca,'FontSize',20)

subplot(1,2,2)
histogram(chain(burn_in+1:end,2),50,'edgecolor','k','Normalization','pdf'); hold on;
plot(params.parmu0(2),0,'*', 'Markersize', 15,'LineWidth',3);
Tylocs = linspace(lobounds(2),upbounds(2),200);
plot(Tylocs,normpdf(Tylocs, prior_mean(2), prior_var(2)  ),'LineWidth',5);
% xlim([h_min, h_max])
xlabel('h')
legend('Posterior Samples','True','Prior','location','northwest')
set(gca,'FontSize',20)

fig = figure();
nlags =30;
[ACF_C,lags,bounds] = autocorr(chain(1:end,1), nlags, 0);
[ACF_K,lags,bounds]= autocorr(chain(1:end,2), nlags, 0);
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
cov_matrix  = cov(chain);
corr_matrix = corr(chain);

fwds = [];
for i = burn_in+1:N
    fwds(end+1,:) = ABfun(data.tdata, chain(i,:));
end

fwd_mean = mean(fwds);
fwd_stDev  = std(fwds);
    
    
% plot data and model prediction
fig = figure();
errorbar(data.tdata, data.ydata, data.std, 'ob', 'LineWidth', 1);hold on;
plot(data.tdata,fwd_mean+fwd_stDev, 'k-', 'LineWidth', 4);
plot(data.tdata,fwd_mean, 'r--', 'LineWidth', 3);
plot(data.tdata,fwd_mean-fwd_stDev, 'g', 'LineWidth', 2);

% xlim([25,350])
xlabel('Stress, T')
ylabel('Strain, e')
legend('Data','Model (upper)','Model (mean)', 'Model (lower)',...
    'location','northwest')
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
saveas(fig,"P2_Model_Prediction.eps",'epsc');