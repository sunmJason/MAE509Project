function out = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name)
data.tdata = linspace(10, 66, 15); % x
data.ydata = [96.1 80.12 67.66 57.96 50.90 44.84 39.75 36.16 ...
    33.31 31.15 29.28 27.88 27.18 26.40 25.86]; % T
data.std = [0.2 0.5 0.8 0.45 0.32 0.15 0.7 0.65 ...
    0.54 0.48 0.84 0.56 0.74 0.36 0.75];

k0 = prior_mean;rss = 0.5; kopt=k0;
n=length(data.tdata); p = length(k0);
if use_lsq == 1
    k0   = [-18 0];
    [kopt, rss] = fminsearch(@ABss,k0,[],data);
    xjac = jacob(@ABfun, data.tdata, kopt); 
    mse = rss/(n-p); % mean squared error
    % J'*J is numerically singular, so increase the smallest singular value
    [u,s,v]=svd(xjac); s=diag(s);s(s<1e-3)=1e-3;
    % covariance matrix of the lsq estimates, used as initial proposal
    cmat = v'*diag(1./s.^2)*v*mse;
    qcov = cmat.*2.38^2./p;
else
    qcov = diag(prior_std.^2);
end

% create input arguments for the dramrun function

model.ssfun    = @ABss2;
model.priorfun = @ABprior;

params.par0    = k0; % initial parameter values
params.n       = n;    % number of observations
params.sigma2  = rss;  % prior for error variance sigma^2
params.n0      = 1;    % prior accuracy for sigma^2
params.bounds  = [-40 -5; 0 5];
	
params.parmu0   =  prior_mean;%params.bounds(1,:)+(params.bounds(2,:)-params.bounds(1,:)).*rand(1,2); %;            % prior mean of theta
params.parsig0  = prior_std;            % prior std of theta

options.nsimu    = nsimu;               % size of the chain
options.adaptint = adaptint;            % adaptation interval
options.drscale  = drscale;
options.qcov     = qcov;      % initial proposal covariance 
% options.qcov     = [1 0.9;0.9 1];

% run the chain
[results,chain, s2chain] = dramrun(model,data,params,options);

fig=figure();
subplot(2,1,1)
plot(chain(:,1),'.');ylabel('k_1');
title(sprintf('%s chain. Accepted %.1f%%',name,results.accepted*100))
%title(sprintf('\\tau = %.1f',tau(1)));
subplot(2,1,2)
plot(chain(:,2),'.');ylabel('k_2');
fig.Position = [100 100 900 700];
%title(sprintf('\\tau = %.1f',tau(2)));
saveas(fig,strcat('Figures/Chain_',name,".png"),'png');
N = nsimu;
burn_in = round(N*0.2);
lobounds = params.bounds(1,:);
upbounds =params.bounds(2,:);
prior_mean = params.parmu0;
prior_var = params.parsig0;
fig=figure();
title(strcat('Marginal posterior distributions -',name))
subplot(1,2,1)
histogram(chain(burn_in+1:end,1),50,'edgecolor','k','Normalization','pdf'); hold on;
plot(params.parmu0(1),0,'*', 'Markersize', 15,'LineWidth',3);
Elocs = linspace(lobounds(1),upbounds(1),200);
plot(Elocs,normpdf(Elocs, prior_mean(1), prior_var(1)  ),'LineWidth',5);
% xlim([phi_min, phi_max])
xlabel('phi');
legend('Posterior Samples','True','Prior','location','northwest');
set(gca,'FontSize',20)

subplot(1,2,2)
histogram(chain(burn_in+1:end,2),50,'edgecolor','k','Normalization','pdf'); hold on;
plot(params.parmu0(2),0,'*', 'Markersize', 15,'LineWidth',3);
Tylocs = linspace(lobounds(2),upbounds(2),200);
plot(Tylocs,normpdf(Tylocs, prior_mean(2), prior_var(2)  ),'LineWidth',5);
% xlim([h_min, h_max])
xlabel('h');
legend('Posterior Samples','True','Prior','location','northwest');
set(gca,'FontSize',20);
fig.Position = [100 100 900 700];
saveas(fig,strcat('Figures/Marginal_Distribution_',name,".png"),'png');

fig = figure();
nlags =30;
[ACF_C,lags,bounds] = autocorr(chain(1:end,1), nlags, 0);
[ACF_K,lags,bounds]= autocorr(chain(1:end,2), nlags, 0);
plot(lags,ACF_C,'bo-',lags,ACF_K,'r*-','linewidth',3);
title(strcat('Autocorrelation-',name))
ylabel('Autocorrelation');
xlabel('Lag');
grid on;
legend('C','K');
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
fig.Position = [100 100 900 700];
saveas(fig,strcat('Figures/Autocorrelation_',name,".png"),'png');
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
title(strcat('Forward Uncertainty -',name))
% xlim([25,350])
xlabel('Stress, T')
ylabel('Strain, e')
legend('Data','Model (upper)','Model (mean)', 'Model (lower)',...
    'location','northwest')
% set(gca,'FontSize',24)
set(gca,'FontSize',20)
fig.Position = [100 100 900 700];
saveas(fig,strcat('Figures/Forward Uncertainty_',name,".png"),'png');
m1 = mean(chain(burn_in+1:end,:));
s1 = mean(std(chain(burn_in+1:end,:)));
out = [ACF_C(end),ACF_K(end),s1,m1(1),m1(2)];