clc;clear;close all; 
true_res = [-18.41 0.00191];
nsimu= 10000; 
lb= [-40 -5]; 
ub =  [ 0 5];
set(0,'DefaultFigureVisible','off');
% for i=1:2
% name = strcat("MH Random Trial",num2str(i));
% drscale = 0; 
% adaptint = 0 ;
% prior_mean = lb + (ub-lb).*rand(1,2); 
% prior_std =  [10 0.1] + ([60 10]-[10 0.1]).*rand(1,2);
% use_lsq=0;
% out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% median_autocorr_phi_mh = median(out(:,1)); 
% median_autocorr_h_mh = median(out(:,2)); 
% median_prior_std_mh = median(out(:,3));
% median_rmse_mh = median(rmse);
% best_autocorr_phi_mh = min(out(:,1)); 
% best_autocorr_h_mh = min(out(:,2)); 
% best_prior_std_mh = min(out(:,3));
% best_rmse_mh = min(rmse);
% 
% for i=1:10
% name = strcat("DR Random Trial",num2str(i));
% drscale = 10; 
% adaptint = 0 ;
% prior_mean = lb + (ub-lb).*rand(1,2); 
% prior_std =  [10 0.1] + ([60 10]-[10 0.1]).*rand(1,2);
% use_lsq=0;
% out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% 
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% median_autocorr_phi_dr= median(out(:,1)); 
% median_autocorr_h_dr = median(out(:,2)); 
% median_prior_std_dr = median(out(:,3));
% median_rmse_dr = median(rmse);
% best_autocorr_phi_dr = min(out(:,1)); 
% best_autocorr_h_dr = min(out(:,2)); 
% best_prior_std_dr = min(out(:,3));
% best_rmse_dr = min(rmse);
% 
% for i=1:10
% name = strcat("AM Random Trial",num2str(i));
% drscale = 0; 
% adaptint = 100 ;
% prior_mean = lb + (ub-lb).*rand(1,2); 
% prior_std =  [10 0.1] + ([60 10]-[10 0.1]).*rand(1,2);
% use_lsq=0;
% out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% median_autocorr_phi_am = median(out(:,1)); 
% median_autocorr_h_am = median(out(:,2)); 
% median_prior_std_am = median(out(:,3));
% median_rmse_am = median(rmse);
% best_autocorr_phi_am = min(out(:,1)); 
% best_autocorr_h_am = min(out(:,2)); 
% best_prior_std_am = min(out(:,3));
% best_rmse_am = min(rmse);
% 
% 
% for i=1:10
% name = strcat("DRAM Random Trial",num2str(i));
% drscale = 10; 
% adaptint = 100 ;
% prior_mean = lb + (ub-lb).*rand(1,2); 
% prior_std =  [10 0.1] + ([60 10]-[10 0.1]).*rand(1,2);
% use_lsq=0;
% out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% 
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% median_autocorr_phi_dram = median(out(:,1)); 
% median_autocorr_h_dram  = median(out(:,2)); 
% median_prior_std_dram  = median(out(:,3));
% median_rmse_dram  = median(rmse);
% best_autocorr_phi_dram  = min(out(:,1)); 
% best_autocorr_h_dram  = min(out(:,2)); 
% best_prior_std_dram  = min(out(:,3));
% best_rmse_dram  = min(rmse);
% %% Least Squares Results
% 
% for i=1
% name = strcat("MH Least Squares",num2str(i));
% drscale = 0; 
% adaptint = 0 ;
% prior_mean = [-18 0]; 
% prior_std =  [35 5];
% use_lsq=1;
% out = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% best_autocorr_phi_mh_ols = out(:,1); 
% best_autocorr_h_mh_ols   = out(:,2); 
% best_prior_std_mh_ols   = out(:,3);
% best_rmse_mh_ols   = rmse;
% 
% for i=1
% name = strcat("DR Least Squares",num2str(i));
% drscale = 10; 
% adaptint = 0 ;
% prior_mean = [-18 0]; 
% prior_std =  [35 5];
% use_lsq=1;
% out = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% best_autocorr_phi_dr_ols = out(:,1); 
% best_autocorr_h_dr_ols   = out(:,2); 
% best_prior_std_dr_ols   = out(:,3);
% best_rmse_dr_ols   = rmse;
% 
% for i=1
% name = strcat("AM Least Squares",num2str(i));
% drscale = 0; 
% adaptint = 100 ;
% prior_mean = [-18 0]; 
% prior_std =  [35 5];
% use_lsq=1;
% out = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% 
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% best_autocorr_phi_am_ols = out(:,1); 
% best_autocorr_h_am_ols   = out(:,2); 
% best_prior_std_am_ols   = out(:,3);
% best_rmse_am_ols   = rmse;
% 
% for i=1
% name = strcat("DRAM Least Squares",num2str(i));
% drscale = 10; 
% adaptint = 100 ;
% prior_mean = [-18 0]; 
% prior_std =  [35 5];
% use_lsq=1;
% out = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
% end
% 
% rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
% best_autocorr_phi_dram_ols = out(:,1); 
% best_autocorr_h_dram_ols   = out(:,2); 
% best_prior_std_dram_ols   = out(:,3);
% best_rmse_dram_ols   = rmse;
% 
% 
% x = 1:4; 
% y = [mean([median_autocorr_phi_mh,median_autocorr_h_mh]),mean([best_autocorr_phi_mh,best_autocorr_h_mh]),mean([best_autocorr_phi_mh_ols,best_autocorr_h_mh_ols]);...
%     mean([median_autocorr_phi_dr,median_autocorr_h_dr]),mean([best_autocorr_phi_dr,best_autocorr_h_dr]),mean([best_autocorr_phi_dr_ols,best_autocorr_h_dr_ols]);...
%     mean([median_autocorr_phi_am,median_autocorr_h_am]),mean([best_autocorr_phi_am,best_autocorr_h_am]),mean([best_autocorr_phi_am_ols,best_autocorr_h_am_ols]);...
%     mean([median_autocorr_phi_dram,median_autocorr_h_dram]),mean([best_autocorr_phi_dram,best_autocorr_h_dram]),mean([best_autocorr_phi_dram_ols,best_autocorr_h_dram_ols])];
% fig = figure();
% b= bar(x,y);
% for i=1:3
%     xtips1 = b(i).XEndPoints;
%     ytips1 = b(i).YEndPoints;
%     labels1 = string(b(i).YData);
%     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
% set(gca,'FontSize',18)
% fig.Position = [100 100 1100 900];
% xticklabels(["Metropolis-Hastings" , "Delayed Rejection", "Adaptation", "Delayed Rejection and Adaptation"])
% legend(["Random -Median", "Random -Best", "Least Squares"]);
% ylabel("Mean Autocorrelation of \phi and h")
% saveas(fig,strcat('Figures/Autocorrelation_BarGraph',".png"),'png');
% 
% x = 1:4; 
% y = [median_prior_std_mh,best_prior_std_mh,best_prior_std_mh_ols;...
%      median_prior_std_dr,best_prior_std_dr,best_prior_std_dr_ols;...
%      median_prior_std_am,best_prior_std_am,best_prior_std_am_ols;...
%      median_prior_std_dram,best_prior_std_dram,best_prior_std_dram_ols];
% 
% fig = figure();
% b= bar(x,y);
% for i=1:3
%     xtips1 = b(i).XEndPoints;
%     ytips1 = b(i).YEndPoints;
%     labels1 = string(b(i).YData);
%     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
% set(gca,'FontSize',18)
% fig.Position = [100 100 1100 900];
% ylabel("Posterior Standard Deviation")
% xticklabels(["Metropolis-Hastings" , "Delayed Rejection", "Adaptation", "Delayed Rejection and Adaptation"])
% legend(["Random -Median", "Random -Best", "Least Squares"]);
% saveas(fig,strcat('Figures/Posterior_std_BarGraph',".png"),'png');
% 
% x = 1:4; 
% y = [median_rmse_mh,best_rmse_mh,best_rmse_mh_ols;...
%      median_rmse_dr,best_rmse_dr,best_rmse_dr_ols;...
%      median_rmse_am,best_rmse_am,best_rmse_am_ols;...
%      median_rmse_dram,best_rmse_dram,best_rmse_dram_ols];
% fig = figure();
% b= bar(x,y);
% for i=1:3
%     xtips1 = b(i).XEndPoints;
%     ytips1 = b(i).YEndPoints;
%     labels1 = string(b(i).YData);
%     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
% set(gca,'FontSize',18)
% fig.Position = [100 100 1100 900];
% ylabel("RMSE")
% xticklabels(["Metropolis-Hastings" , "Delayed Rejection", "Adaptation", "Delayed Rejection and Adaptation"])
% legend(["Random -Median", "Random -Best", "Least Squares"]);
% saveas(fig,strcat('Figures/Rmse_BarGraph',".png"),'png');

%% Parametric Testing 
dr = [2 4 8 16];
for i=1:4
    name = strcat("DR = ",num2str(dr(i)));
    drscale = dr(i); 
    adaptint = 0 ;
    prior_mean = [-18 0]; 
    prior_std =  [35 5];
    use_lsq=0;
    out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
end
rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
fig=figure();
plot(dr,0.5*(out(:,1)+out(:,2)),'-*')
xlabel("Delayed Rejection Scale")
ylabel("Mean Autocorrelation of \phi and h")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/DR_Scatter_Autocorrelation',".png"),'png');

fig=figure();
plot(dr,out(:,3),'-*')
xlabel("Delayed Rejection Scale")
ylabel("Posterior Standard Deviation")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/DR_Scatter_prior_std',".png"),'png');

fig=figure();
plot(dr,rmse,'-*')
xlabel("Delayed Rejection Scale")
ylabel("RMSE")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/DR_Scatter_RMSE',".png"),'png');

am = [10 50 100 500];
for i=1:4
name = strcat("AM = ",num2str(am(i)));
drscale = 0; 
adaptint = am(i) ;
prior_mean = lb + (ub-lb).*rand(1,2); 
prior_std =  [10 0.1] + ([60 10]-[10 0.1]).*rand(1,2);
use_lsq=0;
out(i,:) = MC_run(nsimu,drscale,adaptint,prior_mean,prior_std,use_lsq,name); 
end
rmse = sqrt(mean((out(:,4:5)-true_res).^2,2));
fig=figure();
plot(dr,0.5*(out(:,1)+out(:,2)),'-*')
xlabel("Iterations till Adaptation")
ylabel("Mean Autocorrelation of \phi and h")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/AM_Scatter_Autocorrelation',".png"),'png');

fig=figure();
plot(dr,out(:,3),'-*')
xlabel("Iterations till Adaptation")
ylabel("Posterior Standard Deviation")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/AM_Scatter_prior_std',".png"),'png');

fig=figure();
plot(dr,rmse,'-*')
xlabel("Iterations till Adaptation")
ylabel("RMSE")
set(gca,'FontSize',18)
fig.Position = [100 100 1100 900];
saveas(fig,strcat('Figures/AM_Scatter_RMSE',".png"),'png');