clc;clear;


mean1 = optimizableVariable('mean1',[-40,1]);
var1 = optimizableVariable('var1',[0,40]);
mean2 = optimizableVariable('mean2',[0,20]);
var2 = optimizableVariable('var2',[0,10]);
vars = [mean1,var1,mean2,var2];

results = bayesopt(@(vars)AutoOpt_manaswin(vars),vars,'MaxObjectiveEvaluations',50')
