clc;clear;


mean1 = optimizableVariable('mean1',[0,30]);
var1 = optimizableVariable('var1',[0,20]);
mean2 = optimizableVariable('mean2',[0,10]);
var2 = optimizableVariable('var2',[0,10]);
vars = [mean1,var1,mean2,var2];

results = bayesopt(@(vars)AutoOpt(vars),vars,'MaxObjectiveEvaluations',50')
