function y=ABss2(x,data)
% evaluate the sum of squares
% test case "A<->B"

global SS_count
SS_count = SS_count+1;

tdata = data.tdata;
ydata = data.ydata;

s = ABfun(tdata,x);
y = -sum((ydata-s).^2./data.std.^2)*0.5;