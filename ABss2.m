function loglike=ABss2(x,data)
% evaluate the sum of squares
% test case "A<->B"
tdata = data.tdata;
ydata = data.ydata;
misfit = (ABfun(tdata,x) - ydata).^2 ./ (data.std.^2) ;
loglike = -sum( 0.5*misfit );
% y = -sum((ydata-s).^2./data.std.^2)*0.5;