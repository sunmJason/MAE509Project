function T=ABfun(x,input)
% model function for "A<->B" example

% k1 = b(1);
% k2 = b(2);
% 
% A0 = 1;
% B0 = 0;
% 
% kk = k2/(k1+k2);
% 
% y = kk + (A0-kk).*exp(-(k1+k2).*x);

phi = input(1);
h = input(2);
a = 0.95; b = 0.95;
L = 70; 
k = 2.37;
Tamb = 21.29;

gamma = sqrt(2*(a+b)*h/a/b/k);
c1 = -phi/k/gamma * (exp(gamma*L)*(h+k*gamma))/...
    ((exp(-gamma*L)*(h-k*gamma) + exp(gamma*L)*(h+k*gamma)));
c2 = phi/k/gamma + c1;

T = c1 * exp(-gamma * x) + c2 * exp(gamma * x) + Tamb;
