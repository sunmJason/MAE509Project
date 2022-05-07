function T = ht(x,param)
k = 2.37; T_amb = 21.29;
a = 0.95; b = 0.95; L = 70;
gmma = sqrt(2*(a+b)*param(2)/(a*b*k));
c1 = -param(1)/(k*gmma)*((exp(gmma*L)*(param(2)+k*gmma))/(exp(-gmma*L)*(param(2)-k*gmma)+exp(gmma*L)*(param(2)+k*gmma)));
c2 = param(1)/(k*gmma)+c1;
T = c1*exp(-gmma*x)+c2*exp(gmma*x)+T_amb;
end