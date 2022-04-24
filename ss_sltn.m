function [T] = ss_sltn(input)
    phi = input(1);
    h = input(2);
    a = 0.95; b = 0.95;
    L = 70; x = 20;
    k = 2.37;
    Tamb = 21.29;
    
    gamma = sqrt(2*(a+b)*h/a/b/k);
    c1 = -phi/k/gamma * (exp(gamma*L)*(h+k*gamma))/((exp(-gamma*L)*(h-k*gamma) + exp(gamma*L)*(h+k*gamma)));
    c2 = phi/k/gamma + c1;
    
    T = c1 * exp(-gamma * x) + c2 * exp(gamma * x) + Tamb;

end