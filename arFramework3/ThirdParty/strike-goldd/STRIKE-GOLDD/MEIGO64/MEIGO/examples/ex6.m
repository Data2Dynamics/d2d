function [fval] = ex6(x)
% Ackley's optimization test function. Problem: many local minima, one
% steep global minimum fval = 0 at x = [0,...,0]

a = 20;
b = 0.2;
c = 2*pi;
dim = length(x);

fval = -a*exp(-b*sqrt(sum(x.^2)/dim)) - exp(sum(cos(c*x))/dim) + a + exp(1);

end

