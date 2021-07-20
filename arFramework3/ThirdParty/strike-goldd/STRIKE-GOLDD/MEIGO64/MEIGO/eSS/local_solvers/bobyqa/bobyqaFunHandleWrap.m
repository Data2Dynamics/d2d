function [ y ] = bobyqaFunHandleWrap(x, calfun)
% This is a wrapper to compute the objective function value. 
% Called from within the mexbobyqa routine. 
% funHandleWrap currently uses a global variable in bobyqa.m, 
% the argument calfun is not used.
%
% Input:
%   x: parameter vector
%   calfun: objective function handle or string representation
%
% Output:
%   y: objective function value

fun = bobyqa();
y = fun(x);

end