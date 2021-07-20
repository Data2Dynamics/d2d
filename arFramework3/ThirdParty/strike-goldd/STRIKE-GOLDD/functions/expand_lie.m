function [exp_e,exp_p] = expand_lie(pMax,p)
%% Expand, as a Taylor series, the paramter to be removed and exp(epsilon)
% Input:    pMax --> Maximum degree of Lie series
%           p --> Parameter
% Output:   exp_e --> Vector with coefficients of Taylor series of exp(epsilon)
%           exp_p --> Vector with coefficients of Taylor series of
%                     parameter
%
%%
    syms A epsilon
    [n,d]   = numden(p);

    temp1   = taylor(1/(1-A)-1, A, 'Order', pMax);
    temp2   = taylor(exp(epsilon)-1, epsilon, 'Order', pMax);
    
    exp_e   = children(temp2);
    exp_p   = children(temp1);
    if d==1
        exp_p   = subs(exp_p,A,1-1/p);
    else
        exp_p   = subs(exp_p,A,1-p);
    end
    
end

