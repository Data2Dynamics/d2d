% arReset([silent])
% 
% Reset parameters ar.p=ar.pTrue and updates the merrit
% function accordingly by calling arCalcMerit
%
%   silent  A logical variable that indicates whether output at the command line 
%           should be prevented [false]
%               false: no output
%               true:  arCalcMerit prints to the command line
% 
% See also arCalcMerit, arEvaluate

function arReset(silent)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if(nargin==0)
    silent = false;
end

if(isfield(ar, 'pTrue'))
    ar.p = ar.pTrue;
    
    if(silent)
        arCalcMerit(false);
    else
        arCalcMerit;
    end
else
    ar.pTrue = ar.p;
end