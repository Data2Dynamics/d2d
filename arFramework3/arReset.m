% Reset parameter settings to simulated values
%
% arReset(silent)

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
        arChi2(false);
    else
        arChi2;
    end
else
    ar.pTrue = ar.p;
end