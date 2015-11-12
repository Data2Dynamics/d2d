% Removes some fields from ar to save memory

function arIn = arCompress

global ar

if(nargout > 0)
    arIn = ar;
end

% we can throw away these fields since arSimu initializes them again 
% before running the c code

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu = [];
        ar.model(m).condition(c).svFineSimu = [];
        ar.model(m).condition(c).sxFineSimu = [];
        ar.model(m).condition(c).szFineSimu = [];
    end
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = [];
            ar.model(m).data(d).systdFineSimu = [];
        end
    end
end

% Remove vector with last simulated parameters, so that future simulations
% do not skip simulation
if isfield( ar, 'pLastSimulated' )
    ar = rmfield( ar, 'pLastSimulated' );
end