% arIn = arCompress([dumpsyms])
% 
% Removes some fields from ar to save memory. This function is called by
% arSave.
% 
%   dumpsyms    [0]
%               1: some symbolic variables are also removed
% 
% Examples
% arSimu(true,true,true)              % calculate sensitivities in fine grid
% ar.model(1).condition(1).sxFineSimu % show sensitivities
% arCompress                          % compression removes them
% ar.model(1).condition(1).sxFineSimu % check removal
% 
% See also arSave

function arIn = arCompress( dumpsyms )

global ar

if(nargout > 0)
    arIn = ar;
end

if ~exist( 'dumpsyms', 'var' )
    dumpsyms = 0;
end

% we can throw away these fields since arSimu initializes them again 
% before running the c code

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu = [];
        ar.model(m).condition(c).svFineSimu = [];
        ar.model(m).condition(c).sxFineSimu = [];
        ar.model(m).condition(c).szFineSimu = [];
        if ( dumpsyms )
            ar.model(m).condition(c).sym = [];
        end
    end
    if ( isfield( ar.model(m), 'ss_condition' ) )
        for c=1:length(ar.model(m).ss_condition)
            ar.model(m).ss_condition(c).suFineSimu = [];
            ar.model(m).ss_condition(c).svFineSimu = [];
            ar.model(m).ss_condition(c).sxFineSimu = [];
            ar.model(m).ss_condition(c).szFineSimu = [];
            if ( dumpsyms )
                ar.model(m).ss_condition(c).sym = [];
            end
        end    
    end
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = [];
            ar.model(m).data(d).systdFineSimu = [];
            if ( dumpsyms )
                ar.model(m).sym = [];
            end
        end
    end
end

% Invalidate cache so simulations do not get skipped
arCheckCache(1);