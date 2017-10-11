% Perform Local Parameter Sensitivity Analysis
%
% Arguments:
%    obj          - Function which computes 'objective'
%    lpsaPars     - Parameters to compute sensitivities of
%
% Optional arguments: (these should be provided in pairs e.g., 'perturbation', 0.5)
%    perturbation - Perturbation to test (default = 0.5)
%    ids          - How many different objectives obj can return (default: 1)
%
%   obj is a function that takes 1 argument (an index, which iterates from
%   1 to Nfuncs).
%
%   Returns matrix of sensitivities

function S = LPSA( obj, lpsaPars, varargin )
    global ar;
   
    switches = { 'perturbation', 'ids' };
    extraArgs = [ 1, 1, 1, 1 ];
    description = { ...
        {'', 'Specified perturbation'}, ...
        {'', 'Specified call ids'} };
    
    opts = argSwitch( switches, extraArgs, description, 1, varargin );
    
    ids = opts.ids_args;
    if isempty(ids)
        ids = 1;
    end    
    perturbation = opts.perturbation_args;
    if isempty(perturbation)
        perturbation = 0.5;
    end
    
    %% LPSA
    lpsaFunc = @(ids)arrayfun(obj, ids);

    % Get reference simulation
    arSimu(false, true, true);
    ref = lpsaFunc(ids) + 0;
    
    pRefLin = ar.p;
    pRefLin(ar.qLog10==1) = 10.^pRefLin(ar.qLog10==1);

    arPush;
    S = zeros(numel(lpsaPars), numel(ids));
    for a = 1 : numel( lpsaPars )
        arPush;
        curPar = lpsaPars(a);
        
        % Reduce parameter by 50%
        if ar.qLog10(curPar)
            ar.p(curPar) = ar.p(curPar) + log10(perturbation);
            pLin = 10^ar.p(curPar);
        else
            ar.p(curPar) = ar.p(curPar) * perturbation;
            pLin = ar.p(curPar);
        end
        arSimu(false, true, true);

        perturb = lpsaFunc(1 : length(ids)) + 0;
        S(a,:) = ((perturb-ref)./ref) ./ ((pLin-pRefLin(curPar))./(pRefLin(curPar)));

        arPop;
    end
    arPop;
    
end