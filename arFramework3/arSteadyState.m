% Add steady state pre-simulation
%
% Usage:
%   arSteadyState((ar), model, conditionSS, conditionAffected, (states), (tstart) )
%
%   Adds a pre-equilibration to a number of conditions. ConditionSS refers
%   to the condition used for equilibration to steady state (source).
%   Events beloning to the steady state simulation condition are ignored
%   when equilibrating. Input functions will behave as they were for the
%   source condition. In many cases, it is therefore desirable to use a
%   steady state condition with no inputs.
%
% For the target conditions, note that any event occuring immediately at
% the start of the simulation will be overwritten by the event that
% initializes this condition to steady state. If an event immediately 
% upon start is required (for example, concentration = 2*steadystate), 
% consider adding time points before the desired event.
%
% Advanced use:
%   By default, all values of a steady state equilibration are copied into 
%   the destination initial condition. In some cases, this is undesirable.
%   If you wish to omit states; supply a cell array of state names that are
%   to be omitted from the equilibration.
%
%   The optional argument tstart refers to the starting timepoint of the 
%   steady state simulation.
%
% Hint: If you wish to see how the conditions are set up, you can use 
% the command 'arShowDataConditionStructure'
%
% Hint: To clear arSteadyState's changes to the model, invoke
% "arClearEvents". Note that this also clears other events.
%
% Notes:
%    arSteadyState invokes arLink
%    arSteadyState is currently only supported for deterministic
%    simulations
%    Steady state simulations will not check whether a stable steady state
%    exists. If this is not the case, pre-equilibration will fail and
%    an error will be thrown

function ar = arSteadyState( varargin )

    global ar;
        
    if ( nargin < 3 )
        error( 'Function needs at least three arguments.' );
    end
    if ( isstruct(varargin{1}) )
        ar = varargin{1};
        varargin = varargin(2:end);
    end
    
    logCall( 'arSteadyState', varargin{:} );
    
    m       = varargin{1};
    cSS     = varargin{2};
    cTarget = varargin{3};
    if ( length(varargin) > 3 )
        varargin = varargin(4:end);
    end
    
    omissions = [];
    if iscell( varargin{1} )
        for a = 1 : length( varargin{1} )
            state = find(strcmp(ar.model.x, varargin{1}{a}));
            if ( isempty( state ) )
                error( 'Cannot find state %s in the supplied model', varargin{1}{a} );
            else
                omissions = union( omissions, state );
            end
        end
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            varargin = {};
        end
    end
    
    tstart = 0;
    if ( length(varargin)>0 )
        try
            tstart = varargin{1};
        catch
            error( 'Invalid argument passed for tstart' );
        end
    end
    
    if ( ischar(cTarget) )
        if ( strcmp( cTarget, 'all' ) )
            cTarget = 1 : length(ar.model(m).condition);
        else
            error( 'Invalid target condition. Specify either a number or ''all''.' );
        end
    end
        
    nStates = numel( ar.model(m).x );    
    
    % Set up the steady state condition
    % Note that the explicit manual copy is intentional since in R2013
    % structure copies tend to be shallow copies.
    origin                  = ar.model(m).condition(cSS);
    ss_condition.src        = cSS;
    ss_condition.fkt        = origin.fkt;
    ss_condition.checkstr   = origin.checkstr;
    ss_condition.p          = origin.p;
    ss_condition.pLink      = origin.pLink;
    ss_condition.status     = 0;
    ss_condition.dLink      = [];
    ss_condition.pNum       = zeros(size(origin.pNum));
    ss_condition.qLog10     = zeros(size(origin.pNum));
    ss_condition.uNum       = zeros(size(origin.uNum));
    ss_condition.vNum       = zeros(size(origin.vNum));    
    ss_condition.suNum      = zeros(size(origin.suNum));
    ss_condition.svNum      = zeros(size(origin.svNum));
    ss_condition.dvdxNum    = zeros(size(origin.dvdxNum));
    ss_condition.dvduNum    = zeros(size(origin.dvduNum));
    ss_condition.dvdpNum    = zeros(size(origin.dvdpNum));
    ss_condition.svNum      = zeros(size(origin.svNum));    
    ss_condition.y_atol     = origin.y_atol;
    ss_condition.tstart     = tstart;
    ss_condition.tstop      = 1;
    ss_condition.tFine      = [tstart,inf];
    ss_condition.t          = [tstart,inf];
    ss_condition.tEq        = nan;
    ss_condition.dLink      = [];
    ss_condition.tEvents    = [];
    ss_condition.modx_A     = [];
    ss_condition.modx_B     = [];
    ss_condition.modsx_A    = [];
    ss_condition.modsx_B    = [];
    ss_condition.has_tExp   = 0;
    ss_condition.qEvents    = 0;
    ss_condition.uFineSimu  = zeros(2, size(origin.uFineSimu,2));
    ss_condition.vFineSimu  = zeros(2, size(origin.vFineSimu,2));
    ss_condition.xFineSimu  = zeros(2, size(origin.xFineSimu,2));
    ss_condition.zFineSimu  = zeros(2, size(origin.zFineSimu,2));
    ss_condition.suFineSimu = [];
    ss_condition.svFineSimu = [];
    ss_condition.sxFineSimu = [];
    ss_condition.szFineSimu = [];
    
    ss_condition.qMS        = 0;
    ss_condition.start      = 0;
    ss_condition.stop       = 0;
    ss_condition.stop_data  = 0;
    ss_condition.dxdt       = zeros(size(origin.dxdt));
    ss_condition.ddxdtdp    = zeros(size(origin.ddxdtdp));
    ss_condition.ssLink     = cTarget;
    
    % Which states to map to the target condition (default = all)
    ssStates                = ones(1,nStates);
    ssStates(omissions)     = 0;
    ss_condition.ssStates   = ssStates;
    ss_condition.ssIgnore   = omissions;
    
    if ( isfield( ar.model(m), 'ss_condition' ) )
        ar.model.ss_condition(end+1) = ss_condition;
    else
        ar.model.ss_condition = ss_condition;
    end
    insertionPoint = length(ar.model.ss_condition);
    
    h = waitbar(0, 'Linking up steady state');
    % Link up the target conditions
    for a = 1 : length( cTarget )
        waitbar(a/length(cTarget), h, sprintf('Linking up steady state %d/%d', a, length(cTarget)));
        
        % Map the parameters from the ss condition to the target condition
        fromP   = ar.model(m).condition(cSS).p;
        toP     = ar.model(m).condition(cTarget(a)).p;
        map     = mapStrings( fromP, toP );
        
        % Certain parameters may be unmapped (do not exist in target condition).
        ar.model(m).condition(cTarget(a)).ssUnmapped = [];
        L = find( isnan(map) );
        if ~isempty( L )
            warning( 'Certain parameters in target condition not present in steady state reference!' );
            warning( 'The following sensitivities will *not* take the equilibration into account:' );
            arFprintf(1, '%s\n', ar.model(m).condition(cTarget(a)).p{L} );
            ar.model(m).condition(cTarget(a)).ssUnmapped = L;
            map(L) = 1;
        end
    
        nPars   = numel( toP );
        ar.model(m).condition(cTarget(a)).ssParLink = map;
        ar.model(m).condition(cTarget(a)).ssLink = insertionPoint;
        
        vals = zeros(1,nStates);
        sens = zeros(1,nStates,nPars);
        
        ar = arAddEvent(ar, m, cTarget(a), ...
            ar.model(m).condition(cTarget(a)).tstart, ...
            [1:nStates], vals, vals, sens, sens, false );
    end
    close(h);
    
    ar.ss_conditions = true;
    
    % The event addition requires linking (silent link)
    arLink(true);
    
    % Show steady state count (to prevent users from forgetting arClearEvents)
    cnt = 0;
    for a = 1 : length( ar.model )
        cnt = cnt + length(ar.model(m).ss_condition);
    end
    arFprintf(1, 'Number of steady state equilibrations: %d', cnt );
    
    % Show any errors
    arPrintSteadyState(m, 2);
    
    % Invalidate cache so simulations do not get skipped
    arCheckCache(1);
end

function ID = mapStrings( str1, str2 )
    ID = nan(length(str2), 1);
    for a = 1 : length( str2 )
        for b = 1 : length( str1 )
            if strcmp( str1(b), str2(a) )
                if isnan( ID(a) )
                    ID(a) = b;
                else
                    error( 'Duplicate parameter name in ar.model.condition.p' );
                end
            end
        end
    end
end

function logCall( fn, varargin )
    global ar;
    
    if ~isfield(ar, 'eventLog')
        ar.eventLog = {};
    end
    call = [fn '('];
    if length( varargin ) > 0
        call = sprintf('%s%s', call, mat2str(varargin{1}) );
    end
    for a = 2 : length( varargin )
        call = sprintf('%s, %s', call, mat2str(varargin{a}) );
    end
    call = [call ')'];
    ar.eventLog{length(ar.eventLog)+1} = call;
end