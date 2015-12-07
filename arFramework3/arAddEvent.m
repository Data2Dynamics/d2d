% Add solver re-initialization point
%
% Usage:
%   arAddEvent((ar), model, condition, time point )
%
% Adds an event to a condition. Event is added for a specific condition at
% a specific time point. If you wish to see how the conditions are set up,
% you can use the command 'arShowDataConditionStructure'
%
% If a change at the event is required, also specify a state name, and a
% change. Event changes are of the form Ax+B; where x represents the old
% state variable before the event.
%
%  arAddEvent((ar), model, condition, time point, (state name), (A), (B))
%
% Advanced use:
% If one wishes to set the entire state vector and sensitivity vector, it 
% is also possible to supply a (1 x nStates) vector for A and B; and a
% (1 x nStates x nParameter) matrix for sA and sB. Note that no checks are
% made whether the supplied sensitivity matrix is consistent with the state
% change!
%
% arAddEvent((ar), model, condition, time point, (state name), (A), (B), (sA), (sB))
%
% Note:
%    arAddEvent invokes arLink; if this is undesirable, pass an extra
%    argument with false past the steady state arguments

function ar = arAddEvent( varargin )

    if ( nargin < 2 )
        error( 'Function needs at least two arguments.' );
    end
       
    global ar;
    if ( isstruct(varargin{1}) )
        ar = varargin{1};
        varargin = varargin(2:end);
    end

    % If we were called from a high level event function, since the higher
    % level function will already be in the command log and this one doesn't 
    % need to be added.
    s = dbstack(1);
    if ( (length(s)==0) || ( (~strcmp( s(1).file, 'arSteadyState.m' )) && (~strcmp( s(1).file, 'arFindInputs.m' )) ) )
        logCall( 'arAddEvent', varargin{:} );
    end
    
    m = varargin{1};
    c = varargin{2};
    t = varargin{3};

    nStates = numel( ar.model.x );
    nPars   = numel( ar.model(m).condition(c).p );    
    
    if ( length( varargin ) > 3 )
        if ischar( varargin{4} )
            state = find( strcmp(ar.model(m).x, varargin{4}) == 1 );
        else
            state = varargin{4};
        end
    end
    
    stateChanges = 0;
    if ( length( varargin ) > 4 )
        stateChanges = 1;
        A = varargin{5};
        if ( length( varargin ) > 5 )
            B = varargin{6};
        else
            B = 0;
        end
    end
    
    manualSensitivity = 0;
    if ( length( varargin ) > 6 )
        manualSensitivity = 1;
        sA = varargin{7};
        sB = varargin{8};
        try
            if ( sum( size(sA) == [1, nStates, nPars] ) < 3 )
                error( 'Sensitivity matrix A has the wrong dimensions. Should be 1 x nStates x nPars' );
            end
            if ( sum( size(sB) == [1, nStates, nPars] ) < 3 )
                error( 'Sensitivity matrix B has the wrong dimensions. Should be 1 x nStates x nPars' );
            end
        catch
            error( 'Sensitivity matrix has the wrong dimensions. Should be 1 x nStates x nPars' );
        end
    end
       
    doLink = true;
    if ( length( varargin ) > 8 )
        doLink = varargin{9};
    end    
    
    % Insert event
    [ar.model(m).condition(c).tEvents, ~, i] = union( ar.model(m).condition(c).tEvents, t );
    
    % Did we actually add a new event, or just modify a modifier?
    if (length(i) > 0)
        newEvent = 1;
    else
        newEvent = 0;
    end
    
    % Determine where the event was inserted
    I = find( ar.model(m).condition(c).tEvents > t, 1 )-1;
    if isempty(I)
        I = length( ar.model(m).condition(c).tEvents );
    end
    
    % Insertion
    modx_A_ins = ones(1, nStates);
    modx_B_ins = zeros(1, nStates);
    modsx_A_ins = ones(1, nStates, nPars);
    modsx_B_ins = zeros(1, nStates, nPars);    
    
    % Merge mod matrices if required
    if ( newEvent )
        ar.model(m).condition(c).modt     = union( ar.model(m).condition(c).modt, t );
        ar.model(m).condition(c).modx_A   = [ ar.model(m).condition(c).modx_A(1:I-1,:) ; modx_A_ins ; ar.model(m).condition(c).modx_A(I:end,:) ];
        ar.model(m).condition(c).modx_B   = [ ar.model(m).condition(c).modx_B(1:I-1,:) ; modx_B_ins ; ar.model(m).condition(c).modx_B(I:end,:) ];
        ar.model(m).condition(c).modsx_A  = [ ar.model(m).condition(c).modsx_A(1:I-1,:,:) ; modsx_A_ins ; ar.model(m).condition(c).modsx_A(I:end,:,:) ];
        ar.model(m).condition(c).modsx_B  = [ ar.model(m).condition(c).modsx_B(1:I-1,:,:) ; modsx_B_ins ; ar.model(m).condition(c).modsx_B(I:end,:,:) ];        
    end
    
    % Ready to add the event
    if ( stateChanges )
        ar.model(m).condition(c).modx_A(I,state) = A;
        ar.model(m).condition(c).modx_B(I,state) = B;
    
        if ( manualSensitivity )
            ar.model(m).condition(c).modsx_A(I,:,:) = sA;
            ar.model(m).condition(c).modsx_B(I,:,:) = sB;
        else
            % Only the multiplicative part of the sensitivity factors in
            ar.model(m).condition(c).modsx_A(I,state,:) = A;
        end
    end
    
    % Activate event system for this condition
    ar.model(m).condition(c).qEvents = 1;
    
    % Activate event handling in general
    ar.config.useEvents = 1;
    
    % Link (silent)!
    if (doLink)
        arLink(true);
    end
    
    % Invalidate cache so simulations do not get skipped
    arCheckCache(1);
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