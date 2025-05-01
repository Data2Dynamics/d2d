% ar = arAddEvent([ar], model, condition, timepoints, [statename], [A], [B],  [sA], [sB])
%
% Add solver re-initialization point
%
%   ar            ar struct, if not specified the global one is modified
%   model         index of model
%   conditions    indicies of conditions in which the event is set, 'all' -> all conditions
%   timepoints    timepoints when a event is added
%   statename     string, neccessary when the concentration of a state is changed
%   A             change value of state with "statename" to Ax+B where x represents
%   B             the old state variable before the event.
%   sA            sensitivity matrix
%   sB            sensitivity matrix
%
% Output variable
%   ar            ar struct, if only local ar should be changed
%
% Adds an event to a condition. Event is added for a specific condition at
% a specific time point. If you wish to see how the conditions are set up,
% you can use the command 'arShowDataConditionStructure'. Specifying 'all'
% for the conditions sets the event for all model conditions.
%
% If a change at the event is required, also specify a state name, and a
% change. Event changes are of the form Ax+B; where x represents the old
% state variable before the event.
%
% Advanced use:
%  If one wishes to set the entire state vector and sensitivity vector, it
%  is also possible to supply a (1 x nStates) vector for A and B; and a
%  (1 x nStates x nParameter) matrix for sA and sB. Note that no checks are
%  made whether the supplied sensitivity matrix is consistent with the state
%  change!
%
% Example usages:
%   arAddEvent(ar, model, condition, timepoint )
%   arAddEvent(model, condition, timepoint )
%   arAddEvent(model, condition, timepoint , statename, A, B)
%   arAddEvent(model, condition, timepoint , statename, A, B, sA, sB)
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
if isempty(varargin{1})
    error('First variable has to be an ar struct of a model index m.')
end

% If we were called from a high level event function, since the higher
% level function will already be in the command log and this one doesn't
% need to be added.
s = dbstack(1);
if ( (length(s)==0) || ( (~strcmp( s(1).file, 'arSteadyState.m' )) && (~strcmp( s(1).file, 'arFindInputs.m' )) ) )
    logCall( 'arAddEvent', varargin{:} );
end

m       = varargin{1};
condis  = varargin{2};
tL      = varargin{3};

if ( ischar(condis) && strcmpi( condis, 'all' ) )
    condis = 1 : length( ar.model(m).condition );
end

if isempty( condis )
    error( 'Attempted to call arAddEvent with an empty list for the conditions to apply the event to' );
end

nStates = numel( ar.model(m).x );

for jC = 1 : numel( condis )
    c = condis(jC);
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
        %
        if ~isnumeric(A)
            A = str2num(A);
        end
        % if ~isnumeric(B)
        %     if ~isempty(str2num(B)) 
        %         B = str2num(B);
        %     end
        % end
        % %        
        % if isnumeric(A)
        %     A = num2str(A);
        % end
        % if isnumeric(B)
        %     B = num2str(B);
        % end
    end
    
    manualSensitivity = 0;
    if ( length( varargin ) > 6 )
        manualSensitivity = 1;
        sA = varargin{7};
        sB = varargin{8};
        try
            [sAsize_1, sAsize_States, sAsize_Pars]  = size(sA);
            [sBsize_1, sBsize_States, sBsize_Pars]  = size(sB);
            
            if ( sum( [sAsize_1, sAsize_States, sAsize_Pars] == [1, nStates, nPars] ) < 3 )
                error( 'Sensitivity matrix A has the wrong dimensions. Should be 1 x nStates x nPars' );
            end
            if ( sum( [sBsize_1, sBsize_States, sBsize_Pars] == [1, nStates, nPars] ) < 3 )
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
    
    for nt = 1 : length( tL )
        t = tL(nt);
        
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
        
        % Add the point to the tExtra list so that the plots look good
        % (otherwise the plot interpolation might jump over the event
        % point which looks weird and worries users)
        if ( isfield( ar.model(m).condition(c), 'dLink' ) )
            for jd = ar.model(m).condition(c).dLink
                if isfield( ar.model(m).data(jd), 'tExtra' )
                    ar.model(m).data(jd).tExtra = union(ar.model(m).data(jd).tExtra, [t, t - 1e-6 * t]);
                else
                    ar.model(m).data(jd).tExtra = [t, t - 1e-6 * t];
                end
            end
        end
    end
end

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