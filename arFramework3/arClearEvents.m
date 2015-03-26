% Clears all steady state pre-simulations and events

function ar = arClearEvents( varargin )

    global ar;
        
    if ( nargin > 0 )
        if ( isstruct(varargin{1}) )
            ar = varargin{1};
            varargin = varargin(2:end);
        end
    end

    % Clear all events
    h = waitbar(0, 'Clearing events');
    for m = 1 : length( ar.model )
        if isfield(ar.model(m), 'condition' )
            for c = 1 : length( ar.model(m).condition )
                waitbar(0.9*c/length( ar.model(m).condition ), h, sprintf('Clearing events for model %d [condition %d/%d]', m, c, length( ar.model(m).condition )));
                ar.model(m).condition(c).modt = [];
                ar.model(m).condition(c).modx_A = [];
                ar.model(m).condition(c).modx_B = [];
                ar.model(m).condition(c).modx_B = [];
                ar.model(m).condition(c).modsx_B = [];
                ar.model(m).condition(c).tEvents = [];
            end
        end
        if isfield(ar.model(m), 'data' )
            for d = 1 : length( ar.model(m).data )
                ar.model(m).data(d).tEvents = [];
            end
        end   
    end
    
    % Wipe steady state conditions
    h = waitbar(0.9, h, 'Removing steady state conditions');
    if isfield(ar.model, 'ss_condition')
        ar.model = rmfield(ar.model,'ss_condition');
    end
    
    ar.ss_conditions = 0;
    h = waitbar(0.95, h, 'Linking model');
    % The event removal requires linking (silent link)
    arLink(true);
    
    close(h);
end