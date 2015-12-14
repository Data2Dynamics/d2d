% Clears all steady state pre-simulations and events

function arClearEvents( varargin )

    global ar;
    global arOutputLevel;

    % Clear all events
    if ( arOutputLevel > 1 )
        h = waitbar(0, 'Clearing events');
    end
    for m = 1 : length( ar.model )
        if isfield(ar.model(m), 'condition' )
            for c = 1 : length( ar.model(m).condition )
                if ( arOutputLevel > 1 )
                    waitbar(0.9*c/length( ar.model(m).condition ), h, sprintf('Clearing events for model %d [condition %d/%d]', m, c, length( ar.model(m).condition )));
                end
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
    if ( arOutputLevel > 1 )
        h = waitbar(0.9, h, 'Removing steady state conditions');
    end
    if isfield(ar.model, 'ss_condition')
        ar.model = rmfield(ar.model,'ss_condition');
    end
    
    % Clear the event log
    if ( isfield( ar, 'eventLog' ) )
        ar = rmfield(ar,'eventLog');
    end
    
    ar.ss_conditions = 0;
    if ( arOutputLevel > 1 )
        h = waitbar(0.95, h, 'Linking model');
    end
    % The event removal requires linking (silent link)
    arLink(true);
    
    if ( arOutputLevel > 1 )
        close(h);
    end
    
    % Invalidate cache so simulations do not get skipped
    arCheckCache(1);
end