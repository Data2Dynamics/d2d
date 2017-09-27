% Advanced function, use with caution
%
% arRestrictSimulation restricts simulation to those conditions which are
% actually being used by fitted data. Other conditions are *not* simulated.
%
% Undo the setting with arRestrictSimulations(1)

function arRestrictSimulations( varargin )
    global ar;
    
    enable = 0;
    if (nargin > 0)
        enable = varargin{1};
    end
    
    if ( length( varargin ) < 2 )
        % Make a list of active conditions
        for m = 1 : length( ar.model )
            conditions(m).ac = [];
            for d = 1 : length( ar.model(m).data )
                fit = max(ar.model(m).data(d).qFit);
                if ( fit | enable )
                    conditions(m).ac = union( conditions(m).ac, ar.model(m).data(d).cLink );
                end
            end
            disp( sprintf( 'Model %d - Active conditions %d/%d', m, length(conditions(m).ac), length(ar.model(m).condition) ) );
        end
    else
        for m = 1 : length( ar.model )
            conditions(m).ac = [];
        end
        if ( varargin{1} > length( ar.model ) )
            error( 'Incorrect model ID' );
        else
            m = varargin{1};
        end
        if ( min( varargin{2} ) < 1 )
            error( 'Incorrect condition ID (smaller than 1)' );
        end
        if ( max( varargin{2} ) > length( ar.model(m).condition ) )
            error( 'Incorrect condition ID (exceeds limit)' );
        end
        
        conditions(varargin{1}).ac = varargin{2};
    end
    
    populateThreads( 'threads', 'condition', conditions );
    
    % When only non-dynamic parameters are fitted, some models do not
    % simulate the system anymore. Make sure we simulate it here to make
    % sure that the simulation doesn't contain old simulation results
    if ( enable )
        arSimu(false); arCalcMerit(false);
    end
end

function populateThreads( thread_fieldname, condition_fieldname, conditions )
    global ar;
    
    % populate threads
    ar.config.(thread_fieldname) = [];
    ar.config.(thread_fieldname)(1).id = 0;
    ar.config.(thread_fieldname)(1).n = 0;
    ar.config.(thread_fieldname)(1).nd = 0;
    ar.config.(thread_fieldname)(1).ms = int32([]);
    ar.config.(thread_fieldname)(1).cs = int32([]);
    ithread = 1;
    ar.config.nTasks = 0;
    for m = 1:length(ar.model)
        if (isfield(ar.model(m), condition_fieldname))
            for ac = 1 : length( conditions(m).ac )
                if(length(ar.config.(thread_fieldname))<ithread)
                    ar.config.(thread_fieldname)(ithread).id = ithread-1;
                    ar.config.(thread_fieldname)(ithread).n = 0;
                    ar.config.(thread_fieldname)(ithread).nd = 0;
                    ar.config.(thread_fieldname)(ithread).ms = int32([]);
                    ar.config.(thread_fieldname)(ithread).cs = int32([]);
                end
                c = conditions(m).ac(ac);
                ar.config.nTasks = ar.config.nTasks + 1;
                ar.config.(thread_fieldname)(ithread).n = ...
                    ar.config.(thread_fieldname)(ithread).n + 1;
                ar.config.(thread_fieldname)(ithread).nd = ...
                    ar.config.(thread_fieldname)(ithread).nd + ...
                    length(ar.model(m).(condition_fieldname)(c).dLink);        
                ar.config.(thread_fieldname)(ithread).ms(end+1) = int32(m-1);
                ar.config.(thread_fieldname)(ithread).cs(end+1) = int32(c-1);
                ithread = ithread + 1;
                if(ithread>ar.config.nParallel)
                    ithread = 1;
                end
            end
        end
    end
end