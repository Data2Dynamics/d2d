% Function to clear priors in the model
%
% Without arguments, the function clears non-bound priors off of all model parameters
%
%   Optional arguments:
%       'ssclear'       - Clears steady state priors
%       'condclear'     - Clears intercondition constraints
%       'all'           - Clears all forms of priors

function arClearPriors( varargin )

    global ar;
    
    switches = { 'all', 'steadystate', 'condclear' };
    extraArgs = [ 0, 0, 0 ];
    description = { ...
    {'', ''}, ...
    {'', ''}, ...
    {'', ''}, ...
    };

    [opts, ~] = argSwitch( switches, extraArgs, description, 0, varargin );

    ssclear = 0;
    parclear = 1;
    condclear = 0;
    if ( opts.steadystate )
        ssclear = 1;
        parclear = 0;
    end
    if ( opts.condclear )
        condclear = 1;
        parclear = 0;
    end
    if (opts.all)
        ssclear = 1;
        parclear = 1;
        condclear = 1;
    end
    
    cleared = {};
    if ( parclear )
        cleared{end+1} = 'parameter';
        for jp = 1 : numel( ar.pLabel )
            arSetPrior( ar.pLabel{jp}, 0, 0, 1 );
        end
    end
    if ( ssclear )
        cleared{end+1} = 'steady state';
        for jm = 1 : numel( ar.model )
            for jc = 1 : numel( ar.model(jm).condition )
                ar.model(jm).condition(jc).qSteadyState = 0 * ar.model(jm).condition(jc).qSteadyState;
            end
        end
    end
    if ( condclear )
        cleared{end+1} = 'intercondition';
        arClearConditionConstraints;
    end
    list = sprintf( '%s, ', cleared{:} );
    fprintf( 'Cleared %s constraints.\n', list(1:end-2) );
end