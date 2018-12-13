% arToggle( parameter mask, value, [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])
%
%   mask         Parameter mask
%   value		 value of the  parameter
%   qFit         0=fixed, 1=fitted, 2=constant          [unchanged]
%   qLog10       0=normal, 1=log10 parameter values     [unchanged]
%   lb           lower parameter bound                  [unchanged]
%   ub           upper parameter bound                  [unchanged]
%   type         0=box prior, 1=normal prior            [unchanged]
%   meanp		 mean of normal prior                   [unchanged]
%   stdp         standard deviation of normal prior     [unchanged]
%
% Function to quicky and non-destructively try new parameters without losing 
% the old values and parameter fitting settings.
%
% Toggle temporarily stores the old parameters. When toggling parameters,
% new parameters are substituted, but the original parameters (prior to the
% first toggle) are stored. Toggling maintains the original parameter
% values in a global struct called arToggle.
%
% To restore the original values, simply call arToggle without setting
% parameter values. Once the original values are restored, the toggle
% backup is cleared. Manually clearing the toggle cache can be done by
% invoking arToggle( 'clear' ). This destroys the parameter backup and a
% next call to arToggle will store the current ones as the backup.
%
% Examples
%   arToggle( '_act_', 0, 2, 0 );
%       Sets all parameters with _act_ in the name to 0 and disables fitting
%   arToggle( '_act_', 1, 2, 0 );
%       Sets those same parameters to 1.
%   arToggle( '_act_' )
%       Returns the parameters to their original values prior to the first
%       toggle.
%
%   arToggle( 'all' ) or arToggle( )
%       Restores all parameter values to their pre-toggled values.
%
%   arToggle( '_act_', 0, 2, 0 );
%       Sets all parameters with _act_ in the name to 0 and disables fitting
%   arToggle( )
%       Restores parameters and clears toggle cache
%   arSetPars( '_act_', 0, 2, 0 );
%   arToggle( '_act_', 1, 2, 0 );
%   arToggle( )
%       Restores parameters to 0 (note the cache clear and setPars).
%
% Hint: arToggle has been superseded by arPush/arPop.
%
% See also arPush, arPop, arLoadPars, arSetPars and arSetParsPattern
%

function arToggle( ID, varargin )

    global ar;
    global arToggle;

    % First time it's used, store 0s.
    if (isempty(arToggle))
        arToggle.used       = zeros(size(ar.p));
        arToggle.checkstr   = ar.checkstr;
    else
        if ( ~strcmp( ar.checkstr, arToggle.checkstr ) )
            warning( 'Model different from last use of arToggle. Resetting cache...' );
            arToggle.used       = zeros(size(ar.p));
            arToggle.checkstr   = ar.checkstr;            
        end
    end
    
    if ( nargin < 1 )
        ID = 'all';
    end
    if ( ~isnumeric(ID) )
        if ( strcmp( ID, 'all' ) )
            ID = 1 : length(ar.p);
        elseif ( strcmp( ID, 'clear' ) )
            % Invalidate the complete cache and return
            arToggle.used = 0*arToggle.used;
            return;
        else
            ID = arFindPar(ID);
        end
    end
    
    % If there was an old one in the cache; restore it first
    restoredValues = 0;
    for a = 1 : length( ID )
        if ( arToggle.used(ID(a)) == 1 )
            restoredValues          = 1;
            ar.qLog10(ID(a))        = arToggle.qLog10(ID(a));
            ar.p(ID(a))             = arToggle.p(ID(a));
            ar.qFit(ID(a))          = arToggle.qFit(ID(a));
            ar.ub(ID(a))            = arToggle.ub(ID(a));
            ar.lb(ID(a))            = arToggle.lb(ID(a));
            ar.type(ID(a))          = arToggle.type(ID(a));
            ar.mean(ID(a))          = arToggle.mean(ID(a));
            ar.std(ID(a))           = arToggle.std(ID(a));
            arToggle.used(ID(a))    = 0;
        end
    end
    
    if ( length(varargin) > 0 )
        % Store old values
        for a = 1 : length( ID )
            arToggle.used(ID(a))    = 1;
            arToggle.qLog10(ID(a))  = ar.qLog10(ID(a));
            arToggle.p(ID(a))       = ar.p(ID(a));
            arToggle.qFit(ID(a))    = ar.qFit(ID(a));
            arToggle.ub(ID(a))      = ar.ub(ID(a));
            arToggle.lb(ID(a))      = ar.lb(ID(a));
            arToggle.type(ID(a))    = ar.type(ID(a));
            arToggle.mean(ID(a))    = ar.mean(ID(a));
            arToggle.std(ID(a))     = ar.std(ID(a));
        end

        % Set the new parameter values
        for a = 1 : length( ID )
            try
                arSetPars( ar.pLabel{ID(a)}, varargin{:} );
            catch ME
                warning('Failed to set %s: %s', ar.pLabel{ID(a)}, ME.message);
            end
        end
    else
        if ( restoredValues == 0 )
            warning( 'There were no values to restore from the toggle cache.' );
        end
    end
end
