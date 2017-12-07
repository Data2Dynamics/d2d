% function arAddCustomResidual( name, fn, overwrite )
% Arguments:
%    name           - Name of the residual
%    fn             - Function handle of the residual
%    overwrite      - Should an existing function handle with the same name
%                     be overwritten?
%
% Adds a custom residual. Custom residuals must be specified with and name, function handle.
%

function arAddCustomResidual( name, fn, overwrite )
    global ar;
    
    % If it doesn't exist, create the structure
    if ( ~isfield( ar.config, 'user_residual_fun' ) || isempty( ar.config.user_residual_fun ) )
        ar.config.user_residual_fun         = struct;
        ar.config.user_residual_fun.fn      = {};
        ar.config.user_residual_fun.qFit    = [];
        ar.config.user_residual_fun.name    = {};
    end
    
    % If we have an older version of the custom residual function, convert
    % it to the newer version
    if ~isstruct( ar.config.user_residual_fun )
        if isa(ar.config.user_residual_fun, 'function_handle')
            warning( 'Old style custom residual function. Converting to new style.' );
            fnHandle = ar.config.user_residual_fun;
            ar.config = rmfield( ar.config, 'user_residual_fun' );
            arAddCustomResidual( 'Legacy', fnHandle, 0 );
        end
    end
    
    if ( nargin < 2 )
        warning( 'No function was specified for addition' );
        return;
    end
    
    % Typically we only want one custom residual of one type
    if ( nargin < 3 )
        overwrite = 1;
    end
    
    idx = find( ismember( ar.config.user_residual_fun.name, name ) );
    if ( ~isempty(idx) )
        if ( overwrite )
            warning( 'Overwriting old residual with name %s', name );
            if ( numel(idx) > 1 )
                warning( 'Multiple residuals with name %s. Overwriting first one.', name )
                idx = idx(1);
            end
        else
            warning( 'Custom residual with name %s already exists. If you wish to overwrite the residual with that name, please set overwrite flag to 1', name );
            idx = numel(ar.config.user_residual_fun.fn) + 1;
        end
    else
        idx = numel(ar.config.user_residual_fun.fn) + 1;
    end
    
    ar.config.user_residual_fun.fn{idx}     = fn;
    ar.config.user_residual_fun.qFit(idx)   = 1;
    ar.config.user_residual_fun.name{idx}   = name;
    
    % The objective function changed, so invalidate the cache!
    arCheckCache(1);
end