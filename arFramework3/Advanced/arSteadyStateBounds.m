% This function can be used to constrain steady state solutions.
%
% It applies a quadratic penalty when model states or derived variables leave a certain range.
%       (x-ub)^2        if x > ub
%       (x-lb)^2        if x < lb
%       0               otherwise
%
% Usage:
%   arSteadyStateBounds(m, c, xl, xu, zl, zu, logx, logz, weightx, weightz, showConstraints)
%
%   m               - Model index
%   c               - Condition index
%   xl              - Lower bound states
%   xu              - Upper bound states
%   zl              - Lower bound derived variables
%   zu              - Upper bound derived variables
%   logx            - Treat x in logarithmic space
%   logz            - Treat z in logarithmic space
%   weightx         - Weight(s) for the penalty on the states (inverse of SD)
%   weightz         - Weight(s) for the penalty on the derived variables
%   showConstraints - Show when the constraints are active

function arSteadyStateBounds(m, c, xl, xu, zl, zu, logx, logz, weightx, weightz, showConstraints)

    global ar;
    if ( ( numel( xl ) ~= numel( ar.model(m).x ) ) || ( numel( xu ) ~= numel( ar.model(m).x ) ) )
        error( 'Bound vectors need to be the number of states in length!' );
    end
    if ( ( numel( zl ) ~= numel( ar.model(m).z ) ) || ( numel( zu ) ~= numel( ar.model(m).z ) ) )
        error( 'Bound vectors need to be the number of derived variables in length!' );
    end
    if ( sum( xl > xu ) )
        sta = sprintf( '%s ', ar.model.x{xl>xu} );
        error( 'Inconsistent bounds for state variables: %s', sta );
    end    
    if ( sum( zl > zu ) )
        der = sprintf( '%s ', ar.model.z{zl>zu} );
        error( 'Inconsistent bounds for derived variables: %s', der );
    end
    if ~exist( 'logx', 'var' )
        logx = ones( size( xl ) );
    end
    if ~exist( 'logy', 'var' )
        logz = ones( size( zl ) );
    end
    if ~exist( 'weightx', 'var' )
        weightx = 10 * ones( size(xl) );
    end
    if ~exist( 'weightz', 'var' )
        weightz = 10 * ones( size(zl) );
    end
    if ~exist( 'showConstraints', 'var' )
        showConstraints = 0;
    end
    
    serializeVector = @(x)sprintf('%g ', x);
    fprintf( 'arSteadyStateBounds(%d, %d, [%s], [%s], [%s], [%s], [%s], [%s], [%s], [%s], %d)\n', m, c, serializeVector(xl), serializeVector(xu), serializeVector(zl), serializeVector(zu), serializeVector(logx), serializeVector(logz), serializeVector(weightx), serializeVector(weightz), showConstraints );   
    
    % Enforce logical
    logx = logx == 1;
    logz = logz == 1;
    
    if ( xl(logx) <= 0 )
        error( 'Cannot specify an lower bound of zero for x when operating in log mode' );
    end
    if ( zl(logz) <= 0 )
        error( 'Cannot specify an lower bound of zero for z when operating in log mode' );
    end    

    % Transform the bounds which have to be log-trafo'd
    xl(logx) = log10(xl(logx));
    xu(logx) = log10(xu(logx));
    zl(logz) = log10(zl(logz));
    zu(logz) = log10(zu(logz));
    
    res_fun = @()residual_concentrationConstraintsL2(m, c, xl, xu, zl, zu, weightx, weightz, logx, logz);
    arAddCustomResidual( 'SteadyStateBounds', res_fun, 1 );
    ar.config.show_ss_constraints = showConstraints;
end

% This function places a soft bound on concentrations
function [res_user, sres_user, res_type] = residual_concentrationConstraintsL2(m, c, xl, xu, zl, zu, wx, wz, logx, logz)

    global ar
    
    x_active    = ( ~( isnan( xl ) | isnan( xu ) ) );
    z_active    = ( ~( isnan( zl ) | isnan( zu ) ) );
    nx          = sum(x_active);
    nz          = sum(z_active);
    logx        = logx(x_active);
    logz        = logz(z_active);
    
    np          = size( ar.model(m).ss_condition(c).sxFineSimu, 3 );
    pLink       = ar.model(m).ss_condition(c).pLink;
    xss         = ar.model(m).ss_condition(c).xFineSimu(end,x_active);
    zss         = ar.model(m).ss_condition(c).zFineSimu(end,z_active);
    sxss        = squeeze(ar.model(m).ss_condition(c).sxFineSimu(end,x_active,:));
    szss        = squeeze(ar.model(m).ss_condition(c).szFineSimu(end,z_active,:));
    log10s      = ar.qFit(ar.model(m).ss_condition(c).pLink)==1;
    
    % Transform the state and derived variable sensitivities in case they
    % are specified in log10 parameters
    sxss(:,log10s) = sxss(:,log10s) .* repmat( ar.model(m).ss_condition(c).pNum(log10s) * log(10), nx, 1 );
    if ( ~isempty( zl ) )
        szss(:,log10s) = szss(:,log10s) .* repmat( ar.model(m).ss_condition(c).pNum(log10s) * log(10), nz, 1 );
    end
    
	% Transform the sentivities of the ones that have to be penalized in
	% log (note that this has to be done before the states are transformed,
	% since we need the untransformed states for this).
    %   dlog10(y(p))/dp = (1/(y(p)*log(10))) dy(p)/dp
    sxss(logx, :) = sxss(logx, :) .* repmat( 1 ./ (xss(logx) * log(10)), np, 1 ).';
    if ( ~isempty( zl ) )
        szss(logz, :) = szss(logz, :) .* repmat( 1 ./ (zss(logz) * log(10)), np, 1 ).';    
    end
    
    % Transform the states if they are to be penalized in log10
    xss(logx) = log10(xss(logx));
    if ( ~isempty( zl ) )
        zss(logz) = log10(zss(logz));
    end
    
    % Determine which bounds are active
    xl          = xl(x_active);
    xu          = xu(x_active);
    zl          = zl(z_active);
    zu          = zu(z_active);
    wx          = wx(x_active);
    wz          = wz(z_active);
    
    % Compute the residual
    xlower      = (xss < xl) .* (xl - xss) .* wx;
    xupper      = (xss > xu) .* (xss - xu) .* wx;
    if ( ~isempty( zl ) )
        zlower      = (zss < zl) .* (zl - zss) .* wz;
        zupper      = (zss > zu) .* (zss - zu) .* wz;
    end

    if ( ~isempty( zl ) )
        res_user    = [ xlower, xupper, zlower, zupper ];
    else
        res_user    = [ xlower, xupper ];
    end
    
    tot = sum( res_user );
    if ( isnan( tot ) || isinf( tot ) )
        xs = isinf(xlower) | isnan(xupper);
        zs = isinf(zlower) | isnan(zupper);
        ix = find(x_active); iz = find(z_active);
        xstr = sprintf( '%s ', ar.model.x{ix(xs)} );
        zstr = sprintf( '%s ', ar.model.z{iz(zs)} );
        fprintf( 'Infinity found in residuals: %s, %s', xstr, zstr );
    end
    
    if ( ar.config.show_ss_constraints )
        xs = (xss < xl) | (xss > xu);
        zs = (zss < zl) | (zss > zu);
        ix = find(x_active); iz = find(z_active);
        xstr = sprintf( '%s ', ar.model.x{ix(xs)} );
        zstr = sprintf( '%s ', ar.model.z{iz(zs)} );
        fprintf( 'Active bounds: %s, %s\n', xstr, zstr );
    end
    
    % Compute the sensitivities
    sxlower     = - repmat( (xss < xl) .* wx, np, 1 ).' .* sxss;
    sxupper     =   repmat( (xss > xu) .* wx, np, 1 ).' .* sxss;
    if ( ~isempty( zl ) )
        szlower     = - repmat( (zss < zl) .* wz, np, 1 ).' .* szss;
        szupper     =   repmat( (zss > zu) .* wz, np, 1 ).' .* szss;    
    end
    
    % Sres with respect to all inner parameters
    if ( ~isempty( zl ) )
        assembled   = [ sxlower; sxupper; szlower; szupper ];
    else
        assembled   = [ sxlower; sxupper ];
    end
    
    res_type    = ones(size(res_user)); % Treat the constraint like data
    
    sres_user   = zeros( size(assembled, 1), numel(ar.p) );
       
    sres_user( :, pLink ) = assembled;
end
