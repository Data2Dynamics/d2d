% arFractionalFlux( species, [varargin] )
%
% Plot contributions of the different fluxes to a state.
%
%   species     Refers to a model state (e.g. 'ATP').
%   tmax        Specify 'tmax' followed by value to set maximum time to 
%               plot.
%   condition   Specify 'condition' followed by index to specify which 
%               condition to plot.
%   model       Specify 'model' followed by index to specify which model to
%               plot.
%
function arFractionalFlux( species, varargin )

    global ar;

    cond = 1;
    model = 1;

    args = { 'tmax', 'condition', 'model' };
    extraArgs = [1, 1, 1];
    opts = argSwitch( args, extraArgs, {}, 0, varargin );

    if ~exist( 'species', 'var' )
        error( 'Please specify a model state of interest to study' );
    end
    if opts.model
        model = opts.model_args;
        if ( ~isnumeric( model ) || ( model > numel( ar.model ) ) )
            error( 'Invalid argument passed for model' );
        end
    end
    if opts.condition
        cond = opts.condition_args;
        if ( ~isnumeric( cond ) || ( cond > numel( ar.model(model).condition ) ) )
            error( 'Invalid argument passed for condition' );
        end    
    end
    if opts.tmax
        tLimOld = ar.model(model).tLim(2);
        ar.model(model).tLim(2) = opts.tmax_args;
        arLink;    
    end
    
    % Simulate the system
    arSimu(false,true,true);
    
    % Fetch required data
    state               = ismember( ar.model(model).x, species );
    if sum( state ) == 0
        str = sprintf('\n%s', ar.model(model).x{:} );
        error( 'State does not exist. Pick from %s\n', str );
    end
    stoich              = ar.model(model).N( state, : );
    t                   = ar.model(model).condition(cond).tFine;
    v                   = ar.model(model).condition(cond).vFineSimu;
    v                   = v.*repmat(stoich, size(v,1), 1);
    involvedReactions   = find( abs(stoich)>0 );
    v                   = v(:,involvedReactions);
    maxes               = max( abs( v ) );
    vNames              = ar.model(model).v( involvedReactions );
    [~,I]               = sort(maxes, 'descend');

    variables = {};
    for jR = 1 : numel( involvedReactions )
        variables = union( variables, symvar( ar.model(model).fv{involvedReactions(jR)} ) );
        [~, involvedStates] = intersect( ar.model(model).x, variables );
    end

    ax(1) = subplot(2,2,1);
    for a = 1 : length( I )
        plot( t, v( :, I(a) ), qColor(a) );
        hold on;
    end
    plot( t, sum(v, 2), 'k--' );
    legend( union( sanitize( vNames(I)), 'Net flux', 'stable' ) );
    
    %u = max( median(ar.model(model).condition.vFineSimu)+.4*mad(ar.model(model).condition.vFineSimu) );
    %l = min( median(ar.model(model).condition.vFineSimu)-.4*mad(ar.model(model).condition.vFineSimu) );
    %ylim([l, u]);
    ylabel( sprintf( 'Contribution to %s flux', ar.model(model).x{state} ) );
    xlabel( 'Time' );

    ax(2) = subplot(2,2,2);
    total = max(v,[],2)-min(v,[],2);
    v = v ./ repmat(total, 1, size(v,2));
    for a = 1 : length( I )
        plot( t, v( :, I(a) ), qColor(a) );
        hold on;
    end
    legend( sanitize( vNames(I) ) );
    ylabel( sprintf( 'Fractional %s flux', ar.model(model).x{state} ) );
    xlabel( 'Time' );

    % Plot all concentrations involved in the rxns
    ax(3) = subplot(2,2,3);
    for jS = 1 : numel( involvedStates )
        plot( t, ar.model(model).condition(cond).xFineSimu(:, involvedStates(jS)), qColor(jS) );
        hold on;
    end
    legend( sanitize( ar.model(model).x(involvedStates) ) );
    ylabel( 'Concentration' );
    xlabel( 'Time' );
    title( 'Concentration' );

    ax(4) = subplot(2,2,4);
    plot( t, ar.model(model).condition(cond).xFineSimu(:, state) );
    ylabel( sprintf( 'Concentration %s', ar.model(model).x{state} ) );
    xlabel( 'Time' );
    linkaxes(ax, 'x');
    xlim([0, ar.model(model).tLim(2)]);

    % Link back with old time list
    if ( exist( 'tLimOld', 'var' ) )
        ar.model(model).tLim(2) = tLimOld;
        arLink;
    end
end

function col = qColor( a )
    colors = 'krbmgcykrbmgcy';
    if ( a > 6 )
        app = '--';
    else
        app = '';
    end
    col = [colors(a), app];
end

function str = sanitize( str )
    str = strrep( str, '_', '\_' );
end