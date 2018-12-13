% Plot model equilibration
%   This function plots the model pre-equilibration up to the point where
%   equilibration was successful. It is primarily meant for debugging
%   purposes.
%
%   function arPlotEquilibration( model, ss_condition )
%   
%       model           - Model index
%       ss_condition    - Steady state condition to plot
%

function arPlotEquilibration( model, ss_condition, N )

    global ar;
    tStart = ar.model(model).ss_condition(ss_condition).tstart;
    tEnd = ar.model(model).ss_condition(ss_condition).tEq;
        
    ar.model(model).ss_condition(ss_condition).tstop = tEnd;
    
    % Number of time points to show
    if ( nargin < 3 )
        N = 1000;
    end
    ts = tStart : ( tEnd - tStart ) / ( N - 1 ) : tEnd;
    ar.model(model).ss_condition(ss_condition).tFine = ts;
    
    oldss = arDeepCopy( ar.model(model).ss_condition(ss_condition) );
    
    % Prepare data
    ar.model(model).ss_condition(ss_condition).uFineSimu = zeros(N, size(ar.model(model).ss_condition(ss_condition).uFineSimu,2));
    ar.model(model).ss_condition(ss_condition).vFineSimu = zeros(N, size(ar.model(model).ss_condition(ss_condition).vFineSimu,2));
    ar.model(model).ss_condition(ss_condition).xFineSimu = zeros(N, size(ar.model(model).ss_condition(ss_condition).xFineSimu,2));
    ar.model(model).ss_condition(ss_condition).zFineSimu = zeros(N, size(ar.model(model).ss_condition(ss_condition).zFineSimu,2));
    
    np = numel( ar.model(model).ss_condition(ss_condition).pNum );
    nx = numel( ar.model(model).x );
    nv = numel( ar.model(model).v );
    nz = numel( ar.model(model).z );
    nu = numel( ar.model(model).u );
    
    ar.model(model).ss_condition(ss_condition).suFineSimu = zeros(N, nu, np);
    ar.model(model).ss_condition(ss_condition).svFineSimu = zeros(N, nv, np);
    ar.model(model).ss_condition(ss_condition).sxFineSimu = zeros(N, nx, np);
    ar.model(model).ss_condition(ss_condition).szFineSimu = zeros(N, nz, np);
       
    ar.stop = 0;
    ar.model(model).ss_condition(ss_condition).status = 0;
    ar.model(model).ss_condition(ss_condition).stop = 0;
    
    % propagate parameters
    arCheckCache(1);
    ar.stop = 0;
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).ss_condition)   
            ar.model(m).ss_condition(c).x0_override = []; % remove initial condition overrides which may have been used for rootfinding
            ar.model(m).ss_condition(c).status = 0;
            ar.model(m).ss_condition(c).pNum = ar.p(ar.model(m).ss_condition(c).pLink);
            ar.model(m).ss_condition(c).qLog10 = ar.qLog10(ar.model(m).ss_condition(c).pLink);
            ar.model(m).ss_condition(c).pNum(ar.model(m).ss_condition(c).qLog10 == 1) = ...
                10.^ar.model(m).ss_condition(c).pNum(ar.model(m).ss_condition(c).qLog10 == 1);
            ar.model(m).ss_condition(c).pNum
            ar.model(m).ss_condition(c).start = 0;
            ar.model(m).ss_condition(c).stop = 0;
            ar.model(m).ss_condition(c).stop_data = 0;
        end
    end

    % Simulate the system
    feval(ar.fkt, ar, true, false, 1, false, 'ss_condition', 'ss_threads', 0);
    
    styles = { 'k', 'r', 'b', 'm', 'k--', 'r--', 'b--', 'm--', 'k-.', 'r-.', 'b-.', 'm-.' };
    styles = { 'k', 'r', 'b', 'k--', 'r--', 'b--' };
    nStyles = numel(styles);
    nPlots = ( ceil( nx/nStyles ) + ceil( nz/nStyles ) + ceil( nv/nStyles ) + ceil( nu/nStyles ) );
    NX = ceil(sqrt(nPlots));
    NY = ceil(nPlots/NX);
    
    figure;
    cplot = plotGroup( model, ss_condition, 'uFineSimu', 'u', styles, NX, NY, 1 );
    cplot = plotGroup( model, ss_condition, 'xFineSimu', 'x', styles, NX, NY, cplot );
    cplot = plotGroup( model, ss_condition, 'zFineSimu', 'z', styles, NX, NY, cplot );
    cplot = plotGroup( model, ss_condition, 'vFineSimu', 'v', styles, NX, NY, cplot );
 
    ar.model(model).ss_condition(ss_condition) = oldss;
end

function cplot = plotGroup( model, ss_condition, field, namefield, styles, NX, NY, cplot )
    global ar;
    
    pidx    = 1;    % Number of lines in this plot
    sidx    = 1;    % Number of lines in this simulation field
    names   = ar.model(model).(namefield);
    sims    = ar.model(model).ss_condition(ss_condition).(field);
    t       = ar.model(model).ss_condition(ss_condition).tFine;
    nStyles = numel(styles);
    
    if ( size(sims,2) > 0 )
        for a = 1 : size( sims, 2 )
            if ( pidx == 1 )
                subplot(NX, NY, cplot); hold on;
                title( namefield );
                cplot = cplot + 1;
            end
            plot( t, sims(:, a), styles{pidx} )
            pidx    = pidx + 1;
            if ( pidx > nStyles )
                pidx = 1;
                legend( cellfun( @(x)strrep(x, '_', '\_' ), names(a-nStyles+1:a), 'UniformOutput', false ) );
            end
        end
        legend( cellfun( @(x)strrep(x, '_', '\_' ), names( size( sims, 2 )-pidx+2:size( sims, 2 ) ), 'UniformOutput', false ) );
    end
end