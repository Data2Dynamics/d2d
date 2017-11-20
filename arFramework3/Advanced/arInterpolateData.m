% Function to interpolate data
%
% This function aids in finding a good fit in the early stages of model
% development. It linearly interpolates between points data to force the model 
% to more strongly weight a specific dataset. This data should be removed
% again before any rigorous statistical analysis
%
% You can remove the data again with arClearInterpolatedData.
%
% Usage:
%
%   function arInterpolateData( m, ds, obs, 'tmin', tmin, 'tmax', tmax, 'steps', steps, .... )
%
%    m        - Model index
%    ds       - Data indices, find using arFindData
%    obs      - Observable
%   
%   Optional arguments :
%    tmin       - Minimal time to interpolate (default = -inf)
%    tmax       - Maximal time to interpolate (default = inf)
%    steps      - Number of steps to use (default = 50)
%    plot       - Show interpolation
%    plotOnly   - Only plot the interpolant. Do not add the data
%    modfactor  - Multiply estimated noise (single component noise model)
%                 with a factor (default = 1)   

function arInterpolateData( m, ds, obs, varargin )  

    global ar;

    switches = { 'tmin', 'tmax', 'steps', 'plot', 'plotonly', 'modfactor' };
    extraArgs = [ 1, 1, 1, 0, 0, 1 ];
    description = { ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', 'Only plotting: not adding the data'} ...
    {'', ''} ...
    };
    [opts, ~] = argSwitch( switches, extraArgs, description, 0, 'softmatching', varargin );    
    
    if nargin < 3
        help arInterpolateData;
        error( 'Function needs at least 3 arguments.' );
    end
    
    if ~opts.modfactor
        modfactor = 1;
    else
        modfactor = opts.modfactor_args;
    end
    if ~opts.tmin
        tmin = -inf;
    else
        tmin = opts.tmin_args;
    end
    if ~opts.tmax
        tmax= inf;
    else
        tmax = opts.tmax_args;
    end
    if ~opts.steps
        steps = 50;
    else
        steps = opts.steps_args;
    end
    if ( tmin > tmax )
        error( 'Tmin > Tmax' );
    end
    
    if ~iscell( obs )
        obs = {obs};
    end
    
    figure;
    NX = ceil( sqrt( numel(ds ) ) );
    NY = ceil( numel(ds) / NX );
    for jo = 1 : numel( obs )
        for jd = 1 : numel(ds)
            obsID = find( ismember( ar.model(m).data(ds(jd)).y, obs{jo} ) );
            if isempty( obsID )
                obsies = sprintf( '%s ', ar.model(m).data(ds(jd)).y{:} );
                warning( 'Cannot find observable %s in data with ID %d. Ignoring this one.\nAvailable observables are:\n%s\n', obs{jo}, ds(jd), obsies )
            else
                [t, y, tReq, yN] = generateData( m, ds(jd), obsID, steps, tmin, tmax, modfactor );

                if ( (opts.plot == 1) || (opts.plotonly == 1) )
                    subplot(NX, NY, jd);
                    plot( t, y, '.' ); hold on;
                    plot( tReq, yN, 'r.' );
                    title( strrep( ar.model(m).data(ds(jd)).name, '_', '\_' ) );
                    if ( jd == numel(ds) )
                        legend( { 'Original data', 'Interpolated data' } );
                    end
                    ylabel( strrep( obs{jo}, '_', '\_' ) );
                end
                if ( ~opts.plotonly )
                    addToData( m, ds(jd), obsID, tReq, yN )
                end
            end
        end
    end
    
    % Link the model
    if ( ~opts.plotonly )
        arLink;
    end
end

% Function which adds the data points to the arStruct
function addToData( m, d, obs, tExpNew, yExpNew )
    global ar;
    
    nExps = numel( yExpNew );
    nans = nan( nExps, numel( ar.model(m).data(d).y ) );
    nData = nans;
    nData(:, obs) = yExpNew;

    if ~isfield( ar.model(m).data(d), 'interpolatedData' )
        ar.model(m).data(d).interpolatedData = [];
    end
    ar.model(m).data(d).interpolatedData( size( ar.model(m).data(d).yExp, 1 ) + 1 : size( ar.model(m).data(d).yExp, 1 ) + size( nData, 1 ), : ) = 1;
    ar.model(m).data(d).yExp = [ ar.model(m).data(d).yExp; nData ];
    ar.model(m).data(d).tExp = [ ar.model(m).data(d).tExp; tExpNew.' ];
    ar.model(m).data(d).yExpStd = [ ar.model(m).data(d).yExpStd; nans ];
end

% Function which generates the data points
function [t, y, tReq, yN] = generateData( m, d, obs, steps, tmin, tmax, modfactor )
    global ar; 
    
    t = ar.model(m).data(d).tExp;
    y = ar.model(m).data(d).yExp(:, obs);
    [tUnique, ~] = unique(t);
    
    % Determine rough estimate of the noise on the data
    yRes = [];
    for a = 1 : numel( tUnique )
        yC          = y(t == tUnique(a));
        yMean(a)    = nanmean(yC); %#ok
        if ( ~isnan( yMean(a) ) )
            yRes        = [yRes; yC - yMean(a)]; %#ok
        end
    end
    nans = find( isnan( yMean ) );
    yMean(nans) = [];
    tUnique(nans) = [];  
    stdest = nanstd( yRes ) * modfactor;
    
    % Restrict to desired time range
    tR = t( (t > tmin) & (t < tmax) );  
    tReq = min(tR) : (max(tR)-min(tR))/(steps-1) : max(tR);
    if ( numel( tUnique ) == 1 )
        error( 'Only one datapoint in interpolation range.');
    end
    
    % Make sure the tReq's don't fall on points where there is actually
    % data. Shift those points in between the existing ones.
    tUq = tUnique + 0;
    for a = 1 : numel( tReq )
        idx = find(ismember(tUq, tReq(a)));
        if ~isempty( idx )
            if ( idx == 1 )
                tReq(a) = mean( tUq(idx:idx+1) );
            else
                tReq(a) = mean( tUq(idx-1:idx) );
            end
            tUq = union( tUq, tReq(a) );
        end
    end
    tReq = setdiff( tReq, tUnique );
    
    yN = zeros( numel( tReq ), 1 );
    
    for a = 1 : numel( tReq )
        ID_end = find( tUnique > tReq(a), 1 );
        y2 = yMean(ID_end);
        y1 = yMean(ID_end-1);
        t2 = tUnique(ID_end);
        t1 = tUnique(ID_end-1);
        yN(a) = ( (y2 - y1) / (t2 - t1) ) * ( tReq(a) - t1 ) + y1;
    end

    yN = yN + stdest*randn(size(yN));
end
