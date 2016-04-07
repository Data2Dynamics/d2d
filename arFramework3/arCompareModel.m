% Simple helper function to determine for which datasets two models predict 
% differences
%
% Usage:
%   arCompareModel(ar1, m1, ar2, m2)
%       where ar1 and ar2 refer to ar structs and m1 and m2 to model indices.
%   arCompareModel(m1, m2)
%       where the default ar struct is used, and two models within are
%       compared
%
% Note: If interactivity mode is on, you can click the matrix to see the 
% curves corresponding to the different models. Active interactivity mode 
% by invoking "arInteractivity on".

function arCompareModel(varargin)
    
    figure('units','normalized','outerposition',[0 0 1 1], 'Name', 'Model comparison. Green means M1 is better.');

    % Special switches
    switches    = { 'relerr' };
    extraArgs   = [ 0 ];
    descriptions = {    { 'Using relative residual', '' }, ...
                    };

    if ( length( varargin ) < 2 )
        error('Insufficient arguments for this mode. See help arCompareModel.');
    end
    
    % An ar struct was supplied
    if ( isstruct( varargin{1} ) )
        if ( length( varargin ) < 4 )
            error( 'Insufficient arguments for this mode. See help arCompareModel.' );
        end
        if ( ~isstruct(varargin{3}) )
            error( 'Incorrect argument #3. Should be an ar struct' );
        end
        if ( ~(isnumeric(varargin{2})&&(numel(varargin{2})==1)&&isnumeric(varargin{4})&&(numel(varargin{4})==1)) )
            error( 'Second and fourth argument should be model IDs (numeric)' );
        end
        [ar1, m1, ar2, m2] = deal( varargin{1:4} );
        flags = varargin{5:end};
    else
        global ar; %#ok
        ar1 = ar;
        ar2 = ar;
        [m1, m2] = deal( varargin{1:2} );
        flags = varargin{3:end};
    end
    if ( m1 > length( ar1.model ) )
        error( 'Model ID 1 was specified out of range' );
    end
    if ( m2 > length( ar2.model ) )
        error( 'Model ID 2 was specified out of range' );
    end

    % Process additional flags
    opts = argSwitch( switches, extraArgs, descriptions, flags );
    
    name1 = ar1.model(m1).name;
    name2 = ar2.model(m2).name;

    %states      = intersect( {ar1.model(m1).x{:}, ar1.model(m1).z{:}}, {ar2.model(m2).x{:}, ar2.model(m2).z{:}} );
    observables = intersect( getObsNames(ar1.model(m1)), getObsNames(ar2.model(m2)) );
       
    dataFiles   = intersect( getDataNames(ar1.model(m1)), getDataNames(ar2.model(m2)) );

    changeMatrix = zeros( numel(dataFiles), numel(observables) );
    totalErr     = zeros( numel(dataFiles), numel(observables) );
    for a = 1 : length( dataFiles )
        I1 = find( ismember({ar1.model(m1).data.name}, dataFiles{a}) );
        I2 = find( ismember({ar2.model(m2).data.name}, dataFiles{a}) );

        % Verify that we have the same experimental data
        for c = 1 : length( I1 )
            if ( nanmax( nanmax(ar1.model(m1).data(I1(c)).yExp - ar1.model(m1).data(I2(c)).yExp) ) > 0 )
                warning( 'Problem finding corresponding dataset for %s\n', dataFiles{a} );
            end
            
            % List by observables
            for b = 1 : length( ar1.model(m1).data(I1(c)).fy )
                [~, obsIndex] = intersect(observables, ar1.model(m1).data(I1(c)).y{b});
                
                if( ar1.model(m1).data(I1(c)).qFit(b) == 1 )
                    if ( ~opts.relerr )
                        if( (ar1.config.useFitErrorMatrix == 0 && ar1.config.fiterrors == 1) || ...
                                (ar1.config.useFitErrorMatrix==1 && ar1.config.fiterrors_matrix(m1,I1(c))==1) )
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2err(b);
                        end
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2(b);
                    else
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + sum((ar1.model(m1).data(I1(c)).yExpSimu(:,b)-ar1.model(m1).data(I1(c)).yExp(:,b)).^2);
                    end
                end
            end
            
            for b = 1 : length( ar2.model(m2).data(I2(c)).fy )
                [~, obsIndex] = intersect(observables, ar2.model(m2).data(I2(c)).y{b});
                
                if( ar2.model(m2).data(I2(c)).qFit(b) == 1 ) 
                    if ( ~opts.relerr )
                        if( (ar2.config.useFitErrorMatrix == 0 && ar2.config.fiterrors == 1) || ...
                                (ar2.config.useFitErrorMatrix==1 && ar2.config.fiterrors_matrix(m2,I2(c))==1) )
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2err(b);
                        end
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2(b);
                    else
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - sum((ar2.model(m2).data(I2(c)).yExpSimu(:,b)-ar2.model(m2).data(I2(c)).yExp(:,b)).^2);
                        totalErr( a, obsIndex )     = totalErr( a, obsIndex ) + sum((ar2.model(m2).data(I2(c)).yExpSimu(:,b)-ar2.model(m2).data(I2(c)).yExp(:,b)).^2);
                    end
                end
            end            
        end
    end
    
    if ( opts.relerr )
        changeMatrix = changeMatrix ./ totalErr;
        changeMatrix(isinf(changeMatrix)|isnan(changeMatrix)) = 0;
    end
    
    imagesc( changeMatrix ); hold on;
    for a = 1 : length( dataFiles )
        for b = 1 : length( observables )
            if ( abs( changeMatrix( a, b ) ) > 1e-2 )
                text( b, a, sprintf( '%.1f', changeMatrix( a, b ) ), 'Color', 'white', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
            end
        end
    end
    change = max(max(abs(changeMatrix)));
    
    % If the interactivity system is enabled, register the callbacks
    % and provide arInteractivity with the required data.
    if ( arInteractivity )
        arInteractivity( 'arCompareModel', dataFiles, observables, ar1, m1, ar2, m2, getPlotIDs(ar1.model(m1), dataFiles), getPlotIDs(ar2.model(m2), dataFiles) );
    end
    
    set(gca, 'YTick', 1 : numel(dataFiles) );
    set(gca, 'XTick', 1 : numel(observables) );
    set(gca, 'YTickLabel', dataFiles );
    set(gca, 'XTickLabel', observables );
    set(gca, 'CLim', [-change, change] );
    colormap( redgreencmap(256, 'Interpolation', 'sigmoid') );
    title( sprintf( '%s - %s', strTrafo(name1), strTrafo(name2) ) );
    colorbar;
end

% Find the corresponding data plots for both models being compared
function plotIDs = getPlotIDs(model, dataFiles)
    plotIDs = cell(1, length(dataFiles));
    for a = 1 : length( dataFiles )
        plotIDs{a} = [];
        dataIDs = find( strcmp( {model.data.name}, dataFiles{a} ) );
        for b = 1 : length( model.plot )
            % Is this plot involved with this dataset? => Add it to the
            % list
            if ( max( ismember( model.plot(b).dLink, dataIDs ) ) )
                plotIDs{a} = [ plotIDs{a}, b ];
            end
        end
    end
end

function str = strTrafo( str )
    str = strrep( str, '_', '\_');
end

function names = getObsNames(model)
    names = {};
    for a = 1 : length( model.data )
        names = union( names, model.data(a).y );
    end
end
    
function names = getDataNames(model)
    names = {};
    for a = 1 : length( model.data )
        if ( max( model.data(a).qFit ) == 1 )
            names = union( names, model.data(a).name );
        end
    end
end

function [opts] = argSwitch( switches, extraArgs, description, varargin )

    for a = 1 : length(switches)
        if ( extraArgs(a) == 0 )
            opts.(lower(switches{a})) = 0;
        else
            opts.(lower(switches{a})) = {};
        end
    end
    
    a = 1;
    if ~isempty( varargin )
        while( a <= length( varargin ) )
            if ( max( strcmpi( switches, varargin{a} ) ) == 0 )
                error( 'Invalid switch argument was provided %s', varargin{a} );
            end
            if ( extraArgs( strcmpi( switches, varargin{a} ) ) == 0 )
                opts.(lower(varargin{a})) = 1;
            else
                try
                    opts.(lower(varargin{a})) = varargin{a+1};
                    a = a + 1;
                catch
                    error( 'Did not provide arguments for flag %s', varargin{a} );
                end
            end
            a = a + 1;
        end
    end
    for a = 1 : length( switches )
        if ( extraArgs(a) == 0 )
            fprintf( '%s\n', description{a}{2-opts.(lower(switches{a}))} );
        else
            fprintf( '%s\n', description{a}{1+isempty(opts.(lower(switches{a})))} );
        end
    end
end