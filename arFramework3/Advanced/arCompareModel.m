% Simple helper function to determine for which datasets two models predict 
% differences
%
% Usage:
%   arCompareModel(ar1, m1, ar2, m2, (flags))
%       where ar1 and ar2 refer to ar structs and m1 and m2 to model indices.
%   arCompareModel(m1, m2, (flags))
%       where the default ar struct is used, and two models within are
%       compared
%
% Optional flags:
%   'relerr'    - Use relative errors instead of chi2 change
%   'absrsq'    - Show R^2 per dataset for Model 1 instead of differences
%   'nogrid'    - Don't draw the grid
%
% Note: If interactivity mode is on, you can click the matrix to see the 
% curves corresponding to the different models. Activate interactivity mode 
% by invoking "arInteractivity on".

function arCompareModel(varargin)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true); 
    if length(filenames)>2
       error('Error: Comparison of more than two models is not supported.') 
    end
    for j=1:length(filenames)
        fname = ['./Results/' filenames{j} '/workspace.mat'];
        if(exist(fname,'file'))
            S=load(fname);
            varargin{2*(j-1)+1} = S.ar;
            varargin{2*(j-1)+2} = 1; % only working for model 1
        else
            error('Error: No workspace found in %s',fname) 
        end    
    end
end

    figure('units','normalized','outerposition',[0 0 1 1], 'Name', 'Model comparison. Green means M1 is better.');

    % Special switches
    switches    = { 'relerr', 'absrsq', 'nogrid', 'titles' };
    extraArgs   = [ 0, 0, 0, 1 ];
    descriptions = {    { 'Using relative residual', '' }, ...
                        { 'Showing absolute R squared of M1 only', '' }, ...
                        { 'Suppressing grid', '' }, ...
                        { 'Using custom titles', '' }, ...
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
        varargin(1:4) = []; %#ok
    else
        global ar; %#ok
        ar1 = ar;
        ar2 = ar;
        [m1, m2] = deal( varargin{1:2} );
        varargin(1:2) = [];
    end
    if ( m1 > length( ar1.model ) )
        error( 'Model ID 1 was specified out of range' );
    end
    if ( m2 > length( ar2.model ) )
        error( 'Model ID 2 was specified out of range' );
    end
    
    % Validate ar structs
    ar1 = arInitFields(ar1);
    ar2 = arInitFields(ar2);

    % Process additional flags
    opts = argSwitch( switches, extraArgs, descriptions, 1, varargin );
    
    if ( opts.relerr && opts.absrsq )
        error( 'Incompatible option flags' );
    end
    
    name1 = ar1.model(m1).name;
    name2 = ar2.model(m2).name;

    %states      = intersect( {ar1.model(m1).x{:}, ar1.model(m1).z{:}}, {ar2.model(m2).x{:}, ar2.model(m2).z{:}} );
    observables = intersect( getObsNames(ar1.model(m1)), getObsNames(ar2.model(m2)) );
       
    dataFiles   = intersect( getDataNames(ar1.model(m1)), getDataNames(ar2.model(m2)) );

    changeMatrix = zeros( numel(dataFiles), numel(observables) );
    totalErr     = zeros( numel(dataFiles), numel(observables) );
    if ( opts.absrsq )
        allData = cell( numel(dataFiles), numel(observables) );
    end
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
                    if ( ~opts.absrsq )
                        if ( ~opts.relerr )
                            % Use change in chi2 as measure
                            if( ar1.config.fiterrors==1 || sum(ar1.qFit(ar1.qError==1)==1)>0 ) 
                                changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2err(b);
                            end
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2(b);
                        else
                            % Use (residual1-residual2)/(residual1) as change measure
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + sum((ar1.model(m1).data(I1(c)).yExpSimu(:,b)-ar1.model(m1).data(I1(c)).yExp(:,b)).^2,'omitnan');
                            totalErr( a, obsIndex )     = totalErr( a, obsIndex ) + sum((ar1.model(m1).data(I1(c)).yExpSimu(:,b)-ar1.model(m1).data(I1(c)).yExp(:,b)).^2,'omitnan');
                        end
                    else
                        % Compute R^2 of M1 instead
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + sum((ar1.model(m1).data(I1(c)).yExpSimu(:,b)-ar1.model(m1).data(I1(c)).yExp(:,b)).^2,'omitnan');
                        allData{ a, obsIndex }      = [ allData{ a, obsIndex }; ar1.model(m1).data(I1(c)).yExp(:,b) ];
                    end
                end
            end
            
            if ( ~opts.absrsq )
                % Compare with M2
                for b = 1 : length( ar2.model(m2).data(I2(c)).fy )
                    [~, obsIndex] = intersect(observables, ar2.model(m2).data(I2(c)).y{b});

                    if( ar2.model(m2).data(I2(c)).qFit(b) == 1 ) 
                        if ( ~opts.relerr )
                            if( ar2.config.fiterrors==1 || (ar2.config.fiterrors==0 && sum(ar2.qFit(ar2.qError==1)<2)>0) ) 
                                changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2err(b);
                            end
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2(b);
                        else
                            changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - sum((ar2.model(m2).data(I2(c)).yExpSimu(:,b)-ar2.model(m2).data(I2(c)).yExp(:,b)).^2);
                        end
                    end
                end
            end
        end
    end
    
    filtered = isinf(changeMatrix)|isnan(changeMatrix);
    
    % Use (residual1-residual2)/(residual1) as change measure
    if ( opts.relerr )
        filtered = filtered | (abs(totalErr)<1e-4);
        changeMatrix = changeMatrix ./ totalErr;
        changeMatrix(filtered) = 0;
    end
    
    % Compute R^2 of M1 instead
    if ( opts.absrsq )
        for a = 1 : size( allData, 1 )
            for b = 1 : size( allData, 2 )
                totalErr( a, b ) = sum( ( allData{a,b} - nanmean(allData{a,b}) ).^2 ,'omitnan');
            end
        end
        
        filtered = filtered | (abs(totalErr)<1e-4);
        changeMatrix = 1 - changeMatrix ./ totalErr;
        changeMatrix(filtered) = 0.5;
    end
        
    imagesc( changeMatrix ); hold on;
    for a = 1 : length( dataFiles )
        for b = 1 : length( observables )
            if ( ( abs( changeMatrix( a, b ) ) > 1e-2 ) && ~filtered( a, b ) )
                text( b, a, sprintf( '%.1f', changeMatrix( a, b ) ), 'Color', 'white', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
            end
        end
    end
    change = max(max(abs(changeMatrix)));
    
    % If the interactivity system is enabled, register the callbacks
    % and provide arInteractivity with the required data.
    if ( arInteractivity )
        arInteractivity( 'arCompareModel', dataFiles, observables, ar1, m1, ar2, m2, getPlotIDs(ar1.model(m1), dataFiles), getPlotIDs(ar2.model(m2), dataFiles), opts.titles_args );
    end
    
    set(gca, 'YTick', 1 : numel(dataFiles) );
    set(gca, 'XTick', 1 : numel(observables) );
    set(gca, 'YTickLabel', arNameTrafo(dataFiles) );
    set(gca, 'XTickLabel', arNameTrafo(observables) );
    if ( ~opts.absrsq )
        set(gca, 'CLim', [-change, change] );
        colormap( redgreencmap(256, 'Interpolation', 'sigmoid') );
    else
        set(gca, 'CLim', [0, 1] );
        colormap( flipud(redgreencmap(256, 'Interpolation', 'sigmoid') ) );
    end
    if ( opts.titles )
        title( sprintf( '%s - %s', arNameTrafo(opts.titles_args{1}), arNameTrafo(opts.titles_args{2}) ) );
    else
        title( sprintf( '%s - %s', arNameTrafo(name1), arNameTrafo(name2) ) );
    end
    colorbar;
    
    if ( ~opts.nogrid )
        for a = 1 : numel(dataFiles)
            line( [-0.5, numel(observables)+0.5], [a a]+0.5, 'Color', [0.15, 0.15, 0.15]) ;
        end
        for a = 1 : numel(observables)
            line( [a a]+0.5, [-0.5 numel(dataFiles)+0.5], 'Color', [0.15, 0.15, 0.15]) ;
        end
    end
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
