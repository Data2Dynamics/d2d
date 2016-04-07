% Simple helper function to determine for which datasets two models predict 
% differences
%
% Usage:
%   arCompareModel(ar1, m1, ar2, m2)
%
%   where ar1 and ar2 refer to ar structs and m1 and m2 to model indices.
%
% Note: If interactivity mode is on, you can click the matrix to see the 
% curves corresponding to the different models. Active interactivity mode 
% by invoking "arInteractivity on".

function arCompareModel(ar1, m1, ar2, m2)
        
    name1 = ar1.model(m1).name;
    name2 = ar2.model(m2).name;

    states      = intersect( {ar1.model(m1).x{:}, ar1.model(m1).z{:}}, {ar2.model(m2).x{:}, ar2.model(m2).z{:}} );
    observables = intersect( getObsNames(ar1.model(m1)), getObsNames(ar2.model(m2)) );
       
    dataFiles   = intersect( getDataNames(ar1.model(m1)), getDataNames(ar2.model(m2)) );

    changeMatrix = zeros( numel(dataFiles), numel(observables) );
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
                    if( (ar1.config.useFitErrorMatrix == 0 && ar1.config.fiterrors == 1) || ...
                            (ar1.config.useFitErrorMatrix==1 && ar1.config.fiterrors_matrix(m1,I1(c))==1) )
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2err(b);
                    end
                    changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) + ar1.model(m1).data(I1(c)).chi2(b);
                end
            end
            
            for b = 1 : length( ar2.model(m2).data(I2(c)).fy )
                [~, obsIndex] = intersect(observables, ar2.model(m2).data(I2(c)).y{b});
                
                if( ar2.model(m2).data(I2(c)).qFit(b) == 1 ) 
                    if( (ar2.config.useFitErrorMatrix == 0 && ar2.config.fiterrors == 1) || ...
                            (ar2.config.useFitErrorMatrix==1 && ar2.config.fiterrors_matrix(m2,I2(c))==1) )
                        changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2err(b);
                    end
                    changeMatrix( a, obsIndex ) = changeMatrix( a, obsIndex ) - ar2.model(m2).data(I2(c)).chi2(b);
                end
            end            
        end
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
    
    set(gca, 'YTick', [1 : numel(dataFiles)] );
    set(gca, 'XTick', [1 : numel(observables)] );
    set(gca, 'YTickLabel', dataFiles );
    set(gca, 'XTickLabel', observables );
    set(gca, 'CLim', [-change, change] );
    colormap( redgreencmap(256, 'Interpolation', 'sigmoid') );
    title( sprintf( '%s - %s', strTrafo(name1), strTrafo(name2) ) );
    colorbar;
end

% Find the corresponding data plots for both models being compared
function plotIDs = getPlotIDs(model, dataFiles)
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
