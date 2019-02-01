% SSout = arSimilarityScore( [m], [cond], [thresh] )
%
%       m        Model index                                [1]
%       c        Condition index (or 'all')                 ['all']
%       thresh   R^2 threshold for when to show the group   [0.99]
% 
%       SSout    matrix of SimmilarityScores
% 
%   Computes the R^2 between predictions after performing linear regression 
%   between each pair. High R^2 indicates potential for model reduction. 
%   When multiple conditions are passed, the minimum is returned (states 
%   have to be the same over all conditions for the model reduction to make 
%   sense).
%   Creates a plot (imagesc) of Similarity scores.

function SSout = arSimilarityScore( m, cond, thresh )
    global ar;
    
    if ( nargin < 1 )
        m = 1;
    end
    if ( nargin < 2 )
        cond = 'all';
    end
    if ( nargin < 3 )
        thresh = 0.99;
    end
    
    nx = numel( ar.model.x );
    SS = ones( nx, nx );
    
    if ( strcmp( cond, 'all' ) )
        condis = 1 : numel( ar.model(m).condition );
    else
        condis = cond;
    end
    
    for c = condis
        for jx1 = 1 : nx
            for jx2 = 1 : nx
                x1 = ar.model(m).condition(c).xFineSimu(:,jx1);
                x2 = ar.model(m).condition(c).xFineSimu(:,jx2);
    
                % Use R^2 as metric
                res = x1 - x1\x2*x2;
                SS(jx1, jx2) = min( SS(jx1, jx2), 1 - sum(res.^2) / sum(x1.^2) );
            end
        end
    end
    
    figure;
    imagesc(SS);
    set(gca, 'CLim', [0, 1]);
    title( 'State similarity score (R^2 after linear regression)' );
    
    uq = unique(SS>thresh, 'rows');
    uq = uq( sum(uq, 2) > 1, : );
    
    if ~isempty( uq )
        fprintf( '(Nearly) identical states with R^2 threshold of %g:\n', thresh );
        for a = 1 : size( uq, 1 )
            fprintf( '%s ', ar.model(m).x{ uq(a,:) } );
            fprintf( '\n' );
        end
    end
    colorbar;
    set(gca, 'YTick', 1:numel(ar.model(m).x));
    set(gca, 'YTickLabels', strrep(ar.model(m).x, '_', '\_') );
    set(gca, 'XTick', 1:numel(ar.model(m).x));
    set(gca, 'XTickLabels', strrep(ar.model(m).x, '_', '\_') );    
    set(gca, 'XTickLabelRotation', 45 );
    
    if ( nargout > 1 )
        SSout = SS;
    end