% function textBar(labels, y, [ax], [col], [lb], [ub])
%
% Generates a bar diagram with text labels.
%
%   labels      cell array of labels
%   y           bar height
%   col         color
%   lb          lower bound of the figure
%   ub          upper bound of the figure

function textBar(labels, y, ax, col, lb, ub)

    clamp = @(x)max(min([1,1,1], x), [0,0,0]);
    if ( ~exist( 'ax', 'var' ) || isempty(ax) )
        ax = axes;
    end
    if ( ~exist( 'col', 'var' ) || isempty(col) )
        col = [.2 0 .8];
    end

    if ( size( col, 1 ) == 1 )
        col = repmat( col, numel(y), 1 );
    end
    
    badIdx = isinf(y) | isnan(y);
    
	for i = 1 : numel(y)
        b = barh( ax, i, y(i) ); hold on;
        set( b, 'EdgeColor', col(i,:) );
        set( b, 'FaceColor', clamp(col(i,:) + 0.9) );
        set( b, 'LineWidth', 1.05 );
        if ( exist('ub', 'var' ) )
            line( [lb(i), ub(i)], [i, i], 'LineWidth', 1.05, 'Color', col(i,:) );
            line( [lb(i), lb(i)], [i-.2, i+.2], 'LineWidth', 1.05, 'Color', col(i,:) );
            line( [ub(i), ub(i)], [i-.2, i+.2], 'LineWidth', 1.05, 'Color', col(i,:) );
        end
    end
    
    set( ax, 'YTick', 1 : numel(y) )
    set( ax, 'YTickLabel', labels, 'fontsize', 9 );
    
    ys = y(~badIdx);
    xlims = [ min([0, min(ys) - 0.1 * abs(min(ys))]), max(ys) + 0.1 * abs(max(min(ys))) ];
    if ( exist('ub', 'var' ) )
        xlims = [ min([0, min(lb) - 0.1 * abs(min(lb))]), max(ub) + 0.1 * abs(max(min(ub))) ];
    end
    if ( xlims(2) == xlims(1) )
        xlims(2) = xlims(1)+0.1;
    end
    
    xlim(xlims);
    ylim( [ 0.5, numel(y) + 0.5 ] );
    hold on;
    plot( repmat( xlims(1) + 0.05 * ( xlims(2) - xlims(1) ), 1, sum(badIdx) ), find(badIdx), 'xr', 'MarkerSize', 10 );
    
    set(gca, 'YDir', 'Reverse' );
    grid off; set(gca, 'XGrid', 'on');
end