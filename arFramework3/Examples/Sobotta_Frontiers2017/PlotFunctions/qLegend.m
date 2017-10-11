% Function to move the legend to a position far from the figure axes (smaller chance the plot will be obfuscated)

function h = qLegend( varargin )
    global errorBarMode
    if ( errorBarMode < 3 )
        h = legend( varargin{:}, 'Box', 'off', 'Location', 'Best' );
        k = get( h, 'position' );
        set( h, 'position', [k(1) .2*k(2) k(3:4)] );
    else
        h = legend( varargin{:}, 'Box', 'off' );
        set(h,'Location','SouthEast');
    end
end
