% function arSubplot(m, n, id, [name], [varargin])
%
% A little subplot manager that prevents you from having to keep track of
% the plot handles. This routine works just like subplot, except that it 
% checks if the plot already exists and makes that axis current if so.
% If the specifications of the figure change (m or n) then the figure is
% reset.
%
% If no name is specified, the function uses the name of the matlab
% function that called it. This is so that this function can simply be used
% as a drop in replacement in functions which only deal with a single
% figure.
%
% If m is left empty ([]), then it will automatically look for a nice plot
% placement. m then needs to contain the number of plots desired.
%

function ax = arSubplot(m, n, id, name, varargin)
    global arSubplotMgr;
    
    if ~isstruct( arSubplotMgr )
        arSubplotMgr = struct;
    end
    if ~exist('name', 'var')
        stack = dbstack;
        name = strrep( stack(end).file, '.', '___' );
    end
    if ~exist( 'm', 'var' )
        error( 'Please specify number of plots' );
    end
    if ~exist( 'n', 'var' )
        error( 'Please specify number of plots' );
    end
    if ~exist( 'id', 'var' )
        error( 'Please specify plot index' );
    end
    
    
    % Automatically decide on the layout
    if ( isempty( m ) )
        [m, n] = fd( n );
    end
    
    if ~isfield( arSubplotMgr, name ) || ~isfield( arSubplotMgr.(name), 'fig' ) || ( ~isvalid( arSubplotMgr.(name).fig ) )
        arSubplotMgr.(name) = struct;
        figure;
        set(gcf, 'Name', name);
        newPlot(name, m, n, id, varargin{:} );
    else
        % Make figure current
        figure( arSubplotMgr.(name).fig );
        
        % Check if the specifications changed. In this case, we need to
        % reset the plot
        if ( ( arSubplotMgr.(name).m ~= m ) || ( arSubplotMgr.(name).n ~= n ) )
            clf;
            newPlot( name, m, n, id, varargin{:} );
        end
        
        % Does the axis exist already? => Use it
        if ( ~isnan( arSubplotMgr.(name).axes(id) ) && ishandle( arSubplotMgr.(name).axes(id) ) )
            axes( arSubplotMgr.(name).axes(id) );
        else
            arSubplotMgr.(name).axes(id) = subplot( m, n, id, varargin{:} );
        end
    end
    
    ax = arSubplotMgr.(name).axes(id);
end

function newPlot( name, m, n, id, varargin )
    global arSubplotMgr;
    
    arSubplotMgr.(name).m        = m;
    arSubplotMgr.(name).n        = n;
    arSubplotMgr.(name).fig      = gcf;
    arSubplotMgr.(name).axes     = nan( m * n, 1 );
    arSubplotMgr.(name).axes(id) = subplot( m, n, id, varargin{:} );
end

function [n, m] = fd( num )
    nx = ceil( sqrt( num ) );
    ny = ceil( num / nx );

    % Test adjacent ones for a better 'fit' (full use of subplots)
    % Don't test too far from the square layout, since this would make
    % poor use of screen space.
    mnx = max( [ floor( nx / 1.5 ), 2 ] );
    mxx = ceil( nx * 1.5 );
    testRange = mnx:mxx;
    
    % Sort by the one closest to an equal distribution along X and Y
    [~, i] = sort( abs( testRange - nx ) );
    testRange = testRange(i);
    
    % Try to find a nice 'natural' break
    for nxc = testRange
        if ceil( num / nxc ) == ( num / nxc )
            n = nxc;
            m = num / nxc;
            return
        end
    end
    n = nx;
    m = ny;
end