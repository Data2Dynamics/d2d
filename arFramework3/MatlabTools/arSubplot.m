% function arSubplot(m, n, id, name, varargin)
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

function ax = arSubplot(m, n, id, name, varargin)
    global arSubplotMgr;
    
    if ~isstruct( arSubplotMgr )
        arSubplotMgr = struct;
    end
    if ~exist('name', 'var')
        stack = dbstack;
        name = strrep( stack(end).file, '.', '___' );
    end
    
    if ~isfield( arSubplotMgr, name )
        arSubplotMgr.(name) = struct;
        clf;
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
