% Function used to control some additional interactivity enhancement
% callbacks. Since these may slow down the plots on some systems, they are
% off by default. Activate them by calling "interactivity on".
%
% Note that interactivity has to be on when the plot is created for the
% callbacks to be registered. This is to keep the profile low when
% interactivity is not desired.
%
% Currently supported interactive modes:
%   plePlot =>  Click on the individual parameter curves to change legend
%               entries.


function interactive = arInteractivity( varargin )

    global arInteractivity;
    
    if ( isempty( arInteractivity ) )
        initInteractivity;
    end
    if ( ~ischar( varargin{1} ) )
        error( 'First argument should be a string' );
    end
    
    % No arguments means we want to know whether the interactivity system
    % is operational
    if (nargin == 0)
        interactive = arInteractivity.active;
        return;
    end

    if strcmp( lower(varargin{1}), 'on' )
        arInteractivity.active = 1;
    end
    if strcmp( lower(varargin{1}), 'off' )
        arInteractivity.active = 0;
    end
    
    % PLE interactivity functions
    if strcmp( varargin{1}, 'ple' );
        arInteractivity.ple.legend = varargin{2};
        set(gcf,'WindowButtonDownFcn', @(hObject, eventData)pleFcn2(hObject, eventData) );
    end
    
end

function initInteractivity()
    global arInteractivity;
    
    arInteractivity.active = 0;
end

function pleFcn2(hObject, eventData)
    global arInteractivity;
    
    if ( arInteractivity.active )
        userData = arInteractivity.ple.legend;
        cp      = get(gca, 'CurrentPoint'); 
        try
            name = userData.legends( find( gco == userData.handles ) );
            %text( cp(1,1), cp(1,2), name );
            if ( ~isempty( name ) )
                userData.currentLegend.handles(userData.currentLegend.ID) = gco;
                userData.currentLegend.legends(userData.currentLegend.ID) = name;
                userData.currentLegend.ID = mod( userData.currentLegend.ID, 5 ) + 1;
                legend( gca, userData.currentLegend.handles, userData.currentLegend.legends );
            end
        catch
            error( 'Unknown interactivity error in callback' );
        end
        arInteractivity.ple.legend = userData;
    end
end