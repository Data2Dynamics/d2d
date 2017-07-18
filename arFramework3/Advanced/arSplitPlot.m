% Split off a set of conditions for a separate plot
%
% Sometimes there are simply too many conditions in a single plot. This
% function can split these into separate plots.
%
% Usage:
%   arSplitPlot( (m), name, (title), filter )
%
%   m           -   Model number (default: 1)
%   name        -   Name of the plot (needs to match exactly)
%   title       -   Title of the new plot
%   filter      -   Cell array of filters for the conditions we wish to
%                   extract for the new plot
%
% The format of the filter is as follows
%   { condition name 1, anonymous function 1, condition name 2, anonymous function 2 }
%
% The anonymous function will be evaluated with the condition listed in
% condition name for each dataset in the original source plot. If it
% evaluates to true for all conditions, the data will be copied to the new 
% plot. Otherwise, the data will be ignored.
%
% The conditions are passed as symbolic values, which need to be cast to a
% string in case a string comparison is desired. For numeric values, they
% can simply be compared directly.
%
% Note that arguments between brackets are optional
%
% Example usage:
%   arSplitPlot(1, 'myExperiment', 'nobmp', {'input_bmp2', @(x)x==0} )

function arSplitPlot( varargin )
    
    global ar;

    if ( isnumeric(varargin{1} ) )
        m = varargin{1};
        varargin = varargin(2:end);
    else
        m = 1;
    end
    
    if ( ischar( varargin{1} ) )
        name = varargin{1};
        varargin = varargin(2:end);
    else
        error( 'Please supply a plot name' );
    end

    plotID = arFindPlot(name, 'exact');
    newPlot = ar.model(m).plot(plotID);
    dataSets = ar.model(m).plot(plotID).dLink;
    
    if ( isempty( varargin ) )
        error( 'Please specify conditions on which to filter' );
    end
    
    if ( ischar( varargin{1}  ) )
        outputName = varargin{1};
        varargin = varargin(2:end);
    end
    
    if ( isempty( varargin ) )
        error( 'Please specify conditions on which to filter' );
    end
    
    if ( iscell( varargin{1} ) )
        filterList = varargin{1};
    end
       
    % Filter on condition variables
    ds = [];
    str = ar.model(m).plot(plotID).name;
    for d = 1 : numel( dataSets )
        accept = 1;
        for a = 1 : 2 : numel( filterList )
            if ( d == 1 )
                str = [ str '_' filterList{a} ' -> ' func2str( filterList{a+1} ) ]; %#ok
            end
            condiVal = arGetCondi( m, dataSets(d), filterList{a} );
            if ~isempty( condiVal )
                accept = accept & feval( filterList{a+1}, sym(condiVal) );
            end
        end
        if ( accept )
            ds = [ ds, dataSets(d) ];
        end
    end
    
    if exist( 'outputName', 'var' )
        newPlot.name = outputName;
    else
        newPlot.name = str;
    end
    
    if isempty( ds )
        error( 'No datasets were accepted' );
    end
    
    newPlot.dLink = ds;
    ar.model(m).plot(end + 1) = newPlot;
    ar.model(m).qPlotYs(end + 1) = 1;
    ar.model(m).qPlotXs(end + 1) = 0;
    ar.model(m).qPlotVs(end + 1) = 0;
    
    % Add a field which indicates that this plot is not an 'original' plot
    ar.model(m).plot(end).origin = plotID;
    fprintf( 'Generated new plot %s (%d) with %d conditions\n', newPlot.name, numel(ar.model(m).plot), numel(ds) );
end
