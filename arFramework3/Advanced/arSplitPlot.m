% arSplitPlot( [m], name, [title], filter )
%
% Split off a set of conditions for a separate plot.
%
% Sometimes there are simply too many conditions in a single plot. This
% function can split these into separate plots.
%
%   m           Model index                                             [1]
%   name        Name of the plot (needs to match exactly)
%   title       Title of the new plot                                   [same as input]
%   filter      Cell array of filters for the conditions we wish to
%               extract for the new plot. The format of the filter is as 
%               follows:
%               { condition name 1, anonymous function 1, condition name 2, anonymous function 2 }
%
% This function splits off a set of conditions for a separate plot. The 
% anonymous function will be evaluated with the condition listed in
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
% Example(s):
%   arSplitPlot(1, 'myExperiment', 'nobmp', {'input_bmp2', @(x)x==0} )
%
% It is also possible to specify multiple conditions to one anonymous function.
%   arSplitPlot(1, 'myExperiment', 'nobmp,negdorso', {{'input_bmp2', 'input_dorso'}, @(bmp,dcf)(bmp>0)||(dcf<0)} )
%
% See also arPlot, arMergePlot, arMergePlotMulti
function arSplitPlot( varargin )
    
    global ar;

    if ( nargin < 2 )
        help arSplitPlot;
        error( 'Must supply at least two arguments' );
    end
    
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
    legid = [];
    str = ar.model(m).plot(plotID).name;
    for d = 1 : numel( dataSets )
        accept = 1;
        for a = 1 : 2 : numel( filterList )
            % We want a cell array, even if it is one (single code path)
            if ~iscell( filterList{a} )
                filterList{a} = {filterList{a}};
            end
            
            % Collect variables required for anonymous function
            filterString = '';
            condiVals = {};
            collectedRequired = 1;
            for jf = 1 : numel( filterList{a} )
                filterString = [ filterString, ' ', filterList{a}{jf} ]; %#ok
                condiVal = arGetCondi( m, dataSets(d), filterList{a}{jf} );
                if ~isempty( condiVal )
                    condiVals{end+1} = sym(condiVal); %#ok
                else
                    collectedRequired = 0;
                end
            end
            if ( collectedRequired )
                try
                    accept = accept & feval( filterList{a+1}, condiVals{:} );
                catch
                    if ( nargin( filterList{a+1} ) > numel( condiVals ) )
                        error( 'Anonymous function takes more arguments than number of conditions specified %s -> %s', filterString, func2str( filterList{a+1} ) );
                    end
                    if ( nargin( filterList{a+1} ) < numel( condiVals ) )
                        error( 'Anonymous function takes fewer arguments than number of conditions specified %s -> %s', filterString, func2str( filterList{a+1} ) );
                    end
                    error( 'Unsuccessfull call to anonymous function %s -> %s', filterString, func2str( filterList{a+1} ) );
                end
            else
                warning( 'Missing condition. No filtering will occur for this anonymous function %s -> %s', filterString, func2str( filterList{a+1} ) );
            end
            if ( d == 1 )
                str = [ str '_' filterString ' -> ' func2str( filterList{a+1} ) ]; %#ok
            end
        end
        if ( accept )
            legid = [ legid d ]; %#ok
            ds = [ ds, dataSets(d) ]; %#ok
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
    
    % Update the legend
    if ( numel(legid) > 1 )
        newPlot.condition = ar.model(m).plot(plotID).condition(legid);
    else
        newPlot.condition = {};
    end
    
    ar.model(m).plot(end + 1) = newPlot;
    ar.model(m).qPlotYs(end + 1) = 1;
    ar.model(m).qPlotXs(end + 1) = 0;
    ar.model(m).qPlotVs(end + 1) = 0;
    
    % Add a field which indicates that this plot is not an 'original' plot
    ar.model(m).plot(end).origin = plotID;
    fprintf( 'Generated new plot %s (%d) with %d conditions\n', newPlot.name, numel(ar.model(m).plot), numel(ds) );
end
