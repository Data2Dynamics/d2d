% Function which can be used to disable datasets by name or reference ar
% structure
% 
% Usage:
%   function ar = arDisableData( (ar), data, (verbose) )
%
%   The variable "data" either contains a cell array of strings referring 
%   to the datasets to disable, the keyword 'all' or a reference ar
%   structure.
%
%   Note: arFindData may be used to find the names of datasets more easily
%
% Examples:
%
%   ar = arDisableData( ar, {'myData1', 'myData2'} )
%   Disables myData1 and myData2 in the ar structure.
%
%   ar = arDisableData( ar, arFindData( ar, 'myData', 'names' ) )
%   Disables all datasets containing 'myData' in the name

function tar = arDisableData( varargin )

    global ar;
    if ( isstruct( varargin{1} ) )
        ar = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    end
    from = varargin{1};
    
    if ( length( varargin ) > 1 )
        verbose = varargin{2};
    else
        verbose = 0;
    end
    
    disableAll = 0;
    if isstruct( from )
        % Get the datasets to disable
        dataSets = {};     
        for c = 1 : length( from.model )
            for d = 1 : length( from.model(c).data )
                dataSets{end+1} = from.model(c).data(d).name;
            end
        end
    else
        if strcmp( from, 'all' )
            disableAll = 1;
        end
        dataSets = from;
    end
    
    % Make sure they are not used again
    for c = 1 : length( ar.model )
        for d = 1 : length( ar.model(c).data )
            if ( disableAll || ( max( strcmp( dataSets, ar.model(c).data(d).name ) ) ) )
                if ( verbose )
                    if ( max( ar.model(c).data(d).qFit ) > 0 )
                        disp( sprintf( 'Disabled dataset %s', ar.model(c).data(d).name ) );
                    else
                        disp( sprintf( 'Dataset %s already inactive', ar.model(c).data(d).name ) );
                    end
                end
                ar.model(c).data(d).qFit = 0 * ar.model(c).data(d).qFit;
            end
        end
    end
    
    tar = ar;
end
