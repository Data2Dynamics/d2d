% Finds datasets in the ar structure and returns their data indices 
% as a vector. Also checks for the specific condition parameters.
%
% Usage:
%   arFindData( ar, (model no), name, conditions )
%
% Example:
%   arFindData( ar, (model no), 'mydata' )
%       Returns all condition IDs whose name contains "mydata"
%   arFindData( ar, (model no), 'mydata', 'dose', '100', 'actd', '1' )
%       Returns all condition IDs whose name contains "mydata" and who
%       correspond to dose being set to 100 and actd being set to 1.
%   arFindData( ar, (model no), 'verbose', 'mydata' )
%       Returns all condition IDs whose name contains "mydata" and prints
%       them.
%
% The argument ar is optional. If not specified, the global ar structure is
% used. 
%
% Returns: List of IDs that correspond to the query.

function [olist names m] = arFindDataSet( varargin )

    global ar;
    if ( isstruct( varargin{1} ) )
        ar = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    end
    
    % Did we specify a model number? If not, assume 1
    m = 1;
    if ( isnumeric( varargin{1} ) )
        m = varargin{1};
        if ( m > length(ar.model) )
            error( 'Model %d does not exist in ar structure (only has %d models)', m, length(ar.model) );
        end
        if (length( varargin ) > 1)
            varargin = varargin(2:end);
        else
            error( 'Need to supply data name' );
        end
    end
        
    if ( ischar( varargin{1} ) )
        string = varargin{1};
        varargin = varargin(2:end);
    else
        error( 'Need to supply data name' );
    end
    
    verbose = 0;
    if ( nargin > 2 )
        if ( ischar( varargin{1} ) )
            if ( strcmp( varargin{1}, 'verbose' ) )
                verbose = 1;
                varargin = varargin(2:end);
            end
        end
    end
    
    % Return names instead
    returnNames = false;
    if (length( varargin ) > 0)
        if isnumeric( varargin{1} )
            if (varargin{1} == 2)
                returnNames = true;
            end
            varargin = varargin(2:end);
        end
    end
    
    if ( mod(length(varargin), 2) ~= 0 )
        error( 'Uneven number of condition arguments' );
    end
       
    condList = ones(length(varargin),1);

    % Find datasets with the correct name
    olist    = [];
    names    = {};
    if ~iscell( string )
        for a = 1 : length( ar.model.data )
            if ~isempty( findstr(lower(ar.model.data(a).name), lower(string) ) )
                olist = [ olist a ];
                names = { names{:}, ar.model.data(a).name };
            end
        end
    else
        for b = 1 : length( string )
            for a = 1 : length( ar.model.data )
                if ~isempty( findstr(lower(ar.model.data(a).name), lower(string{b}) ) )
                    olist = [ olist a ];
                    names = { names{:}, ar.model.data(a).name };
                end
            end
        end
    end
    
    % Filter based on condition variables
    filt = ones(size(olist));
    for a = 1 : length( olist )
        for b = 1 : length(ar.model.data(olist(a)).condition)  
            par = ar.model.data(olist(a)).condition(b).parameter;
            val = ar.model.data(olist(a)).condition(b).value;
            for c = 1 : 2 : length(varargin)
                if ( strcmp( par, varargin(c) ) )
                    condList(c) = 0;
                    if ( ~strcmp( val, num2str(varargin{c+1}) ) )
                        filt(a) = 0;
                    end
                end
            end
        end
    end
    
    % Remove the ones that don't pass the condition filter
    olist( filt == 0)   = [];
    names( filt == 0 )  = [];
    
    % Output
    for a = 1 : length( olist )
        str = [];
        for b = 1 : length(ar.model.data(olist(a)).condition)     
            par = ar.model.data(olist(a)).condition(b).parameter;
            val = ar.model.data(olist(a)).condition(b).value;
            str = sprintf( '%s | %s = %s', str, par, val );
        end
        if ( verbose )
            disp( sprintf( '[%d] %s %s', olist(a), ar.model.data(olist(a)).name, str ) );
        end
    end
    
    condList = condList( 1 : 2 : end );
    if ( max( condList ) > 0 )
        for a = 1 : length( condList )
            if ( condList(a) > 0 )
                warning on;
                warning( 'Condition %s = %s was not used as a filter! Please check the inputs.', varargin{2*a-1}, num2str(varargin{2*a}) );
            end
        end
    end
    
    if (returnNames == true)
        olist = names;
    end