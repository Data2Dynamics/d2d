% Finds datasets in the ar structure and returns their data indices 
% as a vector. Also checks for the specific condition parameters.
%
% Usage:
%   arFindData( (ar), (model no), name, conditions )
% or
%   arFindData( (ar), (model no), 'state', statename / number, conditions )
%
% Examples of the first usage mode:
%   arFindData( ar, (model no), 'mydata' )
%       Returns all data IDs whose name contains "mydata"
%   arFindData( ar, (model no), {'mydata', 'potato'} )
%       Returns all data IDs whose name contains "mydata" or "potato"
%   arFindData( ar, (model no), 'mydata', 'dose', '100', 'actd', '1' )
%       Returns all data IDs whose name contains "mydata" and who
%       correspond to dose being set to 100 and actd being set to 1.
%   arFindData( ar, (model no), 'mydata', 'verbose' )
%       Returns all data IDs whose name contains "mydata" and prints
%       them.
%   arFindData( ar, (model no), 'mydata', 'names' )
%       Only return the full names of the found datasets
%
% Examples of the second usage mode:
%    arFindData( ar, (model no), 'state', 'SOCS3' )
%       Returns all data IDs that have the state SOCS3 in their observation
%       function.
%    arFindData( ar, (model no), 'state', 'SOCS3', 'verbose' )
%       Provides some additional debug information
%    arFindData( ar, (model no), 'state', 'SOCS3', 'names', 'il6', 100 )
%       Returns the same by name, and also filters on the value of il6 in
%       the experiment.
%
% Parameters enclosed by brackets are optional.
%
% Returns: List of IDs that correspond to the query and a cell array of the
% data names which match the search criterion. Note that the condition
% specifiers always have to come at the end. After the names and verbose
% flags.

function [olist names m] = arFindData( varargin )

    if nargin == 0
        help arFindData;
        return;
    end

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
        
    observable = [];
    if ( ischar( varargin{1} ) )
        if strcmp( varargin{1}, 'state' )
            if ( length( varargin ) < 2 )
                error('Insufficient arguments');
            end
            observable = varargin{2};
            varargin = varargin(3:end);
        else
            string{1} = varargin{1};
            varargin = varargin(2:end);
        end
    elseif ( iscell( varargin{1} ) )
        string = varargin{1};
        varargin = varargin(2:end);
    else
        error( 'Please supply a data name' );
    end
    
    verbose = 0;
    if ( length(varargin) > 0 )
        if ( ischar( varargin{1} ) )
            if ( strcmp( varargin{1}, 'verbose' ) )
                verbose = 1;
                varargin = varargin(2:end);
            end
        end
    end
    
    if ( ~isempty(observable) )
        string = scanObservables(m, observable, verbose);
    end
    
    % Return names instead
    returnNames = false;
    if (length( varargin ) > 0)
        if ischar( varargin{1} )
            if ( strcmp( varargin{1}, 'names' ) )
                returnNames = true;
                varargin = varargin(2:end);
            end
        end
    end

    if ( mod(length(varargin), 2) ~= 0 )
        error( 'Uneven number of condition arguments' );
    end
       
    condList = ones(length(varargin),1);

    % Find datasets with the correct name
    olist    = [];
    for b = 1 : length( string )
        for a = 1 : length( ar.model(m).data )
            if ~isempty( findstr(lower(ar.model(m).data(a).name), lower(string{b}) ) )
                olist = union( olist, a );
            end
        end
    end
    
    % Filter based on condition variables
    filt = ones(size(olist));
    for a = 1 : length( olist )
        for b = 1 : length(ar.model(m).data(olist(a)).condition)  
            par = ar.model(m).data(olist(a)).condition(b).parameter;
            val = ar.model(m).data(olist(a)).condition(b).value;
            for c = 1 : 2 : length(varargin)
                % Is it in the condition variable list?
                if ( strcmp( par, varargin(c) ) )
                    condList(c) = 0;
                    % Does it pass?
                    if ( str2num(val) ~= varargin{c+1} )
                        filt(a) = 0;
                    end
                else
                    % It's not: is it in the model default variable list?
                    Q = find(strcmp(ar.model.p,varargin(c)));
                    if ~isempty(Q)
                        val = str2num(ar.model.fp{Q});
                        % Is it a value?
                        if ~isempty(val)
                            condList(c) = 0;
                            if ( val ~= varargin{c+1} )
                                filt(a) = 0;
                            end
                        else
                            % It is a string; Does it refer to a parameter?
                            Q = find(strcmp(ar.pLabel, ar.model.p{Q}));
                            if ~isempty(Q)
                                condList(c) = 0;
                                val = ar.p(Q);
                                if ( val ~= varargin{c+1} )
                                    filt(a) = 0;
                                end
                                if (verbose)
                                    disp('Warning: Filtering on open parameter');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
                
    % Remove the ones that don't pass the condition filter
    olist( filt == 0 )   = [];
    
    % Fetch the names
    names = cell( length( olist ), 1 );
    for a = 1 : length( olist )
        names{a} = ar.model(m).data(olist(a)).name;
    end
    names = unique(names);
    
    % Output
    for a = 1 : length( olist )
        str = [];
        for b = 1 : length(ar.model(m).data(olist(a)).condition)     
            par = ar.model(m).data(olist(a)).condition(b).parameter;
            val = ar.model(m).data(olist(a)).condition(b).value;
            str = sprintf( '%s | %s = %s', str, par, val );
        end
        if ( verbose )
            disp( sprintf( '[%d] %s %s', olist(a), ar.model(m).data(olist(a)).name, str ) );
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
end

function sets = scanObservables( m, state, verbose )
    global ar;
    
    % Convert choice to numeric state number and see whether it exists
    if ( ischar( state ) )
        stateID = find(strcmp(ar.model.x,state));
        
        if ( isempty( stateID ) )
            stateID = find( strcmp( ar.model.z, state ) );
            if ( isempty( stateID ) )
                fprintf('Model states:\n');
                fprintf('%s\n', ar.model.x{:});
                fprintf('Model derived variables:\n');
                fprintf('%s\n', ar.model.z{:});
                error( 'Specified state or derived variable does not exist. Check capitalization?' );
            else
                stateString = ar.model.z{stateID};
            end
        else
            stateString = ar.model.x{stateID};
        end
    end
    
    % Find derived variables this state/variable appears in
    matches = [];
    for a = 1 : length( ar.model(m).fz )
        if max( find( strcmp(strsplit(ar.model(m).fz{a}, {'(',')','/','*','^','+','-',' '}), stateString) ) )
            matches = union( matches, a );
        end
        if max( find( strcmp(strsplit(ar.model(m).z{a}, {'(',')','/','*','^','+','-',' '}), stateString) ) )
            matches = union( matches, a );
        end
    end
    
    if ( verbose )
        fprintf( 'State/derived variable %s appears in the following (derived) variable: \n', stateString );
        for a = 1 : length( matches )
            fprintf( '%s = %s\n', ar.model(m).z{ matches(a) }, ar.model(m).fz{ matches(a) } );
        end      
    end
    
    if ( verbose )
        fprintf('\nAppears in datasets:\n' );
    end
    
    % Find dataset observables, these appear in
    dataSets = {};
    for a = 1 : length( ar.model(m).data )
        for b = 1 : length( ar.model(m).data(a).fy )
            matched = 0;
            dataset = ar.model(m).data(a).name;
            for c = 1 : length( matches )
                % Match against observable
                if max(strcmp(strsplit(ar.model(m).data(a).fy{b}, {'(',')','/','*','^','+','-',' '}), ar.model.z{matches(c)} ));
                    matched = 1;
                    if( verbose )
                        fprintf('%s in %s\n', ar.model.z{matches(c)}, dataset);
                    end
                end
            end
            if( max(strcmp(strsplit(ar.model(m).data(a).fy{b}, {'(',')','/','*','^','+','-',' '}), stateString )) );
                matched = 1;
                if (verbose)
                    fprintf('%s in %s\n', stateString, dataset);
                end
            end
            if (matched)
                dataSets{end+1} = ar.model(m).data(a).name;
            end
        end
    end
    sets = dataSets;
    if ( verbose )
        fprintf('\n\nFiltering...\n\n');
    end
end