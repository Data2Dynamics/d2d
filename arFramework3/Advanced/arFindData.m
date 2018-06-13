% Finds datasets in the ar structure and returns their data indices 
% as a vector. Also checks for the specific condition parameters.
%
% Usage:
%   arFindData( (ar), (model no), name, conditions )
% or
%   arFindData( (ar), (model no), 'state', statename / number, conditions )
% or
%   arFindData( (ar), (model no), 'input', cell array of input strings )
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
% Examples of the third usage mode:
%    arFindData( ar, (model no), 'input', {'input_dcf * step1(0, 10, 5)'} )
%       Returns all data IDs that have an input with the function
%       'input_dcf * step(0, 10, 5)'. Note that the input function is not
%       parsed and must appear exactly as in the .def file.
%
% By default arFindData is permissive. It will return all conditions that
% match the criterion and include those that don't match or where
% insufficient information is available whether they match. To reverse this 
% behaviour, add the flag 'conservative'. This will reject any dataset
% where the flag is not set.
%
% Parameters enclosed by brackets are optional.
%
% Returns: List of IDs that correspond to the query and a cell array of the
% data names which match the search criterion. Note that the condition
% specifiers always have to come at the end. After the names and verbose
% flags.

function [olist, names, m] = arFindData( varargin )

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
    
    switches = { 'exact', 'verbose', 'names', 'state', 'conservative', 'input'};
    extraArgs = [ 0, 0, 0, 1, 0, 1 ];
    description = { ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} };
    [opts, varargin] = argSwitch( switches, extraArgs, description, 0, 'softmatching', varargin );
    
    exact = 0;
    if ( opts.exact )
        exact = 1;
    end
    
    returnNames = false;
    if ( opts.names )
        returnNames = true;
    end
            
    verbose = 0;
    if ( opts.verbose )
        verbose = 1;
    end
    
    if ( opts.state )
        observable = opts.state_args;
    end
    
    if ( length( varargin ) == 0 )
        string{1} = '';
    else
        if ( ischar( varargin{1} ) )
            string{1} = varargin{1};
            varargin = varargin(2:end);
        elseif ( iscell( varargin{1} ) )
            string = varargin{1};
            varargin = varargin(2:end);
        else
            string = '';
        end
    end
    
    if ( opts.state )
        hasObservables = scanObservables(m, observable, verbose);
    end
    
    if ( mod(length(varargin), 2) ~= 0 )
        error( 'Uneven number of condition arguments' );
    end
       
    condList = ones(length(varargin),1);

    % Find datasets with the correct name
    olist    = [];
    for b = 1 : length( string )
        for a = 1 : length( ar.model(m).data )
            % Does it have the observable we look for?
            if ( ( opts.state == 0 ) || (ismember( ar.model(m).data(a).name, hasObservables ) ) )
                % Do we match the name exactly or not?
                if ( ~exact )
                    % Does it have the name we look for?
                    if isempty(lower(string{b})) || ~isempty( strfind(lower(ar.model(m).data(a).name), lower(string{b}) ) )
                        olist = union( olist, a );
                    end
                else
                    if strcmp(ar.model(m).data(a).name, string{b} )
                        olist = union( olist, a );
                    end
                end
            end
        end
    end
    
    % Filter based on condition variables
    filt = ones(size(olist));
    for a = 1 : length( olist )
        for c = 1 : 2 : length(varargin)
            checked = false;
            
            % Is it in the condition variable list?
            ID = find( strcmp( ar.model(m).data(olist(a)).pold, varargin(c) ) );
            
            % No? Is it in the condition struct?
            val = '';
            if ( isempty( ID ) )
                nFields = length(ar.model(m).data(olist(a)).condition);
                for jcs = 1 : nFields
                    if strcmp( ar.model(m).data(olist(a)).condition(jcs).parameter, varargin(c) )
                        val = ar.model(m).data(olist(a)).condition(jcs).value;
                    end
                end
            else
                val = ar.model(m).data(olist(a)).fp{ID};
            end
            
            numval = str2num(val); %#ok
            
            if ( ~isempty(val) )
                % Is it a value?
                if ~isempty(numval)
                    condList(c) = 0;
                    % Does it pass? If so ==> OK
                    chk = varargin{c+1};
                    
                    if ( ~isnumeric(chk) )
                        chk = str2num( chk ); %#ok
                    end
                    
                    if ( numval ~= chk )
                        filt(a) = 0;
                    else
                        checked = true;
                    end
                else
                    % It is a string; Does it refer to a parameter?
                    Q = find(strcmp(ar.pLabel, ar.model(m).fp{ID}));
                    if ~isempty(Q)
                        condList(c) = 0;
                        numval = arGetPars(ar.pLabel{Q},0);
                        
                        if ( numval ~= varargin{c+1} )
                            filt(a) = 0;
                        else
                            checked = true;
                        end
                        
                        if (verbose)
                            disp('Warning: Filtering on open parameter');
                        end
                    end
                end
            end
                        
            % Did not match desired condition
            if ( opts.conservative && (checked == false) )
                filt(a) = 0;
            end
        end
        
        % Check whether it has the desired input
        if ( opts.input )
            if ( length( intersect( ar.model(m).data(olist(a)).fu, opts.input_args ) ) ~= length(opts.input_args) )
                filt(a) = 0;
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
            disp( sprintf( '[%d] %s %s', olist(a), ar.model(m).data(olist(a)).name, str ) ); %#ok
        end
    end
    
    condList = condList( 1 : 2 : end );
    if ( max( condList ) > 0 )
        for a = 1 : length( condList )
            if ( condList(a) > 0 )
                warning on; %#ok
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
        stateID = find(strcmp(ar.model(m).x,state));
        
        if ( isempty( stateID ) )
            stateID = find( strcmp( ar.model(m).z, state ) );
            if ( isempty( stateID ) )
                fprintf('Model states:\n');
                fprintf('%s\n', ar.model(m).x{:});
                fprintf('Model derived variables:\n');
                fprintf('%s\n', ar.model(m).z{:});
                error( 'Specified state or derived variable does not exist. Check capitalization?' );
            else
                stateString = ar.model(m).z{stateID};
            end
        else
            stateString = ar.model(m).x{stateID};
        end
    end
    
    % Find derived variables this state/variable appears in
    matches = [];
    for a = 1 : length( ar.model(m).fz )
        if find( strcmp(strsplit(ar.model(m).fz{a}, {'(',')','/','*','^','+','-',' '}), stateString), 1, 'last' )
            matches = union( matches, a );
        end
        if find( strcmp(strsplit(ar.model(m).z{a}, {'(',')','/','*','^','+','-',' '}), stateString), 1, 'last' )
            matches = union( matches, a );
        end
    end
    
    if verbose == 1
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
                if max(strcmp(strsplit(ar.model(m).data(a).fy{b}, {'(',')','/','*','^','+','-',' '}), ar.model(m).z{matches(c)} ))
                    matched = 1;
                    if( verbose )
                        fprintf('%s in %s\n', ar.model(m).z{matches(c)}, dataset);
                    end
                end
            end
            if( max(strcmp(strsplit(ar.model(m).data(a).fy{b}, {'(',')','/','*','^','+','-',' '}), stateString )) )
                matched = 1;
                if (verbose)
                    fprintf('%s in %s\n', stateString, dataset);
                end
            end
            if (matched)
                dataSets{end+1} = ar.model(m).data(a).name; %#ok<AGROW>
            end
        end
    end
    sets = dataSets;
    if ( verbose )
        fprintf('\n\nFiltering...\n\n');
    end
end