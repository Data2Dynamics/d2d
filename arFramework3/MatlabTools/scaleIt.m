%
% Function for merging CSV files
%
% Scales data based on error model.
% Note: This function is currently work in progress
%
% Usage:
%   scaleIt( names, outFileName, options )
%
%   CSV file must have following columns
%     nExpID  = replicate
%     input_  = inputs per condition (default, can be overriden by setting "inputMask")
%     time    = dependent variable must be called t, T, time or Time
%
%   All other columns are assumed data. Scaling factors will be estimated per
%   observable/replicate combination.
%
%   Options:
%     delimiter          Use a custom delimiter for the csv file (default = ',')
%     obsGroups          Specify observations that share a scaling parameter (e.g. same antibody)
%                        You have to specify these as cell array e.g:
%                        { { 'nSTAT3', 'cSTAT3' }, { 'npSTAT3', 'cpSTAT3' } }
%     inputMask          Use a custom mask to denote that something should be 
%                        used as condition variable (default is input_)  
%     depVar             Use different dependent variable (default = 't', 'T', 'time' or 'Time').
%     expID              Use a different column name to hold the experiment
%                        index (default = nExpID). This column will be used
%                        in combination with the observation to indicate the 
%                        scaling factor.
%     restrictObs        Restrict estimation to list of observables
%     ignoreMask         Data columns to ignore. Start element with * to
%                        indicate that an exact match is not required, only
%                        a substring which contains the mask
%     removeMask         Data columns to remove. Start element with * to
%                        indicate that an exact match is not required, only
%                        a substring which contains the mask
%     twocomponent       Use two component error model (warning: poorly
%                        tested so far)
%     logtrafo           Do the scaling in logarithmic space
%     prescale           Prescale all data with a simple factor prior to
%                        any analysis (with logtrafo, this is done before
%                        the trafo). Supply either function or value.
%     rescale            Adjust scaling factors s.t. maximum becomes 1
%     range              Specify range of values to use for independent
%                        variable
%     varanalysis        Show variance analysis
%     appendcolumn       Add a column with a value to the data
%     splitconditions    Provide independent scaling for each individual condition
%     excludeConditions  Cell array of condition filters (example:
%                        {{'input_il6', @(t)t>10}}
%     nofileincrement    Do not assume different files have their own
%                        scaling factor
%     samescale          Do not scale any data
%     firstoutput        Which field should come first? Note: You may need
%                        when specifying a different dependent variable and 
%                        want to use the data in D2D (which requires time to 
%                        be the first column)
%
% To do: Offsets

function out = scaleIt( names, outFile, varargin )

    if ( nargin < 2 )
        help scaleIt;
        error( 'Insufficient arguments' );
    end

    % Load options
    verbose = 0;
    switches = { 'nofileincrement', 'delimiter', 'obsgroups', 'inputmask', 'depvar', 'expid', 'restrictobs', 'ignoremask', 'twocomponent', 'logtrafo', 'rescale', 'range', 'varanalysis', 'samescale', 'prescale', 'appendcolumn', 'splitconditions', 'removemask', 'excludeconditions', 'firstoutput', 'dlfudge' };
    extraArgs = [ 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0 ];
    description = { ...
    {'', 'Do not increment nExpID for different files'} ...
    {'', 'Custom delimiter specified'} ...
    {'', 'Using custom observation/scaling factor pairing'} ...
    {'', 'Using input mask'} ...
    {'', 'Using custom dependent variable name'} ...
    {'', 'Using custom experiment ID column'} ...
    {'', 'Restricting to specific specified observables'} ...
    {'', 'Using ignore mask'} ...,
    {'', 'Using two component error model'} ...
    {'', 'Using lognormal error model'} ...
    {'', 'Rescaling result'} ...
    {'', 'Specified time point range to use'} ...
    {'', 'Variance analysis requested'} ...
    {'', 'Not fitting scales'} ...
    {'', 'Prescaler specified'} ...
    {'', 'Appending a column'} ...
    {'', 'Splitting conditions'} ...
    {'', 'Using removal mask'} ...
    {'', 'Excluding specific conditions'}, ...
    {'', 'Specifying first output'}, ...
    {'', 'Using detection limit fudge (not recommended)'}, ...
    };
    opts = argSwitch( switches, extraArgs, description, verbose, varargin );
    
    timeVars = {'t', 'T', 'time', 'Time'};
    expVar = 'nExpID';
    inputMask = 'input_';
    obsGroups = [];

    delimiter = ',';
    if ( opts.delimiter )
        delimiter = opts.delimiter_args;
    end
    if ( opts.obsgroups )
        obsGroups = opts.obsgroups_args;
    end
    if ( opts.inputmask )
        inputMask = opts.inputmask_args;
    end
    if ( opts.depvar )
        timeVars = {opts.depvar_args};
    end
    if ( opts.expid )
        expVar = opts.expid_args;
    end
    restrictobs = {};
    if ( opts.restrictobs )
        restrictobs = opts.restrictobs_args;
    end
    ignoreMask = {};
    if ( opts.ignoremask )
        ignoreMask = opts.ignoremask_args;
    end
    removeMask = {};
    if ( opts.removemask )
        removeMask = opts.removemask_args;
    end
    if ( opts.prescale )
        if ( isnumeric(opts.prescale_args) )
            if ( numel(opts.prescale_args) == 1 )
                opts.prescale_args=@(x)x*opts.prescale_args;
            else
                error( 'Only single value or function allowed for prescaling' );
            end
        end
    end
    
    % Default mapping
    trafo           = @(x) x;
    invTrafo        = @(x) x;
    dTrafo          = @(x) 1;
    diTrafo         = @(x) 1;
    
    if ( opts.logtrafo )
        trafo       = @(x) log10(x);
        invTrafo    = @(x) 10.^(x);
        % Transform derivatives
        dTrafo      = @(x) 1./(x.*log(10));
        diTrafo     = @(x) 10.^x * log(10);
    end
    
    % Simplest error model
    errModel        = @(fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks) model(fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks, trafo, invTrafo, dTrafo, diTrafo);
    errModelPars    = 0;
    
    if ( opts.twocomponent )
        % Two component error model
        errModel        = @(fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks) flexible_model(fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks, @twocomponent);
        errModelPars    = 2;
    end
     
    %% Collect csv files
    fieldNames = {};
    for jN = 1 : length( names )
        data{jN} = readCSV( names{jN}, delimiter ); %#ok
        fieldNames = union( fieldNames, fieldnames(data{jN}) );
    end
    
    % For datasets which have no experiment numbers, we add this column,
    % but show a warning
    for jN = 1 : length( names )
        fields = fieldnames( data{jN} );
        if ~ismember( fields, expVar )
            sprintf( 'Did not find variable which indicates experiment ID for dataset %s. Adding it.\n', names{jN} );
            data{jN}.(expVar) = cell(size(data{jN}.(fields{1}))); %#ok<AGROW>
            data{jN}.(expVar)(:) = {'1'}; %#ok<AGROW>
        end
    end
    fieldNames = union( fieldNames, expVar );
    
    % Match up the headers
    expField = 0;
    for jN = 1 : length( fieldNames )
        % Is this the field indicating the experimental replicate number?
        isExpField = strcmp( fieldNames{jN}, expVar );
        
        out.(fieldNames{jN}) = {};
        for jD = 1 : length( data )
            dataName = data{jD};
            
            
            % Does it have this field?
            if ~isfield( data{jD}, fieldNames{jN} )
                fNames = fieldnames( data{jD} );
                dud = num2cell( NaN( numel( data{jD}.(fNames{1})), 1 ) );
                out.(fieldNames{jN}) = [ out.(fieldNames{jN}); dud ];
            else
                newData = data{jD}.(fieldNames{jN});
                if ( isExpField )
                    newData     = num2cell(cellfun(@(a)plus(a,expField), data{jD}.(fieldNames{jN})));
                    % +1 is added to make sure that we never overlap even if user starts counting 
                    % from 0 or 1 inconsistently in different files
                    if ( ~opts.nofileincrement )
                        if ( isnumeric( cell2mat(data{1}.(expVar)) ) )
                            expField    = expField + max(cell2mat(data{1}.(expVar))) + 2;
                        else
                            expField    = expField + max(str2num(cell2mat(data{1}.(expVar)))) + 2;
                        end
                    end
                end
                out.(fieldNames{jN}) = [ out.(fieldNames{jN}); newData ];
            end
        end
    end
    
    % Merge in the conditions
    if ( opts.splitconditions )
        input_fields = {fieldNames{cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))}, expVar};
        for ji = 1 : numel( input_fields )
            Q(:, ji) = cell2mat(out.(input_fields{ji}));
        end
        uq = unique(Q, 'rows');
        F = zeros(size(Q,1), 1);
        for ju = 1 : size(uq,1)
            F( min(Q==repmat(uq(ju,:), size(Q,1), 1), [], 2) ) = ju;
        end
        out.(expVar) = num2cell(F);
    end   
       
    % Put NaN's in the empty fields
    filt = struct;
    for jN = 1 : length( fieldNames )               
        % Filter out columns with text content
        Q = ~cellfun(@isnumeric, out.(fieldNames{jN}));
        if ( sum(Q) > 0 )
            filt.(fieldNames{jN}) = out.(fieldNames{jN});
            out = rmfield( out, fieldNames{jN} );
            fprintf('Filtered field [%s] because of non numeric values\n', fieldNames{jN});
        else
            % Filter out columns with no content
            Q = cellfun(@isnan, out.(fieldNames{jN}));
            if ( sum(Q) == length(Q) )
                fprintf( 'Filtered field [%s] since it has no content\n', fieldNames{jN} )
                filt.(fieldNames{jN}) = out.(fieldNames{jN});
                out = rmfield( out, fieldNames{jN} );
            end
        end
    end
       
    % Find the time variable
    timeVar = fieldNames( ismember(fieldNames, timeVars) );
    if ( numel( timeVar ) > 1 )
        fprintf( 'Found independent variable column headers: \n');
        fprintf( '%s\n', timeVar{:} );
        error( 'Multiple columns whose header matches the independent variable' );
    end
    if ( numel( timeVar ) == 0 )
        error( 'Did not find independent variable %s', timeVars{:} );
    end
    
    % Filter based on time range
    if ( opts.range )
        lims = opts.range_args;
        if numel(lims)>2
            error( 'Invalid range vector specified. Should be of the form [lb, ub]' );
        else
            f = @(x)(x<lims(1))||(x>lims(2));
            removeMaskTime = find( cellfun(f, out.(timeVar{1})) );
            fprintf( 'Removing %d data points based on specified time range\n', sum(removeMaskTime) );
            
            for a = 1 : length( fieldNames )
                if ( isfield( out, fieldNames{a} ) )
                    out.(fieldNames{a})(removeMaskTime) = [];
                end
            end
        end
    end
    if ( opts.excludeconditions )
        for jec = 1 : numel( opts.excludeconditions_args )
            rejectionFilter = opts.excludeconditions_args{1};
            if ( isfield( out, rejectionFilter{1} ) )
                filterList = cellfun(rejectionFilter{2}, out.(rejectionFilter{1}));
            end
            filNames = fieldnames( out );
            for jf = 1 : numel( filNames )
                out.(filNames{jf})(filterList) = [];
                if ( isempty( find( ~cellfun(@isnan, out.(filNames{jf})) ) ) )
                    warning( 'Filtered all data for %s', filNames{jf} );
                    out = rmfield( out, filNames{jf} );
                end
            end
        end
    end
    
    if ( isempty(fieldnames(out)) )
        error( 'No more data to scale. Did you filter it all?' );
    end   
    
    % Dump out things in the ignore mask
    if ( ~isempty( ignoreMask ) || ~isempty( removeMask ) )
        searchpat = union( ignoreMask, removeMask, 'stable' );
        fprintf('Filtering following fields:\n');
        K = fieldnames(out);
        for a = 1 : length( K )
            if ~iscell( searchpat )
                searchpat = {searchpat};
            end
            for b = 1 : length( searchpat )
                if ( searchpat{b}(1) == '*' )
                    filterActive = ~isempty( strfind( K{a}, filterField(searchpat{b}(2:end)) ) );
                else
                    filterActive = strcmp( K(a), searchpat{b} );
                end
                if filterActive
                    if ( isfield( out, K{a} ) )
                        % Was it an ignore pattern?
                        if ( sum(ismember(searchpat{b}, ignoreMask)) )
                            prt = sprintf('Ignored: %s\n', K{a});
                            ignore.(K{a}) = out.(K{a});
                        else
                            prt = sprintf('Removed: %s\n', K{a});
                        end
                        fprintf( prt );
                        out = rmfield( out, K{a} );
                        fieldNames = setdiff( fieldNames, K{a} );
                    end
                end
            end
        end
    end    
    
    %fieldNames = fieldnames(out);
    % Uncomment to test with control case
    % out = unitTest()
    
    % Estimate correct scaling factors
    %if ( opts.prescaleProblem )
    %    q = fieldnames( out );
    %    maxi = 0;
    %    for jq = 1 : length( q )
    %        maxi = max( [maxi, max( max( cell2mat( out.(q{jq}) ) ) ) ] );
    %    end
    %    for jq = 1 : length( q )
    %        for ji = 1 : length( out.(q{jq}) )
    %            out.(q{jq}){ji} = out.(q{jq}){ji} / maxi;
    %        end
    %    end
    %end
    
    % Set default inputs for unset inputs
    inputs = fieldNames(cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false)));
    for a = 1 : numel( inputs )
        defaultValue = 0;
        nans = cellfun(@isnan, out.(inputs{a}));
        if ( sum(nans) > 0 )
            warning( 'Using default value %g for input %s\n', defaultValue, inputs{a} );
            out.(inputs{a})(nans) = {defaultValue};
        end
    end

    % Set values below the detection limit to DL/2
    if ( opts.dlfudge )
        % Find zeroes
        fNames = fieldnames( out );
        variablesOfInterest = setdiff(fNames(cell2mat(cellfun(@(a)isempty(findstr(a, inputMask)), fNames, 'UniformOutput', false))), timeVar);
        
        for i = 1 : numel( variablesOfInterest )
            predictionBDL = cellfun(@(x)x==0, out.(variablesOfInterest{i}));
            values = out.(variablesOfInterest{i});
            minimal = min( [ values{ ~predictionBDL } ] );           
            out.(variablesOfInterest{i})(predictionBDL) = {minimal / 2};
        end
    end
    
    % Don't scale individual 
    if ( opts.samescale )
        oldExpVar = out.(expVar);
        out.(expVar)(:) = {1};
        [ out, dataFields, fieldNames ] = estimateScaling( errModel, errModelPars, out, expVar, timeVars, inputMask, obsGroups, restrictobs, trafo, invTrafo, opts.rescale, opts.prescale_args );
        out.(expVar) = oldExpVar;
    else
        [ out, dataFields, fieldNames ] = estimateScaling( errModel, errModelPars, out, expVar, timeVars, inputMask, obsGroups, restrictobs, trafo, invTrafo, opts.rescale, opts.prescale_args );
    end
    
    % Append a column?
    if ( opts.appendcolumn )
        out.(opts.appendcolumn_args{1}) = out.(expVar);
        out.(opts.appendcolumn_args{1})(:) = {opts.appendcolumn_args{2}};
    end
       
    % Find unique conditions (unrelated to time this time)
    groupNames = union(fieldNames{ismember(fieldNames, expVar)}, {fieldNames{cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))}});
    %groupNames = {fieldNames{ismember(fieldNames, expVar)}};
    for jC = 1 : length(groupNames)
        groups( jC, : ) = out.(groupNames{jC}); %#ok<AGROW>
    end
    [groupData, ~, groupIDs] = unique( cell2mat( groups ).', 'rows' );
    
    % Plot results
    colors = colmap;    
    nX = floor(sqrt(length(dataFields)));
    nY = ceil( length(dataFields)/ nX );    
    for jD = 1 : length( dataFields )
        subplot(nX, nY, jD);
        for a = 1 : length( unique(groupIDs) )
            ID          = groupIDs == a;
            time        = cell2mat( out.(fieldNames{1})(ID) );
            y           = out.([dataFields{jD} '_mean'])(ID);
            yU          = out.([dataFields{jD} '_lb'])(ID);
            yL          = out.([dataFields{jD} '_ub'])(ID);
            yRep        = out.([dataFields{jD} '_scaled'])(ID);
            [time,I]    = sort(time);
            plot( time, y(I), 'Color', colors(a,:)  ); hold on;
            plot( time, yU(I), '--', 'Color', fadeColor(colors(a,:)), 'LineWidth', 0.5  ); hold on;
            plot( time, yL(I), '--', 'Color', fadeColor(colors(a,:)), 'LineWidth', 0.5  ); hold on;
            plot( time, yRep(I), '.', 'Color', colors(a,:) ); hold on;
        end
        title(dataFields{jD});
    end
    
    if (opts.varanalysis)
        analyseVariances( out, opts.varanalysis_args, trafo );
    end
        
    % Columns that did not take part in the estimation
    stringNames = fieldnames( filt );
    
    % Names of scaled columns
    newNames = fieldnames( out );
    
    % Find the ones that have been rescaled => these must be replaced with
    % their scaled variants
    for a = 1 : length( fieldNames )
        if ismember( [fieldNames{a} '_scaled'], newNames )
            fieldNames{a} = [fieldNames{a} '_scaled'];
        end
    end
    
    % Add the columns that were ignored first back
    if ( exist( 'ignore', 'var' ) )
        startNames = fieldnames( ignore );
        fieldNames = { fieldNames{:} startNames{:} }; %#ok
        for jN = 1 : numel( startNames )
            out.(startNames{jN}) = ignore.(startNames{jN});
        end
    end
    
    if ( opts.firstoutput )
        tColumn = ismember( fieldNames, opts.firstoutput_args{1} );
        fieldNames = { fieldNames{tColumn} fieldNames{~tColumn} };
    end
    
    fid = fopen( outFile, 'w' );   
    fprintf( fid, '%s', fieldNames{1} );
    for jN = 2 : length( fieldNames )
        fprintf( fid, ',%s', fieldNames{jN} );
    end
    for jN = 1 : length( stringNames )
        fprintf( fid, ',%s', stringNames{jN} );
    end
    fprintf( fid, '\n' );
    for jD = 1 : length( out.(fieldNames{1}) )
        fprintf( fid, '%d', out.(fieldNames{1}){jD} );
        for jN = 2 : length( fieldNames )
            outArray = out.(fieldNames{jN});
            if iscell( outArray )
                if ( isnumeric( outArray{jD} ) )
                    fprintf( fid, ',%d', outArray{jD} );
                else
                    fprintf( fid, ',%s', outArray{jD} );
                end
            else
                fprintf( fid, ',%d', outArray(jD) );
            end
        end
        for jN = 1 : length( stringNames )
            fprintf( fid, ',%s', filt.(stringNames{jN}){jD} );
        end
        fprintf( fid, '\n' );
    end
    
    fclose(fid);
end

function analyseVariances(out, desiredObs, trafo)
      
    for jO = 1 : length( desiredObs )
        reps = out.([desiredObs{jO} '_scaled']);
        means = out.([desiredObs{jO} '_mean']);
        lb = out.([desiredObs{jO} '_lb']);
        ub = out.([desiredObs{jO} '_ub']);
    
        uMeans = unique( means );
        sdu = zeros( size( uMeans ) );
        for jM = 1 : length( uMeans )
            sdu(jM) = std(trafo(reps(uMeans(jM)==means)));
        end
        
        figure;
        subplot(1,2,1);
        if ( isfield( out, 'nExpID' ) )
            [cs, ~, c] = unique(cell2mat(out.nExpID));
            cols = 'rkbmgyp';
            for a = 1 : length( cs )
                cid = (a-1)-length(cols)*floor(a/length(cols)) + 1;
                plot( trafo(means(c==a)), trafo(reps(c==a))-trafo(means(c==a)), [cols( cid ), '.'] ); hold on;
            end
        else
            plot( trafo(means), trafo(reps)-trafo(means), '.' ); hold on;
        end
        ylabel( [ 'Residual of ', escapeChars(desiredObs{jO}) ] );
        xlabel( [ 'Estimated mean ' escapeChars(desiredObs{jO}) ] );        
        
        subplot(1,2,2);
        plot( means, reps, '.' ); hold on;
        ylabel( [ 'Replicate ' escapeChars( desiredObs{jO} ) ] );
        xlabel( [ 'Estimated mean ' escapeChars(desiredObs{jO}) ] );
    end
end

function str = escapeChars(str)
    str = strrep(str, '_', '\_');
end

function colors = colmap()
    colors = [ ...
            27, 158, 119; ...
            217, 95, 2; ...
            117, 112, 179; ...
            231, 41, 138; ...
            102, 166, 30; ...
            230, 171, 2; ...
            166, 118, 29; ...
            102, 102, 102 ];
        
    colors = repmat( colors, 30, 1 ) / max(max(colors));
end

function col = fadeColor( col )
    col = col * 0.3 + 0.7;
end

function out = unitTest() %#ok
    out.time    = {1,2,3,  1,2,3,  1,2,3        1,2,3,  1,2,3,  1,2,3}.';
    out.nExpID  = {1,1,1,  2,2,2,  3,3,3        1,1,1,  2,2,2,  3,3,3}.';
    out.input_L = {1,1,1,  1,1,1,  1,1,1,       2,2,2,  2,2,2,  3,3,3}.'; 
    out.y       = {1,2,3,  2,4,6,  10,20,30     3,6,9,  6,12,18 30,60,90}.';
end

function [ out, dataFields, fieldNames ] = estimateScaling( errModel, errModelPars, out, expVar, timeVars, inputMask, obsGroups, obsList, trafo, invTrafo, rescale, prescale )
    verbose = 1;
    fieldNames = fieldnames(out);
    
    % Sort by fill
    fill = zeros(1, length(fieldNames));
    for jN = 1 : length( fieldNames )
        fill(jN) = sum(~cellfun(@isnan, out.(fieldNames{jN})));
    end
    fill(ismember(fieldNames, expVar))=1e30;
    timeVar = false(size(fieldNames));
    for a = 1 : length( timeVars )
        timeVar = timeVar | ismember(fieldNames, timeVars{a});
    end
    fill(timeVar)=inf;
        
    [~, fill] = sort( fill, 'descend' );
    fieldNames = fieldNames(fill);
    
    % Assemble condition matrix
    clear data;
    for jD = 1 : length(fieldNames)
        data( jD, : ) = out.(fieldNames{jD}); %#ok
    end
        
    % Find unique conditions (note that a different time point (field 1) is
    % considered a different condition in this setting)
    conditionFields = union( fieldNames{1}, fieldNames(cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))) );
    for jC = 1 : length(conditionFields)
        conds( jC, : ) = out.(conditionFields{jC}); %#ok<AGROW>
    end
    [~, ~, jcondition] = unique( cell2mat( conds ).', 'rows' );
    
    for a = unique(jcondition).'
        if ( sum(jcondition==a) < 2 )
            warnLine = '';
            ID = find(jcondition==a,1);
            for jF = 1 :length( conditionFields )
                warnLine = sprintf( '%s\n%s: %d', warnLine, conditionFields{jF}, out.(conditionFields{jF}){ID} ); %#ok
            end
            warning( 'Found condition with very few data points:\n%s\n', warnLine );
            %pause;
        end
    end    
    
    % Determine data fields
    dataFields  = setdiff( fieldNames, union( conditionFields, expVar ) );
        
    if ~isempty( obsList )
        dataFields = intersect( dataFields, obsList );
    end
        
    % Assign data fields to common groups
    if ( isempty( obsGroups ) )
        disp( 'No groups specified.' );
        for jG = 1 : length( dataFields )
            obsGroups{jG} = jG;
        end
    else
        remainingDatasets = dataFields;
        for jG = 1 : length( obsGroups )
            currentGroups = [];
            for jD = 1 : length( obsGroups{jG} )
                if ( ~ismember( obsGroups{jG}{jD}, remainingDatasets ) > 0 )
                    error( sprintf( 'Observable %s does not exist or is already in a group', obsGroups{jG}{jD} ) ); %#ok
                end
                currentGroup = find( ismember( dataFields, obsGroups{jG}{jD} ) );
                if ~isempty( currentGroup )
                    currentGroups = [ currentGroups, currentGroup ]; %#ok
                else
                    warning( 'Observable %s in grouping not found in datasets. Did you misspell it?', obsGroups{jG} );
                end
            end
            obsGroups{jG} = currentGroups;
            remainingDatasets = setdiff( remainingDatasets, dataFields(obsGroups{jG}) );
        end
        if (~isempty( remainingDatasets ) )
            for jD = 1 : length( remainingDatasets )
                obsGroups{end+1} = find( ismember( dataFields, remainingDatasets{jD} ) ); %#ok
                if ( verbose )
                    fprintf( '%s was not in a group, adding as separate group.\n', remainingDatasets{jD} );
                end
            end
        end
    end
    nGroups = length( obsGroups );
    
    expVar = fieldNames(cell2mat(cellfun(@(a)strcmp(a, expVar), fieldNames, 'UniformOutput', false)));
    nObs  = length( dataFields );

    % Set up a unique scaling factor for each group
    scales = cell(1, nObs);
    current = 0;
    for jG = 1 : nGroups
        for jD = obsGroups{jG}
            scales{jD} = cell2mat(out.(expVar{1})) + current;
        end
        current = max([current; scales{jD}]);
    end
    
    Nrows = numel(out.(dataFields{1})); 
    for jG = 1 : nGroups
        curScale = []; curData = []; curConditions = []; obsjD = [];
        for jD = obsGroups{jG}
            fprintf('Processing obs %s\n', dataFields{jD} );
            nData  = cell2mat(out.(dataFields{jD}));                % Process single observation
            if ~isempty( prescale )
                nData  = prescale( nData );                         % Apply forced prescaling function
            end
            Q{jD}           = find(~isnan( nData ));                %#ok Fetch points that have data
            curData         = [curData; nData(Q{jD})];              %#ok
            curConditions   = [curConditions, jcondition(Q{jD}).']; %#ok Fetch condition IDs this observable appears in
            obsjD           = [obsjD; jD*ones(size(Q{jD}))];        %#ok Observable => This is paired with the conditions of observables in the same group, to ensure they get separate means in the estimation process
            curScale        = [curScale; scales{jD}(Q{jD})];        %#ok
        end
        
        % ConditionLinks links the replicates to the "true" IDs (i.e. true conditions, t, cond pairs)
        % ScaleLinks refers to which scale corresponds to which datapoint.
        % Each "condition" needs a 'true' mean
        % Each "gel/observable group combo" needs its own scale
        [scaleTargets, ~, scaleLinks] = unique(curScale);                                     % ScaleTargets provides a link to the global scales
        [conditionTargets, ~, conditionLinks] = unique([curConditions; obsjD.'].', 'rows');   % ConditionTargets provides a link to global conditions. ConditionTargets(conditionLinks) equates to condition IDs
        conditionTargets = conditionTargets(:,1).';
        
        % The result of the estimation will be means (in the order
        % specified in conditionTargets) and scaling (in the order
        % specified in scaleTargets).
        checkDims( data, conds, scaleTargets, conditionTargets );
        [means, mlb, mub, scalings, offsets] = estimateObs( errModel, errModelPars, curData, scaleTargets, scaleLinks, conditionTargets, conditionLinks, rescale, trafo ); 
        replicates = trafo( (curData - offsets) ./ scalings(scaleLinks) );
        
        % To get where these are globally we have to map them back
        % means(conditionLinks) and scalings(scaleLinks) will sort them
        % back to the order they were in when they entered this function.
        loc = 0;
        for jD = obsGroups{jG}
            % Number of unique points in this dataset
            nData = length( Q{jD} );
            locs = 1 + loc : nData + loc;
            
            % Returned quantities are indexed by condition ID
            out.([dataFields{jD} '_mean'])          = NaN( Nrows, 1 );
            out.([dataFields{jD} '_mean'])(Q{jD})   = invTrafo( means(conditionLinks(locs)) );

            % Returned quantities are indexed by condition ID
            out.([dataFields{jD} '_lb'])            = NaN( Nrows, 1 );
            out.([dataFields{jD} '_lb'])(Q{jD})     = invTrafo( mlb(conditionLinks(locs)) );          

            % Returned quantities are indexed by condition ID
            out.([dataFields{jD} '_ub'])            = NaN( Nrows, 1 );
            out.([dataFields{jD} '_ub'])(Q{jD})     = invTrafo( mub(conditionLinks(locs)) );

            out.([dataFields{jD} '_scale'])         = NaN( Nrows, 1 );
            out.([dataFields{jD} '_scale'])(Q{jD})  = scalings(scaleLinks(locs));

            out.([dataFields{jD} '_scaled'])        = NaN( Nrows, 1 );
            out.([dataFields{jD} '_scaled'])(Q{jD}) = invTrafo( replicates(locs) );
            
            loc = loc + nData;
        end
    end
    
    if (verbose)
        fprintf( '\n\nEstimated scalings for following observables:\n');
        for jG = 1 : length( obsGroups )
            fprintf( 'Scale group %d: ', jG )
            for jD = 1 : length( obsGroups{jG} )
                fprintf( '%s ', dataFields{obsGroups{jG}(jD)} );
            end
            fprintf( '\n' );
        end
    end
end

function checkDims( data, condition, scaleTargets, conditionTargets )
    nPars = length( scaleTargets ) + length( conditionTargets );
    nRes  = length(data);

    if ( nPars > nRes )
        warning( 'More parameters than equations for' );
        condition(:,conditionTargets)
    end
end

function [means, mlb, mub, scalings, offsets] = estimateObs( eModel, errModelPars, data, scaleTargets, scaleLinks, conditionTargets, conditionLinks, rescale, trafo )

    % Parameters are [ means_for_each_condition , scalings ]
    Nmeans = size(conditionTargets,2);
    means  = zeros(1,Nmeans);
    
    % Obtain initial means from data
    for a = 1 : length( conditionTargets )
        loc = find(conditionLinks==a);
        means(a) = trafo( data(loc(1)) );
    end

    Nscalings = numel(scaleTargets)-1;
    Noffsets = 0;
    offsetLinks = ones( size( conditionLinks ) );
    offsets = 1e-9*ones(Noffsets, 1);
    
    lb = [ -1e12 * ones(Nmeans,1); zeros(Nscalings,1); 1e-15*ones(Noffsets,1); 1e-15*ones(errModelPars,1) ];
    if ( rescale )
        fixedScale = nanmax(nanmax(data));
        fprintf( ' <Rescale: %.5g>\n', fixedScale );
    else
        fixedScale = 1;
    end
    
    % Set initial scalings   
    scalings = fixedScale * ones(Nscalings,1);
    initPar  = [ means.'; scalings; offsets; 1e-8 * ones( errModelPars, 1) ];
    
    fprintf( 'Number of free scaling factors %d\n', Nscalings );    
    
    options    = optimset('TolFun', 1e-7, 'TolX', 1e-6, 'MaxIter', 1e4, 'MaxFunevals', 1e5, 'Display', 'Off', 'Jacobian', 'On', 'DerivativeCheck', 'Off' ); % 
    tol        = 1e-7;
    errModel   = @(pars) eModel(fixedScale, pars, data, scaleLinks, length(conditionTargets), length(scalings), Noffsets, conditionLinks, offsetLinks );
    
    if exist('debug', 'var')
        initPar = randn(size(initPar));
        [~,J]=errModel(initPar);
        Jfin=findiff(errModel, initPar);
    end
    
    [p, rn, r, ~, ~, ~, J]  = lsqnonlin( errModel, initPar, lb, 1e125*ones(size(initPar)), options );
    rLast = inf;
    round = 1;
    while( ((rLast - rn)/rn > tol) )
        rLast = rn;
        [p, rn, r, ~, ~, ~, J]   = lsqnonlin( errModel, p, lb, 1e125*ones(size(initPar)), options );
        fprintf('Round %d: Difference last and current optimization: %g\n', round, rLast-rn);
        round = round + 1;
    end
    ci = nlparci(p,r,'Jacobian',J,'alpha',0.05);
    
    fprintf( 'Final objective: %g\n', sum(r.^2) );
    
    means    = p(1:numel(means));
    mlb      = ci(1:numel(means),1);
    mub      = ci(1:numel(means),2);
    scalings = [fixedScale; p(numel(conditionTargets)+1:end - errModelPars - Noffsets)];
    
    if ( Noffsets > 0 )
        offsets  = p(end-errModelPars-Noffsets+1:end-errModelPars);
    else
        offsets  = 0;
    end
end

function [res, J] = model( fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks, trafo, invTrafo, dTrafo, diTrafo )
    scales = [fixedScale; pars(Ncond+1:Ncond+Nscale)];
    if ( Noffsets > 0 )
        res = trafo( scales(scaleLinks) .* invTrafo( pars(conditionLinks) ) + pars(offsetLinks+Ncond+Nscale) ) - trafo( data );  
    else
        res = trafo( scales(scaleLinks) .* invTrafo( pars(conditionLinks) ) ) - trafo( data );  
    end
    
    % Jacobian if desired
    if ( nargout > 1 )
        J = zeros(numel(pars), numel(res));
        
        if ( Noffsets > 0 )
            dTraf = dTrafo( scales(scaleLinks) .* invTrafo( pars(conditionLinks) ) + pars(offsetLinks+Ncond+Nscale) );
        else
            dTraf = dTrafo( scales(scaleLinks) .* invTrafo( pars(conditionLinks) ) );
        end
        % Mean derivatives
        for jC = 1 : Ncond
            J(jC, 1:numel(scaleLinks))                  = (conditionLinks==jC) .* dTraf .* ( diTrafo(pars(conditionLinks)) .* scales(scaleLinks) );
        end
        % Scale derivatives
        for jS = 2 : Nscale+1
            J(jS + Ncond - 1, 1:numel(scaleLinks))      = (scaleLinks==jS) .* dTraf .* invTrafo( pars(conditionLinks) );
        end
        % Offset derivatives
        for jO = 1 : Noffsets
            J(jO + Ncond + Nscale, 1:numel(scaleLinks)) = (offsetLinks==jO) .* dTraf;
        end
        J = J.';
    end
end

% Does not deal with trafos yet
function [res, J] = flexible_model( fixedScale, pars, data, scaleLinks, Ncond, Nscale, Noffsets, conditionLinks, offsetLinks, errorModel )
    fix     = 50;
    scales  = [fixedScale; pars(Ncond+1:Ncond+Nscale)];
    %sigma   = errorModel(pars(Ncond+Nscale+1:end), trafo( pars(conditionLinks) ), trafo( scales(scaleLinks) ));
    [sigma, dsigds, dsigdu, dsigdsig] = errorModel(pars(Ncond+Nscale+1:end), pars(conditionLinks), scales(scaleLinks) );
    if ( max( (2 * log( sigma ) + fix) < 1 ) > 0 ) 
        error( 'Sigma term negative' );
    end

    res = [ ( scales(scaleLinks) .* pars(conditionLinks) - data ) ./ sigma; sqrt( 2 * log( sigma ) + fix ) ];
    
    % Jacobian if desired
    if ( nargout > 1 )
        J = zeros(numel(pars), numel(res));
        
        for jC = 1 : Ncond
            J(jC, 1:numel(scaleLinks))              = (conditionLinks==jC) .* (scales(scaleLinks) ./ sigma + res(1:numel(scaleLinks)) .* (-1./sigma) .* dsigdu);
            J(jC, numel(scaleLinks)+1:end)          = (conditionLinks==jC) .* (1./(sigma .* sqrt(2 * log( sigma ) + fix))) .* dsigdu;
        end
        for jS = 2 : Nscale+1
            J(jS + Ncond - 1, 1:numel(scaleLinks))      = (scaleLinks==jS) .* ( pars(conditionLinks) ./ sigma + res(1:numel(scaleLinks)) .* (-1./sigma) .* dsigds );
            J(jS + Ncond - 1, numel(scaleLinks)+1:end)  = (scaleLinks==jS) .* (1./(sigma .* sqrt(2 * log( sigma ) + fix))) .* dsigds;
        end      

        for l = 1 : size( dsigdsig, 2 )
            pIDs = Ncond+Nscale+l;
            J(pIDs, 1:numel(scaleLinks))     = res(1:numel(scaleLinks)) .* (-1./sigma) .* dsigdsig(:,l);
            J(pIDs, numel(scaleLinks)+1:end) = (1./(sigma .* sqrt(2 * log( sigma ) + fix))) .* dsigdsig(:,l);
        end
        
        J = J.';
    end
end

function [sigma, dsigds, dsigdu, dsigdsig] = twocomponent( noisepars, mus, ~ ) % Last arg scales is not needed
    sigma    = sqrt( noisepars(1)*noisepars(1) + noisepars(2)*noisepars(2)*mus.*mus );
    dsigds   = zeros( numel(sigma), 1 );
    dsigdu   = mus .* noisepars(2)*noisepars(2) ./ sigma;
    
    dsigdsig = [ones(numel(mus),1)*noisepars(1)./sigma, noisepars(2).*mus.*mus./sigma];
end

function data = readCSV( filename, delimiter )
    
    if isempty( strfind( filename, '.' ) )
        filename = [filename '.csv'];
    end
    
    try
        fid = fopen(filename, 'r');
    catch
        error( 'Failed to open %s', filename );      
    end
    
    % Fetch header
    try
        C = textscan(fid, '%s\n',1,'Delimiter','');
    catch
        error( 'Couldn''t find file %s', filename );
    end
    
    % Grab header items
    headers = textscan(C{1}{1}, '%q', 'Delimiter', delimiter);
    headers = headers{1};
    
    % Fetch the rest
    i = 1;
    C = textscan(fid,'%s\n',1,'Delimiter','');
    while(~isempty(C{1}))
        C = textscan(C{1}{1},'%q','Delimiter',delimiter); C = C{1};

        %dataBlock{i, :} = cell(1,size(headers, 1)); %#ok
        for j=1:length(headers)
            if(j>length(C))
                dataBlock(i, j) = {NaN}; %#ok
            else
                if isempty(C{j})
                    dataBlock(i, j) = {NaN}; %#ok
                else
                    val = str2num(C{j}); %#ok
                    if ~isempty( val )
                        dataBlock(i, j) = {val}; %#ok
                    else
                        dataBlock(i, j) = C(j); %#ok
                    end
                end
            end
        end
        C = textscan(fid,'%s\n',1,'Delimiter','');
        i = i + 1;
    end
    
    for a = 1 : length( headers )
        if (~isempty(headers{a}))
            data.(filterField(headers{a})) = dataBlock(:,a);
        end
    end
    fclose(fid);
end

function str = filterField( str )
    str = strrep( str, '+', '_plus_' );
    str = strrep( str, ' ', '_' );
    str = strrep( str, '-', '_minus_' );
end
