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
%     ignoreMask         Add mask for which data columns to ignore
%     twocomponent       Use two component error model (warning: poorly
%                        tested so far)
%     logtrafo           Do the scaling in logarithmic space
%
% To do: Log trafo scaling, offsets and two component error model scaling

function scaleIt( names, outFile, varargin )

    if ( nargin < 2 )
        help scaleIt;
        error( 'Insufficient arguments' );
    end

    % Load options
    verbose = 0;
    switches = { 'delimiter', 'obsgroups', 'inputmask', 'depvar', 'expid', 'restrictobs', 'ignoremask', 'twocomponent', 'logtrafo' };
    extraArgs = [ 1, 1, 1, 1, 1, 1, 1, 0, 0 ];
    description = { ...
    {'', 'Custom delimiter specified'} ...
    {'', 'Using custom observation/scaling factor pairing'} ...
    {'', 'Using input mask'} ...
    {'', 'Using custom dependent variable name'} ...
    {'', 'Using custom experiment ID column'} ...
    {'', 'Restricting to specific specified observables'} ...
    {'', 'Using ignore mask'} ...,
    {'', 'Using two component error model'} ...
    {'', 'Using lognormal error model'} };
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
    ignoreMask = [];
    if ( opts.ignoremask )
        ignoreMask = opts.ignoremask_args;
    end
    
    % Default mapping
    trafo           = @(x) x;
    invTrafo        = @(x) x;
    
    if ( opts.logtrafo )
        trafo       = @(x) log10(x);
        invTrafo    = @(x) 10.^(x);
    end
    
    % Simplest error model
    errModel        = @(pars, data, scaleLinks, Ncond, Nscale, conditionLinks) model(pars, data, scaleLinks, Ncond, Nscale, conditionLinks, trafo);
    errModelPars    = 0;
    
    if ( opts.twocomponent )
        % Two component error model
        errModel        = @(pars, data, scaleLinks, Ncond, Nscale, conditionLinks) flexible_model(pars, data, scaleLinks, Ncond, Nscale, conditionLinks, @twocomponent, trafo);
        errModelPars    = 2;
    end
     
    % Collect csv files
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
                    expField    = expField + max(cell2mat(data{jD}.(fieldNames{jN}))); % + 1;
                end
                out.(fieldNames{jN}) = [ out.(fieldNames{jN}); newData ];
            end
        end
    end

    % Dump out things in the ignore mask
    if ~isempty( ignoreMask )
        K = fieldnames(out);
        for a = 1 : length( K )
            if ~iscell( ignoreMask )
                ignoreMask = {ignoreMask};
            end
            for b = 1 : length( ignoreMask )
                if ~isempty( strfind( K{a}, ignoreMask{b} ) )
                    if ( isfield( out, K{a} ) )
                        out = rmfield( out, K{a} );
                        fieldNames = setdiff( fieldNames, K{a} );
                    end
                end
            end
        end
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
    %fieldNames = fieldnames(out);
    % Uncomment to test with control case
    % out = unitTest()
    
    % Estimate correct scaling factors
    [ out, dataFields, fieldNames ] = estimateScaling( errModel, errModelPars, out, expVar, timeVars, inputMask, obsGroups, restrictobs, trafo, invTrafo );
       
    % Find unique conditions (unrelated to time this time)
    groupNames = union(fieldNames{ismember(fieldNames, expVar)}, {fieldNames{cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))}});
    %groupNames = {fieldNames{ismember(fieldNames, expVar)}};
    for jC = 1 : length(groupNames)
        groups( jC, : ) = out.(groupNames{jC}); %#ok<AGROW>
    end
    [~, ~, groupIDs] = unique( cell2mat( groups ).', 'rows' );
    
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
    
    fid = fopen( outFile, 'w' );
    fprintf( fid, '%s ', fieldNames{1} );
    for jN = 2 : length( fieldNames )
        fprintf( fid, ', %s', fieldNames{jN} );
    end
    for jN = 1 : length( stringNames )
        fprintf( fid, ', %s', stringNames{jN} );
    end
    fprintf( fid, '\n' );
    for jD = 1 : length( out.(fieldNames{1}) )
        fprintf( fid, '%d ', out.(fieldNames{1}){jD} );
        for jN = 2 : length( fieldNames )
            outArray = out.(fieldNames{jN});
            if iscell( outArray )
                fprintf( fid, ', %d', outArray{jD} );
            else
                fprintf( fid, ', %d', outArray(jD) );
            end
        end
        for jN = 1 : length( stringNames )
            fprintf( fid, ', %s', filt.(stringNames{jN}){jD} );
        end
        fprintf( fid, '\n' );
    end
    
    fclose(fid);
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

function [ out, dataFields, fieldNames ] = estimateScaling( errModel, errModelPars, out, expVar, timeVars, inputMask, obsGroups, obsList, trafo, invTrafo )
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
        timeVar = timeVar |  ismember(fieldNames, timeVars{a});
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
            ID = find( sum(jcondition==a) );
            for jF = 1 :length( conditionFields )
                warnLine = sprintf( '%s\n%s: %d', warnLine, conditionFields{jF}, out.(conditionFields{jF}){ID} ); %#ok
            end
            warning( 'Found condition with very few data points:\n%s\n', warnLine );
            pause;
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
            fprintf('Processing obs group %s', dataFields{jD});
            nData           = cell2mat(out.(dataFields{jD}));        %    Process single observation
            Q{jD}           = find(~isnan( nData ));                 %#ok Fetch points that have data
            curData         = [curData; nData(Q{jD})];               %#ok
            curConditions   = [curConditions, jcondition(Q{jD}).'];  %#ok Fetch condition IDs this observable appears in
            obsjD           = [obsjD; jD*ones(size(Q{jD}))];         %#ok Observable => This is paired with the conditions of observables in the same group, to ensure they get separate means in the estimation process
            curScale        = [curScale; scales{jD}(Q{jD})];         %#ok
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
        [means, mlb, mub, scalings] = estimateObs( errModel, errModelPars, curData, scaleTargets, scaleLinks, conditionTargets, conditionLinks ); 
        replicates = trafo(curData ./ scalings(scaleLinks));
        
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

function [means, mlb, mub, scalings] = estimateObs( eModel, errModelPars, data, scaleTargets, scaleLinks, conditionTargets, conditionLinks )

    % Parameters are [ means_for_each_condition , scalings ]
    means = zeros(size(conditionTargets));
    
    % Obtain initial means from data
    for a = 1 : length( conditionTargets )
        loc = find(conditionLinks==a);
        means(a) = data(loc(1));
    end
    
    % Set initial scalings
    scalings = ones(numel(scaleTargets)-1,1);
    initPar  = [ means.'; scalings; 1e10 * ones( errModelPars, 1) ];
    
    lb = zeros(size(initPar));
    lb( end - errModelPars : end ) = 1e-9;
    
    options                 = optimset('TolFun', 0, 'TolX', 1e-9, 'MaxIter', 1e4, 'MaxFunevals', 1e5, 'Display', 'Iter' );
    errModel                = @(pars) eModel(pars, data, scaleLinks, length(conditionTargets), length(scalings), conditionLinks);
    p                       = lsqnonlin( errModel, initPar, lb, 1e25*ones(size(initPar)), options );
    [p, ~, r, ~, ~, ~, J]   = lsqnonlin( errModel, p, lb, 1e25*ones(size(initPar)), options );
    ci                      = nlparci(p,r,'Jacobian',J,'alpha',0.05);
    
    fprintf( 'Final objective: %g\n', sum(r.^2) );
    
    means    = p(1:numel(means));
    mlb      = ci(:,1);
    mub      = ci(:,2);
    scalings = [1; p(numel(conditionTargets)+1:end - errModelPars)];
end

function res = model( pars, data, scaleLinks, Ncond, Nscale, conditionLinks, trafo )
    scales = [1; pars(Ncond+1:Ncond+Nscale)];
    res = trafo( scales(scaleLinks) .* pars(conditionLinks) ) - trafo( data );
end

function res = flexible_model( pars, data, scaleLinks, Ncond, Nscale, conditionLinks, errorModel, trafo )
    fix     = 1e10;
    scales  = [1; pars(Ncond+1:Ncond+Nscale)];
    sigma   = errorModel(pars(Ncond+Nscale+1:end), trafo( pars(conditionLinks) ), trafo( scales(scaleLinks) ));
    
    if ( max( sigma*fix < 1 ) > 0 ) 
        error( 'Sigma term negative :(' );
    end
    
    res = [ ( scales(scaleLinks) .* pars(conditionLinks) - data ) / (sigma / fix), sqrt( 2 * log( sigma * fix ) ) ];
end

function sigma = twocomponent( noisepars, mus, scales )
    sigma = sqrt( noisepars(1)*noisepars(1) + noisepars(2)*noisepars(2)*mus.*mus );
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
    C = textscan(fid, '%s\n',1,'Delimiter','');
    
    
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