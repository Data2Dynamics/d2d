%
% Work in progress function for merging CSV files
%
% Scales data based on minimization of linear model
%
% Usage:
%   mergeCSV( names, outFileName, delimiter )
%
%   CSV file must have following columns
%     nExpID  = replicate
%     input_  = inputs per condition
%     time    = dependent variable must be called t, T, time or Time
%
% All other columns are assumed data. Scaling factors will be estimated per
% observable/replicate combination.

function mergeCSV( names, outFile, delimiter )
    if ( nargin < 3 )
        delimiter = ',';
    end
    
    timeVars = {'t', 'T', 'time', 'Time'};
    expVar = 'nExpID';
    inputMask = 'input_';

    % Collect csv files
    fieldNames = {};
    for jN = 1 : length( names )
        data{jN} = readCSV( names{jN}, delimiter ); %#ok
        fieldNames = union( fieldNames, fieldnames(data{jN}) );
    end
    
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
                    expField    = expField + max(cell2mat(data{jD}.(fieldNames{jN}))) + 1;
                end
                out.(fieldNames{jN}) = [ out.(fieldNames{jN}); newData ];
            end
        end
    end
    
    % Put NaN's in the empty fields
    for jN = 1 : length( fieldNames )               
        % Filter out columns with no content
        Q = cellfun(@isnan, out.(fieldNames{jN}));
        if ( sum(Q) == length(Q) )
            out = rmfield( out, fieldNames{jN} );
        end
    end
    
    % Uncomment to test with control case
    % out = unitTest()
    
    % Estimate correct scaling factors
    [ out, dataFields, fieldNames ] = estimateScaling( out, expVar, timeVars, inputMask );
       
    % Find unique conditions (unrelated to time this time)
    groupNames = union(fieldNames{ismember(fieldNames, expVar)}, {fieldNames{cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))}});
    for ( jC = 1 : length(groupNames) )
        groups( jC, : ) = out.(groupNames{jC});
    end
    [~, ~, groupIDs] = unique( cell2mat( groups ).', 'rows' );
    
    % Plot results
    colors = 'rgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmc';
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
            plot( time, y(I), colors(a)  ); hold on;
            plot( time, yU(I), [colors(a), '--'], 'LineWidth', 0.5  ); hold on;
            plot( time, yL(I), [colors(a), '--'], 'LineWidth', 0.5  ); hold on;
            plot( time, yRep(I), [colors(a) '.'] ); hold on;
        end
        title(dataFields{jD});
    end
        
    newNames = fieldnames(out);
    
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
        fprintf( fid, '\n' );
    end
    
    fclose(fid);
end

function out = unitTest()
    out.time    = {1,2,3,  1,2,3,  1,2,3        1,2,3,  1,2,3,  1,2,3}.';
    out.nExpID  = {1,1,1,  2,2,2,  3,3,3        1,1,1,  2,2,2,  3,3,3}.';
    out.input_L = {1,1,1,  1,1,1,  1,1,1,       2,2,2,  2,2,2,  3,3,3}.'; 
    out.y       = {1,2,3,  2,4,6,  10,20,30     3,6,9,  6,12,18 30,60,90}.';
end

function [ out, dataFields, fieldNames ] = estimateScaling( out, expVar, timeVars, inputMask )
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
    timeVar = find(timeVar);
    fill(timeVar)=inf;
        
    [~, fill] = sort( fill, 'descend' );
    fieldNames = fieldNames(fill);
    
    % Assemble condition matrix
    clear data;
    for ( jD = 1 : length(fieldNames) )
        data( jD, : ) = out.(fieldNames{jD});
    end
    
    % Find unique conditions (note that a different time point (field 1) is
    % considered a different condition in this setting)
    conditionFields = union( fieldNames{1}, {fieldNames{cell2mat(cellfun(@(a)~isempty(findstr(a, inputMask)), fieldNames, 'UniformOutput', false))}} );
    for ( jC = 1 : length(conditionFields) )
        conds( jC, : ) = out.(conditionFields{jC});
    end
    [~, ~, jcondition] = unique( cell2mat( conds ).', 'rows' );
    
    expVar = {fieldNames{cell2mat(cellfun(@(a)strcmp(a, expVar), fieldNames, 'UniformOutput', false))}};
    if (isempty(expVar))
        error( 'Did not find variable which indicates experiment ID' );
    end
    scales = cell2mat(out.(expVar{1}));
    
    dataFields  = setdiff( fieldNames, union( conditionFields, expVar ) );
    Nrows = numel(out.(dataFields{1})); 
    for jD = 1 : length(dataFields)
        curData         = cell2mat(out.(dataFields{jD}));     % Process single observation
        Q               = find(~isnan( curData ));            % Fetch points that have data
        curData         = curData(Q);
        curConditions   = jcondition(Q).';                    % Fetch condition IDs this observable appears in
        curScale        = scales(Q);
        
        % ConditionLinks links the replicates to the "true" IDs (i.e. true conditions, t, cond pairs)
        % ScaleLinks refers to which scale corresponds to which datapoint.
        % Each "condition" needs a 'true' mean
        % Each "gel/observable combo" needs its own scale
        [scaleTargets, ~, scaleLinks] = unique(curScale);                   % ScaleTargets provides a link to the global scales
        [conditionTargets, ~, conditionLinks] = unique(curConditions);      % ConditionTargets provides a link to global conditions. ConditionTargets(conditionLinks) equates to condition IDs
        
        % The result of the estimation will be means (in the order
        % specified in conditionTargets) and scaling (in the order
        % specified in scaleTargets).
        checkDims( data, conds, scaleTargets, conditionTargets );
        [means, mlb, mub, scalings] = estimateObs( curData, scaleTargets, scaleLinks, conditionTargets, conditionLinks ); 
        replicates = curData .* scalings(scaleLinks);
        
        % To get where these are globally we have to map them back
        % means(conditionLinks) and scalings(scaleLinks) will sort them
        % back to the order they were in when they entered this function.

        % Returned quantities are indexed by condition ID
        out.([dataFields{jD} '_mean'])      = NaN( Nrows, 1 );
        out.([dataFields{jD} '_mean'])(Q)   = means(conditionLinks);
        
        % Returned quantities are indexed by condition ID
        out.([dataFields{jD} '_lb'])        = NaN( Nrows, 1 );
        out.([dataFields{jD} '_lb'])(Q)     = mlb(conditionLinks);          
        
        % Returned quantities are indexed by condition ID
        out.([dataFields{jD} '_ub'])        = NaN( Nrows, 1 );
        out.([dataFields{jD} '_ub'])(Q)     = mub(conditionLinks);        

        out.([dataFields{jD} '_scale'])     = NaN( Nrows, 1 );
        out.([dataFields{jD} '_scale'])(Q)  = scalings(scaleLinks);

        out.([dataFields{jD} '_scaled'])    = NaN( Nrows, 1 );
        out.([dataFields{jD} '_scaled'])(Q) = replicates;
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

function [means, mlb, mub, scalings] = estimateObs( data, scaleTargets, scaleLinks, conditionTargets, conditionLinks )

    % Parameters are [ means_for_each_condition , scalings ]
    means = zeros(size(conditionTargets));
    
    % Obtain initial means from data
    for a = 1 : length( conditionTargets )
        loc = find(conditionLinks==a);
        means(a) = data(loc(1));
    end
    
    % Set initial scalings
    scalings = ones(numel(scaleTargets)-1,1);
    initPar = [ means.'; scalings ];
    
    options = optimset('TolFun', 1e-12, 'TolX', 1e-9, 'MaxIter', 1e4, 'MaxFunevals', 1e5, 'Display', 'Iter' );
    [p, ~, r, ~, ~, ~, J] = lsqnonlin( @(pars)model(pars, data, scaleLinks, length(conditionTargets), conditionLinks), initPar, zeros(size(initPar)), 1e12*ones(size(initPar)), options );
    ci = nlparci(p,r,'Jacobian',J);
    
    means = p(1:numel(means));
    mlb = ci(:,1);
    mub = ci(:,2);
    scalings = [1; p(numel(conditionTargets)+1:end)];
    
end

function res = model( pars, data, scaleLinks, Ncond, conditionLinks )
    scales = [1; pars(Ncond+1:end)];

    res = pars(conditionLinks) - scales(scaleLinks) .* data;
end

function data = readCSV( filename, delimiter )
    
    if isempty( strfind( filename, '.' ) )
        filename = [filename '.csv'];
    end
    
    fid = fopen(filename, 'r');
    
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
                dataBlock(i, j) = {NaN};
            else
                val = str2num(C{j});
                if ~isempty( val )
                    dataBlock(i, j) = {val};
                else
                    dataBlock(i, j) = {NaN};
                end
            end
        end
        C = textscan(fid,'%s\n',1,'Delimiter','');
        i = i + 1;
    end
    
    for a = 1 : length( headers )
        data.(filterField(headers{a})) = dataBlock(:,a);
    end
    fclose(fid);
end

function str = filterField( str )
    str = strrep( str, '+', '_plus_' );
    str = strrep( str, ' ', '_' );
    str = strrep( str, '-', '_minus_' );
end