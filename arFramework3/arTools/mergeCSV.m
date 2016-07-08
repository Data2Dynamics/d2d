function mergeCSV( names, outFile, delimiter )
    if ( nargin < 3 )
        delimiter = ',';
    end
    
    expVar = 'nExpID';

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
    
    fieldNames = fieldnames(out);
    
    % Sort by fill
    fill = zeros(1, length(fieldNames));
    for jN = 1 : length( fieldNames )
        fill(jN) = sum(~cellfun(@isnan, out.(fieldNames{jN})));
    end
    fill(ismember(fieldNames, expVar))=1e30;
    fill(ismember(fieldNames, 'Time'))=inf;
    fill(ismember(fieldNames, 'time'))=inf;
    fill(ismember(fieldNames, 'T'))=inf;
    fill(ismember(fieldNames, 't'))=inf;
    
    [~, fill] = sort( fill, 'descend' );
    fieldNames = fieldNames(fill);
    
    fid = fopen( outFile, 'w' );
    fprintf( fid, '%s ', fieldNames{1} );
    for jN = 2 : length( fieldNames )
        fprintf( fid, ', %s', fieldNames{jN} );
    end
    fprintf( fid, '\n' );
    for jD = 1 : length( out.(fieldNames{1}) )
        fprintf( fid, '%d ', out.(fieldNames{1}){jD} );
        for jN = 2 : length( fieldNames )
            fprintf( fid, ', %d', out.(fieldNames{jN}){jD} );
        end
        fprintf( fid, '\n' );
    end
    
    fclose(fid);
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
        if ( ~isempty( headers{a} ) )
            data.(headers{a}) = dataBlock(:,a);
        end
    end
    fclose(fid);
end